import os
import pickle
import toolz
import pandas as pd
import scanpy as sc
from pybustools.butterfly import CUHistogram
from sctools import adata_merge
from crukiopy_release.gcsutils import download_from_gs
from crukiopy_release.kallisto import gs_kallisto_nextflow_single_sample
import pathlib

GS = 'gs://cruk_01_kallisto_nextflow_trimmed_polya'
SAMPLE_METADATA = '/home/mstrasse/CRUK-code/crukio-pipeline/sample_metadata.csv'

# hacky:
# SAMPLE_METADATA = pathlib.Path(__file__).parent.parent / 'crukio-pipeline/sample_metadata.csv'

diagnosis_fix_dict = {
    'E26C': 'M', 
    'E19D': 'NS',
}

def fix_diagnosis(adata):
    adata.obs.diagnosis = adata.obs.diagnosis.astype(str)
    for s, d in diagnosis_fix_dict.items():
        adata.obs.loc[adata.obs.samplename==s, 'diagnosis'] = d
        
    adata.obs['sample_diagnosis'] = adata.obs.samplename.astype(str) + '_' + adata.obs.diagnosis.astype(str)
    return adata

def load_cnv(patient):
    """
    loading the precomputed CNV profiles from disk
    """
    CNV_folder = '/home/mstrasse/TB4/CNV_per_patient'
    with open(f'{CNV_folder}/{patient}_cnv.pkl', 'rb') as fh:
        CNV, denoised_CNV = pickle.load(fh)
    return CNV, denoised_CNV


def read_metadata():
    """
    load the sample metadata from disk
    """
    return pd.read_csv(SAMPLE_METADATA).dropna(axis='index') # drops na rows


def patient_to_samples(patient):
    """
    translate a patient to all patients sample-ids
    """
    metadata = read_metadata()
    return list(metadata.query('patient==@patient')['samplename'].values)


def samplename2foldername(samplename):
    """
    convert samplename to foldername (which is used to store fastq,bam and kallisto data)
    """
    metadata = read_metadata()
    metadata_dod = read_DOD_metadata()

    if samplename in metadata.samplename.values:  # bugs if not done on values
        foldername = metadata.set_index('samplename').loc[samplename].foldername
    elif samplename in metadata_dod.samplename.values:
        foldername = metadata_dod.set_index('samplename').loc[samplename].foldername
    else:
        raise ValueError(f'unknown samplename {samplename} ')
    assert isinstance(foldername, str)
    return foldername

def _cellranger_adata_from_gs(samplename, nmolecules_prefilter=10, gslocation=None):
    """
    create adata from the cellranger quantified output.
    Careful, as we store the cellranger stuff in nearline, this incurs cost!
    """
    CELLRANGER_BUCKET = 'gs://cruk-data-nearline'
    assert isinstance(samplename, str)
    if gslocation is None:
        gslocation = CELLRANGER_BUCKET.strip('/')
    else:
        gslocation = gslocation.strip('/')

    foldername = samplename2foldername(samplename)

    gsloc = gslocation + '/' + foldername  + '/raw_feature_bc_matrix'
    tmp_target_name = f'/home/mstrasse/TB4/tmp/cellranger_download/{samplename}'
    os.mkdir(tmp_target_name)
    download_from_gs(gsloc, tmp_target_name)
    local_folder = tmp_target_name
    A = sc.read_10x_mtx(local_folder + '/raw_feature_bc_matrix/')
    A.obs['samplename'] = samplename
    sc.pp.filter_cells(A, min_counts=nmolecules_prefilter)

    return A


def cellranger_adata_from_gs(samplenames, nmolecules_prefilter=10, gslocation=None):
    """
    create adata from the cellranger quantified output.
    Careful, as we store the cellranger stuff in nearline, this incurs cost!
    """
    adatas = []
    for s in samplenames:
        A = _cellranger_adata_from_gs(s, nmolecules_prefilter, gslocation)
        adatas.append(A)

    adatas = adata_merge(adatas, security_check=False)


def adata_from_gs(samplenames, nmolecules_prefilter=10, gslocation=None):
    """
    download the kallisto-data from gcs, and turn into adata
    :param nmolecules_prefilter: to keep mem size small, each sample will be prefiltered to CBs with >threshold UMIs
    """
    assert isinstance(samplenames, list)
    if gslocation is None:
        gslocation = GS.strip('/')
    else:
        gslocation = gslocation.strip('/')

    adatas = []

    for s in samplenames:
        foldername = samplename2foldername(s)
        print('downloading', s)
        A = gs_kallisto_nextflow_single_sample(gslocation+'/'+foldername, TEMPDIR='/home/mstrasse/TB4/tmp')
        A.obs['samplename'] = s
        sc.pp.filter_cells(A, min_counts=nmolecules_prefilter)
        adatas.append(A)
    adatas = adata_merge(adatas, security_check=False)

    return adatas


def CU_from_samplename(samplename, per_cell=False, DB_folder = '/home/mstrasse/TB4/butterfly_all'):
    """
    loading precomputed butterfly histograms for a sample
    """
    foldername = samplename2foldername(samplename)
    pkl = f'{DB_folder}/{foldername}.pkl' if not per_cell else f'{DB_folder}/{foldername}_CB.pkl'
    if not os.path.exists(pkl):
        print(f"{pkl} doenst exist")
        raise ValueError('no CU computed for', samplename)
    with open(pkl, 'rb') as fh:
        CU = pickle.load(fh)
    CU = toolz.valmap(lambda x: x if isinstance(x, CUHistogram) else CUHistogram(x), CU)
    return CU


def download_bus(samplename, gslocation=None, postfix_busoutput='/kallisto/sort_bus/bus_output'):
    """
    downloads the busfolder from GCS.

    :param gslocation: optional; which GS-bucket to use (defautls to gs://cruk_01_kallisto_nextflow_trimmed_polya). Can be useful for alternative bus-quntifications

    :returns: the temporary local folder which it was downloaded to
    """
    if gslocation is None:
        gslocation = GS.strip('/')
    else:
        gslocation = gslocation.strip('/')

    tmp_target_name = f'/home/mstrasse/TB4/tmp/bus_download/{samplename}'
    os.mkdir(tmp_target_name)
    gsloc = gslocation + '/' + samplename2foldername(samplename)  + postfix_busoutput
    download_from_gs(gsloc, tmp_target_name)
    local_folder = tmp_target_name
    return local_folder

def download_R1_fastqs(samplename, local_folder):

    assert pathlib.Path(local_folder).exists()
    foldername = samplename2foldername(samplename)
    gsloc = f'gs://cruk-fastq/{foldername}/**/*_R1_*.fastq.gz'
    download_from_gs(gsloc, local_folder)
    return local_folder


def download_R2_fastqs(samplename, local_folder):

    assert pathlib.Path(local_folder).exists()
    foldername = samplename2foldername(samplename)
    gsloc = f'gs://cruk-fastq/{foldername}/**/*_R2_*.fastq.gz'
    download_from_gs(gsloc, local_folder)
    return local_folder


def create_samplefile(sample_dict:dict, outfile:str):
    """
    sample_dict: sampleD -> gs://location
    gs_location: where to foind the kallisto output in gcs
    """
    with open(outfile, 'w') as fh:
        for sID, gs in sample_dict.items():
            fh.write(f'{sID} {gs}\n')
