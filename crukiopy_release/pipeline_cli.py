import gc
import scanpy as sc
import pandas as pd
import fire
from sctools import adata_merge, pipeline
from sctools.celltypes.celltype_annotation import build_marker_dict, create_meta_markers
from sctools.celltypes.celltype_colors import celltype_colors
from crukiopy_release.kallisto import gs_kallisto_nextflow_single_sample
import scHCLpy.adata
from crukiopy_release.gcsutils import hash_from_adata
import multiprocessing


"""
functions used in the nextflow pipeline to turn kallisto quantifications
into AnnData objects (scanpy)
exposes the essentials as CLI

python pipeline-cli.py download ...
python pipeline-cli.py concat ...
python pipeline-cli.py download ...
"""


def download_kallisto_counts(sampleID: str, gslocation: str, metadatafile: str, outputfile: str, TEMPDIR: str, min_counts:int=1000, min_genes:int=0, percentmito:float=1):
    """
    download kallisto quantification, turn into adata, add metadata,
    do some minimal preprocessing (cell detection, doublet detection) save to disk
    """
    # annotating metadata
    sample_metadata = pd.read_csv(metadatafile)
    current_metadata = sample_metadata.query('foldername==@sampleID')
    assert len(current_metadata) == 1, f"found more/less than one entry in metadata for sample {sampleID}: #entries {len(current_metadata)}"

    # downloading
    adata = gs_kallisto_nextflow_single_sample(f'{gslocation}', TEMPDIR=TEMPDIR, verbose=True)


    # for col in ['foldername', 'samplename', 'patient', 'tissue', 'procedure', 'diagnosis', 'treatment']:
    for col in current_metadata.columns:
        # dirty: values[0] since is still a dataframe (1row)
        adata.obs[col] = current_metadata[col].values[0]

    # prefilter the cells on QC metrics
    # -> makes the later merges etc much smaller
    sc.pp.filter_cells(adata, min_counts=min_counts)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    ix_mito = adata.obs.percent_mito < percentmito
    adata._inplace_subset_obs(ix_mito)
    
    if len(adata) == 0: # if no cells remain after filter, annotate_doublets would produce an error
        print(f'WARNING: {sampleID} had no cells left with umi>{min_counts}')
        adata.obs['doublet_score'] = []
    else:
        adata = pipeline.annotate_doublets(adata, groupby='samplename')

    adata.write_h5ad(outputfile)


def concat_adatas(outputfile: str, *files):
    """
    merge all adata files into a big one.
    """
    adata = None
    for f in files:
        _tmp = sc.read_h5ad(f)
        assert 'samplename' in _tmp.obs.columns, "samplename not set"
        print(f'concat: {f}')
        if not adata:
            adata = _tmp
        else:
            # adata = adata_merge([adata, _tmp])
            adata = adata_merge([adata, _tmp], security_check=False)

            gc.collect()
    # fix the CB index
    adata.obs['CB'] = adata.obs.index.map(lambda x: x.split('-')[0])
    adata.obs.index = adata.obs['CB'] + '-' + adata.obs['samplename'].astype(str)
    adata.obs_names_make_unique()  # just to be on the save side
    adata.write_h5ad(outputfile)


def processing(adatafile: str, outputfile: str, percentmito: float, min_counts=1000, min_genes=500):
    adata = sc.read_h5ad(adatafile)

    the_hash = hash_from_adata(adata)  # hash to remember which samples went into the adata
    adata = adata_processing(adata, percentmito, min_counts, min_genes)
    adata.write_h5ad(outputfile)

    hashfile = outputfile + '.hash'
    with open(hashfile, 'w') as fh:
        fh.write(f'{the_hash}\n')


def _annotation(adata):
    # HCL cell type annotations
    dfadata, _ = scHCLpy.adata.scHCL_adata(adata, verbose=True, n_cores=multiprocessing.cpu_count())

    # get rid of those if the are already present
    for f in ['hcl_score', 'hcl_celltype', 'hcl_refined', 'CLid', 'CL_name']:
        if f in adata.obs.columns:
            adata.obs.drop(f, axis=1, inplace=True)

    adata.obs = adata.obs.merge(dfadata, left_index=True, right_index=True)

    # TODO: add coarse celltypes (based on leiden or nobatch_leiden)

    return adata


def adata_processing(adata, percentmito, min_counts=1000, min_genes=500):
    """
    processing a merged adata:
    """
    
#     sname = "_".join(sorted(adata.obs['samplename'].unique()))
#     adata.obs.to_csv(f'/tmp/{sname}.csv')

    # filter cells with little genes (careful this might remove erys and B-cells
    sc.pp.filter_cells(adata, min_genes=min_genes)
    gc.collect()

    # filter doublets
    ix_no_doublet = adata.obs['doublet_score'] < 0.2
    adata._inplace_subset_obs(ix_no_doublet)
    gc.collect()

    # processing and batch corrction
    adata = pipeline.standard_processing(adata, MITO_CUTOFF=percentmito, detect_doublets=False, UMI_CUTOFF=min_counts)

    adata.obs['sample_diagnosis'] = adata.obs.apply(lambda row: row['samplename'] + '_' + str(row['diagnosis']), axis=1)  #str(diagnosis) in case its NaN

    # # HCL cell type annotations
    adata = _annotation(adata)

    # only keep a few annotations
    ann = ['samplename', 'n_genes', 'n_molecules', 'doublet_score',
           'percent_mito', 'leiden', 'celltypes_auto',
           'diagnosis', 'phase', 'sample_diagnosis', 'patient', 'treatment', 'procedure',
           'hcl_refined', 'hcl_celltype', 'hcl_score',
           'CLid',	'CL_name',
           'nobatch_leiden', 'nobatch_louvain']
    _tmp = pipeline.export_for_cellxgene(adata, [_ for _ in ann if _ in adata.obs.columns])

    return _tmp


def _old_celltype_calling(adata):
    """
    old cell tye calling based on a dictinoary of a few markers and some
    crude scoring based on classification accuracy
    """
    # celltypes
    tissue = adata.obs.tissue.values
    assert all(tissue == tissue[0]), "different tissues merged; not clear how to do celltype annotation"
    tissue = tissue[0]

    print(tissue)
    celltype_marker_dict = build_marker_dict(adata, tissue=tissue)
    df_cluster_celltype, cluster_celltype_mapping = create_meta_markers(adata, celltype_marker_dict, cluster_name='leiden',auc_cutoff=0.85 , SIMPLE=True)
    adata.uns["marker_based_celltype_colors"] = [celltype_colors[_] if _ in celltype_colors  else 'grey' for _ in adata.obs.marker_based_celltype.astype('category').cat.categories.values ]
    adata.obs['celltypes_auto'] = adata.obs['marker_based_celltype']
    return adata


if __name__ == '__main__':
    fire.Fire({
        'download': download_kallisto_counts,
        'concat': concat_adatas,
        'processing': processing
    })
