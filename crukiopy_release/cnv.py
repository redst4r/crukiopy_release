from cnvpy.my_inferCNV import preprocess, inferCNV, denoising
from cnvpy.plotting import diagnosis_cmap
import pickle
import numpy as np
from sctools import adata_merge


def do_cnv(adata, patient, outfolder):
    NS = 'NS'
    ix_sub = adata.obs.query('diagnosis==@NS or patient==@patient').index
    Q = adata[ix_sub].copy()
    Qlog = preprocess(Q)

    ref_groups = Qlog.obs.query('diagnosis == @NS').sample_diagnosis.unique().astype(str)
    print(ref_groups)

    CNV = inferCNV()
    CNV.infer(Qlog, ref_field='sample_diagnosis', ref_groups=ref_groups)
    print('clustering')
    CNV.cluster()
    denoised_CNV = denoising(CNV)
    denoised_CNV.cluster()


    with open(f'{outfolder}/{patient}_cnv.pkl', 'wb') as fh:
        pickle.dump((CNV, denoised_CNV), fh)
    import gc
    del Q, Qlog

    gc.collect()



def do_cnv_against_reference(adata, reference_adata, outname):
    """
    run infercnv again the reference adata
    """
    assert np.all(reference_adata.obs.diagnosis.isin(['NE','NS']))

    adata.obs['cnv_group'] = 'query'
    reference_adata.obs['cnv_group'] = 'reference'

    Q = adata_merge([adata, reference_adata],security_check=False)
    Qlog = preprocess(Q)


    CNV = inferCNV()
    CNV.infer(Qlog, ref_field='cnv_group', ref_groups=['reference'])
    import gc
    del Q, Qlog
    gc.collect()
    print('clustering')
    CNV.cluster()
    denoised_CNV = denoising(CNV)
    denoised_CNV.cluster()


    with open(outname, 'wb') as fh:
        pickle.dump((CNV, denoised_CNV), fh)

    del denoised_CNV, CNV
    gc.collect()
