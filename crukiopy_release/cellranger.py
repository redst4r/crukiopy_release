import tempfile
import scanpy as sc
import os
from sctools import adata_merge, annotate_qc_metrics


def gs_to_adata(gs_locations:list, samle_names:list, targeth5name, cellranger_version=3):
    """
    converts all the locations on google cloud into a single! adata, tagging
    each sample with the given sample name
    """
    assert targeth5name.endswith('.h5ad'), 'use .h5ad as file extension'
    adatas = []
    for samplename, gs in zip(samle_names, gs_locations):

        with tempfile.NamedTemporaryFile(suffix='.h5') as fh:
            tmp_target_name = fh.name
            h5name = 'filtered_feature_bc_matrix.h5' if cellranger_version == 3 else 'filtered_gene_bc_matrices_h5.h5'
            cmd = f'gsutil cp {gs}/{h5name} {tmp_target_name}'
            os.system(cmd)
            tmp_ = sc.read_10x_h5(tmp_target_name)
            tmp_.obs['samplename'] = samplename
            tmp_.var_names_make_unique()
            adatas.append(tmp_)

    
    # Problem: they all have different sets of genes.
    # gather a list of all
    # and fill the missing ones as seros
#     all_genes = set()
#     for a in adatas:
#         all_genes = all_genes | set(a.var.index)
#     all_genes = sorted(list(all_genes))
    
#     all_genes_df = pd.concat([a.var for a in adatas]).drop_duplicates()
    
#     actually: easier: just take the ones present in all!
    shared_genes = set(adatas[0].var.index)
    for a in adatas:
        shared_genes = shared_genes & set(a.var.index)
    shared_genes = sorted(list(shared_genes))
    adatas = [a[:, shared_genes].copy() for a in adatas] # we get a View if we dont copy
    
    # TODO enable the check; but 10x dataformat is SOOO messy, just skip that
    Q = adata_merge(adatas, security_check=False)
    Q.raw = Q.copy()
    Q = annotate_qc_metrics(Q)
    # Q = annotate_gene_metadata(Q)

    Q.write_h5ad(targeth5name)
