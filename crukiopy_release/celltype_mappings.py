"""
Code to categorize the cell type calls into broader classes, like Marcophages, Bcells, ...
"""

MACROPHAGES = [
    'CL:0000235', # macro
    'CL:0000775', # neutrophil
    'CL:0000451', #dendritic
    'CL:0000784', # plasmacorticoid dendtriic cell
    'CL:0000094', #granolucyte
]

TCELLS = [
    'CL:0000084',
    'CL:0000623', # NK
    'CL:0000542' # lympocyte
]

BCELLS = [
    'CL:0000236',
    'CL:0000946' # AB secreting cell
]

FIBROBLASTS = [
    'CL:0000057',
    'CL:0000192', # smooth muscle
    'CL:0000186', # myofibroblast
    'CL:0000737', # striated mucle
    'CL:0002481', # myoid cell
    'CL:0002410', #	pancreatic stellate cell
    'CL:0000077', #meosthelial
]
ENDOTHELIAL = [
    'CL:0000115'
]
MAST = [
    'CL:0000097',
]

OTHER = [
    'CL:0008001', # HSC precurso
    'CL:0000540', #neuron
    'CL:0002573', # schwann cell , neuron
    'CL:0000232', # ery
    'CL:0000351', # trophoblast
    'CL:0000646' #basal cell
]

EPITHELIAL_SQUAM = [
    'CL:0002252', #epithelial cell of esophagus
]
EPITHELIAL_COLUMNAR = [
    'CL:0002180', # mucous cell of stomach
    'CL:0000584', #'enterocyte'
    'CL:1000272', #'lung secretory cell'
    'CL:0000160', #	goblet cell
    'CL:0000082', #	epithelial cell of lung
    'CL:0000155', #	peptic cell
    'CL:0002633', #	respiratory basal cell
    'CL:0017000', #	pulmonary ionocyte
    'CL:0002079', #	pancreatic ductal cell
    'CL:0002326', #	luminal epithelial cell of mammary gland
    'CL:0000182', #	hepatocyte
    'CL:0002181', #	mucus neck cell of gastric gland
    'CL:0002178', #	epithelial cell of stomach
    'CL:0000162', #	parietal cell
    'CL:0000166', #	chromaffin cell
    'CL:0000173', #	pancreatic D cell
    'CL:0000171', #	pancreatic A cell
    'CL:0005019', #	pancreatic epsilon cell
    'CL:0002275', #	pancreatic PP cell
    'CL:0002145', #	ciliated columnar cell of tracheobronchial tree
    'CL:0002063', # type II pneumocyte  # techincally those are squamous!??
    'CL:0002062', # type I pneumocyte # techincally those are squamous!??
    'CL:0000158', #CLub cell
]

STROMA = ENDOTHELIAL+ FIBROBLASTS
IMMUNE = TCELLS + BCELLS + MACROPHAGES + MAST
EPITHELIAL = EPITHELIAL_COLUMNAR + EPITHELIAL_SQUAM

d_tissue = {
    'stromal': STROMA,
    'immune': IMMUNE,
    'epithelial_columnar': EPITHELIAL_COLUMNAR,
    'epithelial_squamous': EPITHELIAL_SQUAM,
    'other': OTHER
}

d_celltype = {
    'Macrophages': MACROPHAGES,
    'Fibroblasts': FIBROBLASTS,
    'Endothelial cells': ENDOTHELIAL,
    'Mast cells': MAST,
    'Bcells': BCELLS,
    'Tcells': TCELLS,
#     'Columnar epithelial cells': EPITHELIAL_COLUMNAR,
#     'Squamous epithelial cells': EPITHELIAL_SQUAM,
    'Epithelial cells': EPITHELIAL,
    'Other': OTHER
}


def _dict_argmax(d):
    return max(d, key=d.get)


def annotate_celltype_per_cluster(adata, cluster_field, coarse_dict):
    """
    adjust the cell types per cluster via majority vote

    :param cluster_field: attribute in adata.obs which to group cells by
    :param coarse_dict:  a mapping from coarse_celltype -> list(CL_ids)
    """
    cluster_to_ct = {}
    for cl, df in adata.obs.groupby(cluster_field):

        a = {name: df.CLid.isin(cl_list).sum() for name, cl_list in coarse_dict.items()}
        argmax_tissue = _dict_argmax(a)
        cluster_to_ct[cl] = argmax_tissue

    cts = adata.obs[cluster_field].apply(lambda x: cluster_to_ct[x])
    return cts


def annotate_epi_stroma_immune(adata, cluster_field):
    """
    annotate the broader layer of cells: Done by cluster, not by individual cells
    """
    tissue_split = annotate_celltype_per_cluster(adata, cluster_field, d_tissue)
    return tissue_split

def annotate_coarse_celltype(adata, cluster_field):
    """
    annotate the coarse celltype: Done by cluster, not by individual cells
    """
    ct_split = annotate_celltype_per_cluster(adata, cluster_field, d_celltype)
    return ct_split
