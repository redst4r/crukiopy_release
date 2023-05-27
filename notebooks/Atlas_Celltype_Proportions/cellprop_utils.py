import scanpy as sc
from tidyverse.dummy import *
from sctools import pipeline
import sys
sys.path.append('/home/mstrasse/CRUK-code/')
from sctools.de import gene_expression_to_flat_df_NEW, scanpy_DE_to_dataframe_fast, get_de_genes
from crukiopy_release.datatools import read_metadata
from sctools.scplotting import recolor, kneeplot_split
from crukiopy_release.colormaps import celltype_order_coarse_celltype, color_dict_diagnosis
import gc


def plot_patient_diagnosis(adata, patient):
    """
    plots only the cells of the patient, color by diagnosis
    """
    plt.figure()
    plt.scatter(adata.obsm['X_umap'][:,0], adata.obsm['X_umap'][:,1], s=0.5, color='grey',alpha=0.1)
    ix_patient = adata.obs.patient==patient
    for d in adata[ix_patient].obs.diagnosis.unique():
        ix_diag = adata.obs.diagnosis == d
        ix = np.logical_and(ix_patient, ix_diag)
        
        color_ix = np.where(adata.obs.diagnosis.cat.categories.values == d)[0][0]
        color = adata.uns['diagnosis_colors'][color_ix]
        plt.scatter(adata.obsm['X_umap'][ix,0], adata.obsm['X_umap'][ix,1], s=5, c=color, label=d)
    plt.legend()
import matplotlib.pyplot as plt
from scipy.stats import beta
def diagnosis_vs_cluster_plot(A, cluster, absolute_numbers=False):
    
    # stupid: the colums are categirocal! which messes up things. Hence cast them to str
    crosstab = pd.crosstab(A.obs.leiden==cluster, A.obs.diagnosis.astype(str))# [['N','N(stomach)','M','D','T']].T
    
    # fill in potentially missing diagnoses:
    for d in ['NE','NS','M','D','T']:
        if not d in crosstab.columns:
            crosstab[d] = 0
    crosstab = crosstab[['NE','NS','M','D','T']].T
    
    if not absolute_numbers:
        q = crosstab.div(crosstab.sum(axis=1), axis=0).T  # nomrlaize, such that within a diagnosis, it sums to 1
    else:
        q = crosstab.T

    _F = pd.DataFrame(q.values, index= q.index,columns=q.columns.astype(str)).reset_index().melt(id_vars='leiden')

    _F = _F.replace({'leiden': {True:cluster, False:'other'}})
    _F = _F.rename({'leiden': 'Cluster'}, axis=1)
    _F['diagnosis'] = pd.Categorical(_F['diagnosis'], ['NE','NS','M','D','T'])
    _F['Cluster'] = pd.Categorical(_F['Cluster'], ['other', cluster])

    
    # confidence intervals
#     crosstab = pd.crosstab(A.obs.leiden==cluster, A.obs.diagnosis)[['N','N(stomach)','M','D','T']].T
    df_errorbar = []

    for name, row in crosstab.iterrows():
        BETA = beta(a=row[True], b=row[False])
#         plt.plot(BETA.pdf(np.linspace(0,1,100)), label=name)
        q5,q50,q95 = BETA.ppf([0.01,0.5, 0.99])
        df_errorbar.append({
            'q5': q5,
            'q50': q50,
            'q95': q95,
            'name': name,
        })
#     plt.legend()
    df_errorbar= pd.DataFrame(df_errorbar)
#     p = None
    p = pn.ggplot() + pn.geom_col( pn.aes(x='diagnosis', y='value', fill='Cluster'), data=_F) + \
        pn.geom_errorbar(mapping=pn.aes(x='name', ymin='q5', ymax='q95'), data=df_errorbar, ) + \
        pn.labs(x='Diagnosis', y='relative cell frequency', title=f"Fraction of cells in cluster {cluster} across diagnosis") +\
        pn.theme(figure_size=(4,3))
    return p, df_errorbar, crosstab