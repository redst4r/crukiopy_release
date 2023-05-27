import pandas as pd
from anndata import AnnData


def load_proteomics():
    """
    load proteomics data into an AnnData
    """
    ecm_metadata = '/home/mstrasse/CRUK-code/ECM/2021_0510_tN_M_D_T_ECM_condition_setup_Sui_v2.xlsx'
    ecm_expression = '/home/mstrasse/CRUK-code/ECM/2021_0503_tN_M_D_T_ECM_protein_quant_pivot_Sui_v1.csv'
    scrnaseq_meta = '/home/mstrasse/CRUK-code/crukio-pipeline/sample_metadata.csv'
    adata_ecm = _load_proteomics(ecm_metadata, ecm_expression, scrnaseq_meta)
    return adata_ecm

def _load_proteomics(ecm_metadata, ecm_expression, scrnaseq_meta):

    df_meta_ecm = pd.read_excel(ecm_metadata)
    df_data = pd.read_csv(ecm_expression)

    df_data = df_data[df_data['PG.Genes'].apply(lambda x: isinstance(x, str))]

    metadata = ['PG.Genes', 'PG.UniProtIds', 'PG.ProteinNames', 'PG.ProteinDescriptions', 'PG.BiologicalProcess', 'PG.CellularComponent', 'PG.MolecularFunction', ]
    adata = AnnData(
        X=df_data.drop(metadata, axis=1).T.values,
        obs=pd.DataFrame(df_data.drop(metadata, axis=1).T.index, columns=['index']).set_index('index'),
        var=df_data[metadata].set_index('PG.Genes'),
    )

    # adding the ECM metadata
    adata.obs = adata.obs.reset_index().merge(df_meta_ecm, left_on='index', right_on='Label', how='left')

    # renaming the samples to canoncial CRUK notation
    def rename_cruk(x):
        if x.startswith('EN0'):
            return x
        return x.upper()[:-1].replace('_', '')
    adata.obs['samplename'] = adata.obs['CRUK Annotation'].apply(rename_cruk)

    # adding hte CRUK metadata
    df_meta = pd.read_csv(scrnaseq_meta)
    adata.obs = adata.obs.reset_index().merge(df_meta, on='samplename', how='left').set_index('index')
    adata.obs.loc[adata.obs.index.map(lambda x: x.startswith('EN0')), 'diagnosis'] = 'N(true)'
    adata.obs.loc[adata.obs.index.map(lambda x: x.startswith('EN0')), 'patient'] = ['EN01', 'EN01', 'EN02', 'EN02', 'EN03', 'EN03']
    adata.obs.loc[adata.obs.index.map(lambda x: x.startswith('EN0')), 'procedure'] = 'resection'  # according to Philippe

    return adata
