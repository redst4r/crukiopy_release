import pandas as pd
import scanpy as sc
import numpy as np
import pathlib
import tqdm

coarse_CT_dict = {'B cell': ['B cell'],
 'Thelper':['CD4+ T cell',  'CD4+ T cell PD1+'],
 'Treg': ['CD4+ Treg'],
 'cytoTcell': ['CD8+ T cell', 'CD8+ T cell PD1+'],
 'Neutrophil': ['Neutrophil'],
 'DC': ['DC'],
 'Macrophage': ['M1 Macrophage', 'M2 Macrophage'],
 'Plasma': ['Plasma'],

 'SquamousEpi': ['Squamou p63+ EGFRhi', 'Squamous Annexin A1+', 'Squamous p63+'],

 'GastricEpi': [
    'Chief',
    'Foveloar',  'Foveloar Ki67+ p53+', 'Foveloar p53+',
    'Parietal',
    'Neck', 'Neck Ki67+ p53+', 'Neck p53+',
    'Neuroendocrine',
  ],
 'IntestinalEpi': [
     'Goblet', 'Goblet Ki67+ p53+', 'Goblet p53+','Paneth'
     ],

 'OtherEpi': [
     'Epithelial', 'Epithelial CD73hi', 'Epithelial CK7+',
     'Epithelial CK7+ p53+', 'Epithelial HLADR+',
     'Epithelial Ki67+ p53+', 'Epithelial MUC1+ Ki67+',
     'Epithelial p53+', 'Epithelial pH2AX+',
    ],

  'Nerve': ['Nerve'],

  'Lymph': ['Lymphatic', 'Lymphatic CD73+'],

  'Endothelial': ['Endothelial', 'Endothelial CD36hi', 'Endothelial aSMAhi'],
  'SMC': ['Smooth Muscle'],
  'Stromal': ['Stroma', 'Stroma CD73+'],
}

even_coarser_dict = {
    'Immune': ['B cell',
               'CD4+ T cell',  'CD4+ T cell PD1+','CD4+ Treg', 'CD8+ T cell', 'CD8+ T cell PD1+',
               'Neutrophil', 'DC', 'M1 Macrophage', 'M2 Macrophage','Plasma'],

 'SquamousEpi': ['Squamou p63+ EGFRhi', 'Squamous Annexin A1+', 'Squamous p63+'],

 'GastricEpi': [
    'Chief',
    'Foveloar',  'Foveloar Ki67+ p53+', 'Foveloar p53+',
    'Parietal',
    'Neck', 'Neck Ki67+ p53+', 'Neck p53+',
    'Neuroendocrine',
  ],
 'IntestinalEpi': [   'Goblet', 'Goblet Ki67+ p53+', 'Goblet p53+',
    'Paneth'],

 'OtherEpi': [
     'Epithelial', 'Epithelial CD73hi', 'Epithelial CK7+',
     'Epithelial CK7+ p53+', 'Epithelial HLADR+',
     'Epithelial Ki67+ p53+', 'Epithelial MUC1+ Ki67+',
     'Epithelial p53+', 'Epithelial pH2AX+',
    ],

  'Nerve': ['Nerve'],

  'Stroma': ['Lymphatic', 'Lymphatic CD73+', 'Endothelial', 'Endothelial CD36hi', 'Endothelial aSMAhi', 'Smooth Muscle', 'Stroma', 'Stroma CD73+'],
 }

def invert_dict(d):
    new_dict = {}
    for k, v_list in d.items():
        for v in v_list:
            assert v not in new_dict
            new_dict[v] = k
    return new_dict

INV_coarse_CT_dict = invert_dict(coarse_CT_dict)
INV_coarse_CT_dict2 = invert_dict(even_coarser_dict)

matt_diagnosis_dict = [
    #rediagnosisMatt    smaplename    #proc      # initial diag        # initial short
    ['T', 'E08_reg001', 'Biopsy', 	'Tumor + Dysplasia', 'T'],
    ['NE','E08_reg002', 'Biopsy', 	'Normal esophagus', 'NE'],
    ['NE','E08_reg003', 'Biopsy', 	'Normal esophagus', 'NE'],
    ['T', 'E08_reg004', 'Biopsy', 	'Tumor + Metaplasia', 'T'],
    ['T', 'E08_reg005', 'Biopsy', 	'Tumor + Metaplasia', 'T'],
    ['M', 'E08_reg006', 'Biopsy', 	'Metaplasia', 'M'],
    ['M', 'E08_reg007', 'Biopsy', 	'Metaplasia', 'M'],
    ['T', 'E11_reg001', 'Biopsy', 	'Tumor + Metaplasia', 'T'],
    ['T', 'E11_reg002', 'Biopsy', 	'Tumor + Metaplasia', 'T'],
    ['T', 'E11_reg003', 'Biopsy', 	'Metaplasia + Tumor', 'T'],
    ['D', 'E11_reg004', 'Biopsy', 	'Metaplasia', 'M'],
    ['?', 'E11_reg005', 'Biopsy', 	'Normal esophagus; Tumor?', 'NE/T'],
    ['?', 'E11_reg006', 'Biopsy', 	'Normal esophagus; Tumor?', 'NE/T'],
    ['T', 'E18_reg001', 'Resection', 	'Metaplasia', 'M'],
    ['T', 'E18_reg002', 'Resection', 	'Metaplasia', 'M'],
    ['T', 'E18_reg003', 'Resection', 	'Tumor', 'T'],
    ['T', 'E18_reg004', 'Resection', 	'Tumor', 'T'],
    ['NE', 'E17_reg001', 'Biopsy', 	'Normal esophagus', 'NE'],
    ['D', 'E17_reg002', 'Biopsy', 	'Metaplasia; Dysplasia?', 'MD'],
    ['M', 'E17_reg003', 'Biopsy', 	'Metaplasia; Dysplasia?', 'MD'],
    ['M', 'E17_reg004', 'Biopsy', 	'Tumor', 'T'],
    ['T', 'E17_reg005', 'Biopsy', 	'Metaplasia', 'M'],
    ['M', 'E17_reg006', 'Biopsy', 	'Metaplasia', 'M'],
    ['NS', 'E19_reg001', 'Resection', 	'Normal esophagus', 'NE'],  #prob NS
    ['T', 'E19_reg002', 'Resection', 	'Normal esophagus', 'NE'],
    ['T', 'E19_reg003', 'Resection', 	'Normal esophagus', 'NE'],
    ['T', 'E19_reg004', 'Resection', 	'Normal esophagus', 'NE'],
    ['T','E12_reg001', 'Resection', 	'Metaplasia + Dysplasia + Tumor', 'T'],
    ['M','E12_reg002', 'Resection', 	'Metaplasia + Dysplasia + Tumor', 'T'],
    ['NS', 'E12_reg003', 'Resection', 	'Metaplasia + Dysplasia + Tumor', 'T'],
    ['NE/NS',  'E12_reg004', 'Resection', 	'Metaplasia + Dysplasia + Tumor', 'T'],
    ['NE/NS', 'E12_reg005', 'Resection', 	'Metaplasia + Dysplasia + Tumor', 'T'],
]
matt_diagnosis_dict = pd.DataFrame(matt_diagnosis_dict, columns=['rediagnosis', 'samplename', 'procedure', 'matt_diagnosis', 'diagnosis'])
matt_diagnosis_dict['patient'] = matt_diagnosis_dict.samplename.apply(lambda x: x.split('_')[0])



def load_sample(samplename, codex_adata_path='/home/mstrasse/TB4/CODEX_adata/'):
    """
    load a single codex sample from disk by its name, e.g. 'E08_reg001'
    annnotates coarse celltype
    """
    adata = sc.read_h5ad(pathlib.Path(codex_adata_path) / f'{samplename}.h5ad' )
    adata.obs['samplename'] = samplename
    adata.obs['cellid'] = adata.obs.index + '_' + samplename
    adata = _codex_postprocess(adata)
    return adata

def _codex_postprocess(adata):
    """
    annotate some coarse cell types
    """
    adata.obs['MS_CL'] = adata.obs['Cell Type'].apply(lambda x: INV_coarse_CT_dict[x])
    adata.obs['MS_CL_coarse'] = adata.obs['Cell Type'].apply(lambda x: INV_coarse_CT_dict2[x])
    adata.obs['MS_CL_coarse'] = pd.Categorical(
                                    adata.obs['MS_CL_coarse'],
                                    categories=['SquamousEpi', 'GastricEpi',
                                                'IntestinalEpi', 'OtherEpi',
                                                'Nerve', 'Stroma', 'Immune'])
    return adata


def iterate_codex_adatas(codex_path="/home/mstrasse/TB4/CODEX_adata/"):
    """
    iterating over all allvaiable codex samples
    """
    for f in tqdm.tqdm(pathlib.Path(codex_path).iterdir()):
        adata = sc.read_h5ad(f)
        sname = f.stem
        adata.obs['samplename'] = sname
        adata.obs['cellid'] = adata.obs.index + '_' + sname
        yield sname, _codex_postprocess(adata)

