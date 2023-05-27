import collections
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap

def rgb(c1, c2, c3):
    return mcolors.rgb2hex([c1/255, c2/255, c3/255])


"""
colors for the cell subtype clusters in cellProportion analysis
"""
_cm_colors = [
    rgb(86, 235, 211),
    rgb(42, 132, 118),
    rgb(168, 230, 103),
    rgb(12, 168, 46),
    rgb(208, 217, 174),
    rgb(112, 142, 48),
    rgb(57, 242, 122),
    rgb(3, 98, 160),
    rgb(166, 182, 249),
    rgb(69, 83, 194),
    rgb(238, 128, 254),
    rgb(149, 37, 186),
    rgb(111, 47, 95)
]
_cm_colors = _cm_colors+ ['grey']
cmap_subtypes = ListedColormap(
    _cm_colors
)


"""
colormap for diagnosis
defaults to grey
"""
color_dict_diagnosis = collections.defaultdict(
    lambda: "#6b6a6a",  #default color
    {
        "NE": "#339933",
        "NS": "#0492C2",
        "M": "#FFCC33",
        "MT": "#FFCC33",
        "MDT": "#FF6600",
        "DT": "#FF6600",
        "D": "#FF6600",
        "T": "#B22222",
    }
)


# cross tissue diagnosis colormap

cross_tissue_diag_dict = collections.defaultdict(lambda: 'grey', {
    'NE': '#051d00',
    'NS': '#00ff23',
    'NL': '#16715c',
    'NC': '#009c78',
    'ME': 'orange',
    'MS': '#ba790d',# 'yellow',
    'MC': '#aa790d',# 'yellow',
    'DL': '#124dcc',
    'DE': '#1d3871',
    'DS': '#12a2de',
    'DC': '#61dfe0',
    'TE': 'red',
    'TL': 'darkred',
    'TS': '#5c20a0',
    'TC': '#d80bd5'
    
})


"""
cell type colormap by Dave
"""
cm = plt.get_cmap('tab20', 20)
NUM_COLORS = 20
colors = [mcolors.rgb2hex(cm(i/NUM_COLORS)) for i in range(NUM_COLORS)]

color_dict_CellLabels = collections.defaultdict(
    lambda: colors[15],
    {
        # immune
        'B_cells':colors[3],
        'NK_cells':colors[12],
        'cd4_Tcells':colors[8],
        'cd8_Tcells':colors[9],
        'mast_cells':colors[6],
        'monocytes_macs_DCs':colors[7],
        'naive_T_cells':colors[13],
        'neutrophils':colors[2],

        # stroma
        'endothelial':colors[18],
        'fibroblasts':colors[1],
        'myofibroblasts':colors[19],

        # epithelial
        'neuroendocrine':colors[5],
        'gi_epithelial':colors[4],
        'hepatoid':colors[11],
        'parietal':colors[16],
        'squamous_epithelial':colors[17],

        # other
        '__':colors[14],
        'NA':colors[15]
    }
)
celltype_order_CellLabels = [
    'mast_cells',
    'monocytes_macs_DCs',
    'neutrophils',
    'B_cells',
    'NK_cells',
    'cd4_Tcells',
    'cd8_Tcells',
    'naive_T_cells',

    'endothelial',
    'fibroblasts',
    'myofibroblasts',

    'gi_epithelial',
    'parietal',
    'neuroendocrine',
    'hepatoid',
    'squamous_epithelial',

    '__',
    'NA'
]


"""
colormap for the coarse_ct from crukiopy_release.celltype_mappings.annotate_coarse_celltype
"""
color_dict_coarse_celltype = collections.defaultdict(
    lambda: '#6b6a6a',
    {
        'Macrophages': colors[7],
        'Fibroblasts': colors[1],
        'Endothelial cells': colors[18],
        'Mast cells': colors[6],
        'Bcells': colors[3],
        'Tcells': colors[8],
        #     'Columnar epithelial cells': EPITHELIAL_COLUMNAR,
        #     'Squamous epithelial cells': EPITHELIAL_SQUAM,
        'Epithelial cells': colors[4],
        'Other': colors[14]
    }
)
celltype_order_coarse_celltype = [
    'Mast cells',
    'Macrophages',
    'Bcells',
    'Tcells',
    'Endothelial cells',
    'Fibroblasts',
    'Epithelial cells',
    'Other'
]
