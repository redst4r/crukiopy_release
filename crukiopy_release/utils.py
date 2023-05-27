import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotnine as pn
import pathlib


def load_metadata():
    raise ValueError('deprecated. use datatools.read_metadata')


def plot_QC_metrics(adata, outdir=None, dpi=300):
    """
    warning: this is very specialized to CRUK samples, assuming the
    sample of origin being stored in "samplename"
    """
    # genes per cell
    fig1 = plt.figure(figsize=(25,5), dpi=dpi)
    plt.subplot(131)
    sns.violinplot(x=adata.obs.samplename.values, y=adata.obs['n_genes'])
    plt.ylabel('# genes')
    plt.xticks(rotation=90)

    plt.subplot(132)
    sns.violinplot(x=adata.obs.samplename.values, y=np.log10(adata.obs['n_genes']))
    plt.ylabel('log10(# genes)')
    plt.xticks(rotation=90)

    plt.subplot(133)
    sns.distplot(np.log10(adata.obs['n_genes']), kde=False, bins=100)
    plt.xlabel('log10(# n_genes)')

    # molecules per cell
    fig2 = plt.figure(figsize=(25,5), dpi=dpi)
    plt.subplot(131)
    sns.violinplot(x=adata.obs.samplename.values, y=adata.obs['n_molecules'])
    plt.ylabel('# molecules')
    plt.xticks(rotation=90)
    plt.subplot(132)
    sns.violinplot(x=adata.obs.samplename.values, y=np.log10(adata.obs['n_molecules']))
    plt.ylabel('log10(# molecules)')
    plt.xticks(rotation=90)

    plt.subplot(133)
    sns.distplot(np.log10(adata.obs['n_molecules']), kde=False, bins=100)
    plt.xlabel('log10(# molecules)')
    # plt.xlim([2.5,3])

    # Mitoconhdiral
    n_samples = len(adata.obs['samplename'].unique())
    fig3 = plt.figure(figsize=(15,5), dpi=dpi)
    plt.subplot(131)
    sns.violinplot(x='samplename', y='percent_mito', data=adata.obs)
    plt.hlines(0.05,xmin=-0.5,xmax=-.5+n_samples, color='r')
    plt.xticks(rotation=90)

    plt.subplot(132)
    sns.violinplot(x='samplename', y='percent_mito', data=adata.obs)
    plt.hlines(0.05,xmin=-0.5,xmax=-0.5+n_samples, color='r')
    plt.ylim([0,0.3])
    plt.xticks(rotation=90)

    plt.subplot(133)
    sns.distplot(adata.obs['percent_mito'], kde=False, bins=100)
    plt.xlabel('percent_mito')

    # fig4 = plt.figure(figsize=(15,5), dpi=dpi)
    # plt.subplot(121)
    # plt.scatter(adata.obs.n_molecules, adata.obs.n_genes, c=adata.obs.percent_mito, s=3, alpha=0.5, cmap=plt.cm.plasma)
    # plt.colorbar()
    # plt.xlabel('n_molecules')
    # plt.ylabel('n_genes')
    # plt.hlines(10**2.5, 0,5000)
    # plt.vlines(10**2.7, 0,1000)
    #
    # plt.subplot(122)
    # plt.scatter(adata.obs.n_molecules, adata.obs.n_genes, c=adata.obs.percent_mito, s=3, alpha=0.5, cmap=plt.cm.plasma)
    # plt.colorbar()
    # plt.xlim(0,5000)
    # plt.ylim(0,1000)
    # plt.xlabel('n_molecules')
    # plt.ylabel('n_genes')
    # plt.hlines(10**2.5, 0,5000)
    # plt.vlines(10**2.7, 0,1000)
    #
    # fig5 = plt.figure(dpi=dpi)
    # d = {string:integer for integer, string in enumerate(adata.obs['samplename'].unique())}
    # # plt.scatter(adata.obs['n_molecules'], adata.obs['n_genes'], alpha=0.5, s=5, c=[d[_] for _ in adata.obs['samplename']], cmap=plt.cm.Set1)
    #
    # plt.subplot(121)
    # sns.scatterplot('n_molecules','n_genes', hue='samplename', data=adata.obs, alpha=0.5, s=15)
    # plt.subplot(122)
    # sns.scatterplot('n_molecules','n_genes', hue='samplename', data=adata.obs, alpha=0.5, s=15)
    # plt.xlim([0,20000])
    # plt.ylim([0,2000])

    if outdir:
        for fh, fname in zip([fig1, fig2, fig3,
                              #fig4, fig5
                              ],
                             ['QC_ngenes.png', 'QC_nmolecules.png','QC_mito.png',
                              # 'QC_mitoscatter1.png', 'QC_mitoscatter2.png'
                              ]):
            fh.savefig(outdir + '/' + fname)


    return \
    pn.ggplot(adata.obs, pn.aes('n_molecules', 'n_genes' , color='percent_mito')) \
           + pn.scale_x_log10() + pn.scale_y_log10()\
           + pn.geom_point(alpha=0.5) \
           + pn.facet_wrap('~samplename')\
           + pn.theme(figure_size=(5,5) )


def pn_qc_plots(adata):

    # scatter of #UMIs vs #genes
    p1 = pn.ggplot(adata.obs, pn.aes('n_molecules', 'n_genes' , color='percent_mito')) \
       + pn.scale_x_log10() + pn.scale_y_log10()\
       + pn.geom_point(alpha=0.5) \
       + pn.facet_wrap('~samplename')\
       + pn.theme(figure_size=(5,5) )


    geom = pn.geom_violin(alpha=0.5)
    p2 = pn.ggplot(adata.obs,
                   pn.aes(y='n_molecules', x='samplename', fill='diagnosis')) \
        + pn.scale_y_log10()\
        + geom\
        + pn.theme(figure_size=(15,5))

    p3 = pn.ggplot(adata.obs,
                   pn.aes(y='n_genes', x='samplename', fill='diagnosis')) \
        + pn.scale_y_log10()\
        + geom\
        + pn.theme(figure_size=(15,5))

    p4 = pn.ggplot(adata.obs,
                   pn.aes(y='percent_mito', x='samplename', fill='diagnosis')) \
        + pn.scale_y_log10()\
        + geom\
        + pn.theme(figure_size=(15,5))

    # cell epxresing at least one transcript of that gene
    var = adata.var.copy()

    nonzero_cells = (adata.raw.X != 0).sum(0)
    nonzero_fraction = np.array(nonzero_cells).flatten() / adata.shape[0]
    zero_fraction = 1- nonzero_fraction

    # mean nonzero expression
    mean_exp = np.array(adata.raw.X.sum(0) / nonzero_cells).flatten()
    mean_exp = np.array(adata.raw.X.mean(0)).flatten()

    mean_exp[np.isnan(mean_exp)] = 0

    # poisson regime
    lam = np.logspace(-3,3,1000)
    zero_frac_poisson = np.exp(-lam)
    mean_exp_poisson = lam #not totally correct, this should be the nonzero mean
    d = pd.DataFrame({'lambda': lam, 'zero_frac': zero_frac_poisson})

    var['mean_expression']=mean_exp
    var['zero_fraction']=zero_fraction

    p5 = pn.ggplot(var, pn.aes('mean_expression', 'zero_fraction')) + pn.geom_point(alpha=0.1) + pn.scale_x_log10() + pn.geom_point(pn.aes(x='lambda', y='zero_frac'),data=d, color='red')
