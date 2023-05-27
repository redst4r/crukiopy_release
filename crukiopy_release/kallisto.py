import tempfile
import os
from sctools import adata_merge, annotate_qc_metrics
from crukiopy_release.gcsutils import download_from_gs
from sctools.kallisto import _load_kallisto_to_adata


def gs_kallisto_nextflow_single_sample(gs, genecounts=True, verbose=False, TEMPDIR=None):
    """
    download a single kallisto quantification from gs
    expects the file structure on gs to be just as nextflow created it!
    i.e `gs` must contain
    |
    |- kallisto
          |
          |- bustools_counts
          |  |
          |  |- bus_output_eqcount --  tcc.barcodes.txt, tcc.ec.txt tcc.mtx
          |  |
          |  |- bus_output_genecounts -- gene.barcodes.txt, gene.ec.txt gene.mtx
          |
          |- bustools_metrics
               |
               |- bus_output.json

    """
    if TEMPDIR:
        assert os.path.exists(TEMPDIR), f'custom TEMPDIR {TEMPDIR} doesnt exist'

    with tempfile.TemporaryDirectory(dir=TEMPDIR) as tmp_target_name:
        # special care with the gs-path: it must not contain a trailing /
        gs = gs.strip('/')

        if verbose:
            print(f'Downloading from {gs} to {tmp_target_name}')

        # actually, only download specific folders
        # gs_metrics = f'{gs}/kallisto/bustools_metrics'

        if genecounts:
            gs_count = f'{gs}/kallisto/bustools_counts/bus_output_genecount'
        else:
            gs_count = f'{gs}/kallisto/bustools_counts/bus_output_eqcount'

        # download_from_gs(gs_metrics, tmp_target_name)
        download_from_gs(gs_count, tmp_target_name)

        if verbose:
            print('DONE Downloading')

        metricsfile = f'{tmp_target_name}/bustools_metrics/bus_output.json'
        if not os.path.exists(metricsfile):
            metricsfile = None

        if genecounts:
            matrixfile = f'{tmp_target_name}/bus_output_genecount/gene.mtx'
            genefile = f'{tmp_target_name}/bus_output_genecount/gene.genes.txt'
            bcfile = f'{tmp_target_name}/bus_output_genecount/gene.barcodes.txt'
        else:
            matrixfile = f'{tmp_target_name}/bus_output_eqcount/tcc.mtx'
            genefile = f'{tmp_target_name}/bus_output_eqcount/tcc.genes.txt'
            bcfile = f'{tmp_target_name}/bus_output_eqcount/tcc.barcodes.txt'

        _tmp = _load_kallisto_to_adata(matrixfile, genefile, bcfile, metricsfile)

        _tmp = _kallisto_post_downloading(_tmp)
        return _tmp


def gs_kallisto_nextflow(gs_locations: dict, targeth5name: str, genecounts=True, verbose=True):
    """
    downloads a set of gs-kallisto quantifications and stores them in a single adata
    gs_locations is a dict of samplename->gs_location
    each cell in the adata gets the correct samplename
    """
    assert targeth5name.endswith('.h5ad'), 'use .h5ad as file extension'
    adatas = []

    if verbose:
        print('downloading samples')
    for samplename, gs in gs_locations.items():
        _tmp = gs_kallisto_nextflow_single_sample(gs, genecounts, verbose)
        _tmp.obs['samplename'] = samplename
        adatas.append(_tmp)

    print('merging adatas')
    adatas = adata_merge(adatas)
    print('done merging adatas')

    adatas.uns['kallisto_sources'] = gs_locations

    print('write_h5ad')
    adatas.write_h5ad(targeth5name)


def _kallisto_post_downloading(adata):

    # adata.var_names_make_unique()
    # kallisto has the ensemble_id ('ENSG00000241950.1') as index
    # whereas cellranger has the gene symbol
    # hence convert this to symbol

    # TODO: THIS SHOULD BE REPLACED BY sctools.kallisto.annotate_gene_symbols
    # print('loading biomart')
    # df_biomart = biomart_query_all()
    #
    # "spliting away the version number"
    # ensebmle_ids = [_.split('.')[0] for _ in adata.var.index]
    # assert len(ensebmle_ids) == len(set(ensebmle_ids)), "ensemble ids not unique"
    # adata.var.index = adata.var.index.map(lambda x: x.split('.')[0])
    # adata.var.index.name = 'gene_ids'
    #
    # # annotate the symbol: problem: there's multiple entries in df_biomart
    # # per gene_ids (diff transcripts)
    # # they should however have the same gene symbol
    #
    # geneid_2_symbol = pd.DataFrame(df_biomart.groupby('ensembl_gene_id').symbol.agg(symbol=lambda x: sorted(list(set(x)))[0]))  # sorted() is important otherwise we might get different members of the set each time
    # adata.var = adata.var.reset_index()
    #
    # # also some genes might not have symbols in biomart
    # ind = []
    # for gid in adata.var['gene_ids']:
    #     if gid in geneid_2_symbol.index:
    #         symbol = geneid_2_symbol.symbol.loc[gid]
    #     else:
    #         symbol = gid
    #     ind.append(symbol)
    # adata.var.index = ind
    # adata.var.index.name = 'index'
    # adata.var_names_make_unique()

    from sctools.kallisto import annotate_gene_symbols
    adata = annotate_gene_symbols(adata)

    # some features cellranger puts in
    adata.var['feature_types'] = 'Gene Expression'
    adata.var['genome'] = 'GRCh38'

    adata.raw = adata.copy()

    print('annotate_qc_metrics')
    adata = annotate_qc_metrics(adata)

    return adata
