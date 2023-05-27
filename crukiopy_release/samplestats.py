import pandas as pd
import tempfile
import multiprocessing as mp
from crukiopy_release.gcsutils import download_from_gs
import os
import json


def get_cruk_10x_metadata(gcs_10x_folder):
    """
    for the gs-bucket, look up the metrics_summary, download and put into a DataFrame
    """
    assert gcs_10x_folder.startswith('gs://')
    assert not gcs_10x_folder.endswith('/')

    print(gcs_10x_folder)
    metrics_file = f'{gcs_10x_folder}/metrics_summary.csv'
#     local_file = download_file_to_tmp(metrics_file)
#     df = pd.read_csv(local_file)

    with tempfile.NamedTemporaryFile() as fp:
        download_from_gs(metrics_file, fp.name)
        try:
            df = pd.read_csv(fp.name)
        except Exception as e:
            print(f'{gcs_10x_folder} not present')

    # some cosmetics on the data, which has ',' as 1000 separator and also % signs
    df = df.apply(lambda x: pd.to_numeric(x.astype(str).str.replace(',','').str.replace('%',''), errors='coerce'))
    df['samplename'] = gcs_10x_folder
    return df


def get_cruk_10x_metadata_parallel(gs_folders, n_cores=16):
    """
    for each gs-bucket, look up the metrics_summary, download and put into a DataFrame
    """
    with mp.Pool(n_cores) as pool:
        results = pool.map(get_cruk_10x_metadata, gs_folders)
    df = pd.concat(results, sort=True)
    return df


def get_kallisto_metadata(gcs_10x_folder):
    """
    for the gs-bucket, look up the metrics_summary, download and put into a DataFrame
    """
    assert gcs_10x_folder.startswith('gs://')
    assert not gcs_10x_folder.endswith('/')

    """
    download the first metrics file, which is the kallisto output ($aligned etc)
    """
    metrics_file = f'{gcs_10x_folder}/kallisto/sort_bus/bus_output/run_info.json'

    with tempfile.NamedTemporaryFile() as fp:
        download_from_gs(metrics_file, fp.name)
        with open(fp.name, 'r') as myfile:
            data = myfile.read()
        try:
            df = pd.DataFrame(json.loads(data), index=[gcs_10x_folder])
        except Exception as e:
            print(e)
            print(gcs_10x_folder)
            return None
    df = df.reset_index().rename({'index': 'foldername'}, axis=1)


    """
    here's aslo the output of bustools inspect
    """
    inspect_file = f'{gcs_10x_folder}/kallisto/bustools_metrics/bus_output.json'
    with tempfile.NamedTemporaryFile() as fp:
        download_from_gs(inspect_file, fp.name)
        with open(fp.name, 'r') as myfile:
            data = myfile.read()
        try:
            df_inspect = pd.DataFrame(json.loads(data), index=[gcs_10x_folder])
        except Exception as e:
            print(e)
            print(gcs_10x_folder)
            return None
    df_inspect = df_inspect.reset_index().rename({'index': 'foldername'}, axis=1)

    return df.merge(df_inspect, on='foldername')



def get_kallisto_metadata_parallel(gs_folders, n_cores=16):
    """
    for each gs-bucket, look up the metrics_summary, download and put into a DataFrame
    """
    with mp.Pool(n_cores) as pool:
        results = pool.map(get_kallisto_metadata, gs_folders)
    df = pd.concat(results, sort=True)
    return df
