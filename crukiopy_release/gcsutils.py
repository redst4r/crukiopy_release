import subprocess
import sys
import anndata
import hashlib
import json


def download_from_gs(gs_url, target_folder):
    """
    Download from google cloud storage. This can be a single file, or an entire folder

    :param target_folder: folder where the download gets saved in
    """
    ret = subprocess.run(
        ['ionice', '-c', '2', '-n', '7', "gsutil", "-m", "cp", '-r', gs_url, target_folder],
        check=True
    )
    if ret.returncode != 0:
        print("Child was terminated by signal", ret, file=sys.stderr)


def _samplenames_to_hash(samplenames_list:list):
    """
    turn a samplelist (not neccesarily unique) into a order invariant MD5 hash
    """
    samplenames = sorted(list(set(samplenames_list)))
    md5 = hashlib.md5(json.dumps(samplenames, sort_keys=True).encode('utf-8')).hexdigest()
    return md5

def hash_from_samplefile(samplefile:str):
    """
    compute a hash for that particular samplefile.
    This has should uniquely identify which samples go into the generated adata.

    Note that the python hash()function doesnt work, it is randomized every session
    """
    samplenames = []
    with open(samplefile, 'r') as fh:
        for l in fh:
            if l.strip() == '':  # skip empty liens
                continue
            sname = l.strip().split()[0]
            samplenames.append(sname)

    return _samplenames_to_hash(samplenames)

def hash_from_adata(adata):
    """
    compute a hash for that adata.
    This has should uniquely identify which samples went into the adata,
    identified by .foldername.
    """
    s = adata.obs.foldername.unique()
    return _samplenames_to_hash(s)
