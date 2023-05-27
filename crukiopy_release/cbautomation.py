import subprocess as sb
import logging
import pathlib
import re
import tempfile
import toolz
import shutil
import os
import hashlib
from crukiopy_release.gcsutils import download_from_gs
import scanpy as sc


# standard config with the usual batch and nonbatch corrected coordinates
cb_config_CRUK = {
    'coords': [
        {"file": "umap_coords.tsv", "shortLabel": "UMAP"},
        {"file": "pca_coords.tsv", "shortLabel": "pca"},
        {"file": "umap_nobatch_coords.tsv", "shortLabel": "UMAP_nobatch"},
        {"file": "pca_original_coords.tsv", "shortLabel": "pca_nobatch" }
    ]
}

logging.basicConfig(level=logging.WARN)

def _check_clean_name(name):
    """
    dataset names cant contain funky characters, like /
    """
    return bool(re.match('^[a-zA-Z0-9\\-\\_]+$', name))

def _cleanup_description(desc):
    """
    cant contain any funny characters either
    """
    return desc.replace("'", r"\'").replace('"', r'\"')


class CB_Collection():
    """
    A collection of datasets in Cellbrowser
    """
    def __init__(self, name: str, description: str, subsample_for_testing=False, hideDownload=False):

        name = name.replace(' ', '_')
        assert _check_clean_name(name), f"{name} contains forbidden characters (e.g whitespaces)" # no withspaces or other funny stuff
        self.name = name
        self.description = _cleanup_description(description)
        self.datasets = []
        self.subsample_for_testing = subsample_for_testing
        self.hideDownload= hideDownload

    def _add_dataset(self, d):
        self.datasets.append(d)

    def do_import(self, target_dir, postprocess_fn=toolz.functoolz.identity):
        """
        first step of the generation of the cellbrowser

        :param postprocess_fn: function f(h5_filename). must take a path to a single h5file and return the path to a single h5file (which can be modified)
        """
        assert self.check_existance(), "certain datasets cant be found on disk!"

        collection_dir = f'{target_dir}/{self.name}' # where the content of this collection goes to
        for d in self.datasets:
            d.do_import(collection_dir, postprocess_fn)
        self._write_cellbrowser_conf(collection_dir)
        self._write_desc_conf(collection_dir)


    def do_build(self, rootdir, target_dir, html_outputdir):
        """
        second step of the generation of the cellbrowser
        """
        cmd = f'cd {target_dir}; CBDATAROOT="{rootdir}" cbBuild -r -o {html_outputdir}'
        print('running', cmd)
        sb.run(cmd, shell=True, check=True)

        # call build on each dataset recursively
        for d in self.datasets:
            if isinstance(d, CB_Collection):
                d.do_build(rootdir, f'{target_dir}/{d.name}', html_outputdir)

    def pretty_print(self, indent=0):
        """
        pretty print the collection
        """
        tabs = "".join([' ']*indent)
        print(f'{tabs}Title {self.name}\n{tabs}Datasets')
        for d in self.datasets:
            d.pretty_print(indent=indent+4)

        print(f'{tabs}---')

    def _write_cellbrowser_conf(self, target_dir):
        """write the cellbrowser.conf for this collection"""
        with open(f'{target_dir}/cellbrowser.conf', 'w') as fh:
            fh.writelines([
                f'name="{self.name}"\n',
                f'shortLabel="{self.name}"\n'
            ])

    def _write_desc_conf(self, target_dir):
        with open(f'{target_dir}/desc.conf', 'w') as fh:
            fh.writelines([
                f'title="{self.name}"\n',
                f'abstract="{self.description}"\n',
                f'hideDownload="true"\n' if self.hideDownload else f'hideDownload="false"'
            ])

    def check_existance(self):
        """
        check that each dataset in this collection (h5ad file) exists
        """
        flags = []
        for ds in self.datasets:
            f = ds.check_existance()
            flags.append(f)
            if not f:
                logging.info(f'{ds.name} doesnt exist')
                if isinstance(ds, CB_Dataset):
                    logging.info(f'missing {ds.h5ad}')
        return all(flags)

N_SUBSAMPLE = 100
class CB_Dataset():
    """
    A single dataset (i.e. a single h5ad object)
    """
    def __init__(self, name:str, description:str, h5ad:str, config:dict, tmpdir: str, subsample_for_testing=False, hideDownload=False):
        name = name.replace(' ', '_')
        assert _check_clean_name(name), f"{name} contains forbidden characters (e.g whitespaces)" # no withspaces or other funny stuff
        self.name = name
#         self.name = name.replace(' ', '_')

        self.description = _cleanup_description(description)
        self.h5ad = h5ad
        self.config = config
        self.tmpdir = tmpdir
        self.subsample_for_testing = subsample_for_testing
        self.hideDownload = hideDownload
    def is_remote(self):
        """
        returns true if the h5ad is a path to a gs bucket
        """
        return self.h5ad.startswith('gs://')

    def do_import(self, target_dir, postprocess_fn):
        """
        calling cbImportScanpy on that dataset

        we create a copy of the .h5ad file in a temporary folder (either local copy or gs:// download
        and call the postprocess_dn on that copy (wont modify the original) and finally build the cbImport

        :param postprocess_fn: function f(h5_filename). must take a path to a single h5file and return the path to a single h5file (which can be modified)
        """

        name = self.name
        fname = self.h5ad.split('/')[-1]
        dataset_dir = f'{target_dir}/{name}'

        job_completed_file = f'{dataset_dir}/completed.txt'

        if os.path.exists(job_completed_file):
            print(f'already done {name} {self.h5ad}')
            return
        else:
            print(f'building {name} {self.h5ad}')

        with tempfile.TemporaryDirectory(prefix='cellbrowser_download',dir=self.tmpdir) as tmpfolder:

            if self.is_remote():
                #download to tempdir
                logging.info(f'downloading {self.h5ad} into {tmpfolder}')
                download_from_gs(self.h5ad, tmpfolder)
                tmp_h5ad_path = tmpfolder + '/' + fname
            else:
                #local copy to tempdir
                logging.info(f'copying {self.h5ad} into {tmpfolder}')

                if False:
                    # make a symlink instead of a copy
                    # TODO this is risky, the postprocessing could modify the symlinked file and hence change the original
                    p = pathlib.Path(self.h5ad)
                    tmp_h5ad_path = tmpfolder + '/' + p.name
                    p.symlink_to(tmp_h5ad_path)

                else:  # create a full copy
                    shutil.copy(self.h5ad, tmpfolder)
                    tmp_h5ad_path = tmpfolder + '/' + fname

            if self.subsample_for_testing:
                print('subsampling')
                a = sc.read_h5ad(tmp_h5ad_path)
                import numpy as np
                ix = np.random.choice(a.obs.index, np.minimum(N_SUBSAMPLE, a.shape[0]))
                b = a[ix]
                b.write_h5ad(tmp_h5ad_path)

            # renaming .var!
            # there's some issues: when a column is named gene_ids
            # cellbrowser starts pulling in EnsemblIDs and thigns get messed up
            # hence: if gene_ids is present rename it to ensembl_gene_id
            a = sc.read_h5ad(tmp_h5ad_path)
            if 'gene_ids' in a.var:
                logging.info('renaming gene_ids!!')

                a.var.rename({'gene_ids': 'ensembl_gene_id'}, axis=1, inplace=True)
                a.write_h5ad(tmp_h5ad_path)

            logging.info('Calling postprocessing')
            new_h5_filename = postprocess_fn(tmp_h5ad_path)

            cmd1 = f'cbImportScanpy -i "{new_h5_filename}" -o "{dataset_dir}" -n "{name}" --clusterField leiden --matrixFormat=mtx'
            logging.info(f'running cbImport: {cmd1}')
            sb.run(cmd1, shell=True, check=True)

        with open(job_completed_file, 'w') as fh:
            fh.write('done')

        # customize some of the files
        self._write_cellbrowser_conf(dataset_dir)
        self._write_desc_conf(dataset_dir)

    def pretty_print(self, indent=0):
        tabs = "".join([' ']*indent)
        print(f'{tabs}Title {self.name}')
        print(f'{tabs}H5AD: {self.h5ad}')
        print(f'{tabs}---')

    def check_existance(self):
        """
        check if the h5ad file exists
        """
        if self.is_remote():
            return gs_exists(self.h5ad)
        else:
            return pathlib.Path(self.h5ad).exists()

    def _write_desc_conf(self, directory):
        content = f"""
abstract = '{self.description}'
methods = 'please fill in method'
hideDownload ="{self.hideDownload}"
"""
        with open(f'{directory}/desc.conf', 'w') as fh:
            fh.write(content)

    def _write_cellbrowser_conf(self, directory):
        cbd = self.config.copy()
        cbd['name'] = self.name
        cbd['shortLabel'] = self.name
        content =  _cbconfig_to_str(cbd)

        with open(f'{directory}/cellbrowser.conf', 'w') as fh:
            fh.write(content)


def _recursive_tostr(x):
    """
    to convert dicts to strings
    """
    if isinstance(x, str):
        return "'"+x+"'"
    if isinstance(x, int):
        return f"{x}"

    if isinstance(x, list):
        s = [_recursive_tostr(_) for _ in x]
        s = ','.join(s)
        return f'[{s}]'
    if isinstance(x, dict):
        s = []
        for k,v in x.items():
            assert isinstance(k, str)
            if isinstance(v, str):
                s.append(f'"{k}":"{v}"')
            elif isinstance(v,  int):
                s.append(f'"{k}":{v}')
            else:
                raise ValueError(f'{k}:{v}, value is not a string or int')
        s = ','.join(s)
        s = '{' + s + '}'
        return s

    raise ValueError('_str conversion didnt work')

def _cbconfig_to_str(cb_config):
    """
    create a string/file version of the CB-config dict
    """
    required_fields = ['name', 'shortLabel']
    for f in required_fields:
        assert f in cb_config

    defaults = {
        'exprMatrix': 'matrix.mtx.gz',
        'meta':'meta.tsv',
        'geneIdType':'auto',
        'defColorField':'leiden',
        'labelField':'leiden',
        'enumFields':['leiden', 'nobatch_leiden'],
        'markers':{"file": "markers.tsv", "shortLabel":"Cluster Markers"},  # todo: check that it exists
        'labelField':'leiden',
        'quickGenesFile': 'quickGenes.tsv',
        'coords': [{ "file": "umap_coords.tsv", "shortLabel": "UMAP"}, {"file": "pca_coords.tsv", "shortLabel": "pca"}]  
    }
    defaults.update(cb_config)

    the_string = []
    for k,v in defaults.items():
        vstring = _recursive_tostr(v)
        the_string.append(f'{k}={vstring}')

    return '\n'.join(the_string)


def parse_dataset(the_dict, tmpdir, subsample_for_testing=False, hideDownload=False):
    """
    turn a dataset dictionary entry into a Dataset object

    :param tmpdir: wher to put temporary files when downloading/processing
    """
    logging.info('parsing dataset: %s', the_dict['title'])
    return CB_Dataset(
        name=the_dict['title'],
        description=the_dict['description'],
        h5ad=the_dict['h5ad'],
        config=the_dict['config'],
        tmpdir=tmpdir,
        subsample_for_testing=subsample_for_testing,
        hideDownload=hideDownload,
    )


def parse_collection(the_dict, tmpdir, subsample_for_testing=False, hideDownload=False):
    """
    main function to turn a dictionary into at CB_Collection.
    """
    logging.info('parse collection: %s', the_dict['title'])

    C = CB_Collection(the_dict['title'], description=the_dict['description'], subsample_for_testing=subsample_for_testing, hideDownload=hideDownload)
    for d in the_dict['collection']:
        logging.info(d['title'])
        if 'collection' in d: # the item itself is another collection
            d_parsed = parse_collection(d, tmpdir, subsample_for_testing, hideDownload=hideDownload)
        else:
            d_parsed = parse_dataset(d, tmpdir, subsample_for_testing, hideDownload=hideDownload)
        C._add_dataset(d_parsed)
    return C


def bucket_to_list(bucketname: str):
    '''
    Return bucket's contents to python list of strings.
    We also slice off the bucket name on each line,
    in case we need to search many buckets for one file.
    '''
    if not bucketname.endswith('/'):
        bucketname = bucketname + '/'
        print('warning, bucketname should have a trailing /')

    return sb.run(
        ['gsutil','ls','-r', bucketname + '**'],
        shell=False,
        text=True,
        stdout=sb.PIPE,
        check=True
    ).stdout.replace(bucketname, "").splitlines()

def gs_exists(blobname: str):
    '''
    check if a file exists
    '''
    matches = sb.run(
        ['gsutil','ls', blobname],
        shell=False,
        text=True,
        stdout=sb.PIPE,
        check=True
    ).stdout.splitlines()
    return len(matches) > 0


"""
Stuff to turn the cellbrowser folder into a webpage that can be hosted out of a google bucket
"""

html_404 = """
<!DOCTYPE html>
<html>
  <head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>Page Not Found</title>

    <style media="screen">
      body { background: #ECEFF1; color: rgba(0,0,0,0.87); font-family: Roboto, Helvetica, Arial, sans-serif; margin: 0; padding: 0; }
      #message { background: white; max-width: 360px; margin: 100px auto 16px; padding: 32px 24px 16px; border-radius: 3px; }
      #message h3 { color: #888; font-weight: normal; font-size: 16px; margin: 16px 0 12px; }
      #message h2 { color: #ffa100; font-weight: bold; font-size: 16px; margin: 0 0 8px; }
      #message h1 { font-size: 22px; font-weight: 300; color: rgba(0,0,0,0.6); margin: 0 0 16px;}
      #message p { line-height: 140%; margin: 16px 0 24px; font-size: 14px; }
      #message a { display: block; text-align: center; background: #039be5; text-transform: uppercase; text-decoration: none; color: white; padding: 16px; border-radius: 4px; }
      #message, #message a { box-shadow: 0 1px 3px rgba(0,0,0,0.12), 0 1px 2px rgba(0,0,0,0.24); }
      #load { color: rgba(0,0,0,0.4); text-align: center; font-size: 13px; }
      @media (max-width: 600px) {
        body, #message { margin-top: 0; background: white; box-shadow: none; }
        body { border-top: 16px solid #ffa100; }
      }
    </style>
  </head>
  <body>
    <div id="message">
      <h2>404</h2>
      <h1>Page Not Found</h1>
      <p>The specified file was not found on this website. Please check the URL for mistakes and try again.</p>
      <h3>Why am I seeing this?</h3>
      <p>This page was generated by the Firebase Command-Line Interface. To modify it, edit the <code>404.html</code> file in your project's configured <code>public</code> directory.</p>
    </div>
  </body>
</html>
"""

html_index = """
<!DOCTYPE html>
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<title></title>

<meta name="viewport" content="width=device-width, initial-scale=1.0">

</head>

<style>

body {

	background-image: url('grunge_patterns.jpg');
	background-attachment: fixed;
	color: #333;
}

.box {
	border-radius: 3px;
	background: rgba(101, 101, 101, 0.7); margin: auto; padding: 12px;
}

.lightbox {
	zoom: 1.5;
	position: fixed;
	top: 0;
	left: 0;
	width: 100%;
	height: 100%;
	background: rgba(10, 10, 10, 0.8);
	text-align: center;
	margin: auto;

}

div.horizontal {
	display: flex;
	justify-content: center;
	height: 100%;
}

div.vertical {
	display: flex;
	flex-direction: column;
	justify-content: center;
	width: 100%;
}

::-webkit-input-placeholder {
   color: #955;
   text-align: center;
}

::-moz-placeholder {
   color: #955;
   text-align: center;
}

:-ms-input-placeholder {
   color: #955;
   text-align: center;
}

</style>

<body>
	<div id="loginbox" class="lightbox" >
		<div class="horizontal">
			<div class="vertical">
				<div class="box">
					<input style="margin: 16px; text-align: center;" id="password" type="password" placeholder="password" /> <br />
					<button id="loginbutton" type="button">Enter</button>
					<p id="wrongPassword" style="display: none">wrong password</p>
				</div>
			</div>
		</div>
	</div>

	<script type="text/javascript" src="https://code.jquery.com/jquery-1.12.0.min.js"></script>
	<script type="text/javascript" src="https://rawcdn.githack.com/chrisveness/crypto/7067ee62f18c76dd4a9d372a00e647205460b62b/sha1.js"></script>
	<script type="text/javascript">
	"use strict";


	function loadPage(pwd) {

		var hash= pwd;
		hash= Sha1.hash(pwd);
		hash = hash.replace('#','');
		var url= hash + "/index.html";

		$.ajax({
			url : url,
			dataType : "html",

			success : function(data) {
				window.location= url;
			},

			error : function(xhr, ajaxOptions, thrownError) {
				parent.location.hash= hash;
				//$("#wrongPassword").show();
				$("#password").attr("placeholder","wrong password");
				$("#password").val("");
			}
		});
	}

	$("#loginbutton").on("click", function() {
		loadPage($("#password").val());
	});

	$("#password").keypress(function(e) {
		if (e.which == 13) {
			loadPage($("#password").val());
		}
	});

	$("#password").focus();

	</script>
	<!-- The core Firebase JS SDK is always required and must be listed first -->
	<script src="/__/firebase/7.10.0/firebase-app.js"></script>

	<!-- TODO: Add SDKs for Firebase products that you want to use
     		https://firebase.google.com/docs/web/setup#available-libraries -->

	<!-- Initialize Firebase -->
	<script src="/__/firebase/init.js"></script>

</body>
</html>
"""


def _create_html_scaffold(outdir:str):
    """
    create the index and 404 webpages in the specified directory
    """
    with open(pathlib.Path(outdir) / 'index.html', 'w') as fh:
        fh.write(html_index)
    with open(pathlib.Path(outdir) / '404.html', 'w') as fh:
        fh.write(html_404)

def cellbrowser_to_webpage(cellbrowser_folder:str, out_folder: str, passwd):
    """
    main function to turn a built cellbrowser folder (resulting from .do_build())
    into a hostable webpage with some password protection

    :param cellbrowser_folder: folder where .do_build() put all its files
    :param out_folder: webpage output folder. must exist already
    :param passwd: password to hide the webpage behind
    """
    assert pathlib.Path(out_folder).exists()
    _create_html_scaffold(out_folder)

    encrypted_cellbrowser_foldername = hashlib.sha1(passwd.encode()).hexdigest()

    target = pathlib.Path(out_folder) / encrypted_cellbrowser_foldername
    print(f'copying to {str(target)}')
    # the contents in cellbrowser_folder/* will end up in the hashfolder
    shutil.copytree(cellbrowser_folder, target)
