# JWST_raider

JWST_raider is a python script CLI tool which wraps over the `Mast.Observations` module of astropy query module to index and download **public uncalibrated** JWST images. Currently only supports NIRCAM/IMAGING and MIRI/Imaging dataset search.

## Dependency
JWST_raider depends mainly on the `astropy`,`tqdm`, and  `multiprocessing` packages. `Multiprocessing` is used when the user requests the download to be performed in parallel.

## Usage

JWST_raider has 2 main functions, `index` and `fetch`. First, the index function queries the MAST archive for the JWST frames and creates an index file. The index file is an ASCII CSV file which contains the JWST images.  The index function currently supports two query modes, `sky` and `pid`. For `sky`, the `sky` function performs a cone search while `pid` performs  a search based on proposal id. For `sky` mode, filtering with proposal ID is also supported. Second, the `fetch` function downloads the images and groups them into folders according to the instrument and the Filters.  While separating the two functions may seem cumbersome to many, it allows the user to slice out selected frames and check the expected number of downloads before running the main download scripts. This is useful for downloading data from projects which have several fields which may or not have the same target names.

To display the command helps, type `python jwst_raider.py --help`, `python jwst_raider.py index --help`, and `python jwst_raider.py fetch --help`
 
## Usage
### Indexing with PID Mode
Here, we will use data from the PRIMER team (PID: 1837) which have data in the GOODS-North and UDS field as an example. Indexing in `pid` mode can be called with the following command

``python jwst_raider.py index --type=pid -pid=1837``

This will search all files associated with PID:1837 and write the all frames (NIRCAM and MIRI) information to the file `./index_pidm_inst_all_proj_1837_obs.csv`. If the user is interested solely NIRCAM images, the flag `--instrum` can be set to `--instrum=1` which will write only NIRCAM images to the index file. Currently, `0=all`, `1=NIRCAM/IMAGING` and `2=MIRI/IMAGING`. An output file directory can also be specified with `-outdir=dir` and a suffix can be appended to the index filename with the flag `--suffix`. For example,

``python jwst_raider.py index --type=pid -pid=1837 -i=1 -o=filelist -s=20231212`` 

will search for NIRCAM images of PID:1837 and write the index file to the `filelist` directory with the suffix 20231212.

### Indexing with sky Mode
Again, we will use data from the PRIMER team (PID: 1837) specifically in the GOODS-North. For example, indexing in `sky` mode can be called with the following command

``python jwst_raider.py sky --type=sky --ra=150.13 --dec=2.36 --rad=30 -pid=1837 -o=filelist -s=20231212_cosmos``

will perform a cone search for NIRCAM images centred on RA=150.13 and Dec=2.36 with a radius of 30 arcminutes that are associated with PID:1837 and write the index file to the `filelist` directory with suffix "_20231212_cosmos". 

The longest operation is the fetching of image data products from the MAST servers. Depending on the number of observations, this can take up to several minutes. In the meantime, sit back, relax, and grab a beer. 

### Downloading
Downloading of the files listed in an index file can be performed with the following command.

`python jwst_raider.py fetch --index=index_pidm_inst_1_proj_1837_obs.csv`

The script will download the images listed in `index_pidm_inst_1_proj_1837_obs.csv` to the directory `./data` or the user can point to an output directory using the flag `--outdir=PATH`. By default, the files are separated into folders according to instrument and filters. During the download, the script also checks if the file has already been downloaded. By default, the script skips files already downloaded but can be forced to re-download the images using the flag `--refetch` (no values). Downloads can be performed in parallel by setting the keyword `--ncpu=` according to the number of threads you want.





