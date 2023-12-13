import sys
import os
import glob
from time import sleep
from multiprocessing import Pool
from astropy.table import Table, vstack
import astropy.units as u
from argparse import ArgumentParser
from tqdm import tqdm
from astroquery.mast import Observations


# main index functions
def get_dataprod(obstable, calib_level):
    """
    Obtains the dataproduct of the JWST observation
    Parameters
    ----------
        - obstable : Astropy Table
        - calib_level : calibration level to obtain currently fixed 
          to 1 "_uncal" files
    Returns
    ----------
        - outtable : Astropy Table
            - astropy table of the frames
    """
    dp_table = Observations.get_product_list(obstable)
    outtable = dp_table[(dp_table['productType']=='SCIENCE') & 
                        (dp_table['calib_level']==calib_level)]
    
    
    for i in obstable.keys():
        outtable[i]=obstable[i]
    sleep(0.1) # added so to not spam MAST server. Seriously......
    return outtable


# grab all public observations
def grab_obs_projid(args):
    """
    Obtains the observations table of the from MAST
    based on the proposal ID
    Parameters
    ----------
        - args : Dict
            -- user supplied dictionary
    Returns
    ----------
        - outtable : Astropy Table
            - astropy table of the observations
    """
    prop_id=args.prop_id
    all_obs = Observations.query_criteria(proposal_id=prop_id, 
                                          intentType='science',
                                          dataRights='public',
                                          provenance_name='CALJWST')
    all_obs.remove_columns(['target_classification','s_region'])
    
    return all_obs

def grab_obs_position(args):
    """
    Obtains the observations table of the from MAST
    using the coordinates position
    Parameters
    ----------
        - args : Dict
            -- user supplied dictionary
    Returns
    ----------
        - outtable : Astropy Table
            - astropy table of the observations
    """
    ra, dec, rad=args.ra, args.dec, args.rad/60  
    prop_id=args.prop_id
    all_obs=Observations.query_object("{} {}".format(ra,dec), 
                                      radius=rad*u.deg)
    
    # keep only JWST
    all_obs=all_obs[(all_obs['provenance_name']=='CALJWST') &
                    (all_obs['dataRights']=='PUBLIC') &
                    (all_obs['intentType']=='science')]
    
    if prop_id!='-1':
        all_obs=all_obs[all_obs['proposal_id']==str(prop_id)] 


    # do a bit of clean up
    all_obs.remove_columns(['target_classification','s_region'])

    return all_obs


def pick_instrument(args, in_table):
    """"
    command to parse the user supplied instrument keys to the table 
    instrument keyword
    Parameters
    ----------
        - key : int
            -- user supplied instrument keyword
    Returns
    ----------
        - keyword on MAST website
    """
    if args.instrum==0:
        # all observations
        print('pick all instrument')
        out_table=in_table
    elif args.instrum==1:
        # pick NIRCAM images
        out_table=in_table[(in_table['instrument_name']=='NIRCAM/IMAGE')]
    elif args.instrum==2:
        # pick MIRI images
        out_table=in_table[(in_table['instrument_name']=='MIRI/IMAGE')]
    else:
        print('Unknown Instrument')
        sys.exit(1)
    return out_table
    

def grab_obs(args):
    """
    Wrapper to choose the observations query
    Parameters : 
    ------------
        - args : Dictionary
            -- user supplied dictionary
    Returns
        - out_cat : Astropy Table
            -- Astropy Table of the JWST 
             observations
    """

    if args.type=='pid':
        out_cat=grab_obs_projid(args)
    elif args.type=='sky':
        out_cat=grab_obs_position(args)
    else:
        raise AssertionError('Unknown task to grab observation')
    return out_cat

def parse_output_fname(args):
    """
    function parse the output index file name
    Parameters : 
    ------------
        - args : Dictionary
            -- user supplied dictionary
    Returns
        - out_filename : str
            -- output filename

    """
    if args.suffix is not None:
        suffix='_'+str(args.suffix)
    else:
        suffix=''
    
    if args.type=='pid':
        fn_string='index_pidm_inst_{}_proj_{}_obs{}.csv'
        out_filename=fn_string.format(args.instrum,args.prop_id, suffix)

    elif args.type=='sky':        
        fn_string='index_skym_inst_{}_proj_{}_obs{}.csv'
        if args.prop_id!='-1':
            out_filename=fn_string.format(args.instrum,args.prop_id, suffix)
        else:
            out_filename=fn_string.format(args.instrum,'all', suffix)

    else:
        raise AssertionError('Unknown task to grab observation')
        
    return out_filename


# main sections
def main_index(args):
    """
    Main body of the index command. searches for the images based on the given 
    on the given parameters
    Parameters
    ----------
        - args : Dictionary
    Returns
    ----------
        - None.
    """
    # grab public observations
    print('Grabbing Observations')
    pub_obs=grab_obs(args)#[0:5]
    # do some clean up
    
    # pick out the intruments you want
    pub_obs=pick_instrument(args, pub_obs)
    
    # grab files, fix calib level to 1 for now
    print('Fetching frames IDs')
    frame_table_list=[get_dataprod(i,1) for i in tqdm(pub_obs)]
    combined_frames_table=vstack(frame_table_list)

    # write output index table
    os.makedirs(args.outdir, exist_ok=True)
    outpath=os.path.join(args.outdir, parse_output_fname(args))
    print('Writing index files to {}'.format(outpath))
    combined_frames_table.write(outpath, format='ascii.csv', overwrite=True)

    
#####  main download fuctions

def download_file(file, outdir, cache):
    """
    command to download the image
    instrument keyword
    Parameters
    ----------
        - file : str
            - filename with MAST key from the 
              "getdataproduc command"
        - outdir : str
            - output directory 
        - cache : Boolean
            - cache flag.
                -- True : the files are checked first already downloaded
                -- False: Force redownload
    Returns
    ----------
        - flag : integer 
            -- interger flag indicating download success
                0=fail
                1=success
    """
    fname=os.path.basename(file)
    try:
        Observations.download_file(file,local_path=os.path.join(outdir, fname),
                                   cache=cache)
        flag=1
    except:
        flag=0
    return flag


def parse_instrum(key):
    """
    command to parse the user supplied instrument keys to the table 
    instrument keyword
    Parameters
    ----------
        - key : int
            -- user supplied instrument keyword
    Returns
    ----------
        - keyword on MAST website
    """
    if key=='NIRCAM/IMAGE':
        out='NRC_IMA'
    elif key=='MIRI/IMAGE':
        out='MIRI_IMA'
    else:
        raise KeyError('Unknown Instrument in the Table')
    return out

def make_path(row,args):
    """
    command to create the paths of each file
    Parameters
    ----------
        - row : list
            -- list of directory paths to  create
        - args : Dict
            -- user supplied parameters
    Returns
    ----------
        - None
    """
    base_dir=args.outdir
    inst_dir=parse_instrum(row['instrument_name'])
    filt_dir=row['filters']
    out_dir=os.path.join(base_dir, inst_dir, filt_dir, 'uncal')
    return out_dir


def make_outputdir(dirlist):
    """
    command to create directory
    Parameters
    ----------
        - dirlist : list
            -- list of directory paths to  create
    Returns
    ----------
        - None
    """
    for i in dirlist: os.makedirs(i, exist_ok=True)

    
def prep_payload(table, args):
    """
    command to prepare the parametrs for the down load function.
    Parameters
    ----------
        - table : Astropy Table objects
            -- Astropy table containing the ObsID of frames to download
    Returns
    ----------
        - dw_results: list
            list containing the flag reporting successful/failed downloads
    """
    package=[]
    for i in table:
        package.append((str(i['dataURI']), str(i['outdir']), bool(i['cache'])))
        
    return package


def run_download_single(package):
    """
    command to perform download in sequence.
    Parameters
    ----------
        - package : tuple
            -- tuple containing the parameters which are to be
            supplied to the download file function
    Returns
    ----------
        - dw_results: list
            list containing the flag reporting successful/failed downloads
    """
    dw_results=[download_file(*i) for i in package]
    return dw_results


def run_download_parallel(package, ncpu=1):
    """
    command to perform download in parallel.
    Parameters
    ----------
        - package : tuple
            -- tuple containing the parameters which are to be
            supplied to the download file function
        - ncpu : int
            -- number of threads to spawn download 
    Returns
    ----------
        - dw_results[0]: list
            list containing the flag reporting successful/failed downloads
    """
    with Pool(processes=ncpu) as pool:
        dw_results=[pool.starmap(download_file, package)]
    return dw_results[0]


def run_download(package, args):
    """
    Wrapper of the download command.
    Parameters
    ----------
        - package : tuple
            -- tuple containing the parameters which are to be
            supplied to the download file function
        - args : Dictionary
            -- user supplied parameters
    Returns
    ----------
        - out : list
            list containing the flag reporting successful/failed downloads
    """
    if args.ncpu==1:
        out=run_download_single(package)
    elif args.ncpu>1:
        out=run_download_parallel(package, ncpu=args.ncpu)
    else:
        raise AssertionError('Unknown download option')
    return out

    
def main_download(args):
    """
    Main body of the fetch command. Performs the download based on the given arguments
    Parameters
    ----------
        - args : Dictionary
    Returns
    ----------
        - None.
    """
        
    print('Running File Downloading')

    # read frame tables
    frame_table=Table.read(args.index)
    frame_table['cache']= not args.refetch
    frame_table['outdir']=[make_path(i, args) for i in frame_table]
    
    # make output directory
    make_outputdir(set(frame_table['outdir']))
    print('Found {} frames in {}'.format(len(frame_table), ','.join(set(frame_table['filters']))))
    
    # prepare download packages
    dw_payload=prep_payload(frame_table, args)
    
    # do downloads, save download results
    dw_success=run_download(dw_payload, args)

    # summerize it
    frame_table['dw_success']=dw_success
    print('Successfully downloaded {}/{} frames'.format(len(frame_table[frame_table['dw_success']==1]),len(frame_table)))
    
    sys.exit(0)
    
# Overall Main function


def main(args):
    """
    Main body of the code. Chooses they search type based on the "type" arguement
    Parameters
    ----------
        - args : Dictionary
    Returns
    ----------
        - None
    
    """
    
    
    if args.mode=='index':
        main_index(args)
    elif args.mode=='fetch':
        main_download(args)
    else:
        print('Unknown command, killing myself')
        sys.exit(1)

if __name__ == '__main__':
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(dest="mode",help='sub-command help')

    # create the parser for the "index" command
    parser_a = subparsers.add_parser('index', help='index help')
    parser_a.add_argument('-t','--type', type=str, default='None', choices=['sky','pid'], help='Observation searcher')
    parser_a.add_argument('--ra', type=float, default=0, help='Right Ascension to search in degrees. Default=0 Deg')
    parser_a.add_argument('--dec', type=float, default=0, help='Declination to search in degrees. Default=0 Deg' )
    parser_a.add_argument('-r','--rad',  type=float, default=1, help='Radius of positional search in arcmin. Default=1 arcmin')
    parser_a.add_argument('-p','--prop_id', type=str, default='-1', help='JWST Proposal ID')
    parser_a.add_argument('-i','--instrum', type=int, default=0, choices=[0,1,2], help='JWST instrument to pickout\n0=All\n1=NIRCAM/IMAGE\n2=MIRI/IMAGE')
    parser_a.add_argument('-o','--outdir', type=str, default='./', help='Output directory string. Default="./"')
    parser_a.add_argument('-s','--suffix', type=str, default=None, help='Suffix to the index file')


    # create the parser for the "fetch" command
    parser_b = subparsers.add_parser('fetch', help='fetch help')
    parser_b.add_argument('-i','--index', type=str, default=None,  help='Frame index files')
    parser_b.add_argument('-o','--outdir', type=str, default='./jwst_data', help='Output directory string. Default="./"')
    parser_b.add_argument('-r','--refetch', action='store_true', help='Forces redownload of the file')
    parser_b.add_argument('-n','--ncpu', type=int, default=1, help='Number of CPUs to run download. Default=1')

    args=parser.parse_args()

    if len(sys.argv)==1:
        print(parser.print_help())
    else:
        if (len(sys.argv)==2) & (args.mode=='index'):
            print(parser_a.print_help())
        elif (len(sys.argv)==2) & (args.mode=='fetch'):
            print(parser_b.print_help())
        else:
            print(args)
            main(args)
    