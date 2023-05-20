#!/usr/bin/env python

__doc__ = """Cortical RSFC Hub Profile Script: This script will take the outputs
from Identify_Hubs.py and create Hub Profiles and a Hub Profile correlation matrix
that can be used with clustering algorithms to create hub categories. 

Cleaned/generalized for sharing May, 2023 - DVD
Please cite Demeter, 2023 [1] if used for analyses.
"""
__references__ = """References
----------
[1] Demeter, D.V., Gordon, E.M., Nugiel, T., Garza, A., Larguinho, T.L., & Church, J.A. (2023). 
Resting-state cortical hubs in youth organize into four categories. Cell Reports, 42(5), 112521. 
https://doi.org/10.1016/j.celrep.2023.112521
[2] https://doi.org/10.5281/zenodo.7814714
"""
__version__ = "0.1.0"

import argparse,datetime,glob,logging,os,subprocess,sys,time
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
## suppress matplotlib font warnings and junk
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)

start_time=time.time()
def main(argv=sys.argv):
    arg_parser = argparse.ArgumentParser(prog='Create_Hub_Profiles.py',
                                         description=__doc__,
                                         formatter_class=argparse.RawDescriptionHelpFormatter,
                                         epilog=__references__,
                                         usage='%(prog)s outdir [OPTIONS]')
    # Check for arguments. #
    if len(sys.argv[1:])==0:
        print('\nArguments required. Use -h option to print FULL usage.\n')
    arg_parser.add_argument('outdir', type=os.path.abspath,
                            help='Path to output directory. (This should be the same dir '
                                 'you used for previous step outputs) /final_hub_indices/ '
                                 'folder should be under this main path.'
                            )
    arg_parser.add_argument('-n', metavar='', action='store', type=str, required=True,
                            help='Name for final hub profile corr mat.',
                            dest='name'
                            )
    arg_parser.add_argument('-orderlist', metavar='', action='store', type=os.path.abspath, required=False,
                            help='Optional subID order list .txt for hub profile correlation '
                                 'matrix. One subID per line. Group IDs to create matrix like '
                                 'Figure S3 (B).',
                            dest='orderlist'
                            )
    arg_parser.add_argument('-v','--version', action='version', version='%(prog)s: ' + __version__)
    args = arg_parser.parse_args()
    # Setting up logger #
    logging.basicConfig(level=logging.DEBUG, format='%(message)s')
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    logger.debug('\n--------------------- setup info ---------------------------------')
    # check input dirs
    hub_indices_dir = os.path.join(args.outdir,'final_hub_indices')
    if os.path.isdir(hub_indices_dir):
        logger.debug(f'hub indices: {hub_indices_dir}')
    else:
        sys.exit(f'hub indices not found! {hub_indices_dir}')
    dist_cens_zmat_dir = os.path.join(args.outdir,'final_csv_outputs')
    if os.path.isdir(dist_cens_zmat_dir):
        logger.debug(f'zmat dir: {dist_cens_zmat_dir}')
    else:
        sys.exit(f'zmat dir not found! {dist_cens_zmat_dir}')
    conn_profiles_dir = os.path.join(args.outdir,'final_conn_profiles')
    if os.path.isdir(conn_profiles_dir):
        pass
    else:
        os.makedirs(conn_profiles_dir)
    # Group name
    logger.debug(f'Group/File Name: {args.name}')
    ## quick numbers check:
    ind_glob = glob.glob(os.path.join(hub_indices_dir,'*'))
    z_glob = glob.glob(os.path.join(dist_cens_zmat_dir,'*'))
    logger.debug(f'Hub Indices files found: {len(ind_glob)}')
    logger.debug(f'zmat files found: {len(z_glob)} ')
    logger.debug('------------------------- end ------------------------------------\n')
    here = os.path.dirname(os.path.abspath(__file__))

    # Zero indexed network indices for averaging - ugly but works
    aud = [9,63,64,65,66,67,68,69,76,101,103,159,170,223,226,229,231,232,238,243,267,268,328,329]
    co = [20,21,26,27,33,39,62,70,71,75,80,81,83,100,102,104,110,111,146,152,179,180,184,186,187,
        191,195,197,218,222,233,234,237,244,245,247,248,273,316,317]
    cp = [11,88,92,172,253]
    dmn = [0,3,5,24,25,43,93,113,115,116,125,126,144,145,149,150,151,153,155,156,161,164,183,185,
        199,219,224,256,258,277,278,289,314,315,320,321,322,323,324,325,330]
    da = [40,41,42,48,50,51,54,73,86,87,90,91,94,99,105,106,109,112,154,188,198,202,207,210,235,
        249,251,252,261,265,270,274]
    fp = [6,8,23,77,95,107,108,147,148,166,167,169,181,239,259,260,271,272,275,276,318,319,326,327]
    none = [10,17,18,72,114,117,118,119,120,121,122,123,124,127,128,132,133,134,141,143,158,171,177,
            178,279,280,281,282,283,284,285,286,287,288,290,291,295,296,299,300,301,302,303,304,305,311,313]
    ret = [12,13,129,142,173,293,294,312]
    sal = [28,82,182,246]
    ssmh = [1,29,30,31,32,34,35,36,37,44,45,46,47,49,53,55,56,57,162,189,190,192,193,194,200,201,203,204,205,
            206,208,209,212,213,214,215,216,269]
    ssmm = [2,38,52,58,163,196,211,217]
    va = [22,59,60,61,74,78,79,84,85,157,160,220,221,225,227,228,230,236,240,241,242,331,332]
    vis = [4,7,14,15,16,19,89,96,97,98,130,131,135,136,137,138,139,140,165,168,174,175,176,250,254,255,257,
        262,263,264,266,292,297,298,306,307,308,309,310]

    def create_profiles(sub,zmat_file,hub_indices):
        # make hub connectivity profiles & save to .txt
        conn_profiles = []
        for hi in hub_indices:
            # get all connectivity values for that hub - as list for easy indexing
            hi = hi - 1
            all_hub_conns = zmat.loc[hi, :].values.tolist()
            corr_removed = 0
            # make list and check for the self correlation
            ## SLOPPY but neurotic
            aud_conns = [all_hub_conns[i] for i in aud] # 24 items
            co_conns = [all_hub_conns[i] for i in co] # 40 items
            cp_conns = [all_hub_conns[i] for i in cp] # 5 items
            dmn_conns = [all_hub_conns[i] for i in dmn] # 41 items
            da_conns = [all_hub_conns[i] for i in da] # 32 items
            fp_conns = [all_hub_conns[i] for i in fp] # 24 items
            none_conns = [all_hub_conns[i] for i in none] # 47 items
            ret_conns = [all_hub_conns[i] for i in ret] # 8 items
            sal_conns = [all_hub_conns[i] for i in sal] # 4 items
            ssmh_conns = [all_hub_conns[i] for i in ssmh] # 38 items
            ssmm_conns = [all_hub_conns[i] for i in ssmm] # 8 items
            va_conns = [all_hub_conns[i] for i in va] # 23 items
            vis_conns = [all_hub_conns[i] for i in vis] # 39 items
            net_conns = [aud_conns,co_conns,cp_conns,dmn_conns,da_conns,
                        fp_conns,none_conns,ret_conns,sal_conns,ssmh_conns,
                        ssmm_conns,va_conns,vis_conns]
            for nc in net_conns:
                if 1.0 in nc:
                    corr_removed = corr_removed + 1
                    nc.remove(1.0)
                else:
                    pass
            if corr_removed == 1:
                pass
            else:
                sys.exit('Something went wrong removing the self-corr value. Exiting...')
            avg_conns_list = [np.mean(aud_conns),np.mean(co_conns),np.mean(cp_conns),
                            np.mean(dmn_conns),np.mean(da_conns),np.mean(fp_conns),
                            np.mean(none_conns),np.mean(ret_conns),np.mean(sal_conns),
                            np.mean(ssmh_conns),np.mean(ssmm_conns),np.mean(va_conns),
                            np.mean(vis_conns)]
            avg_conns_list = [str(round(num, 4)) for num in avg_conns_list]
            profile_list = [sub + '_' + str(hi + 1)] + avg_conns_list
            profile = ' '.join(profile_list)
            conn_profiles.append(profile)   
        ## write hub conn profiles to file for THIS SUBJECT
        cp_outfile_path = os.path.join(conn_profiles_dir, sub + '_HUB_CONN_PROFILES.txt')
        f = open(cp_outfile_path, 'w')
        for hcp in conn_profiles:
            hcp_line = (str(hcp),'\n')
            f.write(' '.join(hcp_line))
        f.close()
        
    def make_corr(IDorder):
        # Create and save hub profile matrix. 
        hcprofiles_dict = {}
        for sub in IDorder:
            hc_profile_path = os.path.join(conn_profiles_dir, sub + '_HUB_CONN_PROFILES.txt')
            with open(hc_profile_path, 'r') as f:
                for i in f.readlines():
                    hcp_vect = i.strip()
                    hcp_vect = [i for i in hcp_vect.split()]
                    the_key = hcp_vect[0]
                    # remove key string & convert to floats
                    del hcp_vect[0]
                    hcp_vect = [float(i) for i in hcp_vect]
                    hcprofiles_dict[the_key] = hcp_vect
        ## put lists into a dataframe and correlate for louvain
        all_hubprofs_df = pd.DataFrame(hcprofiles_dict)
        prof_corr = all_hubprofs_df.corr(method='pearson')
        ## Save prof_corr as .csv
        logger.debug(f'Saving hub profile correlation matrix .csv for clustering...')
        out_mat_path = os.path.join(args.outdir,args.name + '_Hub_Profile_CorrMat.csv')
        np.savetxt(out_mat_path,prof_corr,fmt='%f',delimiter=',')
        # Draw the heatmap with the mask and correct aspect ratio (not very useful at this stage)
        logger.debug(f'Saving hub profile correlation plot (not really useful at this stage)...')
        sns.heatmap(prof_corr, vmin=-1, vmax=1, center= 0, cmap= 'coolwarm')
        plt.savefig(os.path.join(args.outdir,f'{args.name}_Hub_Profile_Correlations.png'), dpi=1200, format='png', bbox_inches='tight')

    ### SCRIPT ENTRY ###
    ## MAKE subject list
    sub_list = []
    for x in ind_glob:
        the_file = os.path.basename(x)
        subid = the_file.replace('_HUB_INDICES.txt','')
        sub_list.append(subid)
    # I. Make hub profiles
    for sub in sub_list:
        logger.debug(f'Working on {sub}...')
        # get zmat
        zmat_file = os.path.join(dist_cens_zmat_dir,
                                sub + '_DIST_CENSORED_ZMAT.csv')
        zmat = pd.read_csv(zmat_file, header=None)
        # get hubs                   
        hub_file = os.path.join(hub_indices_dir,
                                sub + '_HUB_INDICES.txt') 
        with open(hub_file, 'r') as f:
            hub_indices = [i.strip() for i in f.readlines()]
        # make ints
        hub_indices = [int(i) for i in hub_indices]
        ## create hub profiles
        create_profiles(sub,zmat_file,hub_indices)
        # quick check profile files made (match to expected participants)
        all_conn_profiles = glob.glob(os.path.join(conn_profiles_dir,'*'))
    logger.debug(f'Hub profiles made for {len(all_conn_profiles)} participants.\n')
    # II. Make profile corr mat
    if args.orderlist:
        try:
            with open(args.orderlist, 'r') as file:
                IDorder = file.readlines()
                IDorder = [line.strip() for line in IDorder]
                IDorder = [i for i in IDorder if i]
        except FileNotFoundError:
            sys.exit(f'ID order list file not found. Check path: {args.orderlist}')
    else:
        IDorder = sub_list
    make_corr(IDorder)

    full_runtime = time.time() - start_time
    print('\nFull Script Runtime: ' + str(datetime.timedelta(seconds=round(full_runtime))))
if __name__ == '__main__':
    sys.exit(main())
