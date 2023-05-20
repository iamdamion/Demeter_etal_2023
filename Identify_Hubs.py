#!/usr/bin/env python

__doc__ = """Cortical RSFC Hubs Script: This script will identify and categorize 
cortical hubs, using resting state timeseries extracted from the Gordon 333
parcel set. Script to replicate previously published findings. 

Notes:
1. Fully processed 333 parcel timeseries should be motion censored ahead of time. 
    a. They should also be in parcel number order (row 1 = timeseries of parcel01, etc)
2. Distance censoring matrix that's provided works for the Conte69_fs_LR atlas. 
    b. If another atlas was used or this is being run in subject space, a new distance
       censoring matrix should be created.  

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

import argparse,bct,datetime,glob,logging,os,random,shutil,subprocess,sys,time
import networkx as nx
import nibabel as nb
import numpy as np
from scipy import stats

start_time=time.time()
def main(argv=sys.argv):
    arg_parser = argparse.ArgumentParser(prog='Identify_Hubs.py',
                                         description=__doc__,
                                         formatter_class=argparse.RawDescriptionHelpFormatter,
                                         epilog=__references__,
                                         usage='%(prog)s inputlist [OPTIONS]')
    # Check for arguments. #
    if len(sys.argv[1:])==0:
        print('\nArguments required. Use -h option to print FULL usage.\n')

    arg_parser.add_argument('inputlist', type=os.path.abspath,
                            help='Path to .txt file containing your participant list. List should '
                                 'include space separated data, with each line containting: '
                                 '<participantID> <FULLY processed 333 parcel timeseries>'
                            )
    arg_parser.add_argument('-attempts', action='store', type=int, required=False,
                            default=1000,
                            help='Number of Infomap iteration attempts. (Default=1,000)',
                            dest='attempts'
                            )
    arg_parser.add_argument('-nocleanup', action='store_true', required=False,
                            default=False,
                            help='Suppress temporary file cleanup. '
                                 '(Default=False)',
                            dest='nocleanup'
                            )
    arg_parser.add_argument('-o', metavar='', action='store', type=os.path.abspath, required=True,
                            help='Path to output directory',
                            dest='outdir'
                            )
    arg_parser.add_argument('-overlay', action='store', type=int, required=False,
                            choices=[1,2,3],
                            default=1,
                            help='Overlay color fo identified hub parcels. Current choices: '
                                 '1: Pink, 2: Crimson, 3: Lime. (Default=Pink)',
                            dest='overlay'
                            )
    arg_parser.add_argument('-v','--version', action='version', version='%(prog)s: ' + __version__)
    args = arg_parser.parse_args()
    # Setting up logger #
    logging.basicConfig(level=logging.DEBUG, format='%(message)s')
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    logger.debug('\n--------------------- setup info ---------------------------------')
    # load & clean input data list
    try:
        with open(args.inputlist, 'r') as file:
            input_data = file.readlines()
            input_data = [line.strip() for line in input_data]
            input_data = [i for i in input_data if i]
    except FileNotFoundError:
        sys.exit(f'Input data file not found. Check path: {args.inputlist}')
    logger.debug(f'IDs and Paths imported: {args.inputlist}')
    # create temporary file structure
    out_subdirs = ['pajek_files','infomap_outputs','pc_outputs','csv_outputs',
                   'final_avg_pc_percs','final_hub_indices','final_hubs_dlabels']
    if os.path.isdir(args.outdir):
        for osub in out_subdirs:
            sdir = os.path.join(args.outdir,osub)
            if os.path.isdir(sdir):
                pass
            else:
                os.makedirs(sdir)
    else:
        sys.exit(f'Output dir not found. Please check: {args.outdir}')
    logger.debug(f'Output Dir: {args.outdir}')
    # handle basic overlay colors
    if args.overlay == 1:
        overlay_color = (1,0.07843,0.57647,1) # pink
        logger.debug(f'Hub Overlay Color: Pink')
    elif args.overlay == 2:
        overlay_color = (0.8627,0.07843,0.23529,1) # crimson
        logger.debug(f'Hub Overlay Color: Crimson')
    elif args.overlay == 3:
        overlay_color = (0,1,0,1) # lime
        logger.debug(f'Hub Overlay Color: Lime')

    logger.debug('------------------------- end ------------------------------------\n')
    ## global vars
    here = os.path.dirname(os.path.abspath(__file__))
    mult_mat_path = os.path.join(here,'COMBINED333_LR_Distance_MULTIPLICATION_MASK.csv')
    gordon_parcel_path = os.path.join(here,'Parcels_LR.dlabel.nii')
    threshold_list = [0.003, 0.004, 0.005, .01, .015, .02, .025, .03, .035, .04, .045, .05]

    def make_zmat(preprocd_ts):
        sub_corr_mat = np.corrcoef(preprocd_ts)
        np.fill_diagonal(sub_corr_mat,0)
        ztrans_mat = np.arctanh(sub_corr_mat)
        np.fill_diagonal(ztrans_mat,1)
        return ztrans_mat
    
    def dist_censor(ztrans_mat):
        mult_mat = np.genfromtxt(mult_mat_path, delimiter=',')
        dist_censored_zmat = ztrans_mat.copy()
        dist_censored_zmat = dist_censored_zmat * mult_mat
        dist_censored_zmat[dist_censored_zmat==-0]=0 # correct -0 for aesthetics
        # save distance censored csvs for later steps
        csv_dir = os.path.join(args.outdir,'csv_outputs')
        if os.path.isdir(csv_dir):
            pass
        else:
            os.makedirs(csv_dir)
        out_mat_path = os.path.join(csv_dir,subid + '_DIST_CENSORED_ZMAT.csv')
        np.savetxt(out_mat_path,dist_censored_zmat,fmt='%f',delimiter=',')
        return dist_censored_zmat

    def run_infomap(subid, ready_zmat):
        ##  Make pajek & run infomap
        infomap_out_dir = os.path.join(args.outdir,'infomap_outputs')
        pajek_out_dir = os.path.join(args.outdir,'pajek_files')
        for f in glob.glob(os.path.join(infomap_out_dir,subid + '*')):
            os.remove(f)
        for thresh in threshold_list:  
            ##  threshold distance censored zmat
            thresh_mat = bct.threshold_proportional(ready_zmat, float(thresh),copy=True)
            check_nan = np.isnan(np.sum(thresh_mat))
            if check_nan == True:
                sys.exit('NaN found. This is not right! Exiting...')
            else:
                pass
            logger.debug(f'  -Creating pajek file: {thresh}')
            nx_graph = nx.from_numpy_matrix(thresh_mat)
            pajek_path = os.path.join(pajek_out_dir,
                                      subid + '_' + str(thresh) + '_upper_mat.net')
            nx.write_pajek(nx_graph,pajek_path)
            ##  subprocess infomap
            output_name = os.path.basename(pajek_path).replace('.net','')
            rand_num = random.randint(1,9999)
            infomap_comm = ' '.join(['infomap',pajek_path,infomap_out_dir,
                                     '--clu', # cluster indices for all nodes
                                     '-2', # two-level (confirmed with EG)
                                     '-s ' + str(rand_num), # picking random number for the seed of the internal random number generator
                                     '-N ' + str(args.attempts), # number of attempts (iterations) for infomap
                                     '--out-name ' + output_name, # Naming for output files (to not overwrite pajek...)
                                     '--silent'
                                     ])
            logger.debug(f'  - Running infomap ({args.attempts} iterations) silently')
            subprocess.call(infomap_comm, shell=True)
        ## clean up affiliation vector and save to file. 
        clu_files = glob.glob(os.path.join(infomap_out_dir,f'{subid}*'))
        logger.debug(f'  -Cleaning up clu files and saving affiliation vectors')
        for cf in clu_files:
            if os.path.isdir(cf):
                pass
            else:
                clu_tup = []
                with open(cf, "r") as clu_file:
                    for line in clu_file:
                        stripped_line = line.strip()
                        if '#' in stripped_line:
                            pass
                        else:
                            sl = stripped_line.split(' ')
                            clu_tup.append((int(sl[0]),int(sl[1])))
                clu_file.close()
                aff_ref_path = cf.replace('_upper_mat.clu','_aff_vect_REF.txt')
                aff_vect_path = cf.replace('_upper_mat.clu','_aff_vect.txt')
                ## order the tuples by the parcel numbers
                s_tups = sorted(clu_tup, key=lambda x: x[0])
                ## write to the affiliation vector REFERENCE first (reference has the parcel number also)
                f = open(aff_ref_path, 'w')
                for c in s_tups:
                    c_line = (str(c[0]),str(c[1]),'\n')
                    f.write(' '.join(c_line))
                f.close()
                ## next write affiliation vector. (communities only to import for PC step)
                f = open(aff_vect_path, 'w')
                for c in s_tups:
                    c_line = (str(c[1]),'\n')
                    f.write(' '.join(c_line))
                f.close() 

    def run_pc_steps(subid):
        ## Step 6,7,8, and 9: Calc PC, get percents, censor out low degree nodes, get average percent, and save as file. 
        logger.debug(f'  -Load distance censored zmat')
        orig_distcensored_mat_path = os.path.join(args.outdir,
                                                  'csv_outputs',
                                                  subid + '_DIST_CENSORED_ZMAT.csv')
        orig_dist_censored_zmat = np.genfromtxt(orig_distcensored_mat_path, delimiter=',')
        pc_outputs_dir = os.path.join(args.outdir,'pc_outputs')
        for thr in threshold_list:
            logger.debug(f'   thresh: {thr} - calc pc, get percents, censor low degree nodes, avg percs, save file...')
            ## re-threshold (and take upper mat only) original dist censored matrix
            pc_thresh_mat = bct.threshold_proportional(orig_dist_censored_zmat, float(thr))
            pc_upper_mat = np.triu(pc_thresh_mat,1)
            ## load that threshold's affiliation vector created from infomap step
            vect_text = os.path.join(args.outdir,
                                     'infomap_outputs',
                                     '_'.join([subid,str(thr),'aff','vect.txt'])
                                    )
            aff_vect = np.loadtxt(vect_text)
            ## calc PC with AND without nans
            PC_nonans = bct.participation_coef(pc_upper_mat, aff_vect, degree='undirected')
            # with nans:
            W = pc_upper_mat
            ci = aff_vect
            degree = 'undirected'
            if degree == 'in':
                W = W.T
            _, ci = np.unique(ci, return_inverse=True)
            ci += 1
            n = len(W)  # number of vertices
            Ko = np.sum(W, axis=1)  # (out) degree
            Gc = np.dot((W != 0), np.diag(ci))  # neighbor community affiliation
            Kc2 = np.zeros((n,))  # community-specific neighbors
            for i in range(1, int(np.max(ci)) + 1):
                Kc2 += np.square(np.sum(W * (Gc == i), axis=1))
            P = np.ones((n,)) - Kc2 / np.square(Ko)
            PC_wnans = P
            ## retain nan mask so we can put them back in when needed
            PC_nan_mask = np.argwhere(np.isnan(PC_wnans))
            ## get PC percs with original (to avoid nan error and following adult paper methods)
            PC_percs = np.asarray([stats.percentileofscore(PC_nonans, i) for i in PC_nonans])
            ## Censor out any low degree (<25th perc) nodes
            degs = bct.degrees_und(pc_upper_mat)
            low = np.percentile(degs, 25)
            deg_mult = np.where(degs < low, 0, 1)
            deg_cens_PC_percs = PC_percs * deg_mult
            ## convert to nans for nanmean function later
            deg_cens_PC_percs[deg_cens_PC_percs == 0] = np.float64('nan')
            # put back nans from PC step (again for nanmean step)
            cens_PC_perc_wnan = deg_cens_PC_percs.copy() # making a copy to be safe/no overwrite
            cens_PC_perc_wnan[PC_nan_mask] = np.float64('nan')
            # Save as text file for later
            pc_perc_path = os.path.join(pc_outputs_dir,f'{subid}_{thr}_CENS_PC_PERC.txt')
            f = open(pc_perc_path, 'w')
            for cpc_perc in cens_PC_perc_wnan:
                cpc_perc_line = (str(cpc_perc),'\n')
                f.write(' '.join(cpc_perc_line))
            f.close()  
        logger.debug(f'   -Complete: Degree censor & PC cals/perc per thresh')

        ## Now get the average across percentiles using nanmean
        logger.debug(f'  -Averaging PC percentiles')
        all_thresh_PCs_list = []
        for thr in threshold_list:
            cpc_perc_path = os.path.join(pc_outputs_dir,f'{subid}_{thr}_CENS_PC_PERC.txt')
            cpc_perc_vect = np.loadtxt(cpc_perc_path)
            all_thresh_PCs_list.append(cpc_perc_vect.tolist())
        ## convert to numpy array & average
        allPC_percs_wnans = np.array(all_thresh_PCs_list)
        avg_PC_perc = np.nanmean(allPC_percs_wnans, axis=0)
        avg_PC_perc = np.nan_to_num(avg_PC_perc) 
        ## Save subject's avg pc percentage vect as .txt for later
        avg_pc_dir = os.path.join(args.outdir,'final_avg_pc_percs')
        avg_pc_path = os.path.join(avg_pc_dir,f'{subid}_FINAL_AVG_PC_PERCENTAGE.txt')
        f = open(avg_pc_path, 'w')
        for av in avg_PC_perc :
            av_line = (str(av),'\n')
            f.write(' '.join(av_line))
        f.close()

        logger.debug(f'   AVG degree censored PC file calculated and saved.')

    def label_hubs(subid):
        # Step 10: Label parcels above 80 as a hub
        # Load average, degree censored, PC percentile file:
        pc_percs_dir = os.path.join(args.outdir,'final_avg_pc_percs')
        avg_perc_path = os.path.join(pc_percs_dir,f'{subid}_FINAL_AVG_PC_PERCENTAGE.txt')
        avg_vect = np.loadtxt(avg_perc_path)
        ## added per EG - 5.3.21 - get percentile rank of avg perc values
        final_avg_perc_vect = np.asarray([stats.percentileofscore(avg_vect, i) for i in avg_vect])
        cens_avg_vect = np.where(final_avg_perc_vect < 80, 0, final_avg_perc_vect)
        hub_indices = np.asarray(np.nonzero(cens_avg_vect)[0])
        ## add 1 to adjust for python 0 indexing to get parcel # of hub(s)
        hub_indices = hub_indices + 1
        hub_list = hub_indices.tolist()
        logger.debug(f'  -Number of hubs found: {hub_indices.shape[0]}')
        ## SAVE HUB indices for density map creation
        hub_indices_dir = os.path.join(args.outdir,'final_hub_indices')
        hub_indices_path = os.path.join(hub_indices_dir,f'{subid}_HUB_INDICES.txt')
        f = open(hub_indices_path, 'w')
        for hi in hub_list:
            hi_line = (str(hi),'\n')
            f.write(' '.join(hi_line))
        f.close()
        logger.debug('   -Hub index text file saved')

        ### Make basic dlabel of hubs
        logger.debug('  -BONUS: Auto-create hub dlabel.nii file')
        ## Load template dlabel
        img = nb.load(gordon_parcel_path)
        image_data = img.get_fdata()
        image_data = np.asarray(image_data,dtype=np.int16)  # here I retain the shape of the dataobj
        ## get axes info
        axes = [nb.cifti2.cifti2_axes.from_index_mapping(mim) for mim in img.header.matrix]
        ## load my example header and get the original LabelAxis as a dict....
        cifti_header = img.header
        namedmp = cifti_header.get_axis(0)
        brainmodel = cifti_header.get_axis(1)
        label_dict = namedmp.label[0]
        ## edit dict to highlight hubs
        # iterate over parcels in label dict. If not a hub = half transparent. If hub, change color and label name. 
        hub_number = 1
        for lab in label_dict:
            if lab in hub_list:
                # label and color as a hub
                logger.debug(f'   -found hub {hub_number}')
                label_dict[lab] = (f'    -Subject_Hub_{hub_number}', overlay_color) 
                hub_number = hub_number + 1
            else:
                if lab == 0:
                    pass
                else:
                    parc_name = label_dict[lab][0]
                    tup_list = list(label_dict[lab][1])
                    tup_list[-1] = 0
                    color_tup = tuple(tup_list)
                    label_dict[lab] = (parc_name, color_tup)
        ## make the new dict a LabelAxis object
        updated_hub_axis = nb.cifti2.LabelAxis([""], label_dict)
        ## finally, use new label axis and original brain model and write the new file out
        dlabel_output_dir = os.path.join(args.outdir,'final_hubs_dlabels')
        new_header = nb.cifti2.cifti2_axes.to_header([updated_hub_axis, brainmodel])
        new_cifti = nb.cifti2.Cifti2Image(image_data, new_header)
        dlabel_path = os.path.join(dlabel_output_dir,f'{subid}.dlabel.nii')
        new_cifti.to_filename(dlabel_path)
        logger.debug(f'  -Hub surface dlabel.nii file made: {dlabel_path}')

    ### SCRIPT ENTRY ###
    # Iterate over participant data entries 
    for entry in input_data:
        subid = entry.split(' ')[0]
        ts_path = entry.split(' ')[1]
        preprocd_ts = np.loadtxt(ts_path)
        logger.debug(f'{subid}')
        # I. Create z transformed corr mat
        logger.debug(f' -Creating zmat and distance censoring')
        ztrans_mat = make_zmat(preprocd_ts)
        # II. Distance censor zmat
        dist_censored_zmat = dist_censor(ztrans_mat)
        # III. Run infomap
        logger.debug(f' -Running infomap at all thresholds')
        run_infomap(subid, dist_censored_zmat)
        # IV. Calculate PC, get percs, censor, avg percs, save file
        logger.debug(f' -Running Participation Coefficient steps')
        run_pc_steps(subid)
        # V. Identify/Label hubs, make indices file, create hubs dlabel file
        logger.debug(f' -Identifying hub parcels')
        label_hubs(subid)
    # VI. Clean temp files
    if args.nocleanup == True:
        logger.debug(f' -Temp/working files RETAINED. All steps done.')
    else:
        logger.debug(f' -Cleaning FILES from any non-final output folders...')
        all_out_dirs = sorted(glob.glob(os.path.join(args.outdir,'*')))
        for deldir in all_out_dirs:
            if os.path.isdir(deldir):
                if 'final' in deldir:
                    pass
                else:
                    try:
                        shutil.rmtree(deldir)
                    except FileNotFoundError:
                        logger.debug(f'  -Not found? {deldir}')
        logger.debug(f'  -Temp/working files cleaned.')


    full_runtime = time.time() - start_time
    print('\nFull Script Runtime: ' + str(datetime.timedelta(seconds=round(full_runtime))))
if __name__ == '__main__':
    sys.exit(main())
