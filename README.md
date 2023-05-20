## <p align="center">[**Resting-state cortical hubs in youth organize into four categories**](http://doi.org/10.1016/j.celrep.2023.112521)</p>
### <p align="center">[DOI: 10.1016/j.celrep.2023.112521](http://doi.org/10.1016/j.celrep.2023.112521)</p>


### Authors:
Damion V. Demeter*, Evan M. gordon, Tehila Nugiel, AnnaCarolina Garza, Tyler L. Larguinho, Jessica A. Church
* Correspondence: ddemeter@ucsd.edu

![Graphical Abstract](https://ars.els-cdn.com/content/image/1-s2.0-S2211124723005326-fx1_lrg.jpg)

---
## <p align="center">**Scripts and other resources for replication of this work**</p>   

### Scripts:
- Cortical Hub Identification Script: Identify_Hubs.py
  - This script identifies hub parcels and associated files.
- Hub Density Map Creation Script: Hub_Denisty_Map.py
  - This script creates a hub density map across all participants.
  - As in Figure 1 (A)
- Hub Profiles Script: Create_Hub_Profiles.py
  - This script creates hub profiles and a hub profile correlation matrix that can be used to cluster hub profiles into categories.
  - We recommend using the Louvain algorithm and methods described in this paper to identify hub categories, but other clustering methods can be used with this output.
  - The optional correlation matrix plot isn't as useful until after hub category clustering, but can be used to make a matrix as in Figure S3 (B). 

### Python Requirements: these scripts are written in python 3.9.7.
- Python package requirements:
```
Brain Connectivity Toolbox for Python (bct)
NetworkX (networkx)
NiBabel (nibabel)
Numpy (numpy)
SciPy (scipy) 
```
### Other Requirements:
1. Timeseries should be fully processed, motion corrected, etc (appropriate steps for your chosen processing pipeline). This current script requires that timeseries are created from the [Gordon 333 Parcel set](https://balsa.wustl.edu/2Vm69) and exported to a .txt file. (This script does NOT handle vertex-wise data)
 - .txt file format should be 333 rows by X time/TR of scan, exported to .txt file using the [Gordon 333 Cortical Parcel set](https://balsa.wustl.edu/2Vm69)
 - .txt files can be created from dense timeseries files (.dtseries.nii) by using the "wb_command -cifti-parcellate" command from [Connectome Workbench](https://www.humanconnectome.org/software/workbench-command)
2. A .txt file for your participant list is required that has (space separated) <participant ID> <path to timeseries.txt file>. (see example in this repository for format)
3. Infomap should be installed locally and able to be called from the command line. (Tested with infomap versions 1.9.0 & 2.6.0 - Other versions may need adjusted clu file editing.)  

### Basic Outputs:
1. Identify_Hubs.py
  - 3 folders with your final data.
    - /final_avg_pc_percs/ - Average participation coefficient percentile for each cortical hub (In Gordon 333 set parcel number order)
    - /final_csv_outputs/ - All distance censored zmat files
    - /final_hub_indices/ - Parcel indices (Parcel #) for parcels identified as a hub 
    - /final_hubs_dlabels/ - dlabel.nii files with shaded hub parcels, used for visualization

2. Hub_Density_Map.py
   - /<name>_Gordon333_Hub_Counts.txt - Cumulative count of how many hubs were identified for each parcel
   - /<name>_hubs_density_map.pscalar.nii - Hub density map for viewing in workbench view
 
3. Create_Hub_Profiles.py
   - /final_conn_profiles/ - Hub profiles, per subID, for all identified hubs. (Figure 6 (b))
   - /<name>_Hub_Profile_Correlations.png - Hub profile correlation matrix plot (Not very useful until after profile clustering)


- Beyond this, see the -h (help) argument in the scripts for full details of required arguments. 




