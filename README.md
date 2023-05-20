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
- Hub Density Map Creation Script: Hub_Denisty_Map.py
- Hub Categorization Script: Categorize_Hubs.py

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
- Timeseries should be fully processed, motion corrected, etc (appropriate steps for your chosen processing pipeline). This current script requires that timeseries are created from the [Gordon 333 Parcel set](https://balsa.wustl.edu/2Vm69) and exported to a .txt file. (This script does NOT handle vertex-wise data)
 - .txt file format should be 333 rows by X time/TR of scan, exported to .txt file using the [Gordon 333 Cortical Parcel set](https://balsa.wustl.edu/2Vm69)
 - .txt files can be created from dense timeseries files (.dtseries.nii) by using the "wb_command -cifti-parcellate" command from [Connectome Workbench](https://www.humanconnectome.org/software/workbench-command)
- A .txt file for your participant list is required that has (space separated) <participant ID> <path to timeseries.txt file>. (see example in this repository for format)

### Basic Outputs:
1. Identify_Hubs.py
  - 3 folders with your final data.
    - /final_avg_pc_percs/ - Average participation coefficient percentile for each cortical hub (In Gordon 333 set parcel number order)
    - /final_hub_indices/ - Parcel indicies (Parcel #) for parcels identified as a hub 
    - /final_hubs_dlabels/ - dlabel.nii files with shaded hub parcels, used for visualization

2. Hub_Density_Map.py
 
 
3. Categorize_Hubs.py



- Beyond this, see the -h (help) argument in the scripts for full details of required arguments. 




