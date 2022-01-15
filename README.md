# A spatiotemporal route to understanding metal halide perovskite crystallization

<center><h3>Mansoor Ani Najeeb<sup>1</sup>, Rodolfo Keesey<sup>1</sup>, Margaret Zeile<sup>1</sup>, Zhi Li<sup>2</sup>, Venkateswaran Shekar<sup>3</sup>, Nicholas Leiby<sup>4</sup>, Matthias Zeller<sup>5</sup>, Emory Chan<sup>2</sup>, Joshua Schrier<sup>6</sup>, Alexander J. Norquist<sup>1</sup></h3></center>
<br>

1. Department of Chemistry, Haverford College, 370 Lancaster Avenue, Haverford, Pennsylvania 19041, USA
2. Molecular Foundry, Lawrence Berkeley National Laboratory, 1 Cyclotron Road, Berkeley, California 94720, USA
3. Department of Computer Science, Haverford College, 370 Lancaster Avenue, Haverford, Pennsylvania 19041, USA
4. Two Six Technologies, 901 N. Stuart Street, Arlington, Virginia, 22203, USA
5. Department of Chemistry, Purdue University, West Lafayette IN 47907, USA
6. Department of Chemistry, Fordham University, 441 E. Fordham Road, The Bronx, New York, 10458, USA

This repo contains code, data and jupyter notebook related to RAPID_2.

Binder to run RAPID2.ipynb: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/darkreactions/rapid2/HEAD?labpath=RAPID2.ipynb)


## Repository contents
The following sections indicate the folders which contain code and related data

### Jupyter notebooks

1. RAPID2.ipynb - Notebook containing all visualizations and reaction outcomes

### Raw data
All raw data files are located in the ```data``` folder

1. Diffusion_rate_and_crystal_height_all_CSV- Contains .csv files of extracted solvent heights and crystal growth time
2. cifs - Contains the Crystallographic Information Files
3. images - Contains side vial images of each experiment performed
4. xrd/xy - Contains xy files for XRD data
5. 0042.perovskitedata_RAPID.csv - Escalate generated data file including 8 experimental features (with "_rxn_" as header prefix) and 67 chemical features (with "_feat_" as header prefix). The detailed explanations of these features are listed in "Explanation of Features-Descriptors" section in "Perovskite Dataset Description.pdf". This CSV file is used in visualization and machine learning. 
6. 0042.perovskitedata_RAPID_full.csv - This escalate generated data file contains the same experiments as "0042.perovskitedata_RAPID.csv" but has all 787 features, including additional "_raw_" features describing experiment details (see "Explanation of Features-Descriptors" section in "Perovskite Dataset Description" for the explanations of "_raw_" prefix). The csv file is not used for visualization or machine learning. 
7. image_list.json - Keeps track of all image files in the image folder
8. ml_data.pkl - Python pickle file containing ML results
9. inventory.csv - Chemical inventory data
10. organic_inchikey.csv - Inchi keys and chemical names
11. s_spaces.json - Co-ordinates of state space for each amine

### Scripts in src folder
The following python scripts are used in the RAPID.ipynb notebook to generate visualizations

1. diffusion_coeff - Matlab code for running diffusion model
2. 
