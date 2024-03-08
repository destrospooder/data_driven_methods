## Data-Driven Analysis of Low Reynolds # Pitching Airfoil Direct Numerical Simulations

[Analysis Details (WIP)](https://azielben.quarto.pub/airfoil-dns-analysis/)

Data was sourced from Aaron Towne's [database for reduced-complexity modeling of fluid flows](https://deepblue.lib.umich.edu/data/collections/kk91fk98z). If replication of results is desired, here`s the file tree:

`3_airfoilDNS_refactor` <br/> 
 ┣ `dataset` <br/> 
 ┃ ┣ `airfoilDNS_a25f0p05.h5` <br/> 
 ┃ ┣ `airfoilDNS_grid.h5` <br/> 
 ┃ ┗ `airfoilDNS_parameters.h5` <br/> 
 ┣ `dmd_analysis_figs` <br/> 
 ┃ ┣ `mode_amplitudes.png` <br/> 
 ┃ ┣ `ux_dmd_modes.png` <br/> 
 ┃ ┗ `uy_dmd_modes.png` <br/> 
 ┣ `pod_analysis_figs` <br/> 
 ┃ ┣ `squared_sv.png` <br/> 
 ┃ ┣ `squared_sv_truncated.png` <br/> 
 ┃ ┣ `temporal_amplitudes.png` <br/> 
 ┃ ┣ `ux_reconstruction.png` <br/> 
 ┃ ┣ `ux_spatial_modes.png` <br/> 
 ┃ ┣ `uy_reconstruction.png` <br/> 
 ┃ ┗ `uy_spatial_modes.png` <br/> 
 ┣ `sindy_figs` <br/> 
 ┃ ┣ `comparison.png` <br/> 
 ┃ ┣ `sindy_amps (lambda = 0.0135).png` <br/> 
 ┃ ┗ `sindy_amps (lambda = 0.025).png` <br/> 
 ┣ `svd_store` (empty directory necessary to run `pod_analysis.py`)<br/> 
 ┣ `dmd_analysis.py` <br/> 
 ┣ `pod_analysis.py` <br/> 
 ┣ `sindy.py` <br/> 
 ┗ `var_import.py` <br/> 