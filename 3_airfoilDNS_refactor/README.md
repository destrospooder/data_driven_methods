[Analysis Details (WIP)](https://azielben.quarto.pub/airfoil-dns-analysis/)

If replication of results is desired, here's the file tree:

3_airfoilDNS_refactor
 ┣ dataset
 ┃ ┣ airfoilDNS_a25f0p05.h5
 ┃ ┣ airfoilDNS_grid.h5
 ┃ ┗ airfoilDNS_parameters.h5
 ┣ dmd_analysis_figs
 ┃ ┣ mode_amplitudes.png
 ┃ ┣ ux_dmd_modes.png
 ┃ ┗ uy_dmd_modes.png
 ┣ pod_analysis_figs
 ┃ ┣ squared_sv.png
 ┃ ┣ squared_sv_truncated.png
 ┃ ┣ temporal_amplitudes.png
 ┃ ┣ ux_reconstruction.png
 ┃ ┣ ux_spatial_modes.png
 ┃ ┣ uy_reconstruction.png
 ┃ ┗ uy_spatial_modes.png
 ┣ sindy_figs
 ┃ ┣ comparison.png
 ┃ ┣ sindy_amps (lambda = 0.0135).png
 ┃ ┗ sindy_amps (lambda = 0.025).png
 ┣ svd_store
 ┃ ┣ svdeez_s.npy
 ┃ ┣ svdeez_u.npy
 ┃ ┗ svdeez_vh.npy
 ┣ __pycache__
 ┃ ┣ pod_analysis.cpython-311.pyc
 ┃ ┗ var_import.cpython-311.pyc
 ┣ dmd_analysis.py
 ┣ pod_analysis.py
 ┣ sindy.py
 ┗ var_import.py