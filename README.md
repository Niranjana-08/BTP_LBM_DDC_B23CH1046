# Lattice Boltzmann Simulation of Double-Diffusive MHD Natural Convection

Bachelor's Thesis Project, Department of Chemical Engineering, IIT Jodhpur
Author: V Niranjana (B23CH1046) | Supervisor: Dr. Krunal Madhukar Gangawane

This repository holds the solver codes, simulation output, plotting scripts and figures behind the project report and presentation. The solver is the in-house lattice Boltzmann code of the supervisor's research group, used here with only the changes listed in Section 4.

## 1. Repository structure

```
README.md
code_1_base_solver/
    code1_base_solver.cpp                 Code 1: base solver, one case per compilation
code_2_targeted_Ra_Ha_study/
    code2_targeted_Ra_Ha_study.cpp        Code 2: final nine-run Ra-Ha study
results_csv/
    summary_results_Ra_Ha_study.csv       Nu, Sh and iterations of the nine runs
    grid_csv/                             one CSV per run (columns i, j, T, C, u, v)
post_processing/
    LBM_Plot2_Colab_notebook.ipynb        Google Colab plotting notebook (run on the nine final cases)
    plot_fields.py                        isotherms + concentration/streamlines (script version)
    plot_analysis.py                      stream function, mid-height profiles, Nu/Sh plots
    field_io.py                           helper used by the two scripts
results_plots/
    field_plots/                          9 two-panel figures (isotherms; concentration + streamlines)
    stream_function/                      stream function for Ra = 1e3, 1e4, 1e5 (Ha = 0, 50, 100)
    midheight_profiles/                   T, C and v at mid-height
    transport_rates/                      Nu and Sh versus Ra and versus Ha
    schematics/                           cavity schematic, D2Q9 lattice
```

## 2. Running the solver

Needs a C++ compiler (g++). Code 2 runs the nine cases one after the other until the residual is below 1e-8; allow roughly 1 to 2 hours depending on the machine.

```
cd code_2_targeted_Ra_Ha_study
g++ -O3 -o code2 code2_targeted_Ra_Ha_study.cpp
./code2
```

It creates `output_csv_data/` with `summary_results_Ra_Ha_study.csv` and the nine `grid_Ra<Ra>_Ha<Ha>_N0_Pr0.054.csv` files. Code 1 is compiled the same way, but its parameters (Ra, Ha, Le, N, Pr) are set with `#define` lines at the top of the file and must be edited before compiling.

Summary file columns: `Run, Ra, Ha, N, Pr, Le, avgNu_left, avgSh_left, avgNu_bottom, avgSh_bottom, Converged_Iter`. `avgSh_left` is negative because concentration increases from left to right (the report uses its magnitude). The bottom-wall values are zero because the top and bottom walls are insulated and impermeable.

## 3. Making the plots

Google Colab: put the nine grid CSV files in one Google Drive folder, open `LBM_Plot2_Colab_notebook.ipynb` in Colab (Drive: right-click, Open with, Google Colaboratory), set the input folder (`data_dir`) and the output folder (`output_dir`) in the first cell, and run it. The notebook reads only files starting with `grid_`, takes Ra, Ha, N and Pr from each file name and saves one two-panel figure per run in the output folder.

Local Python (numpy, pandas, matplotlib):

```
cd post_processing
python plot_fields.py ../results_csv/grid_csv ../results_plots/field_plots
python plot_analysis.py ../results_csv/grid_csv ../results_plots/regenerated --summary ../results_csv/summary_results_Ra_Ha_study.csv
```

## 4. Final results (Pr = 0.054, Le = 1, N = 0)

| Ra | Ha = 0 | Ha = 50 | Ha = 100 |
|---|---|---|---|
| 1e3, Nu (left wall) | 1.1134 | 1.0055 | 1.0052 |
| 1e4, Nu (left wall) | 2.0854 | 1.0425 | 1.0085 |
| 1e5, Nu (left wall) | 3.7451 | 2.1488 | 1.2487 |

Sherwood numbers equal the Nusselt numbers to within about 0.1% (see the summary CSV).

## 5. Provenance of the data

The nine grid files in `results_csv/grid_csv/` were generated with a restructured implementation of the same algorithm and equations as Code 2. Code 2 in this repository was checked against it for the case Ra = 1e3, Ha = 50: both gave 24 981 iterations and Nu = 1.00548 (1.005484), with the fields agreeing to CSV rounding. Running Code 2 as described above regenerates the nine cases.

