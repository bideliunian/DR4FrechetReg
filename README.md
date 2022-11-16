# SDR4FrechetReg
The code for the paper 'Dimension Reduction for Fr' then echet Regression'.

### Supporting software requirements

R version 4.0 or later.

### Libraries and dependencies used by the code

R packages to run Frechet SDR methods:

* tictoc v 1.1
* frechet v 0.2.0

R packages to analyze mortality data and reproduce figures:

* ggplot2 v 3.4.0
* dplyr v 1.0.10
* plotly v 4.10.1
* pracma v 2.4.2
* RColorBrewer v 1.1-3
* gridExtra v 2.3
* htmlwidgets v 1.5.4

### Folder Structure: 

* `./DR4FrechetReg/Functions/`  code for all functions used in the paper.
* `./DR4FrechetReg/DistData/`  code to reproduce simulations in Section 6.2 Scenario I for distributions.
* `./DR4FrechetReg/SPDmatData/` code to reproduce simulations in Section 6.3 Scenario II for SPD matrix data.
* `./DR4FrechetReg/SphereData/` code to reproduce simulations in Section S.4.2 Scenario III for Spherical data.
* `./DR4FrechetReg/MortalityData/` code to reproduce data analysis for human mortality data in Section 7.
* `./DR4FrechetReg/StrokeData/` code to reproduce data analysis for stroke data in Section S.6.


### Reproducibility workflow

#### Simulations:

* Under `./DR4FrechetReg/DistData/`:
  
  + Run **`benchmark_dist.R`** to get the benchmark error for all models.
  
  + Run **`CV4kernelbd.R`** to get the optimal kernel type and tuning parameter.
  
  + Run the bash commands `bash dist_model12a.bashrc`, `bash dist_model12b.bashrc`, and `bash dist_model34.bashrc` separately  to send 100 separate PBS jobs at each time. Simulation results will be stored under folder `./DR4FrechetReg/DistData/Results/`.
  
  + Run **`result_read.R`** to create Table 1 and Table 2.
  
  + For order determination comparison in Section S.5 in Supplementary Material, run the bash commands `bash dist_model12_a_order.bashrc` and `bash dist_model34_order.bashrc` separately. Simulation results will be stored under folder `./DR4FrechetReg/DistData/OrderResults/`. Then run **`order_result_read.R`** to create Table S.3.

* Under `./DR4FrechetReg/SPDmatData/`:
  
  + Run **`CV4kernelbd.R`** to get the optimal kernel type and tuning parameter.
  
  + Run the bash commands `bash spd_a.bashrc` and `bash spd_b.bashrc` separately  to send 100 separate PBS jobs at each time. Simulation results will be stored under `./DR4FrechetReg/SPDmatData/Results/`.
  
  + Run the script **`result_read.R`** to create Table 3.

  + Run the bash commands `bash spd_order.bashrc`. Simulation results will be stored under `./DR4FrechetReg/SPDmatData/OrderResults/`. Then run **`order_result_read.R`** to create Table S.4.
  
* Under `./DR4FrechetReg/SphereData/`:
  
  + Run **`CV4kernelbd.R`** to get the optimal kernel type and tuning parameter.
  
  + Run the bash commands `bash sphere_a.bashrc` and `bash sphere_b.bashrc` separately  to send 100 separate PBS jobs at each time. Simulation results will be stored under `./DR4FrechetReg/SphereData/Results/`. Run the script **`result_read.R`** to create Table S.1.
  
  + Run the bash commands `bash sphere_order.bashrc`. Simulation results will be stored under `./DR4FrechetReg/SphereData/OrderResults/`. Then run **`order_result_read.R`** to create Table S.5.

* Under `./DR4FrechetReg/Visualization/`, run R scripts **`plot_dist.R`**, **`plot_spd.R`** and **`plot_sphere.R`** separately to generate Figures 1, 2 and S.2.

For a single run of the experiments, for example, Model I-1-(a), execute `Rscript training_dist_model12a.R Args` in a terminal, where `Args` controls the seed.

#### Data Applications:

* Under `./DR4FrechetReg/MortalityData/`, run **`mortality.R`** to do the analysis for mortality data and generate Figures 1 and 4.

* Under `./DR4FrechetReg/StrokeData/`, run **`stroke.R`** to do the analysis for stroke data and generate Figures S.3.
