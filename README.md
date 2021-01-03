# EcoEvoMC

EcoEvoMC is a simulation model for discrete-time population growth with intra- and interspecific competition strength and carrying capacity dependent on an evolving quantitative trait. The full description of the model is in Fielding & Pantel 2020, Eco-evolutionary feedbacks and the maintenance of metacommunity diversity in a changing environment.

# Citations
Fielding, A.P. and Pantel, J.H. 2020. Eco-evolutionary feedbacks and the maintenance of metacommunity diversity in a changing environment. bioRxiv 2020.06.11.145680; doi: https://doi.org/10.1101/2020.06.11.145680

Fielding, A.P., Pantel, J.H. 2020. Eco-Evolutionary feedbacks and the maintenance of metacommunity diversity in a changing environment. Genes 11, 1433.

# Table of Contents
Code for the simulation, and the initial metacommunities used in the simulations, are found in the directory **code**. MATLAB results files produced by the simulation model (the results used in Fielding & Pantel 2020) are in the directory **results**.

# Installation
Use of this code requires MATLAB and access to the files eeed.m, Kx.m, N_ICs_nch68_k3_ic1.mat, N_ICs_nch85_k3_ic1.mat, and N_ICs_nch15_k3_ic1.mat.

# Usage
Download the files eeed.m, Kx.m,  N_ICs_nch68_k3_ic1.mat, N_ICs_nch85_k3_ic1.mat, and N_ICs_nch15_k3_ic1.mat (in the directory code). Place all in the same directory. The file eeed.m is the main simulation file, so executing that will run all 27 simulation conditions described in Fielding & Pantel 2020. Caution that all 27 simulation conditions are set to run in a loop, so the execution of code may take some time. The code produces a set of 27 results (.mat) files that will be saved to the same directory as the code. The function in the file Kx.m is called during the course of eeed.m. The files N_ICs_nch68_k3_ic1.mat, N_ICs_nch85_k3_ic1.mat, and N_ICs_nch15_k3_ic1.mat are the initial communities for each of the 3 levels of niche width (0.68, 0.85, and 1.5). Users can produce their own initial metacommunities, or use the same that were included in the Fielding & Pantel 2020 (Genes) manuscript.

# Credits
This code was written by Aidan P. Fielding, with contributions by J.H. Pantel.

# License
This code adheres to a GNU GPLv3 license.
