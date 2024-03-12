# PersonalizedStrategyUpdates


Code for paper "Evolution of collective cooperation under personalised strategy updates".

The codes include program for numerical simulations (written by Python 3.8.5) and functions for theoretical calculations (written by MathWorks Matlab R2021a). The details for these codes are as follow:

#### Numerical simulation

- "fixation_probability_personalized_rate.py" is used to numerically obtain the fixation probability of cooperation ($\rho_C$) through benefit-to-cost ratio in Fig. 2a-d in the main text. 

  - Input: 

    - the adjacent matrix of any given network ("sf_100_k6.mat" in this repository)
    - individual update rates
    - an array of benefit-to-cost ratio

  - Output: an array of corresponding fixation probability of cooperation ($\rho_C$)

  - **Demo**: We provide the output file of the simulation results with a scale-free network used in Fig. 2a-c, under the settings of update rates ($\lambda_i=1, 1/k_i, k_i$) respectively. The corresponding file names are:

    - "output_sf_n100_k6_IdenticalRate_lambda_1.mat"
    - "output_sf_n100_k6_IdenticalRate_lambda_1_k.mat"
    - "output_sf_n100_k6_IdenticalRate_lambda_k.mat".

    Each file contains an array of benefit-to-cost ratio (`b_array`), and an array of corresponding $\rho_C$ (`rhoc_array`)，which is exactly the results with the scale-free network in Fig. 2a-c. 
   - **Requirements**: see the list of packages needed in requirements.txt.
   - **Other instructions**: Please set the variable `cpu_cores_num` to the number of CPU cores on the local computer.
     
 To reproduce the plots in Fig. 2a-c, run this file and calculate the theoretical $C^*$ using "getBCratioRateUniIni.m" as below. 


#### Theoretical calculation

- "getBCratioRateUniIni.m" is used to calculate the theoretical critical ratio $C^*$ using Eq. (1) in the main text, where the evolutionary game process starts from a single cooperator placed uniformly at random on the network.
  - Input: 
    - the adjacent matrix of any given network
    - individual update rates
  - Output: Critical benefit-to-cost ratio $C^*$
  - **Demo**: We provide the demo file "demo_bcr_PersonalizedRate.m" to calculate theoretical $C^*$ of the scale-free network in Fig. 2a-c. 
- "bcrRateApprox.m" is used to calculate the approximated results shown in Fig. 5b using Eq. (3)  in the main text. 
  - Input: 
    - the adjacent matrix of any given network
    - individual update rates
  - Output: Approximated critical benefit-to-cost ratio $C^*$ 
  - **Demo**: The approximated $C^*$ of the scale-free network in Fig. 2a-c is provided in file "demo_bcr_PersonalizedRate.m".  The accuracy can be checked with the theoretical critical ratio calculated using Eq. (1) in the same demo.
- "OptUpRat.m" is used to optimize the update rate for each individual to minimize the critical ratio $C^*$  given the adjacent matrix of any network;
  - Input: the adjacent matrix of any given network
  - Output: 
    - An array of the critical ratio $C^*$ (iteration steps * 1)
    - An array of update rates (network size * iteration steps)
  - **Demo**: We provide the demo file "demo_OptUpRat.m" to optimize update rates of a given scale-free network (Fig. 5d, e in the main text), and the output file "output_sf_100_k6_OptimalUpdateRate.mat" to show the results of optimization.  `rate_process` records the iteration process of update rates, and `bcr_array` records the iterations of critical ratio. Note that the optimal rates can be obtained by taking the last column, and the corresponding minimal  $C^*$ is the last element of `bcr_array`.


The other functions are subfunctions needed to run the above three functions, and the function of each file is depicted at the beginning of each file.



Note: to run the code, make sure that all files are in the same folder.

