clear
clc
load("sf_100_k6.mat") 
mAdj=A_sf;  % adjacent matrix of a scale-free network

k_array = sum(mAdj,2);  % nodes' degree
rate_1_ki = 1./k_array;
rate_ki = k_array;
rate_identical = ones(length(mAdj),1);
%% Theoretical critical benefit-to-cost-ratio with different update rates
bcr_identical = getBCratioRateUniIni(mAdj,rate_identical); % % \lambda_i=1
bcr_1_ki = getBCratioRateUniIni(mAdj,rate_1_ki); % \lambda_i=1/k_i
bcr_ki = getBCratioRateUniIni(mAdj,rate_ki);  % \lambda_i=k_i


%% Approximation results 
bcr_identical_approx = bcrRateApprox(mAdj,rate_identical);
bcr_1_ki_approx = bcrRateApprox(mAdj,rate_1_ki);
bcr_ki_approx = bcrRateApprox(mAdj,rate_ki);
