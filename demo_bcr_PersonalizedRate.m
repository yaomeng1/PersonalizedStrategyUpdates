clear
clc
load("sf_n100_k6_idx0.mat") 
mAdj=A_sf;  % adjacent matrix of a scale-free network
n = length(mAdj); % network size
k_array = sum(mAdj,2);  % nodes' degree
rate_1_ki = 1./k_array;
rate_ki = k_array;
rate_identical = ones(length(mAdj),1);
%% Theoretical calculations of critical benefit-to-cost-ratio with different update rates
bcr_identical_theoretical = getBCratioRateUniIni(mAdj,rate_identical); % \lambda_i=1
bcr_1_ki_theoretical = getBCratioRateUniIni(mAdj,rate_1_ki); % \lambda_i=1/k_i
bcr_ki_theoretical = getBCratioRateUniIni(mAdj,rate_ki);  % \lambda_i=k_i

%% Critical benefit-to-cost-ratio obtained by numerical simulations
% load numerical results of fixation probability (rhoc) thourgh
% benefit-to=cost ratio (b_array)
load('output_sf_n100_k6_idx0_IdenticalRate_lambda_1.mat')
b_array_identical = b_array;
rhoc_array_indentical = rhoc_array;
load('output_sf_n100_k6_idx0_PersonalizedRate_lambda_1_k.mat')
b_array_1_k = b_array;
rhoc_array_1_k = rhoc_array;
load('output_sf_n100_k6_idx0_PersonalizedRate_lambda_k.mat')
b_array_k = b_array;
rhoc_array_k = rhoc_array;

bcr_identical_numerical = intersection(b_array_identical,rhoc_array_indentical,1/n); % % \lambda_i=1
bcr_1_ki_numerical = intersection(b_array_1_k,rhoc_array_1_k,1/n); % \lambda_i=1/k_i
bcr_ki_numerical = intersection(b_array_k,rhoc_array_k,1/n);  % \lambda_i=k_i

%% Approximation results 
bcr_identical_approx = bcrRateApprox(mAdj,rate_identical);
bcr_1_ki_approx = bcrRateApprox(mAdj,rate_1_ki);
bcr_ki_approx = bcrRateApprox(mAdj,rate_ki);
