clear
clc
load("sf_100_k6.mat") 
mAdj=A_sf; % adjacent marix of network
n = length(mAdj); % network size
[rate_process,bcr_array] = OptUpRat(mAdj);  % optimize individual update rate given a network structure

rate_process = rate_process./sum(rate_process,1)*n; % normalize the update rate with average rate of 1
optimal_rate = rate_process(:,end);  % optimal update rate 
optimal_bcr = bcr_array(end);  % corresponding critical benefit-to-cost ratio 