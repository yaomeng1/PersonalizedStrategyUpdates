function [rate_process,bcr_array] = OptUpRat(mAdj)
% This function is used to find the optimal update rate for each node to
% minimize the critical ratio C*, following the schedule of Algorithm OptUpRat in
% the main text.
% Input: mAdj: network adjacent matrix 
% Output: rate_process: the iteration process of update rates for each node
% Output: bcr_array: the iteration process of the corresponding critical ratio C*


n = length(mAdj); %network size
epoch = 10000;  % maximum steps of iterations

epsilon = 1e-6; % a small constant to prevent division by zero
eta = 1; % learning rate
alpha = 0.9; % decay rate


theta = zeros(n,epoch); % initialzation according to Algorithm OptUpRat in the main text
bcr_array =  zeros(epoch,1); 
gs_accumulate = 0;   % initialzation of the accumulated squared gradient
for idx=1:epoch-1  
    [bcr, grad] = gradientTheta(mAdj,theta(:,idx)); % compute the gradient of C* with respect to theta_i for each i 
    bcr_array(idx) = bcr;
    gs_accumulate = alpha*gs_accumulate + (1-alpha)*grad.^2; % moving average the squared gradient
    theta(:,idx+1) = theta(:,idx)-eta./sqrt(gs_accumulate+epsilon).*grad; % gradient descent
    if idx>1
        if sum(abs(theta(:,idx+1)-theta(:,idx)))/n<1e-6  % end iteration if the variation of theta_i small enough
        break
        end
    end
end

rate_process = exp(theta(:,1:idx)); % obtain the final optimal update rate lambda_i = exp(theta_i)
bcr_array = bcr_array(1:idx); % obtain the minimal critical ratio C*

end

