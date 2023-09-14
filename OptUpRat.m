function [optimal_rate,optimal_bcr] = OptUpRat(mAdj)
% This function is used to find the optimal update rate for each node to
% minimize the critical ratio C*
% the results returns the optimal rate for each node: optimal_bcr
% and the corresponding critical ratio C*: optimal_bcr


n = length(mAdj);
epoch = 10000;
epsilon = 1e-6;
eta = 1; % learning rate
alpha = 0.9; % decay rate


theta = zeros(n,epoch);
bcr_array =  zeros(epoch,1);
gs_accumulate = 0;   
for idx=1:epoch-1
    [bcr, grad] = jacobianTheta(mAdj,theta(:,idx));
    bcr_array(idx) = bcr;
    gs_accumulate = alpha*gs_accumulate + (1-alpha)*grad.^2;
    theta(:,idx+1) = theta(:,idx)-eta./sqrt(gs_accumulate+epsilon).*grad;
    if idx>1
        if sum(abs(theta(:,idx+1)-theta(:,idx)))/n<1e-6
        break
        end
    end
end

optimal_rate = exp(theta(:,idx));
optimal_bcr = bcr_array(idx);

end

