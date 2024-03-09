function [bcr,bcr_gradient] = gradientTheta(mAdj, thetaArray)
% This function is to calculate the gradient of the critical ratio C* with respect to theta_i
% Input: mAdj: network adjacent matrix  
% Input: thetaArray: this variable is defined to calculate the update rate lambda_i = exp(theta_i)
% Output: bcr: critical benefit-to-cost ratio C*
% Output: bcr_gradient: the gradient of the critical ratio C* with respect to theta_i

%% calculate the coalescence time 
n = length(mAdj); % network size
w = sum(mAdj)';
W = sum(w);
pi = (w/W);
P1 = normalizedLaplacian(mAdj)+speye(n);
P2 = P1 * P1;
P3 = P2 * P1;
rateArray = exp(thetaArray);

C = cartProdRate(mAdj,mAdj,rateArray);
s = reshape(w*w',n^2,1)*sum(rateArray)/n*2;
L = spdiags(sum(C,2),0,length(C),length(C)) - C;

diagcoords = 1:n+1:n^2;
othercoords = sparse(1:n^2);
othercoords(diagcoords)=0;
othercoords = find(othercoords);
s = s(othercoords);
L = L(othercoords,othercoords);
remTime = zeros(n);
remTime(othercoords) = L\s;
%% calculate the critical ratio C* 
eta_1 = sum((P1.*remTime).')*pi;
eta_2 = sum((P2.*remTime).')*pi;
eta_3 = sum((P3.*remTime).')*pi;
bcr = eta_2 / (eta_3 - eta_1);
%% calculate the gradient
bcr_gradient = zeros(n,1);

parfor i=1:n
    gradient_s = reshape(w*w',n^2,1)*2/n;
    gradient_s((i-1)*n+1:i*n) = gradient_s((i-1)*n+1:i*n)-w(i)*w.*remTime(i,:)'+w.*(mAdj(i,:)*remTime)';
    gradient_s(i:n:n^2) = gradient_s(i:n:n^2)-w(i)*w.*remTime(i,:)'+w.*(mAdj(i,:)*remTime)';
    rem_grad_rate_i = zeros(n); % calculate \partial{eta_{ij}}\partial{\lambda_i}
    gradient_s = gradient_s(othercoords);
    rem_grad_rate_i(othercoords) =  L\gradient_s;
    eta_1_partial = sum((P1.*rem_grad_rate_i*exp(thetaArray(i))).')*pi;
    eta_2_partial = sum((P2.*rem_grad_rate_i*exp(thetaArray(i))).')*pi;
    eta_3_partial = sum((P3.*rem_grad_rate_i*exp(thetaArray(i))).')*pi;
    bcr_gradient(i) = 1/(eta_3-eta_1)^2*(eta_2_partial*(eta_3-eta_1)-eta_2*(eta_3_partial-eta_1_partial));
end


end