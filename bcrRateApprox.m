function bcr_approx = bcrRateApprox(mAdj, rateArray)
% This function is used to compute the approximation in Eq. (3) in the main text
% Input: mAdj: a ajacent matrix (n*n) of any given network    
% Input: rate: a column array (n*1) of individual update rates
% Output: bcr_approx: the approximation result of the critical
% benefit-to-cost ratio (above which cooperation is favoured) of any given network structure and any individual
% update rates.

    n = length(mAdj); % network size
    w = sum(mAdj)';  % column array of nodes' degree
    W = sum(w);
    P1 = mAdj./w;   % P10(i,j): probability of one-step random walk from i to j 
    P2=P1*P1;     % P20(i,j): probability of two-step random walk from i to j
    mu2 = mean(w.*w); % second moment of the degree distribution
    
    rateArray = n*rateArray/sum(rateArray);	% normalization of the individual update rates with overall rates of n.
    lambda_ij = rateArray'+rateArray;    
    delta_rate = (rateArray'-rateArray)./(rateArray+rateArray');
    eta_mean = (1/(mu2*n)*sum(w*w'./lambda_ij,"all")-1/n*sum(1./(2*rateArray)))/...
        (1+1/(mu2*n)*sum(w.*delta_rate.*mAdj,"all"));  % the mean value of coalescence time starting from any two different nodes
    delta_remInf_approx=-eta_mean*sum(w.*delta_rate.*mAdj,"all")/W^2; %  $Delta_{\tilde{\eta}^{(\infty)}$ in Equation (S17) in Supplementary Information
    delta_rem2_approx= -eta_mean*(w/W)'*sum(P1.*delta_rate.*P1',2); %  $Delta_{\tilde{\eta}^{(2)}$ in Equation (S18) in Supplementary Information
    delta_rem3_approx= -eta_mean*(w/W)'*sum(P2.*delta_rate.*P1',2); %  $Delta_{\tilde{\eta}^{(3)}$ in Equation (S19) in Supplementary Information
    
    % Calculate Eq. (3) for unweighted networks in the main text (and also Equation (S20) for weighted networks in the Supplementary Information)
    nom1 = 1/(mu2*n)*sum(w*w'./lambda_ij,"all")...
        -1/W*(w'*(1./(2*rateArray))+sum(w.*P1./lambda_ij,"all")); 
    dom1 = (w'*diag(P2)/n)/(mu2*W)*sum(w*w'./lambda_ij,"all")...
        -1/W*sum(w.*(P1+P2)./lambda_ij,"all");
    bcr_approx = (nom1 -delta_rem2_approx+W^2/(mu2*n)*delta_remInf_approx)/...
        (dom1-delta_rem2_approx-delta_rem3_approx+W*(w'*diag(P2)/n)/mu2*delta_remInf_approx); 
end