function bcr = getBCratioRateUniIni( mAdj, rateArray)
% Computes the critical benefit-to-cost ratio C* for favouring cooperation 
% starting from any network and update rates according to Eq. (1) in the main text. 
% Input: mAdj: input network adjacent matrix
% Input: rateArray: column array of individual update rates
% Output: bcr: output critical benefit-to-cost ratio in Eq. (1) in the main text


n = length(mAdj); % network size
w = sum(mAdj);  % degree of nodes
W = sum(w);
pi = (w/W).';     

Rem = findRemeetingTimesRateUniIni(mAdj,rateArray); % compute coalescence time according to Eq. (7) in the main text
Rem = Rem - diag(diag(Rem));

P1 = normalizedLaplacian(mAdj)+speye(n);  % P1(i,j): probability of one-step random walk from i to j 
P2 = P1 * P1;  % P2(i,j): probability of 2-step random walk from i to j 
P3 = P2 * P1;  % P3(i,j): probability of 3-step random walk from i to j 
 
eta_1 = sum((P1.*Rem).')*pi; % eta^{(1)} in Equation (S6) in the supplementary information
eta_2 = sum((P2.*Rem).')*pi; % eta^{(2)} in Equation (S6) in the supplementary information
eta_3 = sum((P3.*Rem).')*pi; % eta^{(3)} in Equation (S6) in the supplementary information

bcr = eta_2 / (eta_3 - eta_1); % calculate the critical ratio according to Eq. (1) in the main text

end

