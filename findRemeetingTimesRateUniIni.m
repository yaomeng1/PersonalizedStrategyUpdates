function remTime = findRemeetingTimesRateUniIni( mAdj,rateArray)
% Compute eta_ij shown in Eq. (1) in the main text
% according to Eq. (7) in the main text
% Input: mAdj: adjacent matrix
% Input: rateArray: individual update rates
% Output: remTime: coalscence time

if ~issparse(mAdj)
    mAdj= sparse(mAdj);
end
mAdj = max(mAdj,mAdj');
w = sum(mAdj);
n = length(mAdj);

% To solve the linear equations Lx = s, generate the coefficient matrix L
C = cartProdRate(mAdj,mAdj,rateArray);
s = reshape(w'*w,n^2,1)*sum(rateArray)/n*2;
L = spdiags(sum(C,2),0,length(C),length(C)) - C;


% solve Lx = s to abtain eta_{ij} (i is not equal to j) 
diagcoords = 1:n+1:n^2;
othercoords = sparse(1:n^2);
othercoords(diagcoords)=0;
othercoords = find(othercoords);
s = s(othercoords);
L = L(othercoords,othercoords);
remTime = zeros(n);
remTime(othercoords) = L\s;
end

