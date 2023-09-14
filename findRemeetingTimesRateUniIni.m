function remTime = findRemeetingTimesRateUniIni( mAdj,rateArray, varargin )
% Computes remeeting times of the simple random walk on a graph
% given by its weighted adjacency matrix mAdj
if ~issparse(mAdj)
    mAdj= sparse(mAdj);
end
mAdj = max(mAdj,mAdj');
w = sum(mAdj);
n = length(mAdj);
C = cartProdRate(mAdj,mAdj,rateArray);
s = reshape(w'*w,n^2,1)*sum(rateArray)/n*2;
L = spdiags(sum(C,2),0,length(C),length(C)) - C;

diagcoords = 1:n+1:n^2;
othercoords = sparse(1:n^2);
othercoords(diagcoords)=0;
othercoords = find(othercoords);
s = s(othercoords);
L = L(othercoords,othercoords);
remTime = zeros(n);
remTime(othercoords) = L\s;
% remTime=remTime+diag(diag(1+(LL*remTime+remTime*LL')/2));
end

