function L = normalizedLaplacian(mAdj)
% Computes the normalized Laplacian of the graph given by the weighted
% adjacency matrix mAdj
if ~issparse(mAdj)
    mAdj = sparse(mAdj);
end

mAdj = max(mAdj,mAdj');

n = length(mAdj);
s = sum(mAdj,2);
[nod1,nod2,w]=find(mAdj);
L = sparse(nod1,nod2,w./s(nod1),n,n) - speye(n);


end

