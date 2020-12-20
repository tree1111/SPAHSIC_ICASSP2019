%--------------------------------------------------------------------------
% This function calculates the affinity matrix
% U: the orthogonal basis vectors of each superpixel
% N: the number of superpixels

% A: the affinity matrix
%--------------------------------------------------------------------------
function A = affinityHSI(U,N)

A = zeros(N,N);
r = size(U{1},2);                           % the dimension of subspace 
for i = 1:N
    for j = i+1:N
        aff = sum(sum((U{i}'*U{j}).^2));
        d_2 = r - aff;
        kernel = exp(-d_2/7);               % 2*segma^2 = 7
        A(i,j) = kernel;
        A(j,i) = kernel;                    % the affinity matrix is symmetric
    end
end 
 
end