%--------------------------------------------------------------------------
% This function calculates the first r principal components for 
% each superpixel.
% HSI: L*MN 2D hyperspectral images
% superclass: desired number of superpixels
% L: scaling factor
% r: subspace dimension

% U: the orthogonal basis vectors of each superpixel
%--------------------------------------------------------------------------
function U = orthvector(HSI,super_class,supernum,r)

U =cell(1,supernum);
super_index = reshape(super_class,1,size(super_class,1) * size(super_class,2));
for i = 1:supernum
    index  = find(super_index == i);
    temp = HSI(:,index);
    [D,T] = size(temp);
    [u,~,~]=svds(temp,D,'L');
    u = u(:,1:r);
    U{i} = u;
end

end