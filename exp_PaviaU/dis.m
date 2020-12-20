%--------------------------------------------------------------------------
% This function caculates the distance between two pixels.
% x: a pixel from hyperspectral images
% y: a pixel from hyperspectral images

% d: the distance between x and y
%--------------------------------------------------------------------------
function d = dis(x,y)

d = sin(acos(dot(x,y)/(norm(x)*norm(y))));

end