% implements the function D*vol = [Dx; Dy]*image
%   -> finite differences!
%
%   input:  I is a 2D volume (not vectorized)
%   output: result are 2x 2D images - gradient in x and gradient in y (not vectorized)
%   Dx = result(:,:,1)
%   Dy = result(:,:,2)

function result = opDx(I, bWithCircularBoundaryConds)
    
    if nargin<2
        bWithCircularBoundaryConds = false;
    end

    % this is better
    if ~bWithCircularBoundaryConds
        % Dx and Dy
        result = cat(3, I(:,[2:end end],:)-I, I([2:end end],:,:)-I);
            
    % with convolution
    else
        result = cat(3, circshift(I,[0 -1])-I, circshift(I,[-1 0])-I);
    end
end