% implements the function D'*v = [Dx' Dy']*v
%   -> finite differences!
%
%   input:  v are 2x 2D images (not vectorized)
%   output: I is a 2D image (not vectorized)

function I = opDtx(v, bWithCircularBoundaryConds)
    
    if nargin<2
        bWithCircularBoundaryConds = false;
    end

    if ~bWithCircularBoundaryConds
        % Dx and Dy
        I = [zeros([size(v,1) 1]) v(:,1:end-1,1)] - [v(:,1:end-1,1) zeros([size(v,1) 1])] + ... 
            [zeros([1 size(v,2)]); v(1:end-1,:,2)] - [v(1:end-1,:,2); zeros([1 size(v,2)])];
        
    else
        I = (circshift(v(:,:,1),[0 1])-v(:,:,1)) + (circshift(v(:,:,2),[1 0])-v(:,:,2));
    end
end