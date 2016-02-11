function x_out = Spatial_Tikhonov(x,y,param)
% same  as tikhonov_it.m, 28.11.2011 isik
% This function computes the tikhonov regularization
%
% F(x) = min ||y - x ||^2 + lambda * ||Delta{x}||^2
%
% Delta is the laplacian operator delta[n] = [1 -2 1]; so symmetric, in
% matrix form Delta^T = Delta.

% the iterative algorithm, prox operator probably leads to the same result.

% x(k+1) = x(k) - mu*grad(F(x)), the descent direction is -grad(F(x)) and
% mu is the step size.
% grad(F(x)) = -y+x+lambda*Delta^t*Delta{x}, then

% x(k+1) = x(k) + mu*y - mu*(I + lambda*Delta^T*Delta){x(k)}.

%%%%%%---INPUTS---%%%%%%

% x = 2D matrix, output from fgp
% y = observation
% dim=param.dimTik;% dim can be 2 or 3D.
% iter=param.NitSpat; %internal iterations for tikhonov (not too much)
% lambda=param.LambdaSpat;%regularization coefficient, 1.
% mu=param.stepsize; %step size, small 0.01
% IND=param.VoxelIdx; %[xi,yi,zi] = nonzero indices to convert 2D x to 4D x_vol
% SIZE=param.Dimension; %the size of the 4D data [v1,v2,v3,t]

% IND = [xi,yi,zi]
% SIZE = [v1,v2,v3,t]

%%%%%---OUTPUTS---%%%%%%
% x_out
%%%%%%%%%%%%%%%%%%%%%%%%%

% Written by Isik, 30.08.2011

% disp('inside tikhonov')


dim=param.dimTik;% dim can be 2 or 3D.
iter=param.NitSpat; %internal iterations for tikhonov (not too much)
lambda=param.LambdaSpat;%regularization coefficient, 1.
mu=param.stepsize; %step size, small 0.01
IND=param.VoxelIdx; %[xi,yi,zi] = nonzero indices to convert 2D x to 4D x_vol
SIZE=param.Dimension; %the size of the 4D data [v1,v2,v3,t]

if (length(SIZE) ~= 4)
    error('SIZE should have 4 dimensions.')
end

x_out = zeros(size(x)); % out

%%%% convert to 4D
x_vol = zeros(SIZE);
y_vol = zeros(SIZE);
% temp = zeros(SIZE);

for i = 1:length(IND(:,1));
    x_vol(IND(i,1),IND(i,2),IND(i,3),:) = x(:,i);
    y_vol(IND(i,1),IND(i,2),IND(i,3),:) = y(:,i);
end

k=0;

if (dim==2)
    % take each plane one by one.
    h=[0 1 0; 1 -4 1; 0 1 0];
    
    H = fft2(h,SIZE(1),SIZE(2));
    while(k<iter)
        for i=1:SIZE(4)
            for j=1:SIZE(3)
                x_vol(:,:,j,i) = (1-mu)*x_vol(:,:,j,i) + mu*y_vol(:,:,j,i)  - mu*lambda*ifft2((conj(H).*H).*fft2(x_vol(:,:,j,i),SIZE(1),SIZE(2)));
            end
        end
    k=k+1;    
    end
end

if (dim==3)
    h = zeros(3,3,3);
    h(:,:,1) = [0 0 0; 0 1 0; 0 0 0];
    h(:,:,2) = [0 1 0; 1 -6 1; 0 1 0];
    h(:,:,3) = h(:,:,1);
    
    H = fftn(h,SIZE(1:3));
    while(k<iter)
        for i=1:SIZE(4)
            x_vol(:,:,:,i) = (1-mu)*x_vol(:,:,:,i) + mu*y_vol(:,:,:,i)  - mu*lambda*ifftn(H.*conj(H).*fftn(x_vol(:,:,:,i),SIZE(1:3)));
        end
    k=k+1; 
    end
end


for i=1:length(IND(:,1));
    x_out(:,i) = x_vol(IND(i,1),IND(i,2),IND(i,3),:);
end

end