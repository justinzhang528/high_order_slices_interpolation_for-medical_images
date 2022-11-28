V = niftiread('003_flair.nii');
Mask = niftiread('003_seg.nii');
numberOfSlices = 48;

whos V

slices = zeros(128,128,2);
for k = 1:numberOfSlices
    slices(:,:,k) = rescale(V(:,:,k));
end

for k = 1:numberOfSlices
    imshow(slices(:,:,k),[])
end

% set the slice heights
z = 1:numberOfSlices; % algorithm works with non-equidistant slices distances

%% slice interpolation
R = 3;
lambda = 1000; % the higher the regularization factor lambda it will search for more global motions but ignores local structures
tau = 100;  % step size of the implicit gradient descent
TOL = 0.04;
maxIter = 1000;
borderSize = 0.1;

[slices_interpolated,z_interpolated,vx,vy] = sliceInterp_spline_intensitySpline(slices,z,R,lambda,tau,TOL,maxIter,borderSize);

s = size(slices_interpolated)
s = s(3)

for k = 1:s
    imshow(slices_interpolated(:,:,k),[])
end
niftiwrite(slices_interpolated,'003_flair_interpolated.nii');