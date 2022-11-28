% Sets up a test problem and interpolates it with the intensity spline
% interpolated registered spline trajectories.


%% define the 'circling ellipse' test problem

% size of the domain
m = 90; n = 120;
% center points
center_m = floor(m/2); center_n = floor(n/2);
% minor and major axis of the ellipsoidal motion
R_y = m/4; R_x = n/4;
% minor and major axis of te ellipses
r_y = m/5; r_x = n/5;
% set up discretization of the circular motion
numberOfSlices = 9;
theta = linspace(3*pi/2,2*pi+3*pi/2,numberOfSlices);
% set up the center points of the circles
centerPoints_y = center_m + R_y * sin(-theta);
centerPoints_x = center_n + R_x * cos(-theta);

% define the coloring of the circles
% pI = @(x) 1;
% pI = @(x) x;
% pI = @(x) x.^3;
pI = @(x) sin(x*pi*2);

% create the test problem
slices = zeros(m,n,numberOfSlices);
for k = 1:numberOfSlices
    slices(:,:,k) = drawCircle(m,n,centerPoints_x(k),centerPoints_y(k),r_x,r_y) * pI(k/numberOfSlices);
end
% set the slice heights
z = 1:numberOfSlices; % algorithm works with non-equidistant slices distances

% set up the a cell array of slices in different contrasts
modSize = 1;
Imultimod = {};
for k = 1:numberOfSlices
    for mm = 1:modSize
        Imultimod{k}(:,:,mm) = drawCircle(m,n,centerPoints_x(k),centerPoints_y(k),r_x,r_y); % here you could add different contrasts
    end
end

%% slice interpolation
R = 3;
lambda = 10000; % the higher the regularization factor lambda it will search for more global motions but ignores local structures
tau = 100;  % step size of the implicit gradient descent
TOL = 0.04;
maxIter = 1000;
borderSize = 0.1;

[slices_interpolated,z_interpolated,vx,vy] = sliceInterp_spline_intensitySpline(slices,z,R,lambda,tau,TOL,maxIter,borderSize);

% Update: For the paper we interpolated MR slices that are captured at
% different inversion recovery times, which results to different contrasts
% of the same amnatomical structures.
% for the circeling circle problem we actually used the single mode
% Imultimod with the proper circles. Otherwise, registration problems
% occure when the circle has the same color as the background.
% [slices_interpolated,z_interpolated,vx,vy] = sliceInterp_spline_intensitySpline_multimod(slices,Imultimod,z,R,lambda,tau,TOL,maxIter,borderSize);

%% analytic solution
theta_analytic = linspace(3*pi/2,2*pi+3*pi/2,(numberOfSlices-1)*R+1);
centerPoints_y_analytic = center_m + R_y * sin(-theta_analytic);
centerPoints_x_analytic = center_n + R_x * cos(-theta_analytic);
slices_interpolated_analytic = zeros(m,n,(numberOfSlices-1)*R+1);
for k = 1:(numberOfSlices-1)*R+1
    slices_interpolated_analytic(:,:,k) = drawCircle(m,n,centerPoints_x_analytic(k),centerPoints_y_analytic(k),r_x,r_y);
end


%% visualize the results

figure(1);clf
interpolatedSliceNumber = 1; subPlotNumber = 1;
for originalSliceNumber = 1:numberOfSlices-1
    for intermediateNumber = 1:R
        subplot(numberOfSlices-1,R+1,subPlotNumber)
        imshow(slices_interpolated(:,:,interpolatedSliceNumber),[-1,1])
%         imshowpair(slices_interpolated_analytic(:,:,interpolatedSliceNumber),slices_interpolated(:,:,interpolatedSliceNumber),'Scaling','joint')
        if intermediateNumber == 1
            title(['original slice ' num2str(originalSliceNumber)])
        end
        interpolatedSliceNumber = interpolatedSliceNumber+1; subPlotNumber = subPlotNumber+1;
    end
    subplot(numberOfSlices-1,R+1,subPlotNumber); subPlotNumber = subPlotNumber+1;
    imshow(slices_interpolated(:,:,interpolatedSliceNumber),[-1,1])
%     imshowpair(slices_interpolated_analytic(:,:,interpolatedSliceNumber), slices_interpolated(:,:,interpolatedSliceNumber),'Scaling','joint')
    title(['original slice ' num2str(originalSliceNumber+1)])
end


%% functions
function I = drawCircle(m,n,a,b,rx,ry)
    I = zeros(m,n);
    for i = 1:m
        for j = 1:n
            if ((i-b)/ry)^2 + ((j-a)/rx)^2 <= 1
                I(i,j) = 1;
            end
        end
    end
end
