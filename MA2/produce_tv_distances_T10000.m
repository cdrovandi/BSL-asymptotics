
% estimate TV distances between approximation and true posterior for 100
% simulated datasets of size 10000.
% change the approximation via the .mat file as required

load('results_adj_T10000_K10_n500.mat');

theta_approx = theta_is;
load('results_exact_T10000.mat')
theta_exact = theta;

len = 0.005;
wid = 0.005;
theta1_grid = 0.4:len:0.8;
theta2_grid = 0:wid:0.4;
[X1,X2] = meshgrid(theta1_grid, theta2_grid);
X1 = X1(:);
X2 = X2(:);

tv_dist = zeros(100,1);
parfor i = 1:100
    %i
    [f_true,~] = ksdensity(theta_exact(:,:,i),[X1 X2]);
    
    if (length(unique(w_is(:,i))) == 1) % resampling already done
        thetar = theta_approx(:,:,i);
    else
        r = randsample(1:length(theta_approx(:,1,i)), 1000, 'true', w_is(:,i));
        thetar = theta_approx(r,:,i);
    end
    
    [f_approx,~] = ksdensity(thetar,[X1 X2]);
    tv_dist(i) = 0.5*sum(abs(f_true - f_approx))*len*wid;
    
end
median(tv_dist)

