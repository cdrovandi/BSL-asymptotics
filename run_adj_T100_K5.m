

load('data_10000_100reps.mat');
y = y(1:100,:);
load('results_is_sl_shrink_T100_K5.mat')
n = 50;

gamma = 0;
L = 5;

parfor i = 1:100
    
    theta = theta_is(:,:,i);
    w = w_is(:,i);
    r = randsample(1:10000, 1000, 'true', w);
    theta_resample  = theta(r,:);
    
    theta_adj(:,:,i) = adjust_gp_acf_function(y(:,i),theta_resample,gamma,L,n);
    
end

theta_is = theta_adj;
w = ones(1000,1)/1000;

save('results_adj_T100_K5.mat','theta_is','w');

