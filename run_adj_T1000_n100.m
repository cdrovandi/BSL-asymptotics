

load('data_1000_100reps.mat');
load('results_is_sl_shrink_T1000_n100.mat')
n = 100;

gamma = 0;
L = 20;

parfor i = 1:100
    
    theta = theta_is(:,:,i);
    w = w_is(:,i);
    r = randsample(1:10000, 1000, 'true', w);
    theta_resample  = theta(r,:);
    
    theta_adj(:,:,i) = adjust_gp_acf_function(y(:,i),theta_resample,gamma,L,n);
    
end

theta_is = theta_adj;
w = ones(1000,1)/1000;

save('results_adj_T1000_n100.mat','theta_is','w');

