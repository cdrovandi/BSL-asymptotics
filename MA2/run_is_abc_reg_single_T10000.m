

load('abc_dist_cov_T10000.mat');
load('results_is_abc_single_T10000.mat');
load('data_10000_100reps.mat')

K = 20;


theta_is = abc_regadj_ma_acf(y, theta_is, ones(1000,100), ssx, cov_abc, K);
    


