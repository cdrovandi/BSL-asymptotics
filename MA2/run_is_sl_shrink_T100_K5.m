function run_is_sl_shrink_T100_K5(r)


%c = parcluster('local')
%c.JobStorageLocation = tempdir;
%p = parpool(c,16);


load('data_10000_100reps.mat');
load('results_sl_shrink_T100_K5.mat');

M=10000;
n=50;
K=5;
gamma = 0;
y = y(1:100,r);

tic;
[theta_is, w_is] = is_sl_ma_acf_warton(y,M,n,mean(theta(:,:,r)),2*cov(theta(:,:,r)),K,gamma);
finaltime = toc;


save(['results_is_sl_shrink_T100_K5_' num2str(r) '.mat'],'theta_is','w_is','finaltime');


%delete(p);

