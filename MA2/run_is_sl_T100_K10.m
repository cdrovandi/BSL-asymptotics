function run_is_sl_T100_K10(r)
% 100 observations and 10 autocovariances

%c = parcluster('local')
%c.JobStorageLocation = tempdir;
%p = parpool(c,16);


load('data_10000_100reps.mat');
load('results_sl_T100_K10.mat');

M=10000;
n=500;
K=10;
y = y(1:100,r);

tic;
[theta_is, w_is] = is_sl_ma_acf(y,M,n,mean(theta(:,:,r)),2*cov(theta(:,:,r)),K);
finaltime = toc;


save(['results_is_sl_T100_K10_' num2str(r) '.mat'],'theta_is','w_is','finaltime');


%delete(p);

