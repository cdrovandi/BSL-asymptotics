function run_is_sl_T10000_n500(r)


%c = parcluster('local')
%c.JobStorageLocation = tempdir;
%p = parpool(c,16);


load('data_10000_100reps.mat');
load('results_sl_T10000_n500.mat');

M=10000;
n=500;
K=20;

tic;
[theta_is, w_is] = is_sl_ma_acf(y(:,r),M,n,mean(theta(:,:,r)),2*cov(theta(:,:,r)),K);
finaltime = toc;


save(['results_is_sl_T10000_n500_' num2str(r) '.mat'],'theta_is','w_is','finaltime');


%delete(p);

