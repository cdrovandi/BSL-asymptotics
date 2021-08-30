function run_is_sl_shrink_T10000_n2000(r)


%c = parcluster('local')
%c.JobStorageLocation = tempdir;
%p = parpool(c,16);


load('data_10000_100reps.mat');
load('results_sl_shrink_T10000_n500.mat');

M=10000;
n=2000;
K=20;
gamma = 0;

tic;
[theta_is, w_is] = is_sl_ma_acf_warton(y(:,r),M,n,mean(theta(:,:,r)),2*cov(theta(:,:,r)),K,gamma);
finaltime = toc;


save(['results_is_sl_shrink_T10000_n2000_' num2str(r) '.mat'],'theta_is','w_is','finaltime');


%delete(p);

