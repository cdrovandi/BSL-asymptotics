function run_is_abc_single_T10000(r)


%c = parcluster('local')
%c.JobStorageLocation = tempdir;
%p = parpool(c,16);


load('data_10000_100reps.mat');
load('results_is_abc_single_T10000_findeps.mat')
load('abc_dist_cov_T10000.mat')

M=1e7;
K=20;


tic;
[theta_is, w_is, ssx_is] = is_abc_ma_acf_single(y(:,r),M,mu_eps{r},2*cov_eps{r},K,eps{r},cov_abc);
finaltime = toc;


save(['results_is_abc_single_T10000_' num2str(r) '.mat'],'theta_is','w_is','ssx_is','finaltime');


%delete(p);

