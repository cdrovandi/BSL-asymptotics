function run_is_abc_single_T10000_findeps(r)


%c = parcluster('local')
%c.JobStorageLocation = tempdir;
%p = parpool(c,16);


load('data_10000_100reps.mat');
load('results_abc_T10000_n500_eps1.mat')
load('abc_dist_cov_T10000.mat')

M=1e7;
K=20;
%eps=1;
ess_target = 1000;


[eps, mu_eps, cov_eps] = is_abc_ma_acf_single_adap(y(:,r),M,mean(theta(:,:,r)),2*cov(theta(:,:,r)),K,cov_abc,ess_target);


save(['results_is_abc_single_T10000_findeps_' num2str(r) '.mat'],'eps','mu_eps','cov_eps');


%delete(p);

