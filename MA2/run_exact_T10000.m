function run_exact_T10000(r)
% exact 100 observations

%c = parcluster('local')
%c.JobStorageLocation = tempdir;
%p = parpool(c,16);


load('data_10000_100reps.mat');
cov_rw = [0.000157244082770195	6.22720196676529e-05;
6.22720196676529e-05	9.27296131650419e-05];

M=10000;
y = y(:,r);


[theta, loglike] = bayes_ma_sparse(y',M,cov_rw,[0.6 0.2]);


save(['results_exact_T10000_' num2str(r) '.mat'],'theta','loglike');


%delete(p);

