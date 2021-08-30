function run_exact_T1000(r)
% exact 100 observations

%c = parcluster('local')
%c.JobStorageLocation = tempdir;
%p = parpool(c,16);


load('data_10000_100reps.mat');
cov_rw = [0.00137612338695324,0.000477923334713194;0.000477923334713194,0.00112616980221712];

M=10000;
y = y(1:1000,r);


[theta, loglike] = bayes_ma_sparse(y',M,cov_rw,[0.6 0.2]);


save(['results_exact_T1000_' num2str(r) '.mat'],'theta','loglike');


%delete(p);

