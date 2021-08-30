function run_sl_shrink_T100_K20(r)


%c = parcluster('local')
%c.JobStorageLocation = tempdir;
%p = parpool(c,16);


load('data_10000_100reps.mat');
cov_rw = [0.0224366757783134,-0.0155043087329265;-0.0155043087329265,0.0283409994965521];

M=10000;
n=50;
K=20;
gamma = 0;
y = y(1:100,r);


[theta, loglike] = bayes_sl_ma_acf_warton(y,M,n,cov_rw,[0.6 0.2],K,gamma);


save(['results_sl_shrink_T100_K20_' num2str(r) '.mat'],'theta','loglike');


%delete(p);

