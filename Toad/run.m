
% run BSL with and without shrinkage on toad model (simulated and real data)


%% run standard BSL on simulated data

load('data_toads_model1.mat')

NaN_Pos = isnan(Y);

n = 500;
M = 50000;
lag = [1, 2, 4, 8];
ssy = summStat_quantiles3(Y, lag)';
cov_rw = [0.0843,0.0072,-0.0003;0.0072,0.0028,-0.0003;-0.0003,-0.0003,0.0077];
start = [1.7, 35, 0.6];
simArgs = struct('ntoads',ntoads,'ndays',ndays,'model',1,'d0',NaN);
sumArgs = struct('lag',lag);
sumstat_fun = 'summStat_quantiles3';

tic;
[theta,loglike] = bayes_sl_toads(ssy,n,M,cov_rw,start,simArgs,sumArgs,NaN_Pos,sumstat_fun);
time = toc;

save('results_simulated_bsl_n500.mat','theta','loglike','time');



%% run standard BSL on real data

load('radio_converted.mat');

NaN_Pos = isnan(Y);

n = 500;
M = 50000;
lag = [1, 2, 4, 8];
ssy = summStat_quantiles3(Y, lag)';
cov_rw = [0.0843,0.0072,-0.0003;0.0072,0.0028,-0.0003;-0.0003,-0.0003,0.0077];
start = [1.7, 35, 0.6];
simArgs = struct('ntoads',ntoads,'ndays',ndays,'model',1,'d0',NaN);
sumArgs = struct('lag',lag);
sumstat_fun = 'summStat_quantiles3';

tic;
[theta,loglike] = bayes_sl_toads(ssy,n,M,cov_rw,start,simArgs,sumArgs,NaN_Pos,sumstat_fun);
time = toc;

save('results_real_bsl_n500.mat','theta','loglike','time');



%% run shrinkage BSL on simulated data

load('data_toads_model1.mat')

NaN_Pos = isnan(Y);

n = 50;
M = 50000;
lag = [1, 2, 4, 8];
ssy = summStat_quantiles3(Y, lag)';
cov_rw = [0.0843,0.0072,-0.0003;0.0072,0.0028,-0.0003;-0.0003,-0.0003,0.0077];
start = [1.7, 35, 0.6];
simArgs = struct('ntoads',ntoads,'ndays',ndays,'model',1,'d0',NaN);
sumArgs = struct('lag',lag);
sumstat_fun = 'summStat_quantiles3';

gamma = 0.1;
tic;
[theta,loglike] = bayes_sl_toads_warton(ssy,n,M,cov_rw,start,simArgs,sumArgs,NaN_Pos,sumstat_fun,gamma);
time = toc;

save('results_simulated_bsl_n50_shrink.mat','theta','loglike','time');


%% run shrinkage BSL on real data

load('radio_converted.mat');

NaN_Pos = isnan(Y);

n = 50;
M = 50000;
lag = [1, 2, 4, 8];
ssy = summStat_quantiles3(Y, lag)';
cov_rw = [0.0843,0.0072,-0.0003;0.0072,0.0028,-0.0003;-0.0003,-0.0003,0.0077];
start = [1.7, 35, 0.6];
simArgs = struct('ntoads',ntoads,'ndays',ndays,'model',1,'d0',NaN);
sumArgs = struct('lag',lag);
sumstat_fun = 'summStat_quantiles3';

gamma = 0.1;
tic;
[theta,loglike] = bayes_sl_toads_warton(ssy,n,M,cov_rw,start,simArgs,sumArgs,NaN_Pos,sumstat_fun,gamma);
time = toc;

save('results_real_bsl_n50_shrink.mat','theta','loglike','time');




