function run_sl_simulated_shrinkage(r)

% run shrinkage BSL on 100 replicate simulated datasets for toad example

parpool(16);

load('data_simulated_replicates.mat')

Y = Xr{r};
NaN_Pos = isnan(Y);

n = 50;
M = 50000;
lag = [1, 2, 4, 8];
ssy = summStat_quantiles3(Y, lag)';
cov_rw = [0.0843,0.0072,-0.0003;0.0072,0.0028,-0.0003;-0.0003,-0.0003,0.0077];
start = [1.7, 35, 0.6];
simArgs = struct('ntoads',ntoads,'ndays',ndays,'model',1,'d0',NaN);
sumArgs = struct('lag',lag);
gamma = 0.1;
sumstat_fun = 'summStat_quantiles3';

tic;
[theta,loglike] = bayes_sl_toads_warton(ssy,n,M,cov_rw,start,simArgs,sumArgs,NaN_Pos,sumstat_fun,gamma);
time = toc;

save(['results_bsl_n50_shrink' num2str(r) '.mat'],'theta','loglike','time');

delete(gcp);




