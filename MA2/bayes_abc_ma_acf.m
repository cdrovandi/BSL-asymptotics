function [theta, loglike] = bayes_abc_ma_acf(y,M,n,cov_rw,start,K,eps,cov_abc)
% MCMC ABC for MA(2) model
% y - observed time series
% M is the total number of MCMC iterations
% n is number of simulated data sets for estimating ABC likelihood
% cov_rw is the covariance matrix used in the random walk MCMC
% start - starting value of MCMC chain
% K - number of autocovariances to include in the summary statistics
% eps - ABC tolerance
% cov_abc - weighting matrix for mahalanobis ABC distance

loglike=zeros(M,1);

theta_curr = start;
theta = zeros(M,2); %Initialising vector of estimates for motility and diffusivity
ssy = summStat_ma2_acf(y, K);


T = length(y); %Total number of time periods including time zero

%Simulating n sets of data and computing their summary statistics
x = simulate_ma2(theta_curr,T,n);
ssx = summStat_ma2_acf(x, K);

log_kernels = log(mvnpdf(ssy,ssx,eps*cov_abc));
loglike_ind_curr = -log(n) + logsumexp(log_kernels);

for i = 1:M
    %i
    theta_prop = mvnrnd(theta_curr,cov_rw); %Proposing a new pair of parameters
    if (sum(theta_prop) < -1 || diff(theta_prop) < -1 || theta_prop(2) < -1 || theta_prop(2) > 1) % to ensure invertibility
        theta(i,:) = theta_curr;
        continue;
    end
    
    %Simulating n sets of data and computing their summary statistics
    x = simulate_ma2(theta_prop,T,n);
    ssx = summStat_ma2_acf(x, K);
    
    % using Gaussian kernel for ABC
    log_kernels = log(mvnpdf(ssy,ssx,eps*cov_abc));
    loglike_ind_prop = -log(n) + logsumexp(log_kernels);

    % accept/reject
    if (exp(loglike_ind_prop - loglike_ind_curr) > rand)
        %fprintf('*** accept ***\n');
        theta_curr = theta_prop;
        loglike_ind_curr = loglike_ind_prop;
    end
    theta(i,:) = theta_curr;
    loglike(i,1)=loglike_ind_curr;   
end
end

