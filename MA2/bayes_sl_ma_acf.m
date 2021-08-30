function [theta, loglike] = bayes_sl_ma_acf(y,M,n,cov_rw,start,K)
% MCMC BSL for MA(2) model
% y - observed time series
% M is the total number of MCMC iterations
% n is number of simulated data sets for estimating ABC likelihood
% cov_rw is the covariance matrix used in the random walk MCMC
% start - starting value of MCMC chain
% K - number of autocovariances to include in the summary statistics

loglike=zeros(M,1);

theta_curr = start;
theta = zeros(M,2); 
ssy = summStat_ma2_acf(y, K);


T = length(y); % length of time series

%Simulating n sets of data and computing their summary statistics
y = simulate_ma2(theta_curr,T,n);
ssx = summStat_ma2_acf(y, K);

%Calculating the mean and covariance of the summary statistics
the_mean = mean(ssx);
the_cov = cov(ssx);
L = chol(the_cov);
logdetA = 2*sum(log(diag(L)));

%Calculating log synthetic likelihood
loglike_ind_curr = -0.5*logdetA - 0.5*(ssy-the_mean)*inv(the_cov)*(ssy-the_mean)';

for i = 1:M
    %i
    theta_prop = mvnrnd(theta_curr,cov_rw); %Proposing a new pair of parameters
    if (sum(theta_prop) < -1 || diff(theta_prop) < -1 || theta_prop(2) < -1 || theta_prop(2) > 1) % to ensure invertibility
        theta(i,:) = theta_curr;
        continue;
    end
    
    %Simulating n data sets, finding summary statistics and then getting
    %the mean and covariance of these summary statistics
    y = simulate_ma2(theta_prop,T,n);
    ssx = summStat_ma2_acf(y, K);
    
    the_mean = mean(ssx);
    the_cov = cov(ssx);
    L = chol(the_cov);
    logdetA = 2*sum(log(diag(L)));
    
    %Calculating log synthetic likelihood
    loglike_ind_prop = -0.5*logdetA - 0.5*(ssy-the_mean)*inv(the_cov)*(ssy-the_mean)';

    %accept/reject
    if (exp(loglike_ind_prop - loglike_ind_curr) > rand)
        %fprintf('*** accept ***\n');
        theta_curr = theta_prop;
        loglike_ind_curr = loglike_ind_prop;
    end
    theta(i,:) = theta_curr;
    loglike(i,1)=loglike_ind_curr;   
end
end

