function [theta, loglike] = bayes_ma_sparse(ssy,M,cov_rw,start)
% MCMC sampling for MA(2) model (using sparse matrix calculations to speed up likelihood calculation)
% ssy - observed time series
% M is the total number of MCMC iterations
% n is number of simulated data sets for estimating ABC likelihood
% cov_rw is the covariance matrix used in the random walk MCMC
% start - starting value of MCMC chain

loglike=zeros(M,1);

theta_curr = start;
theta = zeros(M,2); 

loglike_ind_curr = loglike_ma2_sparse(ssy,theta_curr);

for i = 1:M
    %i
    theta_prop = mvnrnd(theta_curr,cov_rw); %Proposing a new pair of parameters
    if (sum(theta_prop) < -1 || diff(theta_prop) < -1) %Rejecting negative and >1 proportions and going to next iteration
        theta(i,:) = theta_curr;
        loglike(i)=loglike_ind_curr;
        continue;
    end
    
    loglike_ind_prop = loglike_ma2_sparse(ssy,theta_prop);

    % accept/reject
    if (exp(loglike_ind_prop - loglike_ind_curr) > rand)
        %fprintf('*** accept ***\n');
        theta_curr = theta_prop;
        loglike_ind_curr = loglike_ind_prop;
    end
    theta(i,:) = theta_curr;
    loglike(i)=loglike_ind_curr;   
end
end

