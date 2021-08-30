function [theta, loglike] = bayes_ma(ssy,M,cov_rw,start)
% MCMC sampling for MA(2) model
% ssy - observed time series
% M is the total number of MCMC iterations
% n is number of simulated data sets for estimating ABC likelihood
% cov_rw is the covariance matrix used in the random walk MCMC
% start - starting value of MCMC chain

loglike=zeros(M,1);

theta_curr = start;
theta = zeros(M,2); 

ns = length(ssy); % length of time series

the_cov = zeros(ns,ns);
for k = 1:ns
    the_cov(k,k) = 1 + theta_curr(1)^2 + theta_curr(2)^2;
    if (k < ns)
        the_cov(k,k+1) = theta_curr(1) + theta_curr(1)*theta_curr(2);
        the_cov(k+1,k) = the_cov(k,k+1);
    end
    if (k < ns-1)
        the_cov(k,k+2) = theta_curr(2);
        the_cov(k+2,k) = the_cov(k,k+2);
    end
end

loglike_ind_curr = -0.5*log(det(the_cov)) - 0.5*ssy*inv(the_cov)*ssy';

for i = 1:M
    i
    theta_prop = mvnrnd(theta_curr,cov_rw); %Proposing a new pair of parameters
    if (sum(theta_prop) < -1 || diff(theta_prop) < -1) %Rejecting negative and >1 proportions and going to next iteration
        theta(i,:) = theta_curr;
        loglike(i)=loglike_ind_curr;
        continue;
    end
    
    the_cov = zeros(ns,ns);
    for k = 1:ns
        the_cov(k,k) = 1 + theta_prop(1)^2 + theta_prop(2)^2;
        if (k < ns)
            the_cov(k,k+1) = theta_prop(1) + theta_prop(1)*theta_prop(2);
            the_cov(k+1,k) = the_cov(k,k+1);
        end
        if (k < ns-1)
            the_cov(k,k+2) = theta_prop(2);
            the_cov(k+2,k) = the_cov(k,k+2);
        end
    end
    
    loglike_ind_prop = -0.5*log(det(the_cov)) - 0.5*ssy*inv(the_cov)*ssy';

    % accept/reject
    if (exp(loglike_ind_prop - loglike_ind_curr) > rand)
        fprintf('*** accept ***\n');
        theta_curr = theta_prop;
        loglike_ind_curr = loglike_ind_prop;
    end
    theta(i,:) = theta_curr;
    loglike(i)=loglike_ind_curr;   
end
end

