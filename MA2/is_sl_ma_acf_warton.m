function [theta, w] = is_sl_ma_acf_warton(y,M,n,mu_imp,cov_imp,K,gamma)
% Importance Sampling shrinkage BSL for MA(2) model
% Determines ABC tolerance using a target on effective sample size (ESS)
% y - observed time series
% M is the total number of MCMC iterations
% n is number of simulated data sets for estimating ABC likelihood
% mu_imp - mean of independent Gaussian proposal
% cov_imp - covariance of independent Gaussian proposal
% K - number of autocovariances to include in the summary statistics
% gamma - shrinkage parameter for estimating covariance of summaries

theta = zeros(M,2); %Initialising vector of estimates for motility and diffusivity
w = zeros(M,1);
ssy = summStat_ma2_acf(y, K);


T = length(y); %Total number of time periods including time zero


for i = 1:M
    %i
    while (1)
        theta_prop = mvnrnd(mu_imp,cov_imp); %Proposing a new pair of parameters
        if (sum(theta_prop) > -1 && diff(theta_prop) > -1 && theta_prop(2) > -1 && theta_prop(2) < 1)
            theta(i,:) = theta_prop;
            break;
        end
    end
    
    %Simulating n data sets, finding summary statistics and then getting
    %the mean and covariance of these summary statistics
    y = simulate_ma2(theta(i,:),T,n);
    ssx = summStat_ma2_acf(y, K);
    
    the_mean = mean(ssx);
    the_cov = warton(ssx,gamma);
    L = chol(the_cov);
    logdetA = 2*sum(log(diag(L)));
    
    %Calculating synthetic likelihood
    loglike_ind_prop = -0.5*logdetA - 0.5*(ssy-the_mean)*inv(the_cov)*(ssy-the_mean)';
    
    w(i) = loglike_ind_prop - log(mvnpdf(theta(i,:), mu_imp, cov_imp));

end

w = w - max(w);
w = exp(w);
w = w/sum(w);

end

