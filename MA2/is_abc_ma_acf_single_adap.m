function [eps_target, mu_imp_target, cov_imp_target] = is_abc_ma_acf_single_adap(y,M,mu_imp,cov_imp,K,cov_abc,ess_target)
% Importance Sampling ABC for MA(2) model
% Determines ABC tolerance using a target on effective sample size (ESS)
% y - observed time series
% M is the total number of MCMC iterations
% n is number of simulated data sets for estimating ABC likelihood
% mu_imp - mean of independent Gaussian proposal
% cov_imp - covariance of independent Gaussian proposal
% K - number of autocovariances to include in the summary statistics
% cov_abc - weighting matrix for mahalanobis ABC distance
% ess_target - target effective sample size to work out ABC tolerance

theta = zeros(M,2); %Initialising vector of estimates for motility and diffusivity
w = zeros(M,1);
ssy = summStat_ma2_acf(y, K);
ssx_sim = zeros(M,K);


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
    
    %Simulating n data sets, compuing summary statistics 
    y = simulate_ma2(theta(i,:),T,1);
    ssx = summStat_ma2_acf(y, K);
    ssx_sim(i,:) = ssx;

end

eps_target = fzero(@(x)compute_ess(x, ssy, ssx_sim, cov_abc, theta, mu_imp, cov_imp) - ess_target,[0.1 5]);

w = log(mvnpdf(ssy,ssx_sim,eps_target*cov_abc)) - log(mvnpdf(theta, mu_imp, cov_imp));

w = w - max(w);
w = exp(w);
w = w/sum(w);

r = randsample(1:M, ess_target, 'true', w);
theta_target = theta(r,:);
mu_imp_target = mean(theta_target);
cov_imp_target = cov(theta_target);

end



function f = compute_ess(eps, ssy, ssx_sim, cov_abc, theta, mu_imp, cov_imp)

w = log(mvnpdf(ssy,ssx_sim,eps*cov_abc)) - log(mvnpdf(theta, mu_imp, cov_imp));

w = w - max(w);
w = exp(w);
w = w/sum(w);

f = 1/sum(w.^2);

end

