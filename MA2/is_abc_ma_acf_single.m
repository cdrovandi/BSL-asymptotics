function [theta, w, ssx_sim] = is_abc_ma_acf_single(y,M,mu_imp,cov_imp,K,eps,cov_abc)
% Importance Sampling ABC for MA(2) model
% y - observed time series
% M is the total number of MCMC iterations
% n is number of simulated data sets for estimating ABC likelihood
% mu_imp - mean of independent Gaussian proposal
% cov_imp - covariance of independent Gaussian proposal
% K - number of autocovariances to include in the summary statistics
% eps - ABC tolerance
% cov_abc - weighting matrix for mahalanobis ABC distance

theta = zeros(M,2); %Initialising vector of estimates for motility and diffusivity
w = zeros(M,1);
ssy = summStat_ma2_acf(y, K);
ssx_sim = zeros(M,K);


T = length(y); % length of time series


for i = 1:M
    %i
    while (1)
        theta_prop = mvnrnd(mu_imp,cov_imp); %Proposing a new pair of parameters
        if (sum(theta_prop) > -1 && diff(theta_prop) > -1 && theta_prop(2) > -1 && theta_prop(2) < 1)
            theta(i,:) = theta_prop;
            break;
        end
    end
    
    %Simulating n data sets, computing summary statistics 
    y = simulate_ma2(theta(i,:),T,1);
    ssx = summStat_ma2_acf(y, K);
    ssx_sim(i,:) = ssx;
    
    w(i) = log(mvnpdf(ssy,ssx,eps*cov_abc)) - log(mvnpdf(theta(i,:), mu_imp, cov_imp));

end

w = w - max(w);
w = exp(w);
w = w/sum(w);

end

