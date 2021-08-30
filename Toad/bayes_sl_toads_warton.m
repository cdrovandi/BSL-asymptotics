function [theta,loglike] = bayes_sl_toads_warton(ssy,n,M,cov_rw,start,simArgs,sumArgs,NaN_Pos,sumstat_fun,gamma)
% runs BSL shrhinkage on toad model
% ssy - observed summary statistic
% n - number of simulated datasets to estimate synthetic likelihood
% M - total number of MCMC iterations
% cov_rw - covariance matrix for MCMC random walk proposal
% start - initial parameter for MCMC chain
% simArgs - additional arguments required in simulation function
% sumArgs - additional arguments required in summary statistic function
% NaN_Pos - indicator matrix specifiying which data points are missing or not missing (some observations in real data are missing)
% summstat_fun - summary statistic function
% gamma - shrinkage parameter of synthetic likelihood covariance

ntoads = simArgs.ntoads;
ndays = simArgs.ndays;

lag = sumArgs.lag;


theta_curr = start; %Initial guesses for parameters
ns = length(ssy);
theta = zeros(M,3);
loglike = zeros(M,1);

% Simulating n sets of data and taking their summary statistics
ssx = zeros(n,ns);
parfor k = 1:n
    X = simulate_toads2(theta_curr,ntoads,ndays);
    X(NaN_Pos) = NaN;
    ssx(k,:) = feval(sumstat_fun,X,lag);
end

% Calculating the mean and covariance of the summary statistics
the_mean = mean(ssx)';
the_cov = warton(ssx, gamma);

loglike_ind_curr = -0.5*log(det(the_cov)) - 0.5*(ssy-the_mean)'/the_cov*(ssy-the_mean);

for i = 1:M
    
    fprintf('i = %i\n',i)
    
    % Proposing a new pair of parameters
    theta_tilde_curr = para_transformation(theta_curr);
    theta_tilde_prop = mvnrnd(theta_tilde_curr,cov_rw);
    theta_prop = para_back_transformation(theta_tilde_prop);
    prob = jacobian_transformation(theta_tilde_prop) / jacobian_transformation(theta_tilde_curr);
            
    % Simulating n sets of data and taking their summary statistics
    ssx = zeros(n,ns);
    parfor k = 1:n
        X = simulate_toads2(theta_prop,ntoads,ndays);
        X(NaN_Pos) = NaN;
        ssx(k,:) = feval(sumstat_fun,X,lag);
    end
    
    % Calculating the mean and covariance of the summary statistics
    the_mean = mean(ssx)';
    the_cov = warton(ssx, gamma);
	
    loglike_ind_prop = -0.5*log(det(the_cov)) - 0.5*(ssy-the_mean)'/the_cov*(ssy-the_mean);

    % accept/reject step
    if (prob * exp(loglike_ind_prop - loglike_ind_curr) > rand)
        fprintf('*** accept ***\n');
        theta_curr = theta_prop;
        loglike_ind_curr = loglike_ind_prop;
    end
    
    theta(i,:) = theta_curr;
    loglike(i) = loglike_ind_curr;   
    
end


end

