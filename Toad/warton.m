function Sigma_warton = warton(y, gamma)

% compute Warton shrinkage covariance estimator

[~,ns] = size(y);

% sample mean and covariance matrix
muhat = mean(y);
Sigmahat = cov(y);

% diagonal entries of Sigmahat
Dhat = diag(diag(Sigmahat));

Dhat1 = Dhat^(0.5);
Dhat2 = Dhat^(-0.5);

% sample correlation matrix
Chat =  Dhat2 * Sigmahat * Dhat2;

corr_warton = gamma*Chat + (1-gamma)*eye(ns);

Sigma_warton = Dhat1 * corr_warton * Dhat1;

end