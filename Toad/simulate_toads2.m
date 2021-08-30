function [X] = simulate_toads2(theta, ntoads, ndays)
% simulate toads location for each toad and day
% INPUT:
% theta - model parameters, alpha, gamma and p0 (the probability of returning to a previous refuge)
% ntoads - the number of individual toads
% ndays - the number of days for simulation
% model - indicator for model version in the paper Marchand, et al (2017),
%         1 stands for random return
% OUTPUT:
% X - a ndays by ntoads matrix contains toads location for each toad and

alpha = theta(1);
gamma = theta(2);
p0 = theta(3);
X = zeros(ndays, ntoads);

for i = 2:ndays
    % toads that stay at new location
    ind = rand(1,ntoads) >= p0;
    deltax = rndlas(sum(ind),gamma,alpha)'; % distance of move (only generate distance for toads that stay at new location) 
    X(i,ind) = X(i-1,ind) + deltax;
    
    % return to one of the previous refuge sites
    % multiple visits to a refuge site increases the weighting
    ind_refuge = randsample(i-1,ntoads-sum(ind),true)';
    idx = sub2ind(size(X),ind_refuge,find(~ind));
    X(i,find(~ind)) = X(idx);
end


		