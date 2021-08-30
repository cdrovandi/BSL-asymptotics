function f = loglike_ma2_sparse(ssy,theta_curr)

ns = length(ssy);

% nOnes = ones(ns, 1) ;
% the_cov = diag((1 + theta_curr(1)^2 + theta_curr(2)^2) * nOnes, 0) + diag((theta_curr(1) + theta_curr(1)*theta_curr(2))*nOnes(1:ns-1), -1) + ...
%     diag((theta_curr(1) + theta_curr(1)*theta_curr(2))*nOnes(1:ns-1), 1) + ...
%     diag(theta_curr(2)*nOnes(1:ns-2), -2) + diag(theta_curr(2)*nOnes(1:ns-2), 2); 

i = 1:ns;
j = 1:ns;
v = ones(1,ns)*(1 + theta_curr(1)^2 + theta_curr(2)^2);

i = [i 1:(ns-1)];
j = [j 2:ns];
v = [v ones(1,ns-1)*(theta_curr(1) + theta_curr(1)*theta_curr(2))];

i = [i 2:ns];
j = [j 1:(ns-1)];
v = [v ones(1,ns-1)*(theta_curr(1) + theta_curr(1)*theta_curr(2))];

i = [i 1:(ns-2)];
j = [j 3:ns];
v = [v ones(1,ns-2)*theta_curr(2)];

i = [i 3:ns];
j = [j 1:(ns-2)];
v = [v ones(1,ns-2)*theta_curr(2)];

the_cov = sparse(i,j,v);

L = chol(the_cov);
logdetA = 2*sum(log(diag(L)));

f = -0.5*logdetA - 0.5*ssy/the_cov*ssy';

end

