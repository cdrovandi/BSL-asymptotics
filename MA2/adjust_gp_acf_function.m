function theta_adj = adjust_gp_acf_function(y,theta,gamma,L,n)
% BSL adjustment function


T = length(y);

J = 200;
B = 200;

% loglikelihood estimate tuning parameters
%n  = 100;

[M,num_params] = size(theta);

mean_approx = mean(theta);
cov_approx = cov(theta);
std_approx = std(theta);

theta_point =  mean_approx;

theta_test = lhsdesign(B,num_params);

theta_min = theta_point - std_approx;
theta_max = theta_point + std_approx;
theta_min = repmat(theta_min,B,1);
theta_max = repmat(theta_max,B,1);


theta_test = theta_min + (theta_max - theta_min).*theta_test;


for b = 1:B
    
    x = simulate_ma2(theta_test(b,:),T,n);
    ssx = summStat_ma2_acf(x, L);
    
    the_mean = mean(ssx);
    the_cov = warton(ssx, gamma);
    
    means{b} = the_mean;
    icovs{b} = inv(the_cov);
    logdetcovs{b} = log(det(the_cov));

end


grad = zeros(J,num_params);
epsilon_fd = std_approx/1000;

for j = 1:J
    
    x = simulate_ma2(theta_point,T,1);
    ssy = summStat_ma2_acf(x, L);
    
    loglike_test = zeros(B,1);
    for b = 1:B
        loglike_test(b) = -0.5*logdetcovs{b} - 0.5*(ssy-means{b})*icovs{b}*(ssy-means{b})';
    end
    
    
    gp_model = fitrgp(theta_test,loglike_test,'KernelFunction','ardsquaredexponential','BasisFunction','none');

    S = gp_model.KernelInformation.KernelParameters(1:num_params);
    
    D = pdist(theta_test,'seuclidean',S).^2;
    D = squareform(D);
    K = gp_model.KernelInformation.KernelParameters(end)^2*exp(-0.5*D);
    
    k = pdist2(theta_point,theta_test,'seuclidean',S).^2;
    k = gp_model.KernelInformation.KernelParameters(end)^2*exp(-0.5*k);
    
    % predict at point estimate
    %A = inv(K + eye(B)*gp_model.Sigma^2)*loglike_test;
    A = (K + eye(B)*gp_model.Sigma^2)\loglike_test;
    
    for d = 1:num_params
        delta_fd = zeros(1,num_params);
        delta_fd(d) = 1;
        
        theta_prop = theta_point + epsilon_fd.*delta_fd;
        k = pdist2(theta_prop,theta_test,'seuclidean',S).^2;
        k = gp_model.KernelInformation.KernelParameters(end)^2*exp(-0.5*k);
        
        qr = k*A;
        
        theta_prop = theta_point - epsilon_fd.*delta_fd;
        k = pdist2(theta_prop,theta_test,'seuclidean',S).^2;
        k = gp_model.KernelInformation.KernelParameters(end)^2*exp(-0.5*k);
        
        ql = k*A;
        
        grad(j,d) = 0.5*(qr - ql)/epsilon_fd(d);
        
    end
    
end



Omega = cov(grad);
Gamma = cov(theta);
P = Gamma*sqrtm(Omega)*sqrtm(inv(Gamma));
theta_adj = theta;
for i = 1:M
    theta_adj(i,:) = theta_point + (P*(theta(i,:) - theta_point)')';
end


end




