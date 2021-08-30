function run_adjust_quantilesumms_gp_n500(r)

load('data_simulated_replicates.mat')
ntoads = ntoads;  ndays = ndays;

Y = Xr{r};

load(['results_bsl_n50_shrink' num2str(r) '.mat']);

tic;

lag = [1, 2, 4, 8];

ssy = summStat_quantiles3(Y, lag);

fn_simulate_toads = @(theta_prop,ntoads,ndays) simulate_toads2(theta_prop,ntoads,ndays);

ns = length(ssy);

J = 200;
B = 200;

% loglikelihood estimate tuning parameters
n  = 500;
gamma = 0.1;

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


parfor b = 1:B
    ssx = zeros(n,ns);

    for k = 1:n
        X = fn_simulate_toads(theta_test(b,:),ntoads,ndays);
        ssx(k,:) = summStat_quantiles3(X,lag);
    end
    
    the_mean = mean(ssx);
    the_cov = warton(ssx,gamma);
    
    means{b} = the_mean;
    icovs{b} = inv(the_cov);
    logdetcovs{b} = log(det(the_cov));

end


grad = zeros(J,num_params);
epsilon_fd = std_approx/1000;

parfor j = 1:J
    
    Y = fn_simulate_toads(theta_point,ntoads,ndays);
    ssy = summStat_quantiles3(Y,lag);
    
    loglike_test = zeros(B,1);
    for b = 1:B
        loglike_test(b) = -0.5*logdetcovs{b} - 0.5*(ssy-means{b})*icovs{b}*(ssy-means{b})';
    end
    
    
    gp_model = fitrgp(theta_test,loglike_test,'KernelFunction','ardsquaredexponential','BasisFunction','none');
    
    
    S = gp_model.KernelInformation.KernelParameters(1:num_params);
    
    D = pdist(theta_test,'seuclidean',S).^2;
    D = squareform(D);
    K = gp_model.KernelInformation.KernelParameters(end)^2*exp(-0.5*D);
    
    % predict at point estimate
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
P = Gamma*Omega^0.5*Gamma^-0.5;
theta_adj = theta;
for i = 1:M
    theta_adj(i,:) = theta_point + (P*(theta(i,:) - theta_point)')';
end


finaltime = toc;

save(['results_bsl_n50_shrink_adj_n500' num2str(r) '.mat'],'theta_adj','finaltime');

end

