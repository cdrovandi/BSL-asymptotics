
% code for adjusting the shrinkage BSL results for simulated data

load('data_toads_model1.mat')
load('results_simulated_bsl_n50_shrink.mat')

lag = [1, 2, 4, 8];
model = 1;
d0 = NaN;
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






% code for plotting results

load('results_simulated_bsl_n50_shrink.mat')
theta_shrink = theta;
load('results_simulated_bsl_n500.mat')




subaxis('MT',0.02,'MB',0.05,'ML',0.03,'MR',0.02,'PL',0.02,'PR',0.01,'Pt',0.015,'PB',0.05);


subaxis(2,3,1);

[kx,ky] = ksdensity(theta(:,1:2));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));

contour(X,Y,Z,'LineWidth',1.5,'ShowText','Off','LineColor','k');
title('shrinkage \gamma = 0.1');
hold on;


[kx,ky] = ksdensity(theta_shrink(:,1:2));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,'LineWidth',1.5,'ShowText','Off','LineColor','r','LineStyle','--');
xlim([1.45, 1.9]);
ylim([30, 40]);
legend('standard','shrinkage');
xlabel('\alpha','FontSize',16);
ylabel('\gamma','FontSize',16);


subaxis(2,3,2);

[kx,ky] = ksdensity(theta(:,[1 3]));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));

contour(X,Y,Z,'LineWidth',1.5,'ShowText','Off','LineColor','k');
title('shrinkage \gamma = 0.1');
hold on;


[kx,ky] = ksdensity(theta_shrink(:,[1 3]));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,'LineWidth',1.5,'ShowText','Off','LineColor','r','LineStyle','--');
xlim([1.45, 1.9]);
ylim([0.56, 0.65]);
xlabel('\alpha','FontSize',16);
ylabel('p_0','FontSize',16);
legend('standard','shrinkage');


subaxis(2,3,3);

[kx,ky] = ksdensity(theta(:,2:3));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));

contour(X,Y,Z,'LineWidth',1.5,'ShowText','Off','LineColor','k');
title('shrinkage \gamma = 0.1');
hold on;


[kx,ky] = ksdensity(theta_shrink(:,2:3));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,'LineWidth',1.5,'ShowText','Off','LineColor','r','LineStyle','--');
xlim([30, 40]);
ylim([0.56, 0.65]);
legend('standard','shrinkage');
xlabel('\gamma','FontSize',16);
ylabel('p_0','FontSize',16);





subaxis(2,3,4);

[kx,ky] = ksdensity(theta(:,1:2));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));

contour(X,Y,Z,'LineWidth',1.5,'ShowText','Off','LineColor','k');
title('adjusted');
hold on;

[kx,ky] = ksdensity(theta_adj(:,1:2));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,'LineWidth',1.5,'ShowText','Off','LineColor','g','LineStyle','--');
title('adjusted');
xlim([1.45, 1.9]);
ylim([30, 40]);
legend('standard','adjusted');
xlabel('\alpha','FontSize',16);
ylabel('\gamma','FontSize',16);







subaxis(2,3,5);

[kx,ky] = ksdensity(theta(:,[1 3]));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));

contour(X,Y,Z,'LineWidth',1.5,'ShowText','Off','LineColor','k');
hold on;

[kx,ky] = ksdensity(theta_adj(:,[1 3]));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,'LineWidth',1.5,'ShowText','Off','LineColor','g','LineStyle','--');
title('adjusted');
xlim([1.45, 1.9]);
ylim([0.56, 0.65]);
legend('standard','adjusted');
xlabel('\alpha','FontSize',16);
ylabel('p_0','FontSize',16);






subaxis(2,3,6);

[kx,ky] = ksdensity(theta(:,2:3));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));

contour(X,Y,Z,'LineWidth',1.5,'ShowText','Off','LineColor','k');
hold on;

[kx,ky] = ksdensity(theta_adj(:,2:3));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,'LineWidth',1.5,'ShowText','Off','LineColor','g','LineStyle','--');
title('adjusted');
xlim([30, 40]);
ylim([0.56, 0.65]);
legend('standard','adjusted');
xlabel('\gamma','FontSize',16);
ylabel('p_0','FontSize',16);




