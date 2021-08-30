
% same as produce results file but not based on weighted sample
% particularly needed for adjusted BSL results

% load results
load('results_adj_T100_K20_n500.mat')

Ir_05 = zeros(100,1);
Ir_10 = zeros(100,1);
Ir_20 = zeros(100,1);

% obtain mcmc acceptance and coverage rates
parfor i = 1:100
    theta = theta_is(:,:,i);
    
    q025 = quantile(theta(:,1),0.025);
    q05 = quantile(theta(:,1),0.05);
    q10 = quantile(theta(:,1),0.1);
    q90 = quantile(theta(:,1),0.9);
    q95 = quantile(theta(:,1),0.95);
    q975 = quantile(theta(:,1),0.975);
    
    theta1_95(i) = q025<0.6 && q975>0.6;
    theta1_90(i) = q05<0.6 && q95>0.6;
    theta1_80(i) = q10<0.6 && q90>0.6;
    
    q025 = quantile(theta(:,2),0.025);
    q05 = quantile(theta(:,2),0.05);
    q10 = quantile(theta(:,2),0.1);
    q90 = quantile(theta(:,2),0.9);
    q95 = quantile(theta(:,2),0.95);
    q975 = quantile(theta(:,2),0.975);
    
    theta2_95(i) = q025<0.2 && q975>0.2;
    theta2_90(i) = q05<0.2 && q95>0.2;
    theta2_80(i) = q10<0.2 && q90>0.2;
    
    
    [f,~] = ksdensity(theta,theta);
    [f_true,~] = ksdensity(theta,[0.6 0.2]);
    
    
    f_05 = quantile(f,0.05);
    f_10 = quantile(f,0.1);
    f_20 = quantile(f,0.2);
    
    if (f_true > f_05)
        Ir_05(i) = 1;
    end
    
    if (f_true > f_10)
        Ir_10(i) = 1;
    end
    
    if (f_true > f_20)
        Ir_20(i) = 1;
    end
    
  
end



