

% script to obtain the results based on the 100 runs of the approximate
% inference algorithm (change the .mat file as needed)

% load results
load('results_is_sl_T10000_n100.mat')

I_05 = zeros(100,1);
I_10 = zeros(100,1);
I_20 = zeros(100,1);

Ir_05 = zeros(100,1);
Ir_10 = zeros(100,1);
Ir_20 = zeros(100,1);

% obtain mcmc acceptance and coverage rates
parfor i = 1:100
    %i
    theta = theta_is(:,:,i);
    w = w_is(:,i);
    ess(i) = 1/sum(w.^2);
    
    q025 = quantile_weighted(theta(:,1),0.025,w);
    q05 = quantile_weighted(theta(:,1),0.05,w);
    q10 = quantile_weighted(theta(:,1),0.1,w);
    q90 = quantile_weighted(theta(:,1),0.9,w);
    q95 = quantile_weighted(theta(:,1),0.95,w);
    q975 = quantile_weighted(theta(:,1),0.975,w);
    
    theta1_95(i) = q025<0.6 && q975>0.6;
    theta1_90(i) = q05<0.6 && q95>0.6;
    theta1_80(i) = q10<0.6 && q90>0.6;
    
    q025 = quantile_weighted(theta(:,2),0.025,w);
    q05 = quantile_weighted(theta(:,2),0.05,w);
    q10 = quantile_weighted(theta(:,2),0.1,w);
    q90 = quantile_weighted(theta(:,2),0.9,w);
    q95 = quantile_weighted(theta(:,2),0.95,w);
    q975 = quantile_weighted(theta(:,2),0.975,w);
    
    theta2_95(i) = q025<0.2 && q975>0.2;
    theta2_90(i) = q05<0.2 && q95>0.2;
    theta2_80(i) = q10<0.2 && q90>0.2;
    
    % for KDE do some resampling to remove the weights
    r = randsample(1:length(theta(:,1)), 1000, 'true', w);
    thetar = theta(r,:);
    
    [f,~] = ksdensity(thetar,thetar);
    [f_true,~] = ksdensity(thetar,[0.6 0.2]);
    
    
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



