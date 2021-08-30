

% produce results from the replicate simulations

% uncomment code as required
load('results_bsl_n50_shrink_adj_n500.mat')


theta_true = [1.7 35 0.6];

% obtain mcmc acceptance and coverage rates
parfor i = 1:50
   
    theta_rep = theta(:,:,i);
    
    ess(i) = multiESS(theta_rep);
    
    q025 = quantile(theta_rep, 0.025);
    q05 = quantile(theta_rep, 0.05);
    q10 = quantile(theta_rep, 0.1);
    q90 = quantile(theta_rep, 0.9);
    q95 = quantile(theta_rep, 0.95);
    q975 = quantile(theta_rep, 0.975);
    
    I05(i,:) = q025 < theta_true & q975 > theta_true;
    I10(i,:) = q05 < theta_true & q95 > theta_true;
    I20(i,:) = q10 < theta_true & q90 > theta_true;
    
    [f,~] = ksdensity(theta_rep(:, [1 2]),theta_rep(:, [1 2]));
    [f_true,~] = ksdensity(theta_rep(:, [1 2]),[theta_true(1) theta_true(2)]);
    
    f_05 = quantile(f,0.05);
    f_10 = quantile(f,0.1);
    f_20 = quantile(f,0.2);
    
    if (f_true > f_05)
        I05_12(i) = 1;
    end
    
    if (f_true > f_10)
        I10_12(i) = 1;
    end
    
    if (f_true > f_20)
        I20_12(i) = 1;
    end
    
    
    [f,~] = ksdensity(theta_rep(:, [1 3]),theta_rep(:, [1 3]));
    [f_true,~] = ksdensity(theta_rep(:, [1 3]),[theta_true(1) theta_true(3)]);
    
    f_05 = quantile(f,0.05);
    f_10 = quantile(f,0.1);
    f_20 = quantile(f,0.2);
    
    if (f_true > f_05)
        I05_13(i) = 1;
    end
    
    if (f_true > f_10)
        I10_13(i) = 1;
    end
    
    if (f_true > f_20)
        I20_13(i) = 1;
    end
    
    [f,~] = ksdensity(theta_rep(:, [2 3]),theta_rep(:, [2 3]));
    [f_true,~] = ksdensity(theta_rep(:, [2 3]),[theta_true(2) theta_true(3)]);
    
    f_05 = quantile(f,0.05);
    f_10 = quantile(f,0.1);
    f_20 = quantile(f,0.2);
    
    if (f_true > f_05)
        I05_23(i) = 1;
    end
    
    if (f_true > f_10)
        I10_23(i) = 1;
    end
    
    if (f_true > f_20)
        I20_23(i) = 1;
    end
    
end



