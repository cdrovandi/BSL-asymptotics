function theta_is_adj = abc_regadj_ma_acf(y, theta_is, w_is, ssx_is, cov_abc, K)
% run regression adjustment ABC on the results for the 100 simulated
% datasets
inv_cov_abc = inv(cov_abc);

for k = 1:100
    
    ys = y(:,k);
    
    theta = theta_is(:,:,k);
    ssx = ssx_is(:,:,k);
    w = w_is(:,k);
    
    theta_adj = theta;
    
    ssy = summStat_ma2_acf(ys, K);
    
    dist = zeros(length(theta(:,1)),1);
    ssdiff = ssx;
    for i = 1:length(theta(:,1))
        ssdiff(i,:) = ssx(i,:) - ssy;
        dist(i) = sqrt(ssdiff(i,:)*inv_cov_abc*ssdiff(i,:)');
    end
    
    
    the_weights = 0.75*(1 - (dist./max(dist)).^2);
    
    b  = glmfit(ssdiff,theta(:,1),'normal','weights',the_weights);
    theta_adj(:,1) = theta(:,1) - ssdiff*b(2:end);
    
    b  = glmfit(ssdiff,theta(:,2),'normal','weights',the_weights);
    theta_adj(:,2) = theta(:,2) - ssdiff*b(2:end);
    
    theta_is_adj(:,:,k) = theta_adj;
    
    
end

end

