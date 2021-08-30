function s = summStat_quantiles3(X,lag)
% summary statistic function for toad model

nlag = length(lag);

x = cell(nlag,1);
s = [];  
for k = 1:nlag
    l = lag(k);
    x = obsMat2deltax(X, l);
    x = abs(x);
    
    return_ind = x<10;
    x_noret = x(~return_ind);
    
    s = [s sum(return_ind) log(diff(quantile(x_noret, 0:0.1:1))) median(x_noret)];
end

end