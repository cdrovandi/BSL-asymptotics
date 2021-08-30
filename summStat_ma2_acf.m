function ssx = summStat_ma2_acf(x, K)
% summary statistic for the MA(2) example based on autocovariances

[T,n] = size(x);
ssx = zeros(n,K);

for k = 1:K
   ssx(:,k) = sum(x(1:(end-k),:).*x((k+1):end,:));
end
ssx = ssx./T;

end