function [F] = Ftest(data,argvals,fd,df)

n = length(argvals);

% --- Null Hypiothesis is that the gene_t = mu_t
% as data is already centered
%[m,r,n] = size(data);

RSS_H0 = 0;
%for i = 1:r
%d = reshape(data(:,i,:),m,n);     
RSS_H0 = RSS_H0 + nansum(data.^2,2);
%end

% --- Alternative Hypiothesis is that the gene_t = f(t)
yhat    = eval_fd(argvals, fd);
RSS_H1  = 0;
%for i = 1:r
%d = reshape(data(:,i,:),m,n);     
RSS_H1 = RSS_H1 + nansum((data-yhat').^2,2);
%end

% --- F - test adjusted for degrees of freedom
F = ((RSS_H0 - RSS_H1)/(df-1))./(RSS_H1/(n-df));

if(any(RSS_H1<1e-8))
ind = find(RSS_H1<1e-8);     
F(ind) = min(F);
end

