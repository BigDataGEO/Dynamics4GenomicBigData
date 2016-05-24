% ZC number of zero crossings in x
% [n] = zc(x) calculates the number of zero crossings in x

function [n] = zc(x)

x = x(:,~isnan(x(1,:)));

n = zeros(size(x,1),1);
for i = 1:size(x,1)    
    s=sign(x(i,:));
    t=filter([1 1],1,s);
    n(i) = length(find(t==0));
end

