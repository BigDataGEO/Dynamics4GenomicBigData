function [rx] = round2(x)
%% Round to two decimal places
 
rx = round(x.*100)./100;
end