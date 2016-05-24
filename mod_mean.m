function[y] = mod_mean(x)

if(size(x,1)>1) 
    y = mean(x);
else
    y = x ;
end

end
