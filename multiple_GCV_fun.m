function [gcv_m]=multiple_GCV_fun(lambda,B,y,R)

    ngenes = size(y,2);
    for g=1:ngenes
    [gcv(g)]=GCV_fun(lambda,B,y(:,g),R);
    end
    
    gcv_m = nansum(gcv);    
    
end