function [gcv]=GCV_fun(lambda,B,y,R)

    n      = length(y);
    
    M      = B'*B + lambda.*R;
    [Bmatfac, rnkp1] = chol(M);
    MB     = Bmatfac\(Bmatfac'\(B'));
    C      = MB*y;
    %Evaluate S which maps y^hat to y
    F      = B*C; 
    %Compute the 1/n * Sum of squared Errors
    Err    = (y-F)'*(y-F);
    %Compute the 1/n * error degrees of freedom
    %Partial derivative of c w.r.t. beta 
    S        = B*(MB);
    df       = trace(S);
    dfe      = n-df;
    %Compute the generalised cross validation
    gcv      = (Err./n)./((n - df)/n)^2;
    
%     u            = (n.*Err);
%     v            = (dfe.^2);
%     d_dfe        = -2.*dfe*trace(-MB'*R*MB);
%     d_Err_d_c    = n.*(-2.*y'*B+2.*C'*(B'*B));
%     d_c_d_lam    = -M'*R*C;
%     Dgcv         = (d_Err_d_c*d_c_d_lam*v -  u*d_dfe )./ (v).^2; 
    
end