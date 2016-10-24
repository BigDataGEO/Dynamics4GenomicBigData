function [F, F_critic] = Ftest(gene_expression, time_points, fd_smooth_coefficients, degrees_of_freedom)

  n = length(time_points);

  % --- Null Hypiothesis is that the gene_t = mu_t
  % as gene_expression is already centered
  %[m,r,n] = size(gene_expression);

  RSS_H0 = 0;
  %for i = 1:r
  %d = reshape(gene_expression(:,i,:),m,n);     
  RSS_H0 = RSS_H0 + nansum(gene_expression.^2,2);
  %end

  % --- Alternative Hypiothesis is that the gene_t = f(t)
  yhat    = eval_fd(time_points, fd_smooth_coefficients);
  RSS_H1  = 0;
  %for i = 1:r
  %d = reshape(gene_expression(:,i,:),m,n);     
  RSS_H1 = RSS_H1 + nansum((gene_expression-yhat').^2,2);
  %end

  % --- F - test adjusted for degrees of freedom
  F = ((RSS_H0 - RSS_H1)/(degrees_of_freedom-1))./(RSS_H1/(n-degrees_of_freedom));
  
  if(any(RSS_H1<1e-8))
    ind = find(RSS_H1<1e-8);     
    F(ind) = min(F);
  end
  
  % The DRGs will be determined through an upper one-tailed F test.
  % The DRGs will be those whose F statistics are greater than the F_{alpha, numerator_df, denominator_df}
  
  d1 = degrees_of_freedom-1;
  d2 = n-degrees_of_freedom;
  
  % The area below the F curve in (-infty, F_critic) is 0.95. Thus the are in (F_critic, infty) is 0.05.
  F_critic = finv(0.95, d1, d2);
end
 
  