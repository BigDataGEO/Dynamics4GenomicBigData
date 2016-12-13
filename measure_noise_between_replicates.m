% Input
% gene_expression_of_replicates is a R-sized cell array where each element is an MxN matrix of double values representing the gene expression of one replicate. The rows of the matrix are the genes and the columns are the time points. It is assumed that all matrices have the same size. That is to say, it is assumed that all replicates provided as input have the same number of genes and time points.

% Output
% A column (vertical) vector of doubles where the k-th element is a measure of the noise observed in the expression of the k-th gene in the expression matrices provided as input.

% Example
% Two 'replicates', A and B, with three genes and seven time points.

%  A = [2 4 5 9 1 2 3; 1 5 9 4 2 1 5; 2 7 9 4 8 5 9];
%  
%  B = [9 6 1 3 4 5 2; 2 2 1 6 4 3 7; 2 9 6 4 7 1 1];
%  
%  gene_expression_of_replicates = [{A} {B}];
%  
%  noise_per_gene = measure_noise_between_replicates(gene_expression_of_replicates);
%  
%  % Returns
%  
%  noise_per_gene =
%  
%      0.6529
%      0.5580
%      0.3754
%
%  %  The third gene is the most consistent across the two replicates.

function noise_per_gene = measure_noise_between_replicates(gene_expression_of_replicates)

  noise_per_gene = [];
  
  for gene_index = 1:size(gene_expression_of_replicates{1},1)
  
    expression_of_current_gene_across_replicates = [];
    
    for replicate_index=1:length(gene_expression_of_replicates)
    
      gene_expression_of_replicate = gene_expression_of_replicates{replicate_index};
      
      expression_of_current_gene_in_current_replicate = gene_expression_of_replicate(gene_index,:);
      
      expression_of_current_gene_across_replicates = [expression_of_current_gene_across_replicates; expression_of_current_gene_in_current_replicate];
    end
    
    % In the following four lines, the coefficient of variation of each gene across the R replicates.
%      mean_of_current_gene = mean(expression_of_current_gene_across_replicates);    
%      std_of_current_gene = std(expression_of_current_gene_across_replicates);    
%      coefficient_of_variation_of_current_gene = std_of_current_gene./mean_of_current_gene;    
%      noise_of_current_gene = mean(coefficient_of_variation_of_current_gene);
    
    % The following line is an alternative to the previous four lines and measures noise as the
    % average of the gene's standard deviation across the R replicates.
    noise_of_current_gene = mean(std(expression_of_current_gene_across_replicates));
    
    noise_per_gene = [noise_per_gene; noise_of_current_gene];
  
  end

end


