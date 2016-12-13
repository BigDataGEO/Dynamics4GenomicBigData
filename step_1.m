% Input
% Parameter samples is an MxN cell array of strings where each element is the accession number of a
% GEO sample. In other words, each element must start with the prefix 'GSM'. M is the number of
% time points and N is the number of replicates that comprise the experimental condition. If N>1
% then the raw gene expression returned is the average of the N replicates.

% Parameter time_points is a Mx1 cell array of strings where each element is a time point expressed
% with units, e.g., '2 hours'. The dimension of this array must be consistent with that of
% parameter samples.

% Output

function [raw_gene_expression, raw_time_points] = step_1(geoStruct, samples, time_points)
  
  GSE_matrix = double(geoStruct.Data);
  
  raw_gene_expression = [];
  
  for replicate_index = 1:size(samples,2)
    samples_of_current_replicate = samples(:,replicate_index);
    
    [indices_of_samples_in_matrix, not_found] = find_in_cell_array_of_strings(geoStruct.Header.Samples.geo_accession, samples_of_current_replicate);
    
    raw_gene_expression_of_replicate{replicate_index} = GSE_matrix(:,indices_of_samples_in_matrix);
    
    raw_gene_expression = cat(3, raw_gene_expression, raw_gene_expression_of_replicate{replicate_index});
  end
  
  raw_gene_expression = mean(raw_gene_expression,3);
  
  % The following two lines incorporate Vahed's method.
%    mean_per_time_point = mean(raw_gene_expression);
%    raw_gene_expression = raw_gene_expression - repmat(mean_per_time_point, size(raw_gene_expression,1), 1);  

  raw_time_points = ExtractTimePoints(time_points');
  
  raw_time_points = cell2mat(raw_time_points(:,1));
  
end
