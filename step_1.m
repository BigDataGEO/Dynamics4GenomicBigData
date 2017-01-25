% Input
% Parameter samples is an MxN cell array of strings where each element is the accession number of a
% GEO sample. In other words, each element must start with the prefix 'GSM'. M is the number of
% time points and N is the number of replicates that comprise the experimental condition. If N>1
% then the raw gene expression returned is the average of the N replicates. The expression of each
% gene, each time point is averaged across the N replicates. In other words, the N matrices are
% averaged, element-wise.

% Parameter time_points is a Mx1 cell array of strings where each element is a time point expressed
% with units, e.g., '2 hours'. The dimension of this array must be consistent with that of
% parameter samples.

% Output

function [raw_gene_expression, raw_time_points] = step_1(geoStruct, samples, time_points)
  
  GSE_matrix = double(geoStruct.Data);
  
  raw_gene_expression = [];
  
  for time_point_index = 1:size(samples,1)
  
    samples_at_current_time_point = samples(time_point_index,:);
    
    [indices_of_samples_in_matrix, not_found] = find_in_cell_array_of_strings(geoStruct.Header.Samples.geo_accession, samples_at_current_time_point);
    
    raw_gene_expression_at_time_point{time_point_index} = median(GSE_matrix(:,indices_of_samples_in_matrix),2);
  
  end
  
  raw_gene_expression_2 = cell2mat(raw_gene_expression_at_time_point);

  raw_time_points = ExtractTimePoints(time_points');
  
  raw_time_points = cell2mat(raw_time_points(:,1));
  
end

%  A = [2 5 3 6; 1 8 5 9; 1 4 8 9];
%  B = [2 5 2 5; 1 0 3 9; 0 0 1 3];
%  C = [2 1 1 5; 1 8 3 5; 5 0 1 0];
%  
%  Z = [];
%  
%  Z = cat(3, Z, A);
%  
%  Z = cat(3, Z, B);
%  
%  Z = cat(3, Z, C);
%  
%  % Method 1
%  Z = mean(Z,3);
%  
%  % Method 2
%  Z = median(Z,3);
%  
%  % Method 3: Vahed's
%  mean_per_time_point = mean(Z);
%  Z = Z - repmat(mean_per_time_point, size(Z,1), 1);  


% The parameters of this function are cell arrays of strings.
% The function returns the indices in the first array where members of the second array were found.
% Example
% Input
% cell_array_to_search = {'A', 'B', 'ABC', 'C', 'EFGZ', 'D', 'F'}
% cell_array_to_search_for = {'ABC', 'EFGZ', 'HOHO'}
% Output
% indices = [3; 5]
% not_found = ['HOHO']
function [indices, not_found] = find_in_cell_array_of_strings(cell_array_to_search, cell_array_to_search_for)
  indices = [];
  not_found = [];
  not_found_idx = [];
  for i=1:length(cell_array_to_search_for)

    idx = find(strcmp([cell_array_to_search(:)], strtrim(cell_array_to_search_for{i})));
      
    if(isempty(idx))
      not_found = [not_found; {cell_array_to_search_for{i}}];
      not_found_idx = [not_found_idx; i];
    else
      indices = [indices; idx];
    end
  end
end
