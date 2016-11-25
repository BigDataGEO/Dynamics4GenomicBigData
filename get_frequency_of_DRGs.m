% Input:

% list_of_statistically_significant_DRGs is a cell array. Each element is a cell array of size Mx2 where the first column is the probe ids and the second column the corresponding gene name, listing the DRGs of a subject/condition.

% Output:

% frequency_of_DRGs is a Nx2 cell array, where the first column is the gene names (as strings) and the second column is the frequency, as numbers.
% common_probes is a Mx1 cell array with the probe ids of that are DRGs across all the subject/conditions.

function [frequency_of_DRGs, common_probes] =  get_frequency_of_DRGs(list_of_statistically_significant_DRGs)
  
  the_intersection = list_of_statistically_significant_DRGs{1}(:,1);
  
  for k=1:length(list_of_statistically_significant_DRGs)
    the_intersection = intersect(the_intersection, list_of_statistically_significant_DRGs{1}(:,1));
  end
  
  A = list_of_statistically_significant_DRGs{1}(:,1);
  
  B = the_intersection;
  
  common_probes = the_intersection;
  
  intersection_of_probes_and_genes = list_of_statistically_significant_DRGs{1}(find(ismember(A,B)),:);
  
  frequency_of_DRGs = get_frequency_of_each_array_element(intersection_of_probes_and_genes(:,2));
  
  frequency_of_probes = get_frequency_of_each_array_element(intersection_of_probes_and_genes(:,1));
  
end


