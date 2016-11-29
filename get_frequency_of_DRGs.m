% This function finds the intersection of probe ids in the list of DRGs provided.

% This intersection of (DRG) probe ids is returned in variable common_probes.

% The list of gene names associated with each one of the common probes is returned in frequency_of_DRGs along with the number of (common) probes where the gene appears (which could be more than one).

% Example: Suppose the input is composed of the following two matrices.

% | Probe | Gene |		| Probe | Gene |
%     A      G1  		   X        G25
%     B      G52 		   B       G52
%     C      G41   		   D       G12
%     H      G15                   C       G41
%     F      G52                   F       G52

% The common probes are B, C, and F.

% The common DRGs are G52 (frequency 1), G41 (frequency 1), and G52 (frequency 2).

% Input:

% list_of_statistically_significant_DRGs is a cell array. Each element is a cell array of size Mx2 where the first column is the (DRG) probe ids and the second column the corresponding gene name, listing the DRGs of a subject/condition. M is the number of (DRG) probes.

% Output:

% frequency_of_DRGs is a Nx2 cell array, where the first column is the gene names (as strings) and the second column is the frequency, as numbers.
% common_probes is a Mx1 cell array with the probe ids of that are DRGs across all the subject/conditions.

function [frequency_of_DRGs, common_probes] =  get_frequency_of_DRGs(list_of_statistically_significant_DRGs)
  
  intersection_of_probes = list_of_statistically_significant_DRGs{1}(:,1);
  
  for k=1:length(list_of_statistically_significant_DRGs)
    intersection_of_probes = intersect(intersection_of_probes, list_of_statistically_significant_DRGs{k}(:,1));
  end
  
  A = list_of_statistically_significant_DRGs{1}(:,1);
  
  B = intersection_of_probes;
  
  common_probes = intersection_of_probes;
  
  intersection_of_probes_and_genes = list_of_statistically_significant_DRGs{1}(find(ismember(A,B)),:);
  
  frequency_of_DRGs = get_frequency_of_each_array_element(intersection_of_probes_and_genes(:,2));
  
  frequency_of_probes = get_frequency_of_each_array_element(intersection_of_probes_and_genes(:,1));
  
end


