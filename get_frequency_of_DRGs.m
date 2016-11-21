function frequency_of_DRGs =  get_frequency_of_DRGs(list_of_statistically_significant_DRGs)
  
  all_genes_with_repetition = [];
  
  for k=1:size(list_of_statistically_significant_DRGs,2)
%      all_genes_with_repetition = [all_genes_with_repetition; unique(list_of_statistically_significant_DRGs{k})];
    all_genes_with_repetition = [all_genes_with_repetition; list_of_statistically_significant_DRGs{k}];
  end
 
  frequency_of_DRGs = get_frequency_of_each_array_element(all_genes_with_repetition);
  
end

function frequency_per_element = get_frequency_of_each_array_element(the_array)

  [a b c] = unique(the_array);
  d = hist(c,length(a));
  P = [a num2cell(d')];
  
  [B I] = sort(cell2mat(P(:,2)), 'descend');
  
  frequency_per_element = P(I,:);  
end