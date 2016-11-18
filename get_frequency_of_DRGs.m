function frequency_of_DRGs =  get_frequency_of_DRGs(DRGs)
  
  all_genes_with_repetition = [];
  
  n = size(DRGs,2);
  
  for k=1:n
    all_genes_with_repetition = [all_genes_with_repetition; unique(DRGs{k})];
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