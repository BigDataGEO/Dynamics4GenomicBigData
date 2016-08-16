function [p_value] = wilcoxon_signed_rank_test(expression_of_first_subject, expression_of_second_subject)
  
  first_data_set = mean(expression_of_first_subject, 1);
  
  second_data_set = mean(expression_of_second_subject, 1);
  
%    p_value = signrank(first_data_set, second_data_set);
  
  p_value = ranksum(first_data_set', second_data_set');
