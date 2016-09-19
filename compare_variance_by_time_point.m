function variance_p_values = compare_variance_by_time_point(expression_of_first_subject, expression_of_second_subject)
  variance_p_values = [];
  
  number_of_time_points = min(size(expression_of_first_subject,2), size(expression_of_second_subject,2));
  
  for current_time_point=1:number_of_time_points
    test_data = [expression_of_first_subject(:, current_time_point) expression_of_second_subject(:, current_time_point)];
    [p_value, stats] = vartestn(test_data, 'Display','off', 'TestType', 'LeveneAbsolute');
    variance_p_values(current_time_point) = p_value;
  end
  
end
