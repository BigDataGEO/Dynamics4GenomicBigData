function indices = find_strings_in_cell_array(cell_array_of_strings_to_search, cell_array_of_strings_to_find)

  indices = [];
  
  for i=1:length(cell_array_of_strings_to_find)
    indices_of_current_string = find(strcmp(cell_array_of_strings_to_search, cell_array_of_strings_to_find{i}));
    indices = [indices; indices_of_current_string];
  end
end
