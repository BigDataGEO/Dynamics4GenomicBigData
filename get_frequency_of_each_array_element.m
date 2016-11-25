% Receives a cell array of strings of size Nx1 and returns an Nx2 cell array where the first column is the elements in the input array and the second column is the frequency of each element.

%  the_array = [{'a'}; {'s'}; {'a'}; {'e'}; {'q'}];

% Returns 
%  frequency_per_element = 
%  
%      'a'    [2]
%      'e'    [1]
%      'q'    [1]
%      's'    [1]

function frequency_per_element = get_frequency_of_each_array_element(the_array)

  [a b c] = unique(the_array);
  d = hist(c,length(a));
  P = [a num2cell(d')];
  
  [B I] = sort(cell2mat(P(:,2)), 'descend');
  
  frequency_per_element = P(I,:);  
end