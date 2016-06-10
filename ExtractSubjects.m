% Input argument subjectIDParts is either a string or an array of strings.
% Input argument Matrix is a matrix of strings.

% The function returns the row numbers of Matrix whose cells contain all of the strings in subjectIDParts that do not start with '~' and excluding those who contain at least one of the strings in subjectIDParts that start with '~'.


% subjectIDParts = {'BBB_CCC', '~min'}
% Matrix={'AAAAA_BBB_CCC_25.3hrs_DDD_XXXX', 'AAAAA_BBB_CCC_3435.3min_DDD_XXXX', 'AAAAA_BBB_CCC_65.3hrs_DDD_XXXX', 'AAAAA_BBB_CCC_656.3hrs_DDD_XXXX','AAAAA_BBB_CCC_965.3hrs_DDD_XXXX'}'

% Example output:
%  ans =
%  
%       1
%       3
%       4
%       5

function rowsWhereSubjectWasFound = ExtractSubjects(subjectIDParts, Matrix)

  subjectNamePartsToInclude = [];
  subjectNamePartsToExclude = [];
  
  % Function strjoin, used further below, needs to receive an array of strings. Therefore, if the function received only one string then it must be converted to an array of strings containing only the original one.
  if ischar(subjectIDParts)
    subjectIDParts = {subjectIDParts};
  end
  
  for i=1:length(subjectIDParts)
    % If this condition is true, then the string must be excluded.
    if(strcmp(subjectIDParts{i}(1:1),'~'))
      subjectNamePartsToExclude = [subjectNamePartsToExclude, {subjectIDParts{i}(2:length(subjectIDParts{i}))}];
    % Otherwise the string must be included.
    else
      subjectNamePartsToInclude = [subjectNamePartsToInclude, {subjectIDParts{i}}];
    end
  end
  
  rowsWhereSubjectWasFound = ExtractSubjectsWithInclusionsAndExclusions(subjectNamePartsToInclude, subjectNamePartsToExclude, Matrix);
end


function rowsWhereSubjectWasFound = ExtractSubjectsWithInclusionsAndExclusions(subjectNamePartsToInclude, subjectNamePartsToExclude, Matrix)
  rowsWhereSubjectWasFound = [];
  
  [numberOfRows,numberOfColumns] = size(Matrix);
  
  for i=1:numberOfRows
    for j=1:numberOfColumns
      stringToSearch = Matrix{i,j};
      regularExpressionToSearchFor = strjoin(subjectNamePartsToInclude,'.*');     
      startIndex = regexp(stringToSearch, regularExpressionToSearchFor);      
      % If the following condition is true, it means that stringToSearch at least contains the substrings that in subjectNamePartsToInclude.
      if (startIndex~=0)
	if(length(subjectNamePartsToExclude) == 0)
	  rowsWhereSubjectWasFound = cat(2, rowsWhereSubjectWasFound, i);
	else
	  stringToSearchShouldBeIncluded = 1;
	  for k=1:length(subjectNamePartsToExclude)
	    % if 1, then the string should be included. If [] then the string should not be included.
	    found = regexp(stringToSearch, subjectNamePartsToExclude{k});
	    
	    % If this condition is true, then the string does not contain the substring to be excluded. Thus the string should be included in the output list.
	    if isempty(found)
	      % Nothing to do here.
	    else % Otherwise the string contains the substring to be excluded. Thus the string should not be included in the output list.
	      stringToSearchShouldBeIncluded = 0;
	      break;
	    end
	  end
	  if(stringToSearchShouldBeIncluded)
	    rowsWhereSubjectWasFound = cat(2, rowsWhereSubjectWasFound, i);
	  end
	end
      end
    end
  end
  rowsWhereSubjectWasFound = rowsWhereSubjectWasFound.';
end
