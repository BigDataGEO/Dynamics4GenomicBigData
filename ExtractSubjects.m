function rowsWhereSubjectWasFound = ExtractSubjects(subjectIDParts, Matrix)
  rowsWhereSubjectWasFound = [];
  
  % Function strjoin, used further below, needs to receive an array of strings. Therefore, if the function received only one string then it must be converted to an array of strings containing only the original one.
  if ischar(subjectIDParts)
    subjectIDParts = {subjectIDParts};
  end
  
  [numberOfRows,numberOfColumns] = size(Matrix);
  
  for i=1:numberOfRows
    for j=1:numberOfColumns
      stringToSearch = Matrix{i,j};
      regularExpressionToSearchFor = strjoin(subjectIDParts,'.*');
      startIndex = regexp(stringToSearch, regularExpressionToSearchFor);
      if (startIndex~=0)
	rowsWhereSubjectWasFound = cat(2, rowsWhereSubjectWasFound, i);
      end
    end
  end
  rowsWhereSubjectWasFound = rowsWhereSubjectWasFound.';
end
