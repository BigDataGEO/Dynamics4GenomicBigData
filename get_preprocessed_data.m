% Input

% raw_gene_expression : This is the raw gene expression data of the subject as a MxN matrix, where M is the number of genes and N is the number of samples associated to the subject (this number is commonly the same as the number of time points).

% raw_time_points : This is the raw list of time points as they appear on the GSE record. In some GSE records these time points are not necessarily sorted and may have duplicates. This argument is this list of arguments as a row vector (i.e., as an horizontal vector).

% Output

% preprocessed_gene_expression : The preprocessed gene expression of the subject as a MxN matrix, where M is the number of genes and N is the number of time points (sorted and without duplicates).

% preprocessed_time_points : A column vector with the list of time points (sorted and without duplicates, possibly unlike the list of time points passed as argument raw_time_points).

function [preprocessed_gene_expression, preprocessed_time_points] = get_preprocessed_data(raw_gene_expression, raw_time_points)

  Subject = repmat(1,1,length(raw_time_points))';

  %Obtain the log 2 centered gene expression data for each subject.
  gexp = [];

  %Centered Gene Expression
  preprocessed_gene_expression = [];
  preprocessed_time_points  = [];


  %Get Subject
  ind = 1:length(Subject);

  %Get time
  [tmp,ix] = sort(raw_time_points(ind));

  %Check time > -1
  tmp = tmp(tmp>-1);

  %Check for replicates
  tb       = tabulate(tmp);
  preprocessed_time_points  = tb(tb(:,2)>0,1);
  nt       = length(preprocessed_time_points);

  if(max(tb(:,2))>1)
    for j = 1:nt
      %preprocessed_gene_expression time gene
      idx = find(tb(j,1)==tmp);

      if(any(raw_gene_expression(:,ind(ix(idx)))<0))
	gexp(:,j)      = nanmean((raw_gene_expression(:,ind(ix(idx)))),2);
      else
	gexp(:,j)      = nanmean(log2(raw_gene_expression(:,ind(ix(idx)))),2);
      end
    end
  else
    gexp  = raw_gene_expression(:,ind(ix));
  end
  preprocessed_gene_expression = gexp - repmat(nanmean(gexp,2),1,nt);
end

