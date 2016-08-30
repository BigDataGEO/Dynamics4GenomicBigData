function [gexp2, Time, Subject_n_ame] = get_preprocessed_data(Data, Pos, str_ind)

  Subject_n_ame = str_ind;

  Subject = repmat(1,1,length(Pos))';

  %Obtain the log 2 centered gene expression data for each subject.
  gexp = [];

  %Centered Gene Expression
  gexp2 = [];
  Time  = [];


  %Get Subject
  ind = 1:length(Subject);

  %Get time
  [tmp,ix] = sort(Pos(ind));

  %Check time > -1
  tmp = tmp(tmp>-1);

  %Check for replicates
  tb       = tabulate(tmp);
  Time  = tb(tb(:,2)>0,1);
  nt       = length(Time);

  if(max(tb(:,2))>1)
    for j = 1:nt
      %gexp2 time gene
      idx = find(tb(j,1)==tmp);

      if(any(Data(:,ind(ix(idx)))<0))
	gexp(:,j)      = nanmean((Data(:,ind(ix(idx)))),2);
      else
	gexp(:,j)      = nanmean(log2(Data(:,ind(ix(idx)))),2);
      end
    end
  else
    gexp  = Data(:,ind(ix));
  end
  gexp2 = gexp - repmat(nanmean(gexp,2),1,nt);
end

