function [gexp, gexp2, Time, N, n, Subject_n_ame] = get_preprocessed_data(Data, Subject, Pos, str_ind)

Subject_n_ame = str_ind;

%No. of Subjects
N = 1;

% No. of Probe sets
n = size(Data,1);

%Obtain the log 2 centered gene expression data for each subject.
gexp = cell(1,N);

%Centered Gene Expression
gexp2 = cell(1,N);
Time  = cell(1,N);

for i = 1:N
  %Get Subject
  ind = 1:length(Subject);

  %Get time
  [tmp,ix] = sort(Pos(ind));

  %Check time > -1
  tmp = tmp(tmp>-1);

  %Check for replicates
  tb       = tabulate(tmp);
  Time{i}  = tb(tb(:,2)>0,1);
  nt       = length(Time{i});

  if(max(tb(:,2))>1)
%      gexp_wr{i} = zeros(n,max(tb(:,2)),nt);

    for j = 1:nt

      %gexp2 time gene
      idx = find(tb(j,1)==tmp);

      if(any(Data(:,ind(ix(idx)))<0))
%  	gexp_wr{i}(:,1:tb(j,2),j) = (Data(:,ind(ix(idx))));
	gexp{i}(:,j)      = nanmean((Data(:,ind(ix(idx)))),2);
      else
%  	gexp_wr{i}(:,1:tb(j,2),j) = log2(Data(:,ind(ix(idx))));
	gexp{i}(:,j)      = nanmean(log2(Data(:,ind(ix(idx)))),2);
      end
    end
  else
    gexp{i}  = Data(:,ind(ix));
  end
  gexp2{i} = gexp{i} - repmat(nanmean(gexp{i},2),1,nt);
end

