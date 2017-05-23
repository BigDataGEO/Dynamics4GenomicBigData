% Input

% raw_gene_expression : This is the raw gene expression data of the subject as a MxN matrix, where M is the number of genes and N is the number of samples associated to the subject (this number is commonly the same as the number of time points).

% raw_time_points : This is the raw list of time points as they appear on the GSE record. In some GSE records these time points are not necessarily sorted and may have duplicates. This argument is this list of arguments as a row vector (i.e., as an horizontal vector).

% output : A boolean argument indicating whether the results should be output as files. If true, then results are output to directory 'step_2'.

% Output

% gene_expression : The preprocessed gene expression of the subject as a MxN matrix, where M is the number of genes and N is the number of time points (sorted and without duplicates).

% time_points : A column vector with the list of time points (sorted and without duplicates, possibly unlike the list of time points passed as argument raw_time_points).

% smooth_gene_trajectories : This is the expression of genes with smooth trajectories. Also a matrix analogous to gene_expression, although possibly with less rows.

function [gene_expression, standardized_gene_expression, time_points, gene_ID_type, smooth_gene_trajectories] = step_2_from_csv(input_csv_file, output)

  [raw_gene_expression, raw_time_points, gene_ID_type] = step_1_from_csv(input_csv_file);

  [gene_expression, time_points] = get_preprocessed_data(cell2mat(table2cell(raw_gene_expression(:,3:size(raw_gene_expression,2)))), raw_time_points);
  
  smooth_gene_trajectories = Est_Sub_Sel(time_points, gene_expression, 1);
  
  smooth_gene_trajectories = smooth_gene_trajectories{1};
  
  standardized_gene_expression = zscore(gene_expression')';
  
  A = time_points';
  A = A';
  A = strtrim(cellstr(num2str(A))');
  A = strcat('T_', A);
  
  gene_expression = cell2table([strtrim(table2cell(raw_gene_expression(:,1))) strtrim(table2cell(raw_gene_expression(:,2))) num2cell(gene_expression)], 'VariableNames', [{'Probe_ID'} {'Gene_ID'} A]);
  
  standardized_gene_expression = cell2table([strtrim(table2cell(raw_gene_expression(:,1))) strtrim(table2cell(raw_gene_expression(:,2))) num2cell(standardized_gene_expression)], 'VariableNames', [{'Probe_ID'} {'Gene_ID'} A]);
  
  if(output)
  
    global Dynamics4GenomicBigData_HOME;
  
    outputFolder = 'Step_2';
    mkdir(outputFolder);
  
    h=figure('units', 'centimeters', 'position', [0, 0, 30, 24]);

    clear title;

    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 30 24]);
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperSize', [30 24]);
    axisLabelFontSize = 30;
    
    surf(cell2mat(table2cell(gene_expression(:,3:size(gene_expression,2)))), 'FaceColor', 'interp', 'EdgeColor', 'none');

    xlim([1,length(time_points)]);

    set(gca,'XTick',1:length(time_points),'Xticklabel',time_points);
    set(gca,'FontSize',11);

    ylim([1,size(raw_gene_expression,1)]);

    zlim([min(min(cell2mat(table2cell(gene_expression(:,3:size(gene_expression,2)))))),max(max(cell2mat(table2cell(gene_expression(:,3:size(gene_expression,2))))))]);

    xlabel('Time', 'FontSize', axisLabelFontSize);

    ylabel('All genes', 'FontSize', axisLabelFontSize);

    zlabel('Expression', 'FontSize', axisLabelFontSize);

    title(['Expression of all genes'], 'FontSize', axisLabelFontSize);

    hold on;

    hold off;
    
    saveas(gcf, 'Paper_01.png');
    movefile('Paper_01.png', outputFolder);

    matrix_of_files_descs = [{'File name'} {'Description.'}];
    
    matrix_of_files_descs = [matrix_of_files_descs; [{'Paper_01.pdf'} {'Expression of all genes.'}]];
    
    matrix_of_files_descs = [matrix_of_files_descs; [{'gene_expression.csv'} {'Preprocessed gene expression in matrix form. Rows are probes and columns are time points.'}]];
    matrix_of_files_descs = [matrix_of_files_descs; [{'standardized_gene_expression.csv'} {'Standardized gene expression in matrix form. Rows are probes and columns are time points.'}]];
    matrix_of_files_descs = [matrix_of_files_descs; [{'time_points.csv'} {'Preprocessed time points.'}]];
    matrix_of_files_descs = [matrix_of_files_descs; [{'raw_gene_expression.csv'} {'Non-normalized gene expression in matrix form. Rows are probes and columns are time points.'}]];
    matrix_of_files_descs = [matrix_of_files_descs; [{'raw_time_points.csv'} {'Time points as they appear in series matrix.'}]];
    
    writetable(cell2table(matrix_of_files_descs), 'List_and_description_of_output.csv', 'WriteVariableNames', false);

    movefile('List_and_description_of_output.csv', outputFolder);
    
    cd(outputFolder);
    
    
    writetable(gene_expression, ['gene_expression.csv'], 'WriteVariableNames', true);
    
    writetable(standardized_gene_expression, ['standardized_gene_expression.csv'], 'WriteVariableNames', true);
    
    writetable(raw_gene_expression, ['raw_gene_expression.csv'], 'WriteVariableNames', true);
    
    writetable(cell2table(num2cell(raw_time_points)), ['raw_time_points.csv'], 'WriteVariableNames', false);
    
    writetable(cell2table(num2cell(time_points)), ['time_points.csv'], 'WriteVariableNames', false);
    
    cd('..');
    
    close all;
    
  end
  
end

% Input
% Parameter samples is an MxN cell array of strings where each element is the accession number of a
% GEO sample. In other words, each element must start with the prefix 'GSM'. M is the number of
% time points and N is the number of replicates that comprise the experimental condition. If N>1
% then the raw gene expression returned is the average of the N replicates. The expression of each
% gene, each time point is averaged across the N replicates. In other words, the N matrices are
% averaged, element-wise.

% Parameter preprocessed_time_points is a Mx1 cell array of strings where each element is a time point expressed
% with units, e.g., '2 hours'. The dimension of this array must be consistent with that of
% parameter samples.

% Output
function [raw_gene_expression, raw_time_points, gene_ID_type] = step_1_from_csv(input_csv_file)

  global Dynamics4GenomicBigData_HOME;

%    [geoStruct, list_of_genes, gene_ID_type, list_of_probe_ids] = get_geo_data(GEO_number);

  previous_folder = pwd;
  
  cd(Dynamics4GenomicBigData_HOME);

  cd('Input');
  
  data_as_table = readtable(input_csv_file);
  
  cd('..');
  
  cd(previous_folder);
  
  raw_gene_expression = cell2mat(table2cell(data_as_table(:,3:size(data_as_table,2))));
  list_of_probe_ids = table2cell(data_as_table(:,1));
  list_of_genes = table2cell(data_as_table(:,2));
  gene_ID_type = '';
  
  if isnumeric(list_of_genes{1})
    list_of_genes = cellfun(@num2str, list_of_genes, 'UniformOutput', false);
  end
  
  if isnumeric(list_of_probe_ids{1})
    list_of_probe_ids = cellfun(@num2str, list_of_probe_ids, 'UniformOutput', false);
  end
  
  time_labels = data_as_table.Properties.VariableNames(3:size(data_as_table,2))';
  
  raw_time_points = [];
  
  for i=1:length(time_labels)
    parts = strsplit(time_labels{i}, '_');
    raw_time_points = [raw_time_points; str2num(parts{2})];
  end  
  
  raw_gene_expression = cell2table([strtrim(list_of_probe_ids) strtrim(list_of_genes) num2cell(raw_gene_expression)], 'VariableNames', [{'Probe_ID'} {'Gene_ID'} strcat('T', cellfun(@num2str, num2cell(1:length(raw_time_points))', 'UniformOutput', false)')]);
  
end


function [smooth_gene_trajectories, gene_df, gene_ind] = Est_Sub_Sel(time_points, gene_expression, N, npool, plot_pool)

if nargin < 4, npool = 50;  end
if nargin < 5, plot_pool = 0;  end

%  ---------------Set up options for optimisation and cells  ---------------
options = optimset('LargeScale',  'off',  ...
                   'Display',     'off',   ...
                   'Diagnostics', 'off', ...
                   'GradObj',     'off',     ...
                   'Hessian',     'off',    ...               
                   'TolFun',      1e-8, ...
                   'TolX',        1e-8);

gexp3    = cell(N,1);
IQR    = cell(N,1);
SCR   = cell(N,1);
IND = cell(N,1);
lambda = cell(N,1);
df = cell(N,1);
gene_ind = cell(N,1);
gene_df = cell(N,1);
smooth_gene_trajectories = cell(N,1);
CROS  = cell(N,1);

for i = 1:N

    %  ---------------  set up the b-spline basis  ------------------------
    knots    = time_points';
    norder   = 3;
    nbasis   = length(time_points) + norder - 1;
    basisobj = create_bspline_basis([min(time_points) max(time_points)], nbasis, norder, knots);
    R        = eval_penalty(basisobj,2);
    
    sgenes = 1;
    %  --------------- Identify no. of crossings  ---------------
    CROS{i}  = zc(gene_expression);
    for k=1:4
        %  --------------- Obtain genes with k no of crossings  ---------------
        gexp3{i,k} = gene_expression(logical(CROS{i}==k),:);
        %  --------------- Get the IQR of the Genes with k no of corssings ---------------
        IQR{i}   = iqr(gexp3{i,k},2);
        %  --------------- Sort the Genes by IQR Largest-Smallest  ---------------
        [SCR{i},IND{i,k}] = sort(IQR{i},'descend');
        ngenes = length(IND{i,k});%min(50,length(IND{i}));
%          ngenes = min(50,length(IND{i,k}));
        pool_count = 0;
        for g=1:ngenes
            %  --------------- Obtain the lambda for the genes with k no of crossings  ---------------
            lambda{i,k,g} = fminbnd(@GCV_fun, 10.^-6, 10^6, options, eval_basis(time_points,basisobj), gexp3{i,k}(IND{i,k}(g),:)', R);
            %  --------------- Convert the lambda to the corresponding degrees of freedom  ---------------
            df{i,k,g}     = lambda2df(time_points, basisobj, ones(length(time_points),1),  2, lambda{i,k,g});
             %  --------------- If degrees of freedom are greater than no. of crossings  ---------------
              %  ---------------the gene is misclassified and does not belong to this crossing  ---------------
            if (df{i,k,g}>(k+1))
                gene_ind{i,k,pool_count+1} = IND{i,k}(g);
                gene_df{i,k,pool_count+1}  = df{i,k,g};
                pool_count = pool_count +1;
            else
                continue;
            end
            %  --------------- Only want at most npool genes in each crossing  ---------------
            if(pool_count == npool)
                break;
            end
        end
        %  --------------- Obtain the data for the sample  ---------------
        smooth_gene_trajectories{i}(sgenes:(sgenes+pool_count-1),:) = gexp3{i,k}(horzcat(gene_ind{i,k,:}),:);
        if(plot_pool == 1 && pool_count>0)
        figure(i);
        subplot(2,2,k)
        plot(smooth_gene_trajectories{i}(sgenes:(sgenes+pool_count-1),:)','--b')
        hold on;
        hline(0,'r');
        hold off;
        xlim([1,size(smooth_gene_trajectories{i}(sgenes:(sgenes+pool_count-1),:),2)])
        ylim([min(min(smooth_gene_trajectories{i}(sgenes:(sgenes+pool_count-1),:))),max(max(smooth_gene_trajectories{i}(sgenes:(sgenes+pool_count-1),:)))])
        title(['Subject ',num2str(i),' No. genes ',num2str(pool_count),' Max df ',num2str(max(horzcat(gene_df{i,k,:}))),' Min df ',num2str(min(horzcat(gene_df{i,k,:})))])
        end
        sgenes = sgenes + pool_count;
    end
    
end

end


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
  tb = tabulate(tmp);
  preprocessed_time_points  = tb(tb(:,2)>0,1);
  nt = length(preprocessed_time_points);

  if(max(tb(:,2))>1)
    for j = 1:nt
      %preprocessed_gene_expression time gene
      idx = find(tb(j,1)==tmp);

      if(any(raw_gene_expression(:,ind(ix(idx)))<0))
	gexp(:,j) = nanmean((raw_gene_expression(:,ind(ix(idx)))),2);
      else
	gexp(:,j) = nanmean(log2(raw_gene_expression(:,ind(ix(idx)))),2);
      end
    end
  else
    gexp = raw_gene_expression(:,ind(ix));
  end
  preprocessed_gene_expression = gexp - repmat(nanmean(gexp,2),1,nt);
end

% ZC number of zero crossings in x
% [n] = zc(x) calculates the number of zero crossings in x

function [n] = zc(x)

x = x(:,~isnan(x(1,:)));

n = zeros(size(x,1),1);
for i = 1:size(x,1)    
    s=sign(x(i,:));
    t=filter([1 1],1,s);
    n(i) = length(find(t==0));
end
end