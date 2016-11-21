% Input

% raw_gene_expression : This is the raw gene expression data of the subject as a MxN matrix, where M is the number of genes and N is the number of samples associated to the subject (this number is commonly the same as the number of time points).

% raw_time_points : This is the raw list of time points as they appear on the GSE record. In some GSE records these time points are not necessarily sorted and may have duplicates. This argument is this list of arguments as a row vector (i.e., as an horizontal vector).

% output : A boolean argument indicating whether the results should be output as files. If true, then results are output to directory 'step_2'.

% Output

% gene_expression : The preprocessed gene expression of the subject as a MxN matrix, where M is the number of genes and N is the number of time points (sorted and without duplicates).

% time_points : A column vector with the list of time points (sorted and without duplicates, possibly unlike the list of time points passed as argument raw_time_points).

% smooth_gene_trajectories : This is the expression of genes with smooth trajectories. Also a matrix analogous to gene_expression, although possibly with less rows.

function [gene_expression, time_points, smooth_gene_trajectories] = step_2(raw_gene_expression, raw_time_points, output)

  [gene_expression, time_points] = get_preprocessed_data(raw_gene_expression, raw_time_points);
  
  smooth_gene_trajectories = Est_Sub_Sel(time_points,gene_expression,1);
  smooth_gene_trajectories = smooth_gene_trajectories{1};
  
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
    
    surf(gene_expression,'FaceColor','interp','EdgeColor','none');

    xlim([1,length(time_points)]);

    set(gca,'XTick',1:length(time_points),'Xticklabel',time_points);
    set(gca,'FontSize',11);

    ylim([1,size(raw_gene_expression,1)]);

    zlim([min(min(gene_expression)),max(max(gene_expression))]);

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
    
    create_exel_file('List_and_description_of_output.xls', matrix_of_files_descs, 1, [], Dynamics4GenomicBigData_HOME);

    movefile('List_and_description_of_output.xls', outputFolder);
  end
  
end