function [gexp2, Time, subject_name, yCR] = step_2(gene_expression_of_subject_across_time_points, Pos, str_ind, output)

  [gexp2, Time, subject_name] = get_preprocessed_data(gene_expression_of_subject_across_time_points, Pos, str_ind);
  
  yCR = Est_Sub_Sel(Time,gexp2,1);
  yCR = yCR{1};
  
  n = size(gene_expression_of_subject_across_time_points,1);
  
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

    ind= 0;
    
    N=1;

    for sub = 1:N
    
	i = sub;

	surf(gexp2,'FaceColor','interp','EdgeColor','none');

	xlim([Time(1),length(Time)]);

	set(gca,'XTick',1:length(Time),'Xticklabel',Time);
	set(gca,'FontSize',11);

	ylim([1,n]);

	zlim([min(min(gexp2)),max(max(gexp2))]);

	xlabel('Time', 'FontSize', axisLabelFontSize);

	ylabel('Genes', 'FontSize', axisLabelFontSize);

	zlabel('Expression', 'FontSize', axisLabelFontSize);

	title(['All genes'], 'FontSize', axisLabelFontSize);

	hold on;

	hold off;

    end

%      savefig('Paper_01.fig');
%      movefile('Paper_01.fig', outputFolder);
    
    print('Paper_01.pdf','-dpdf');
    movefile('Paper_01.pdf', outputFolder);

    h=figure('units', 'centimeters', 'position', [0, 0, 30, 24]);

    clear title;

    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 30 24]);
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperSize', [30 24]);
    axisLabelFontSize = 30;

    ind= 0 ;

    for sub = 1:N

	surf(yCR,'FaceColor','interp','EdgeColor','none');

	xlim([Time(1),length(Time)]);

	set(gca,'XTick',1:length(Time),'Xticklabel',Time);
	set(gca,'FontSize',11);

	ylim([1,size(yCR,1)]);

	zlim([min(min(yCR)),max(max(yCR))]);

	xlabel('Time', 'FontSize', axisLabelFontSize);

	ylabel('Genes', 'FontSize', axisLabelFontSize);

	zlabel('Expression', 'FontSize', axisLabelFontSize);

	hold on;

	title(['Genes with smooth trajectories with dof<5'], 'FontSize', 20);

	hold off;

    end

    print('Paper_02.pdf','-dpdf');
    movefile('Paper_02.pdf', outputFolder);

    matrix_of_files_descs = [{'File name'} {'Description.'}];
    
    matrix_of_files_descs = [matrix_of_files_descs; [{'Paper_01.pdf'} {'Expression of all genes.'}]];
    matrix_of_files_descs = [matrix_of_files_descs; [{'Paper_02.pdf'} {'Expression of genes with smooth trajectories with degrees of freedom below 5.'}]];
    
    create_exel_file('List_and_description_of_output.xls', matrix_of_files_descs, 1, [], Dynamics4GenomicBigData_HOME);

    movefile('List_and_description_of_output.xls', outputFolder);
  end
  
end