function [gexp, gexp2, Time, N, n, subject_name, yCR] = step_2(gene_expression_of_subject_across_time_points, Subject, Pos, str_ind, outputFig1, outputFig2)

  [gexp, gexp2, Time, N, n, subject_name] = get_preprocessed_data(gene_expression_of_subject_across_time_points, Subject, Pos, str_ind);
  
  if(outputFig1)
    h=figure('units', 'centimeters', 'position', [0, 0, 30, 24]);

    clear title;

    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 30 24]);
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperSize', [30 24]);
    axisLabelFontSize = 30;

    ind= 0;

    for sub = 1:N
    
	i = sub;

	surf(gexp2{i},'FaceColor','interp','EdgeColor','none');

	xlim([Time{sub}(1),length(Time{sub})]);

	set(gca,'XTick',1:length(Time{sub}),'Xticklabel',Time{sub});
	set(gca,'FontSize',11);

	ylim([1,n]);

	zlim([min(min(gexp2{i})),max(max(gexp2{i}))]);

	xlabel('Time', 'FontSize', axisLabelFontSize);

	ylabel('Genes', 'FontSize', axisLabelFontSize);

	zlabel('Expression', 'FontSize', axisLabelFontSize);

	title(['All genes'], 'FontSize', axisLabelFontSize);

	hold on;

	hold off;

    end

    savefig('Paper_01.fig');
    
    print('Paper_01.pdf','-dpdf');
  end
  
  
  yCR = Est_Sub_Sel(Time,gexp2,N);

  if(outputFig2)
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

	surf(yCR{sub},'FaceColor','interp','EdgeColor','none');

	xlim([Time{sub}(1),length(Time{sub})]);

	set(gca,'XTick',1:length(Time{sub}),'Xticklabel',Time{sub});
	set(gca,'FontSize',11);

	ylim([1,size(yCR{sub},1)]);

	zlim([min(min(yCR{sub})),max(max(yCR{sub}))]);

	xlabel('Time', 'FontSize', axisLabelFontSize);

	ylabel('Genes', 'FontSize', axisLabelFontSize);

	zlabel('Expression', 'FontSize', axisLabelFontSize);

	hold on;

	title(['Genes with smooth trajectories with dof<5'], 'FontSize', 20);

	hold off;

    end

    print('Paper_02.pdf','-dpdf');
  end