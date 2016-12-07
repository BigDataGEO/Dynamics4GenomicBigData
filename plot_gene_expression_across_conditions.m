function plot_gene_expression_across_conditions(matrix_ids_of_probes_to_plot, gene_expression)

  h=figure('units', 'centimeters', 'position', [0, 0, 30, 24]);

  clear title;

  set(gcf, 'PaperPositionMode', 'manual');
  set(gcf, 'PaperUnits', 'centimeters');
  set(gcf, 'PaperPosition', [0 0 30 24]);
  set(gcf, 'PaperUnits', 'centimeters');
  set(gcf, 'PaperSize', [30 24]);
  axisLabelFontSize = 30;
  
  plot(gene_expression','-*b');
  
  hold on;

  plot(mean_curve,'o-r','LineWidth',1.5);
  
  xlim([1,length(time_points)]);

  set(gca,'XTick',1:length(time_points),'Xticklabel',time_points);
  set(gca,'FontSize',11);
  
  ylim([min(min(gene_expression))-.05,max(max(gene_expression))+.05]);

  zlim([min(min(gene_expression)),max(max(gene_expression))]);

  xlabel(x_label, 'FontSize', axisLabelFontSize);

  ylabel(y_label, 'FontSize', axisLabelFontSize);

  zlabel(z_label, 'FontSize', axisLabelFontSize);

  title(plot_title, 'FontSize', axisLabelFontSize);
  
  

  for replicate_index = 1:length(gene_expression)
  
    gene_expression_of_current_replicate = gene_expression{replicate_index};
    
    
  
  end

end