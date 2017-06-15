% gene_expression is a matrix where the rows are the genes and the columns are the time points.

% gene_expression_plot(gene_expression, time_points, 'Gene expression', 'Time', 'Genes', 'Expression')

function gene_expression_plot(gene_expression, time_points, plot_title, x_label, y_label, z_label)

  mean_curve = mean(gene_expression, 1);

  h=figure('units', 'centimeters', 'position', [0, 0, 30, 24]);

  clear title;

  set(gcf, 'PaperPositionMode', 'manual');
  set(gcf, 'PaperUnits', 'centimeters');
  set(gcf, 'PaperPosition', [0 0 30 24]);
  set(gcf, 'PaperUnits', 'centimeters');
  set(gcf, 'PaperSize', [30 24]);
  axisLabelFontSize = 30;


%    surf(gene_expression,'FaceColor','interp','EdgeColor','none');
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

  title(strrep(plot_title, '_', '\_'), 'FontSize', axisLabelFontSize);

end


%  T = readtable('gene_expression_sorted_by_F_value.csv');

%  row_index = 1;
%  
%  gene_name = cell2mat(table2cell(T(row_index,3)));
%  
%  gene_expression = cell2mat(table2cell(T(row_index,5:size(T,2))));
%  
%  time_points = [0, 6, 12, 18, 24, 30, 36, 42];
%  
%  gene_expression_plot(gene_expression, time_points, ['Gene ' gene_name], 'Time', 'Expression', '');

%  row_index = 14783;
%  
%  gene_name = cell2mat(table2cell(T(row_index,3)));
%  
%  gene_expression = cell2mat(table2cell(T(row_index,5:size(T,2))));
%  
%  time_points = [0, 6, 12, 18, 24, 30, 36, 42];
%  
%  gene_expression_plot(gene_expression, time_points, ['Gene ' gene_name], 'Time', 'Expression', '');