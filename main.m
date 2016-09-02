set_paths_and_imports;
  
[~,GEO_num] = xlsread('GEO_list.xlsx');

GEO_number = char(GEO_num(1));
Preprocessing_technique = 'Default';

while true
  
  [list_of_genes, raw_gene_expression, raw_time_points, subject_name, condition, gene_ID_type, number_of_top_DRGs_considered] = step_1(GEO_number);

  output_folder = strcat(Dynamics4GenomicBigData_HOME,'Results/',GEO_number,'/',condition);
      
  mkdir(output_folder);
  cd(output_folder);
  
  if ispc() % If running on Windows.
    options = struct('format','doc','showCode',false,'outputDir',output_folder,'stylesheet','document.xsl');
    publish('Paper.m',options);
  elseif isunix() % If running on Unix/Linux.
    options = struct('format','latex','showCode',false,'outputDir',output_folder,'stylesheet','document.xsl', 'imageFormat', 'png');
    publish('Paper.m',options);
    
    close all;
    
    copyfile([Dynamics4GenomicBigData_HOME, '/latex/bibliography.bib'], output_folder);
    copyfile([Dynamics4GenomicBigData_HOME, '/latex/plos2015.bst'], output_folder);
    
    % The following line compiles the .tex file into a .pdf.
    % Two output arguments (x and y) are used simply to prevent the output from being printed onscreen.
    [x, y]=system([Dynamics4GenomicBigData_HOME 'latex/compile.sh ' output_folder]);
  end

  cd(Dynamics4GenomicBigData_HOME);
      
  prompt = '\nWould you like to start a new analysis with a different set of samples (i.e., a different subject/condition)? ([1 "yes", 0 "no"]) ';
  continue_analysis   = input(prompt);
  
  if(~continue_analysis)
    break;
  end
end
