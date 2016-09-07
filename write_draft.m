function write_draft(GEO_number, list_of_genes, raw_gene_expression, raw_time_points, subject_name, condition, gene_ID_type, number_of_top_DRGs_considered)

  global Dynamics4GenomicBigData_HOME;
  
  output_folder = strcat(Dynamics4GenomicBigData_HOME,'Results/',GEO_number,'/',condition);
      
  mkdir(output_folder);
  cd(output_folder);
  
  if ispc() % If running on Windows.
    options = struct('format','doc','showCode',false,'outputDir',output_folder,'stylesheet','document.xsl');
    
    save('Paper.mat');
    
    publish('Paper.m',options);
    
    delete('Paper.mat');
  elseif isunix() % If running on Unix/Linux.
    options = struct('format','latex','showCode',false,'outputDir',output_folder,'stylesheet','document.xsl', 'imageFormat', 'png');

    save('Paper.mat');
    
    publish('Paper.m',options);
    
    close all;
    
    copyfile([Dynamics4GenomicBigData_HOME, '/latex/bibliography.bib'], output_folder);
    copyfile([Dynamics4GenomicBigData_HOME, '/latex/plos2015.bst'], output_folder);
    
    % The following line compiles the .tex file into a .pdf.
    % Two output arguments (x and y) are used simply to prevent the output from being printed onscreen.
    [x, y]=system([Dynamics4GenomicBigData_HOME 'latex/compile.sh ' output_folder]);
    
    delete('Paper.mat');
  end
  
  cd(Dynamics4GenomicBigData_HOME);
end