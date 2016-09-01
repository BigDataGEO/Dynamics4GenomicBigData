set_paths_and_imports;
  
[~,GEO_num] = xlsread('GEO_list.xlsx');

GEO_number = char(GEO_num(1));
Preprocessing_technique = 'Default';

while true

  [Data_GEO,list_of_genes,titles,Info,PInfo,geoStruct] = Obtain_data_from_GEO_website_user(GEO_number,Preprocessing_technique);

  [raw_gene_expression, raw_time_points, str_ind, pr_ind, tb, Subject_name, number_of_top_DRGs_considered, gene_ID_type] = capture_data(GEO_number, Data_GEO,list_of_genes,titles,Info,PInfo,geoStruct);
      
  [~, ~, con] = LCS(char(tb(pr_ind(1),1)),char(tb(pr_ind(end),1)));
  con = strrep(con,' ','_');
  con = strrep(con,'/','_');
  con = strrep(con,'.','_');
  flder = strcat(Dynamics4GenomicBigData_HOME,'Results/',GEO_number,'/',con);
      
  mkdir(flder);
  cd(flder);
  
  if ispc() % If running on Windows.
    options = struct('format','doc','showCode',false,'outputDir',flder,'stylesheet','document.xsl');
    publish('Paper.m',options);
  elseif isunix() % If running on Unix/Linux.
    options = struct('format','latex','showCode',false,'outputDir',flder,'stylesheet','document.xsl', 'imageFormat', 'png');
    publish('Paper.m',options);
    
    close all;
    
    copyfile([Dynamics4GenomicBigData_HOME, '/latex/bibliography.bib'], flder);
    copyfile([Dynamics4GenomicBigData_HOME, '/latex/plos2015.bst'], flder);
    
    % The following line compiles the .tex file into a .pdf.
    % Two output arguments (x and y) are used simply to prevent the output from being printed onscreen.
    [x, y]=system([Dynamics4GenomicBigData_HOME 'latex/compile.sh ' flder]);
  end

  cd(Dynamics4GenomicBigData_HOME);
      
  prompt = '\nWould you like to start a new analysis with a different set of samples (i.e., a different subject/condition)? ([1 "yes", 0 "no"]) ';
  continue_analysis   = input(prompt);
  
  if(continue_analysis)
    break;
  end
end
