function [chartReport, tableReport] = step_7(list_of_genes, list_of_gene_clusters, indices_of_top_DRGs, gene_ID_type)

  global Dynamics4GenomicBigData_HOME;

  currentFolder = pwd;
  cd(Dynamics4GenomicBigData_HOME);
%    py.importlib.import_module('DAVIDWS');
  
  cd(currentFolder);

  list_of_genes_to_annotate = list_of_genes(indices_of_top_DRGs(:));
  
  output_genes(list_of_genes, indices_of_top_DRGs, list_of_gene_clusters, 'Step_7', path);
  
  chartReport = [];
  tableReport = [];
  
%    [chartReport, tableReport] = annotate_genes(list_of_genes_to_annotate, gene_ID_type, true, true);
  
  annotate_genes_and_output_reports(chartReport, tableReport, gene_ID_type, list_of_genes, indices_of_top_DRGs, list_of_gene_clusters, 'Step_7', Dynamics4GenomicBigData_HOME, true, true);
end

function [chartReport, tableReport] = annotate_genes(list_of_genes_to_annotate, gene_ID_type, includeChartReport, includeTableReport)

  inputIds = py.str(strjoin(list_of_genes_to_annotate',', '));
  gene_ID_type = py.str(gene_ID_type);
  listName = py.str('Annotation');
  listType = py.int(0);
  thd = py.float(0.05);
  count = py.int(2);

  speciesPy = py.DAVIDWS.getSpecies(inputIds, gene_ID_type, listName, listType);

  Species = char(speciesPy);

  chartReport = [];
  tableReport = [];
  header = [];

  if not(strcmp(Species, 'None'))
    % The following four lines are in preparation for the chart report.
    if(includeChartReport)
      chartReportPy = py.DAVIDWS.getChartReport(inputIds, gene_ID_type, listName, listType, thd, count);
      chartReport = convertChartReportFromPythonRecordToMatlabMatrix(chartReportPy);
    end

    % The following lines (down to the end of the while loop) are in preparation for the table report.
    % These lines do the query to the DAVID WS.
    % DAVID is limited in the number of genes it can handle in each single call, thus the gene
    % list is sent in (relatively) small chunks.
    if(includeTableReport)
      % This is the maximum number of gene IDs sent in each call to DAVID. If the number of genes
      % to annotate is larger than this, then several calls will be made inside the while loop
      % that follows.
      maximum_number_of_genes_per_ws_call = 1000;
      number_of_genes_processed = 0;
      display(['A total of ' num2str(number_of_genes_processed) ' genes (out of ' num2str(length(list_of_genes_to_annotate)) ') have been processed so far.']);
      while number_of_genes_processed < length(list_of_genes_to_annotate)
	
	start_index = number_of_genes_processed + 1;
	end_index = length(list_of_genes_to_annotate);
	
	if(end_index - start_index > maximum_number_of_genes_per_ws_call)
	  end_index = start_index + maximum_number_of_genes_per_ws_call-1;
	end
	
	number_of_genes_processed = number_of_genes_processed + (end_index - start_index + 1);
	
	smaller_set_of_ids = list_of_genes_to_annotate(start_index:end_index);
	smaller_set_of_ids = py.str(strjoin(smaller_set_of_ids',', '));
	
	tableReportPy = py.DAVIDWS.getTableReport(smaller_set_of_ids, gene_ID_type, listName, listType);
	smalltableReport = convertTableReportFromPythonRecordToMatlabMatrix(tableReportPy);
	
	tableReport = [tableReport; smalltableReport];
	
	% If this is not the last iteration (i.e., if the condition of the while loop still holds and
	% hence DAVID will be called again) then we must wait ten seconds. This is required by DAVID.
	if number_of_genes_processed < length(list_of_genes_to_annotate)
	  pause(10);
	end
	display(['A total of ' num2str(number_of_genes_processed) ' genes (out of ' num2str(length(list_of_genes_to_annotate)) ') have been processed so far.']);
      end
    end
  end
end

function[table_matrix] = convertTableReportFromPythonRecordToMatlabMatrix(tableReportAsPythonList)

  fileFieldSeparator = ',';
  withinFieldSeparator = ';';

  [uselessVariable numberOfTableEntries] = size(tableReportAsPythonList);

  table_matrix = [];

  for i=1:numberOfTableEntries
    tableRecord=tableReportAsPythonList{i};
    
    % The following lines are only to extract the gene ID.
    cP = cell(tableRecord.values);
    cellP = cell(1, numel(cP));
    for n=1:numel(cP)
	strP = char(cP{n});
	cellP(n) = {strP};
    end
    substring_limits = strfind(strP, '"');
    geneID = strP(substring_limits(1)+1:substring_limits(2)-1); % Voila.
    
    record_matrix = [];
    
    [uselessVariable numberOfAnnotationRecords] = size(tableRecord.annotationRecords);

    for j=1:numberOfAnnotationRecords
      annotationRecord = tableRecord.annotationRecords{j};
      category = char(annotationRecord.category);
      all_terms_in_category = '';
      
      [uselessVariable numberOfTerms] = size(annotationRecord.terms);    
      for k=1:numberOfTerms
	term = char(annotationRecord.terms{k});
	term = term(strfind(term, '$')+1:length(term));
	
	if k == 1
	  all_terms_in_category = term;
	else
	  all_terms_in_category = strcat(all_terms_in_category, withinFieldSeparator, ' ', term);
	end
      end
      record_matrix = [record_matrix; [{category} {all_terms_in_category}]];
    end
    table_matrix = [table_matrix; [{geneID} {char(tableRecord.name)} {char(tableRecord.species)} {record_matrix}]];
  end
end


function[matlabMatrix] = convertChartReportFromPythonRecordToMatlabMatrix(chartReportAsPythonList)

  matlabMatrix = [];

  for i=1:length(chartReportAsPythonList)
    pythonRecord = chartReportAsPythonList{i};
    newRow = cellstr([num2str(pythonRecord.EASEBonferroni); num2str(pythonRecord.afdr); num2str(pythonRecord.benjamini); num2str(pythonRecord.bonferroni); {char(pythonRecord.categoryName)}; num2str(pythonRecord.ease); num2str(pythonRecord.fisher); num2str(pythonRecord.foldEnrichment); {char(pythonRecord.geneIds)}; num2str(pythonRecord.id); num2str(pythonRecord.listHits); {char(pythonRecord.listName)}; num2str(pythonRecord.listTotals); num2str(pythonRecord.percent); num2str(pythonRecord.popHits); num2str(pythonRecord.popTotals); num2str(pythonRecord.rfdr); {''}; {char(pythonRecord.termName)}])';
    matlabMatrix = [matlabMatrix; newRow];
  end
  
  labels = {'EASEBonferroni', 'afdr', 'benjamini', 'bonferroni', 'categoryName', 'ease', 'fisher', 'foldEnrichment', 'geneIds', 'id', 'listHits', 'listName', 'listTotals', 'percent', 'popHits', 'popTotals', 'rfdr', 'scores', 'termName'};
  
  matlabMatrix = [labels; matlabMatrix];
end

function write_gene_cluster_to_csv_file(ids_of_genes_in_current_cluster,output_file_name)

  fid = fopen(output_file_name,'w');
  for row = 1:size(ids_of_genes_in_current_cluster,1)
    ids_of_genes_in_current_cluster = ids_of_genes_in_current_cluster;
    fprintf(fid, repmat('%s\t',1,size(ids_of_genes_in_current_cluster,2)-1), ids_of_genes_in_current_cluster{row,1:end-1});
    fprintf(fid, '%s\n', ids_of_genes_in_current_cluster{row,end});
  end
  fclose(fid);
end

function annotate_genes_and_output_reports(chartReport, tableReport, gene_ID_type, list_of_genes, indices_of_top_DRGs, list_of_gene_clusters, output_dir, path, includeChartReport, includeTableReport)

  list_of_genes_to_annotate = list_of_genes(indices_of_top_DRGs(:));
  
  % Now that the chart and table reports have been obtained, we proceed to export them to files,
  % grouped by gene cluster.
  % The output directory is created in the following two lines.
  mkdir(output_dir);
  cd(output_dir);
  write_gene_cluster_to_csv_file(list_of_genes_to_annotate,strcat('All_DRGs', '.txt'));
  for cluster_iteration_ID = 1:length(list_of_gene_clusters)
    % First create the subfolder for the current cluster.
    mkdir(strcat('M',num2str(cluster_iteration_ID)));
    cd(strcat('M',num2str(cluster_iteration_ID)));
    
    % Then output all genes in the current cluster, for reference.
    ids_of_genes_in_current_cluster = list_of_genes(indices_of_top_DRGs(list_of_gene_clusters{cluster_iteration_ID}));
    write_gene_cluster_to_csv_file(ids_of_genes_in_current_cluster,strcat('Genes_in_M', num2str(cluster_iteration_ID), '.txt'));
    
    header = [];
    if includeChartReport && ~isempty(chartReport)
      header = chartReport(1,:);
    end
    
    chart_report_of_current_cluster = header;
    table_report_of_current_cluster = [];

    genes_grouped_by_pathway = [];
    
    table_reports_grouped_by_pathway = [];
    
    for gene_iteration_ID = 1:length(ids_of_genes_in_current_cluster)
      geneID = ids_of_genes_in_current_cluster{gene_iteration_ID};
      
      % This section finds the gene in the big table report.
      if includeTableReport
	if not(isempty(tableReport))
	  index_of_current_gene_in_table_report = find(not(cellfun('isempty', strfind(tableReport(:,1), geneID))));
	  if(not(isempty(index_of_current_gene_in_table_report)))
	    table_report_of_current_cluster = [table_report_of_current_cluster; tableReport(index_of_current_gene_in_table_report,:)];
	    
	    record_of_current_gene_in_table_report = tableReport(index_of_current_gene_in_table_report,:);
	    
	    fields_and_values_of_current_gene_in_table_report = record_of_current_gene_in_table_report{4};
	    
	    index_of_KEGG_field_for_current_gene_in_TR = find(strcmp([fields_and_values_of_current_gene_in_table_report(:,1)], 'KEGG_PATHWAY'));
	    
	    pathway_names = {'NO_PATHWAY:NO_PATHWAY'};
	    if(not(isempty(index_of_KEGG_field_for_current_gene_in_TR)))%i.e., this gene has a pathway associated to it.
	      pathway_names=fields_and_values_of_current_gene_in_table_report(index_of_KEGG_field_for_current_gene_in_TR,2);
	      pathway_names = strsplit(pathway_names{1}, ';');
	    end
	    
	    for pathway_index=1:length(pathway_names)
	      pathway_name = pathway_names{pathway_index};
	      
	      pathway_name_parts = strsplit(pathway_name, ':');
	      pathway_name = pathway_name_parts(1);
	      pathway_description = pathway_name_parts(2);
	      
	      index_of_pathway_record = [];
	      if(size(genes_grouped_by_pathway,1)>0)
		index_of_pathway_record = find(strcmp([genes_grouped_by_pathway(:,1)], pathway_name));
	      end
		
	      if(isempty(index_of_pathway_record))%There is no record with this pathway
		new_pathway = [pathway_name pathway_description {geneID}];
		genes_grouped_by_pathway = [genes_grouped_by_pathway; new_pathway];
	      else%There is already a record with this pathway
		genes_in_this_pathway = genes_grouped_by_pathway(index_of_pathway_record, 3);
		genes_in_this_pathway = [genes_in_this_pathway{:}; {geneID}];
		genes_grouped_by_pathway{index_of_pathway_record, 3} = genes_in_this_pathway;
	      end

	      idx = [];
	      if(size(table_reports_grouped_by_pathway,1)>0)
		idx = find(strcmp([table_reports_grouped_by_pathway(:,1)], pathway_name));
	      end
	      
	      if(isempty(idx))
		new_pathway = [pathway_name pathway_description {record_of_current_gene_in_table_report}];
		table_reports_grouped_by_pathway = [table_reports_grouped_by_pathway; new_pathway];
	      else
		table_reports_in_this_pathway = table_reports_grouped_by_pathway(idx, 3);
		table_reports_in_this_pathway = [table_reports_in_this_pathway{:}; record_of_current_gene_in_table_report];
		table_reports_grouped_by_pathway{idx, 3} = table_reports_in_this_pathway;
	      end
	      
	    end%for pathway_index=1:length(pathway_names)
	  end%if(not(isempty(index_of_current_gene_in_table_report)))
	end%if not(isempty(tableReport))
      end%if includeTableReport
      
      % This section finds the current gene in the big chart report.
      if includeChartReport && ~isempty(chartReport)
	[number_of_rows_in_chart_report annotation_ncolumns] = size(chartReport);
	for current_chart_report_row = 2:number_of_rows_in_chart_report
	  if(~isempty(strfind(lower(chartReport{current_chart_report_row,9}), lower(geneID))))
	    chart_report_of_current_cluster = [chart_report_of_current_cluster; chartReport(current_chart_report_row,:)];
	  end
	end
      end
      
    end%for gene_iteration_ID = 1:length(ids_of_genes_in_current_cluster)
      
    if includeChartReport && ~isempty(chartReport)
    
      chart_report_of_current_cluster = [chart_report_of_current_cluster(:,1:5) chart_report_of_current_cluster(:,19) chart_report_of_current_cluster(:,6:18)];
    
      [number_of_entries_in_chart_report_of_current_cluster uselessVariable] = size(chart_report_of_current_cluster);
      if number_of_entries_in_chart_report_of_current_cluster > 1
  
	fileFieldSeparator=',';
	withinFieldSeparator = ';';
	
	fid = fopen(strcat('Chart_report_of_M', num2str(cluster_iteration_ID),'.csv'), 'w');
		for row=1:number_of_entries_in_chart_report_of_current_cluster
	  fprintf(fid, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', strrep(chart_report_of_current_cluster{row, 1}, fileFieldSeparator, withinFieldSeparator),strrep(chart_report_of_current_cluster{row, 2}, fileFieldSeparator, withinFieldSeparator),strrep(chart_report_of_current_cluster{row, 3}, fileFieldSeparator, withinFieldSeparator),strrep(chart_report_of_current_cluster{row, 4}, fileFieldSeparator, withinFieldSeparator),strrep(chart_report_of_current_cluster{row, 5}, fileFieldSeparator, withinFieldSeparator),strrep(chart_report_of_current_cluster{row, 6}, fileFieldSeparator, withinFieldSeparator),strrep(chart_report_of_current_cluster{row, 7}, fileFieldSeparator, withinFieldSeparator),strrep(chart_report_of_current_cluster{row, 8}, fileFieldSeparator, withinFieldSeparator),strrep(chart_report_of_current_cluster{row, 9}, fileFieldSeparator, withinFieldSeparator),strrep(chart_report_of_current_cluster{row, 10}, fileFieldSeparator, withinFieldSeparator),strrep(chart_report_of_current_cluster{row, 11}, fileFieldSeparator, withinFieldSeparator),strrep(chart_report_of_current_cluster{row, 12}, fileFieldSeparator, withinFieldSeparator),strrep(chart_report_of_current_cluster{row, 13}, fileFieldSeparator, withinFieldSeparator),strrep(chart_report_of_current_cluster{row, 14}, fileFieldSeparator, withinFieldSeparator),strrep(chart_report_of_current_cluster{row, 15}, fileFieldSeparator, withinFieldSeparator),strrep(chart_report_of_current_cluster{row, 16}, fileFieldSeparator, withinFieldSeparator),strrep(chart_report_of_current_cluster{row, 17}, fileFieldSeparator, withinFieldSeparator),strrep(chart_report_of_current_cluster{row, 18}, fileFieldSeparator, withinFieldSeparator),strrep(chart_report_of_current_cluster{row, 19}, fileFieldSeparator, withinFieldSeparator));
	end
	
	fclose(fid);
      end
    end
      
    if includeTableReport && ~isempty(tableReport)
      [number_of_entries_in_table_report_of_current_cluster uselessVariable] = size(table_report_of_current_cluster);
      if number_of_entries_in_table_report_of_current_cluster > 0
	fileID = fopen(strcat('Table_report_of_M', num2str(cluster_iteration_ID),'.csv'),'w');
	fileFieldSeparator=',';
	withinFieldSeparator = ';';
	for index_in_table_report_of_current_cluster=1:number_of_entries_in_table_report_of_current_cluster
	  geneID = table_report_of_current_cluster{index_in_table_report_of_current_cluster,1};
	  datum2 = table_report_of_current_cluster{index_in_table_report_of_current_cluster,2};
	  datum3 = table_report_of_current_cluster{index_in_table_report_of_current_cluster,3};
	  fprintf(fileID, strcat('%s', fileFieldSeparator,' %s', withinFieldSeparator,' %s\n'), geneID, strrep(datum2, fileFieldSeparator, withinFieldSeparator), strrep(datum3, fileFieldSeparator, withinFieldSeparator));
	    
	  records = table_report_of_current_cluster(index_in_table_report_of_current_cluster,4);
	  records = records{1};
	    
	  [number_of_records uselessVariable] = size(records);
	    
	  for index_in_records=1:number_of_records
	    fprintf(fileID, strcat('%s', fileFieldSeparator, '%s\n'), strrep(records{index_in_records,1}, fileFieldSeparator, withinFieldSeparator), strrep(records{index_in_records,2}, fileFieldSeparator, withinFieldSeparator));
	  end
	  fprintf(fileID, strcat(fileFieldSeparator,fileFieldSeparator,fileFieldSeparator,fileFieldSeparator,fileFieldSeparator,'\n'));
	end
	fclose(fileID);
      end%if number_of_entries_in_table_report_of_current_cluster > 0
      
      if size(genes_grouped_by_pathway,1) > 0
	
	for index_of_pathway_record = 1:size(genes_grouped_by_pathway,1)
	  pathway_name = genes_grouped_by_pathway{index_of_pathway_record, 1};
	  pathway_description = genes_grouped_by_pathway{index_of_pathway_record, 2};
	  genes_in_this_pathway = genes_grouped_by_pathway{index_of_pathway_record, 3};
	  
	  table_report_in_this_pathway = [];%equivalent to table_report_of_current_cluster
	  if(size(table_reports_grouped_by_pathway,1)>0)
	    idx = find(strcmp([table_reports_grouped_by_pathway(:,1)], pathway_name));
	    if(not(isempty(idx)))
	      table_report_in_this_pathway = table_reports_grouped_by_pathway{idx, 3};
	    end
	  end
	  
	  mkdir(pathway_name);
	  cd(pathway_name);
	  
	  if(ischar(genes_in_this_pathway))
	    write_gene_cluster_to_csv_file({genes_in_this_pathway},strcat('Genes_in_', pathway_name, '.txt'));
	  else
	    write_gene_cluster_to_csv_file(genes_in_this_pathway,strcat('Genes_in_', pathway_name, '.txt'));
	  end

	  [number_of_entries_in_table_report_in_this_pathway uselessVariable] = size(table_report_in_this_pathway);
	  if number_of_entries_in_table_report_in_this_pathway > 0
	    fileID = fopen(strcat('Table_report_of_', pathway_name,'.csv'),'w');
	    fileFieldSeparator=',';
	    withinFieldSeparator = ';';
	    for index_in_table_report_in_this_pathway=1:number_of_entries_in_table_report_in_this_pathway
	      geneID = table_report_in_this_pathway{index_in_table_report_in_this_pathway,1};
	      datum2 = table_report_in_this_pathway{index_in_table_report_in_this_pathway,2};
	      datum3 = table_report_in_this_pathway{index_in_table_report_in_this_pathway,3};
	      fprintf(fileID, strcat('%s', fileFieldSeparator,' %s', withinFieldSeparator,' %s\n'), geneID, strrep(datum2, fileFieldSeparator, withinFieldSeparator), strrep(datum3, fileFieldSeparator, withinFieldSeparator));
		
	      records = table_report_in_this_pathway(index_in_table_report_in_this_pathway,4);
	      records = records{1};
		
	      [number_of_records uselessVariable] = size(records);
		
	      for index_in_records=1:number_of_records
		fprintf(fileID, strcat('%s', fileFieldSeparator, '%s\n'), strrep(records{index_in_records,1}, fileFieldSeparator, withinFieldSeparator), strrep(records{index_in_records,2}, fileFieldSeparator, withinFieldSeparator));
	      end
	      fprintf(fileID, strcat(fileFieldSeparator,fileFieldSeparator,fileFieldSeparator,fileFieldSeparator,fileFieldSeparator,'\n'));
	    end
	    fclose(fileID);
	  end%if number_of_entries_in_table_report_in_this_pathway > 0
	  cd('..');
	end%for index_of_pathway_record = 1:size(genes_grouped_by_pathway,1)
      end%if size(genes_grouped_by_pathway,1) > 0
    end
    cd('..');
  end
  cd('..');
end



function output_genes(list_of_genes, indices_of_top_DRGs, list_of_gene_clusters, output_dir, path)

  list_of_genes_to_annotate = list_of_genes(indices_of_top_DRGs(:));
  
  % Now that the chart and table reports have been obtained, we proceed to export them to files,
  % grouped by gene cluster.
  % The output directory is created in the following two lines.
  mkdir(output_dir);
  cd(output_dir);
  write_gene_cluster_to_csv_file(list_of_genes_to_annotate,strcat('All_DRGs', '.txt'));
  for cluster_iteration_ID = 1:length(list_of_gene_clusters)
    % First create the subfolder for the current cluster.
    mkdir(strcat('M',num2str(cluster_iteration_ID)));
    cd(strcat('M',num2str(cluster_iteration_ID)));
    
    % Then output all genes in the current cluster, for reference.
    ids_of_genes_in_current_cluster = list_of_genes(indices_of_top_DRGs(list_of_gene_clusters{cluster_iteration_ID}));
    write_gene_cluster_to_csv_file(ids_of_genes_in_current_cluster,strcat('Genes_in_M', num2str(cluster_iteration_ID), '.txt'));
      
    cd('..');
  end
  cd('..');
end
