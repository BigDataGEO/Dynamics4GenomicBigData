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