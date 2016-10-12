function [raw_gene_expression, raw_time_points, condition, subject_name, number_of_top_DRGs, gene_ID_type] = capture_data(geoStruct, raw_gene_expression_data_in_all_samples, titles, metadata_for_all_samples, PInfo)

  %Ask which condition to analyze
  tb = tabulate(titles);
  list_of_sample_record_titles = strcat(cellstr(arrayfun(@num2str, 1:length(tb(:,1)), 'UniformOutput', false))',' : ',tb(:,1));
  display(list_of_sample_record_titles);    
	
  display(sprintf('\n\nAll the samples associated to %s are listed above. You must now specify the samples that you want to include in the analysis.\n', geoStruct.Header.Series.geo_accession));
  display(sprintf('You must do this by entering a term or string that is common ONLY to the desired samples from the list above.\n'));
  prompt = 'Which samples would you like to include in the analysis? (enter the common string) ';
  common_string_identifying_subject_samples = input(prompt);
  indices_of_selected_sample_record_titles = ExtractSubjects(common_string_identifying_subject_samples, list_of_sample_record_titles);
  display(sprintf('\nThe samples that match your search terms are listed below.\n'));
  display(list_of_sample_record_titles(indices_of_selected_sample_record_titles));

  raw_gene_expression = raw_gene_expression_data_in_all_samples(:,indices_of_selected_sample_record_titles);  
  
  display(sprintf('\nThe list below shows the information available for each one of the selected samples.\n\n'));
    
  %Display the Characteristics of the GEO series
  display(strcat(cellstr(arrayfun(@num2str, 1:length({metadata_for_all_samples{:,indices_of_selected_sample_record_titles(1)}}), 'UniformOutput', false))',' : ',{metadata_for_all_samples{:,indices_of_selected_sample_record_titles(1)}}'));

  %Find out where the time is
  prompt = 'Which item in the list corresponds to the time point? (Enter item number or -1 if none): ';
  index_of_time_row_in_sample_metadata = input(prompt);
  
  if index_of_time_row_in_sample_metadata > 0
    % Traditional case: time points can be read from characteristics matrix.
    format bank;
    raw_time_points    = {metadata_for_all_samples{index_of_time_row_in_sample_metadata,indices_of_selected_sample_record_titles}};
    if ~isempty(cell2mat(strfind(raw_time_points,'Baseline')))
      raw_time_points = strrep(raw_time_points, 'Baseline', '0');
    end
    raw_time_points    = cell2mat(cellfun(@(x) str2num(char(regexp(x,'\d+\.?\d*|-\d+\.?\d*|\.?\d*','match'))), raw_time_points, 'UniformOutput', false));

    sane_check = 0;
    while sane_check == 0
      display(raw_time_points');
      prompt = '\nThese are all the time values measured in hours. Are they correct? (Enter 1 for "Yes" or 0 for "No") ';
      sane_check = input(prompt);
      if sane_check ~= 0
	break;
      end
      timePointsMatrix = InputTimePointsManually();
	
      raw_time_points = [];
      timeUnit = 'N/A';
	
      if not(isempty(timePointsMatrix))
	raw_time_points = cell2mat(timePointsMatrix(:,1))';
	timeUnit = timePointsMatrix(1,2);
      end
    end

    %% Find out where the subject is
%      prompt = '\nWhich item in the list corresponds to the subject/condition? (format [1,2,3] or all) ';
%      index_of_subject_in_sample_metadata = input(prompt);
%      [subject_name,~,Subject] = unique({metadata_for_all_samples{index_of_subject_in_sample_metadata,indices_of_selected_sample_record_titles}});
    
    if iscellstr(common_string_identifying_subject_samples)
      subject_name = strjoin(common_string_identifying_subject_samples, '_');
    else
      subject_name = common_string_identifying_subject_samples;
    end
    Subject = repmat(1, 1, size(raw_time_points,2));
    
    [number_of_top_DRGs] = capture_top_number_of_DRGs(geoStruct.Header.Series.geo_accession, size(raw_gene_expression, 1));
    
    [gene_ID_type] = capture_type_of_gene_ID(geoStruct.Header.Series.geo_accession);
    
  else
    % New case where time points must be read from title field or somewhere else.      
    timePointsMatrix = ExtractTimePoints(list_of_sample_record_titles(indices_of_selected_sample_record_titles));
      
    raw_time_points = [];
    timeUnit = 'N/A';
      
    if not(isempty(timePointsMatrix))
      raw_time_points = cell2mat(timePointsMatrix(:,1))';
      timeUnit = timePointsMatrix(1,2);
    end
	
    sane_check = 0;
    while sane_check == 0
      display(raw_time_points');
      prompt = ['These are all the time values measured in ' timeUnit{1} '. Are they correct? (Enter 1 for "Yes" or 0 for "No") '];
      sane_check = input(prompt);
      if sane_check ~= 0
	break;
      end
      timePointsMatrix = InputTimePointsManually();
	
      raw_time_points = [];
      timeUnit = 'N/A';
	
      if not(isempty(timePointsMatrix))
	raw_time_points = cell2mat(timePointsMatrix(:,1))';
	timeUnit = timePointsMatrix(1,2);
      end
    end

    if iscellstr(common_string_identifying_subject_samples)
      subject_name = strjoin(common_string_identifying_subject_samples, '_');
    else
      subject_name = common_string_identifying_subject_samples;
    end
    Subject = repmat(1, 1, size(raw_time_points,2));
    
    [number_of_top_DRGs] = capture_top_number_of_DRGs(geoStruct.Header.Series.geo_accession, size(raw_gene_expression, 1));
    
    [gene_ID_type] = capture_type_of_gene_ID(geoStruct.Header.Series.geo_accession);
  end
  
  [~, ~, condition] = LCS(char(tb(indices_of_selected_sample_record_titles(1),1)),char(tb(indices_of_selected_sample_record_titles(end),1)));
  condition = strrep(condition,' ','_');
  condition = strrep(condition,'/','_');
  condition = strrep(condition,'.','_');
    
end

function [number_of_top_DRGs] = capture_top_number_of_DRGs(GEO_number, total_number_of_genes)
    prompt = ['Dataset ' GEO_number ' contains a total of ' num2str(total_number_of_genes) ' genes.\n\nEnter the number of top DRGs you want to consider in the analysis (or -1 to include them all): '];
    number_of_top_DRGs = input(['\n' prompt]);
    if(~isnumeric(number_of_top_DRGs))
      number_of_top_DRGs = -1;
    end
end


function [gene_ID_type] = capture_type_of_gene_ID(GEO_number)
    idTypes = {
	      'AFFYMETRIX_3PRIME_IVT_ID'
	      'AFFYMETRIX_EXON_GENE_ID'
	      'AFFYMETRIX_SNP_ID'
	      'AGILENT_CHIP_ID'
	      'AGILENT_ID'
	      'AGILENT_OLIGO_ID'
	      'ENSEMBL_GENE_ID'
	      'ENSEMBL_TRANSCRIPT_ID'
	      'ENTREZ_GENE_ID'
	      'FLYBASE_GENE_ID'
	      'FLYBASE_TRANSCRIPT_ID'
	      'GENBANK_ACCESSION'
	      'GENOMIC_GI_ACCESSION'
	      'GENPEPT_ACCESSION'
	      'ILLUMINA_ID'
	      'IPI_ID'
	      'MGI_ID'
	      'PFAM_ID'
	      'PIR_ID'
	      'PROTEIN_GI_ACCESSION'
	      'REFSEQ_GENOMIC'
	      'REFSEQ_MRNA'
	      'REFSEQ_PROTEIN'
	      'REFSEQ_RNA'
	      'RGD_ID'
	      'SGD_ID'
	      'TAIR_ID'
	      'UCSC_GENE_ID'
	      'UNIGENE'
	      'UNIPROT_ACCESSION'
	      'UNIPROT_ID'
	      'UNIREF100_ID'
	      'WORMBASE_GENE_ID'
	      'WORMPEP_ID'
	      'ZFIN_ID'};
    

    fprintf('\n\n');
    
    for indx = 1:size(idTypes,1)
      display([num2str(indx) ': ' idTypes{indx} '']);
    end
    
    prompt = ['Enter the type of gene ID used in the study associated to ' GEO_number ' (e.g., enter 9 for ENSEMBL_GENE_ID): '];
    index_of_gene_ID_type = input(['\n\n' prompt]);
    
    gene_ID_type = idTypes{1};
    if(isnumeric(index_of_gene_ID_type) & index_of_gene_ID_type > 0 & index_of_gene_ID_type <= size(idTypes,1))
      gene_ID_type = idTypes{index_of_gene_ID_type};
    end
end
