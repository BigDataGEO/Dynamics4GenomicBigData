function [Data, Subject, Pos, str_ind, pr_ind, tb, Subject_name, number_of_top_DRGs, gene_ID_type] = capture_data(GEO_number, Data_GEO,gid,titles,Info,PInfo,geoStruct)

  %Ask which condition to analyze
  tb = tabulate(titles);
  dis = strcat(cellstr(arrayfun(@num2str, 1:length(tb(:,1)), 'UniformOutput', false))',' : ',tb(:,1));
  display(dis);    
	
  display(sprintf('\n\nAll the samples associated to %s are listed above. You must now specify the samples that you want to include in the analysis.\n', GEO_number));
  display(sprintf('You must do this by entering a term or string that is common ONLY to the desired samples from the list above.\n'));
  prompt = 'Which samples would you like to include in the analysis? (enter the common string) ';
  str_ind = input(prompt);
  pr_ind = ExtractSubjects(str_ind, dis);
  display(sprintf('\nThe samples that match your search terms are listed below.\n'));
  display(dis(pr_ind));

  %Get Data for that condtion
  Data     = Data_GEO(:,pr_ind);  
  
  display(sprintf('\nThe list below shows the information available for each one of the selected samples.\n\n'));
    
  %Display the Characteristics of the GEO series
  display(strcat(cellstr(arrayfun(@num2str, 1:length({Info{:,pr_ind(1)}}), 'UniformOutput', false))',' : ',{Info{:,pr_ind(1)}}'));

  %Find out where the time is
  prompt = 'Which item in the list corresponds to the time point? (Enter item number or -1 if none): ';
  tm_ind = input(prompt);
  
  if tm_ind > 0
    % Traditional case: time points can be read from characteristics matrix.
    format bank;
    Pos    = {Info{tm_ind,pr_ind}};
    if ~isempty(cell2mat(strfind(Pos,'Baseline')))
      Pos = strrep(Pos, 'Baseline', '0');
    end
    Pos    = cell2mat(cellfun(@(x) str2num(char(regexp(x,'\d+\.?\d*|-\d+\.?\d*|\.?\d*','match'))), Pos, 'UniformOutput', false));

    sane_check = 0;
    while sane_check == 0
      display(Pos');
      prompt = '\nThese are all the time values measured in hours. Are they correct? (Enter 1 for "Yes" or 0 for "No") ';
      sane_check = input(prompt);
      if sane_check ~= 0
	break;
      end
      timePointsMatrix = InputTimePointsManually();
	
      Pos = [];
      timeUnit = 'N/A';
	
      if not(isempty(timePointsMatrix))
	Pos = cell2mat(timePointsMatrix(:,1))';
	timeUnit = timePointsMatrix(1,2);
      end
    end

    %% Find out where the subject is
    prompt = '\nWhich item in the list corresponds to the subject/condition? (format [1,2,3] or all) ';
    su_ind = input(prompt);
    [Subject_name,~,Subject] = unique({Info{su_ind,pr_ind}});
    
    
    [number_of_top_DRGs] = capture_top_number_of_DRGs(GEO_number, size(Data, 1));
    
    [gene_ID_type] = capture_type_of_gene_ID(GEO_number);
    
  else
    % New case where time points must be read from title field or somewhere else.      
    timePointsMatrix = ExtractTimePoints(dis(pr_ind));
      
    Pos = [];
    timeUnit = 'N/A';
      
    if not(isempty(timePointsMatrix))
      Pos = cell2mat(timePointsMatrix(:,1))';
      timeUnit = timePointsMatrix(1,2);
    end
	
    sane_check = 0;
    while sane_check == 0
      display(Pos');
      prompt = ['These are all the time values measured in ' timeUnit{1} '. Are they correct? (Enter 1 for "Yes" or 0 for "No") '];
      sane_check = input(prompt);
      if sane_check ~= 0
	break;
      end
      timePointsMatrix = InputTimePointsManually();
	
      Pos = [];
      timeUnit = 'N/A';
	
      if not(isempty(timePointsMatrix))
	Pos = cell2mat(timePointsMatrix(:,1))';
	timeUnit = timePointsMatrix(1,2);
      end
    end

    if iscellstr(str_ind)
      Subject_name = strjoin(str_ind, '_');
    else
      Subject_name = str_ind;
    end
    Subject = repmat(1, 1, size(Pos,2));
    
    [number_of_top_DRGs] = capture_top_number_of_DRGs(GEO_number, size(Data, 1));
    
    [gene_ID_type] = capture_type_of_gene_ID(GEO_number);
  end
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
