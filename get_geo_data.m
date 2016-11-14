function [geoStruct, list_of_genes, gene_ID_type] = get_geo_data(GEO_number)

  cache_folder_name = 'GEO_cache';
  path_to_cached_file = [cache_folder_name '/' GEO_number '.txt'];
  
  try
    if(exist(cache_folder_name, 'dir') && exist(path_to_cached_file, 'file'))
      geoStruct = geoseriesread(path_to_cached_file);
    else
      mkdir(cache_folder_name);
      geoStruct = getgeodata2(GEO_number, 'ToFile', path_to_cached_file);
    end
    
    list_of_genes = [];
    gene_ID_type = [];
    
    path_to_cached_file_2 = [cache_folder_name '/' GEO_number '.mat'];
    
    if(exist(cache_folder_name, 'dir') && exist(path_to_cached_file_2, 'file'))
    
      load(path_to_cached_file_2, 'list_of_genes', 'gene_ID_type');
      
    else
      list_of_genes = get_list_of_gene_ids(geoStruct);
      gene_ID_type = capture_type_of_gene_ID(geoStruct.Header.Series.geo_accession);
      
      save(path_to_cached_file_2, 'list_of_genes', 'gene_ID_type');
      
    end
    
  catch causeException
    msgID = 'MATLAB:rmpath:DirNotFound1';
    msg = ['Unable to retrieve dataset ''' GEO_number ''' from the Gene Expression Omnibus.'];
    baseException = MException(msgID,msg);
    
    baseException = addCause(baseException,causeException);
    
    throw(baseException);
  end
end

function platform_struct = get_geo_platfom_data(platform_id)
  cache_folder_name = 'GEO_cache';
  path_to_cached_file = [cache_folder_name '/' platform_id '.txt'];
  
  if(exist(cache_folder_name, 'dir') && exist(path_to_cached_file, 'file'))
    platform_struct = geosoftread(path_to_cached_file);
  else
    mkdir(cache_folder_name);
    platform_struct = getgeodata2(platform_id, 'ToFile', path_to_cached_file);
  end
end

function list_of_genes = get_list_of_gene_ids(geoStruct)

  try
  
    row_identifiers_of_gse_matrix = rownames(geoStruct.Data);
    
    platform_id = geoStruct.Header.Series.platform_id;
    
    platform_struct = get_geo_platfom_data(platform_id);
    
    index_of_gpl_column_with_gene_ids = capture_index_of_gpl_column_with_gene_ids(platform_struct);
    
    fprintf('\n');
    display(['Loading gene IDs/names. This can take some time, please wait...']);
    fprintf('\n');
    
    list_of_genes = get_list_of_genes_from_gpl(platform_struct, index_of_gpl_column_with_gene_ids, row_identifiers_of_gse_matrix);
  
  catch causeException
    msgID = 'MATLAB:rmpath:DirNotFound2';
    msg = ['Unable to retrieve genes from GSE record in the Gene Expression Omnibus.'];
    baseException = MException(msgID,msg);
    
    baseException = addCause(baseException,causeException);
    
    throw(baseException);
  end
end

function list_of_genes = get_list_of_genes_from_gpl(platform_struct, index_of_gpl_column_with_gene_ids, row_identifiers_of_gse_matrix)

  row_identifiers_in_platform_record = [];
  if(isnumeric([platform_struct.Data{:, 1}]))
    
    numeric_values=[platform_struct.Data{:, 1}];
    
    for k=1:length(numeric_values)
      row_identifiers_in_platform_record = [row_identifiers_in_platform_record; {num2str(numeric_values(k))}];
    end
  else
    row_identifiers_in_platform_record = platform_struct.Data(:, 1);
  end

  % Under normal circumstances, variable not_found should be EMPTY after the following line.
  [indices, not_found] = find_in_cell_array_of_strings(row_identifiers_in_platform_record, row_identifiers_of_gse_matrix);
  
  if ~isempty(not_found)
    msgID = 'MATLAB:rmpath:DirNotFound3';
    msg = ['There are ' num2str(size(not_found,1)) ' genes that could not be mapped from GSE matrix to the GPL record.'];
    baseException = MException(msgID,msg);
    
    dlmwrite('UnmappedGenes.csv', not_found, '');
    
    throw(baseException);
  end

  list_of_genes = platform_struct.Data(indices, index_of_gpl_column_with_gene_ids);
end

% The parameters of this function are cell arrays of strings.
% The function returns the indices in the first array where members of the second array were found.
% Example
% Input
% cell_array_to_search = {'A', 'B', 'ABC', 'C', 'EFGZ', 'D', 'F'}
% cell_array_to_search_for = {'ABC', 'EFGZ', 'HOHO'}
% Output
% indices = [3; 5]
% not_found = ['HOHO']
function [indices, not_found] = find_in_cell_array_of_strings(cell_array_to_search, cell_array_to_search_for)
  indices = [];
  not_found = [];
  not_found_idx = [];
  for i=1:length(cell_array_to_search_for)

    idx = find(strcmp([cell_array_to_search(:)], strtrim(cell_array_to_search_for{i})));
      
    if(mod(i,1000)==0)
      display(['Scanned ' num2str(i) ' gene IDs/names out of ' num2str(length(cell_array_to_search_for)) '.']);
    end
      
    if(isempty(idx))
      not_found = [not_found; {cell_array_to_search_for{i}}];
      not_found_idx = [not_found_idx; i];
    else
      indices = [indices; idx];
    end
  end
  display(['Scanned ' num2str(i) ' gene IDs/names out of ' num2str(length(cell_array_to_search_for)) '.']);
end

function gene_ID_type = capture_type_of_gene_ID(GEO_number)
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
