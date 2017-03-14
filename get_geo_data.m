function [geoStruct, list_of_genes, gene_ID_type, list_of_probe_ids] = get_geo_data(GEO_number)

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
    list_of_probe_ids = [];
    
    path_to_cached_file_2 = [cache_folder_name '/' GEO_number '.mat'];
    
    if(exist(cache_folder_name, 'dir') && exist(path_to_cached_file_2, 'file'))
    
      load(path_to_cached_file_2, 'list_of_genes', 'gene_ID_type', 'list_of_probe_ids');
      
    else
      [list_of_genes, list_of_probe_ids] = get_list_of_gene_ids(geoStruct);
%        gene_ID_type = capture_type_of_gene_ID(geoStruct.Header.Series.geo_accession);
      
      save(path_to_cached_file_2, 'list_of_genes', 'gene_ID_type', 'list_of_probe_ids');
      
    end
    
    if isnumeric(list_of_genes{1})
      list_of_genes = strtrim(cellstr(num2str(cell2mat(list_of_genes))));
    end
    
    if isnumeric(list_of_probe_ids{1})
      list_of_probe_ids = strtrim(cellstr(num2str(cell2mat(list_of_probe_ids))));
    end
    
  catch causeException
    msgID = 'MATLAB:rmpath:DirNotFound1';
    msg = ['Unable to retrieve dataset ''' GEO_number ''' from the Gene Expression Omnibus.'];
    baseException = MException(msgID,msg);
    
    baseException = addCause(baseException,causeException);
    
    throw(baseException);
  end
end

function [list_of_genes, list_of_probe_ids] = get_list_of_gene_ids(geoStruct)

  try
  
    row_identifiers_of_gse_matrix = rownames(geoStruct.Data);
    
    list_of_probe_ids = row_identifiers_of_gse_matrix;
    
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

function index_of_gpl_column_with_gene_ids = capture_index_of_gpl_column_with_gene_ids(platform_struct)

    prompt = ['Check platform ' platform_struct.Accession '''s record on GEO and enter the index of the column where the gene IDs/names are located: '];
   
    index_of_gpl_column_with_gene_ids = input(['\n' prompt]);
    if(~isnumeric(index_of_gpl_column_with_gene_ids))
      index_of_gpl_column_with_gene_ids = 1;
    elseif(index_of_gpl_column_with_gene_ids < 1 | index_of_gpl_column_with_gene_ids > size(platform_struct.Data,2))
      index_of_gpl_column_with_gene_ids = 1;
    end
end


% This function attempts to find the gene names from the GPL record.

% This is by finding which rows in the GPL record have the same row identifiers as the GSE table
% (the latter are given as input).

% If all the row identifiers from the GSE table are matched as row identifiers in the GPL table,
% then the function returns the values of those rows at the column index provided as argument.

% If not all the row identifiers from the GSE table are matched as row identifiers in the GPL
% table, then the function returns the same list of row identifiers from the GSE table.

% The latter case can ocassionally occur with some platform records. The assumption (i.e., the GSE
% matrix row identifiers can be used as gene names) made by the function in these cases (i.e., the
% GSE matrix rows cannot be found in the GPL record) is necessary in order not to abort the flow of
% program and is usually correct.
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
  [indices, not_found] = find_matrix_rows_in_gpl_record(row_identifiers_in_platform_record, row_identifiers_of_gse_matrix);
  
  list_of_genes = [];
  
  if isempty(not_found)
    % This is the 'good' case. All identifiers from the GSE matrix could be mapped to the GPL
    % record, thus the list of genes obtained from the GPL record is reliable.
    list_of_genes = platform_struct.Data(indices, index_of_gpl_column_with_gene_ids);
  else
    % This is the 'bad' case. Not all identifiers from the GSE matrix could be mapped to the GPL
    % record, thus the list of genes obtained from the GPL record is not reliable.
    % In order not to abort the program, then simply assume that gene identifiers in the GSE matrix
    % can be used as gene names. And display a warning.
    msg = ['There were ' num2str(size(not_found,1)) ' genes that could not be mapped from GSE matrix to the GPL record.'];
    warning(msg);
    list_of_genes = row_identifiers_of_gse_matrix;
  end
end

function [indices, not_found] = find_matrix_rows_in_gpl_record(cell_array_to_search, cell_array_to_search_for)

  indices = [];
  not_found = [];
  not_found_idx = [];
  for i=1:length(cell_array_to_search_for)

    idx = find(strcmp([cell_array_to_search(:)], strtrim(cell_array_to_search_for{i})));
      
    if(mod(i,1000)==0)
      display(['Read ' num2str(i) ' gene IDs/names out of ' num2str(length(cell_array_to_search_for)) '.']);
    end
      
    if(isempty(idx))
      not_found = [not_found; {cell_array_to_search_for{i}}];
      not_found_idx = [not_found_idx; i];
    else
      indices = [indices; idx];
    end
  end
  display(['Read ' num2str(i) ' gene IDs/names out of ' num2str(length(cell_array_to_search_for)) '.']);
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

function [geoStruct,str] =getgeodata2(accessnum,varargin)
% GETGEODATA retrieves Gene Expression Omnibus (GEO) data. 
%
%   GEOSTRUCT = GETGEODATA(ACCESSNUM) searches for the accession number in
%   the Gene Expression Omnibus (GEO) database, and returns a structure
%   containing information for the object. 
%
%   GEOSTRUCT = GETGEODATA(...,'TOFILE',FILENAME) saves the data returned
%   from the database in the file FILENAME.
%
%   Note that currently Sample (GSM), DataSet (GDS), Series (GSE) and
%   Platform (GPL) records are supported. 
%
%   Example:
%   
%          % Get a sample file
%          geoSample = getgeodata('GSM1768')
%
%          % Get a data set and save it to a file
%          geoDataSet = getgeodata('GDS2602','tofile','gds2602.txt')
%
%          % Get a series data matrix
%          geoSeries = getgeodata('GSE11287')
%
%          % Get a Platform record
%          geoSeries = getgeodata('GPL74')
%
%   See http://www.ncbi.nlm.nih.gov/About/disclaimer.html for information
%   about using the GEO database.
%
%   See also GEOSERIESREAD, GEOSOFTREAD, GETGENBANK, GETGENPEPT.


% Copyright 2003-2012 The MathWorks, Inc.


bioinfochecknargin(nargin,1,mfilename);

if ~usejava('jvm')
    error(message('bioinfo:getgeodata:NeedJVM', mfilename));
end
tofile = false;


if nargin > 1
    if rem(nargin,2) == 0
        error(message('bioinfo:getgeodata:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'tofile',''};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname,okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:getgeodata:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:getgeodata:AmbiguousParameterName', pname));
        else
            switch(k)
                case  1  % tofile
                    if ischar(pval)
                        tofile = true;
                        filename = pval;
                    end
            end
        end
    end
end

% convert accessnum to a string if it is a number
if isnumeric(accessnum)
    accessnum = num2str(accessnum);
end

% error if accessnum isn't a string
if ~ischar(accessnum)
    error(message('bioinfo:getgeodata:NotString'))
end

% create the url that is used
% see
%    http://www.ncbi.nlm.nih.gov/entrez/query/static/linking.html
% for more information
% JCR: Changed URL to https, which is not included in MATLAB's native function.
searchurl = sprintf('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?form=text&acc=%s&view=full',accessnum);

% get the html file that is returned as a string
try
    str=urlread(searchurl);
catch theException
    error(message('bioinfo:getgeodata:URLProblem'));
end

% search for text indicating that there weren't any files found
notfound=strfind(str,'No items found');

% string was found, meaning no results were found
if ~isempty(notfound),
    error(message('bioinfo:getgeodata:BadAccessionNumber', accessnum)) ;
end

% GDS datasets
gdsBrowser = strfind(str,'gds_browse.cgi');
if ~isempty(gdsBrowser)
    % We expect the file on the FTP site to be called ACCNUM.soft.gz
    accessnum = upper(accessnum);
    softFileGZName = sprintf('%s.soft.gz',accessnum);
    softFileName = sprintf('%s.soft',accessnum);
    % Make a temporary directory for copying the file
    tempDirName = tempname;
    mkdir(tempDirName);
    % open FTP connection to NCBI and CD to GDS directory and get the file
    ftpConnection = ftp('ftp.ncbi.nih.gov');
    try
        cd(ftpConnection,'pub/geo/DATA/SOFT/GDS'); %#ok<MCCD>
        mget(ftpConnection,softFileGZName,tempDirName);
        close(ftpConnection);
    catch theErr
        % try to close the connection if there is an error
        close(ftpConnection);
        rethrow(theErr);
    end
    % gunzip the file and remove the .gz file
    localCopy = fullfile(tempDirName,softFileGZName);
    gunzip(localCopy);
    delete(localCopy);
    % Read the file
    geoStruct = geosoftread(fullfile(tempDirName,softFileName));
    if tofile
        copyfile(fullfile(tempDirName,softFileName),filename);
    end
    % Clean up
    delete(fullfile(tempDirName,softFileName));
    rmdir(tempDirName)
    return
end

%GSE series
gseFlag = strncmpi(str,'^SERIES',7);
if gseFlag ~= 0
     % We expect the file on the FTP site ftp://ftp.ncbi.nih.gov
     % in pub/geo/DATA/SeriesMatrix/%ACCNUM%
     % named ACCNUM_series_matrix.txt.gz
    accessnum = upper(accessnum);
    seriesFileGZName = sprintf('%s*_series_matrix.txt.gz',accessnum);
    seriesFileName = sprintf('%s_series_matrix.txt',accessnum);
    % Make a temporary directory for copying the file
    tempDirName = tempname;
    mkdir(tempDirName);
    % open FTP connection to NCBI and CD to GDS directory and get the file
    ftpConnection = ftp('ftp.ncbi.nih.gov');

	% JCR: The following lines set the 'pasive' mode that is required to allow downloads
	% through a firewall.
	if exist('pasv')
		pasv(ftpConnection);
	end
    try
        cd(ftpConnection,sprintf('pub/geo/DATA/SeriesMatrix/%s',accessnum)); %#ok<MCCD>
        mget(ftpConnection,seriesFileGZName,tempDirName);
        close(ftpConnection);
    catch theErr
        % try to close the connection if there is an error
        close(ftpConnection);
        eid = theErr.identifier;
        if strcmpi(eid,'MATLAB:ftp:NoSuchDirectory')
            error(message('bioinfo:getgeodata:NoSuchDirectory'));
        else
            rethrow(theErr)
        end
    end
    % gunzip the file and remove the .gz file
    localCopy = fullfile(tempDirName,seriesFileGZName);
    gunzip(localCopy);
    delete(localCopy);
    if ~exist(fullfile(tempDirName,seriesFileName),'file')
        dirInfo = dir([tempDirName filesep '*.txt']);
        if numel(dirInfo) > 1
            error(message('bioinfo:getgeodata:MultipleGSE', tempDirName));
        else
        ftpPath = sprintf('ftp.ncbi.nih.gov/pub/geo/DATA/SeriesMatrix/%s',accessnum);
        error(message('bioinfo:getgeodata:CannotExtractGSE', ftpPath));
        end
    end
    % Read the file
    geoStruct = geoseriesread(fullfile(tempDirName,seriesFileName));
    if tofile
        copyfile(fullfile(tempDirName,seriesFileName),filename);
    end
    % Clean up
    delete(fullfile(tempDirName,seriesFileName));
    rmdir(tempDirName)
    return
end

str = strrep(str,char(0),' ');
geoStruct = geosoftread(str);

%  write out file?
if tofile == true
    writefile = 'Yes';
    % check to see if file already exists
    if exist(filename,'file')
        % use dialog box to display options
        writefile=questdlg(sprintf('The file %s already exists. Do you want to overwrite it?',filename), ...
            '', ...
            'Yes','No','Yes');
    end

    switch writefile,
        case 'Yes',
            if exist(filename,'file')
                fprintf('File %s overwritten.',filename);
            end
            savedata(filename,str);
        case 'No',
            fprintf('File %s not written.',filename);
    end

end

end

function savedata(filename,str)

fid=fopen(filename,'w');

rows = size(str,1);

for rcount=1:rows-1
    fprintf(fid,'%s\n',str(rcount,:));
end
fprintf(fid,'%s',str(rows,:));

fclose(fid);


end
