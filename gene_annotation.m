function[tableReport,chartReport,ClusterReport] = gene_annotation(inputIds,idType)

if nargin<2
    idType = 'AFFYMETRIX_3PRIME_IVT_ID';
end

inputIds = py.str(strjoin(inputIds',', '));
idType = py.str(idType);
listName = py.str('Annotation');
listType = py.int(0);
thd = py.float(0.05);
count = py.int(2);

speciesPy = py.DAVIDWS.getSpecies(inputIds, idType, listName, listType);

Species = char(speciesPy);

chartReport = [];

if not(strcmp(Species, 'None'))
  chartReportPy = py.DAVIDWS.getChartReport(inputIds, idType, listName, listType, thd, count);

  %  chartReport = char(chartReportPy)
  chartReport = convertPythonList(chartReportPy);
end

SummaryReport = [];
tableReport = [];
ClusterReport = [];


function[matlabMatrix] = convertPythonList(pythonList)

  matlabMatrix = [];

  for i=1:length(pythonList)
    pythonRecord = pythonList{i};
    newRow = cellstr([num2str(pythonRecord.EASEBonferroni); num2str(pythonRecord.afdr); num2str(pythonRecord.benjamini); num2str(pythonRecord.bonferroni); {char(pythonRecord.categoryName)}; num2str(pythonRecord.ease); num2str(pythonRecord.fisher); num2str(pythonRecord.foldEnrichment); {char(pythonRecord.geneIds)}; num2str(pythonRecord.id); num2str(pythonRecord.listHits); {char(pythonRecord.listName)}; num2str(pythonRecord.listTotals); num2str(pythonRecord.percent); num2str(pythonRecord.popHits); num2str(pythonRecord.popTotals); num2str(pythonRecord.rfdr); {''}; {char(pythonRecord.termName)}])';
    matlabMatrix = [matlabMatrix; newRow];
  end
  
%    matlabMatrix = array2table(matlabMatrix,...
%      'VariableNames',{'EASEBonferroni', 'afdr', 'benjamini', 'bonferroni', 'categoryName', 'ease', 'fisher', 'foldEnrichment', 'geneIds', 'id', 'listHits', 'listName', 'listTotals', 'percent', 'popHits', 'popTotals', 'rfdr', 'scores', 'termName'});
  
  labels = {'EASEBonferroni', 'afdr', 'benjamini', 'bonferroni', 'categoryName', 'ease', 'fisher', 'foldEnrichment', 'geneIds', 'id', 'listHits', 'listName', 'listTotals', 'percent', 'popHits', 'popTotals', 'rfdr', 'scores', 'termName'};
  
  matlabMatrix = [labels; matlabMatrix];


function[tableReport,chartReport,ClusterReport] = gene_annotation_OLD(inputIds,idType)

if nargin<2
    idType = 'AFFYMETRIX_3PRIME_IVT_ID';
end

%For explantion of the Reports results see: http://david.abcc.ncifcrf.gov/content.jsp?file=functional_annotation.html

obj = DAVIDWebService;
authenticate(obj,'michelle.carey@mcgill.ca'); % replace the email with your registered email address
getConversionTypes(obj);
%inputIds = '1000_at, 1001_at, 1002_f_at, 1003_s_at, 1004_at, 1005_at, 1006_at, 1007_s_at';
%idType = 'AFFYMETRIX_3PRIME_IVT_ID';
listName = 'Annotation';
listType = 0;
inputIds = strjoin(inputIds',', ');
addList(obj, inputIds, idType, listName, listType);
Species = getSpecies(obj);
%find('Homo sapiens'==Species)
%Species = setCurrentSpecies(obj,1);

%% Functional Annotation Summary
SummaryReport = [];%getSummaryReport(obj);

%% For Each Gene what Information does DAVID provide
tableReport = [];%getTableReport(obj);

%% Functional Annotation Chart Report
% Maximum EASE p-value 0.1
% Minimum number of genes for the corresponding term
chartReport = getChartReport(obj,0.05,2);

%% Functional Annotation Clustering Report
%  overlap 4, initialSeed 4, finalSeed 4, linkage .50, kappa .35
ClusterReport = [];%getTermClusterReport(obj,4,4,4,0.50,1);