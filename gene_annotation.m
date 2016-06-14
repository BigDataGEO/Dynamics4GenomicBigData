function[tableReport,chartReport,ClusterReport] = gene_annotation(inputIds,idType)

if nargin<2
    idType = 'AFFYMETRIX_3PRIME_IVT_ID';
end

inputIds = py.str(strjoin(inputIds',', '));
idType = py.str(idType)
listName = py.str('Annotation')
listType = py.int(0)
thd = py.float(0.05)
count = py.int(2)

speciesPy = py.DAVIDWS.getSpecies(inputIds, idType, listName, listType)

Species = char(speciesPy)

chartReportPy = py.DAVIDWS.getChartReport(inputIds, idType, listName, listType, thd, count)

chartReport = char(chartReportPy)

SummaryReport = [];
tableReport = [];
ClusterReport = [];

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