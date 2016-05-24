function [latex] = Generate_Latex_Tabels(par_opts,rowLab,colLab,Caption,Label,filename)
%% Generate Latex Tabels
if nargin<4 rowLab = 1:size(par_opts,1); end
if nargin<4 colLab = cellstr(arrayfun(@num2str, 1:size(par_opts,2), 'UniformOutput', false)); end
if nargin<4 Caption =[]; end
if nargin<5 Label =[]; end 
if nargin<6 filename =[]; end     

% Now use this table as input in our input struct:
input.tableData = par_opts;
 
% Set the row format of the data values (in this example we want to use
% integers only
input.tableRowLabels     = cellstr(arrayfun(@num2str, rowLab, 'UniformOutput', false));
input.tableCloumnHeaders = cellstr(colLab);
input.tableDataRowFormat = {'%.3f',length(rowLab)};
 
% Column alignment ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'c';
 
% Switch table borders on/off:
input.tableBorders = 1; 
 
% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;
 
% LaTex table caption:
input.tableCaption = Caption;
 
% LaTex table label:
input.tableLabel = Label;

% Now call the function to generate LaTex code:
latex = latexTable(input);

if(~isempty(filename))
file = fopen(filename,'w');
fprintf(file,'%s',latex{:});
fclose(file);
end


