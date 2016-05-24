function [success,theMessage] = create_exel_file(fileName,xlsData,sheetName,startRange,path)

if nargin < 3
    sheetName = Sheet1;
    startRange = '';
elseif nargin < 4
    startRange = '';
end

if(ispc)
[success,theMessage] = xlswrite(fileName,xlsData,sheetName,startRange);
else
% Add Java POI Libs to matlab javapath
javaaddpath([path,'poi_library/poi-3.8-20120326.jar']);
javaaddpath([path,'poi_library/poi-ooxml-3.8-20120326.jar']);
javaaddpath([path,'poi_library/poi-ooxml-schemas-3.8-20120326.jar']);
javaaddpath([path,'poi_library/xmlbeans-2.3.0.jar']);
javaaddpath([path,'poi_library/dom4j-1.6.1.jar']);
javaaddpath([path,'poi_library/stax-api-1.0.1.jar']);
[success,theMessage] = xlwrite(fileName, xlsData, sheetName, startRange);
end

end
