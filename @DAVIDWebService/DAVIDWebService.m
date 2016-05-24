function obj = DAVIDWebService

obj.endpoint = 'https://david.ncifcrf.gov/webservice/services/DAVIDWebService/';
%.DAVIDWebServiceHttpSoap12Endpoint/';
obj.wsdl =     'https://david.ncifcrf.gov/webservice/services/DAVIDWebService?wsdl';

obj = class(obj,'DAVIDWebService');

