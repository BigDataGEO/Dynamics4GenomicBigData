function results = load_analysis(GEO_number, condition)

  path_to_results_file = ['Output/' GEO_number '/Conditions/' condition '/' 'Results.mat'];
  
  if(~exist(path_to_results_file, 'file'))
    msgID = 'MATLAB:rmpath:DirNotFound1';
    msg = ['Unable to retrieve analysis results for condition ' condition ' associated to GEO series ' GEO_number '.'];
    baseException = MException(msgID,msg);    
    throw(baseException);    
  else
    load(path_to_results_file, 'results');
  end
end
