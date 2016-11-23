function run_condition(GEO_number, condition, samples, time_points, number_of_top_DRGs)

  global Dynamics4GenomicBigData_HOME;

  try
    fprintf('\n');
    display(['Loading dataset. This can take some time, please wait...']);
    [geoStruct, list_of_genes, gene_ID_type, list_of_probe_ids] = get_geo_data(GEO_number);
  catch
    fprintf('\n');
    display(['Could not retrieve dataset ' GEO_number ' from the Gene Expression Omnibus.']);
    fprintf('\n');
    display(['This is possibly because the GEO refused the FTP connection or because the dataset does not exist.']);
    fprintf('\n');
    display(['Please download manually ' GEO_number '''s matrix to ' pwd '/GEO_cache/' GEO_number '.soft and try again.']);
    return;
  end

  [raw_gene_expression_array, raw_time_points_array] = step_1(geoStruct, samples, time_points);

  fprintf('\n');
  display(['The analysis of condition "' condition '" is starting.']);

  run_pipeline_analysis_on_condition(GEO_number, list_of_genes, raw_gene_expression_array, raw_time_points_array, condition, condition, gene_ID_type, number_of_top_DRGs, list_of_probe_ids, geoStruct);
    
  fprintf('\n');
  display(['The analysis of condition "' condition '" has been completed.']);
    
  fprintf('\n');
  display(['Results have been output to folder ' Dynamics4GenomicBigData_HOME 'Results/' GEO_number '/Conditions/' condition '/']);

end