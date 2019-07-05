%  cd('C:\Workspace\Pipeline')

set_paths_and_imports;

folder_with_simulated_datasets = [Dynamics4GenomicBigData_HOME filesep 'Input'];

cd(folder_with_simulated_datasets);

dataset_file_names = dir('*.csv');

for j = 1:numel(dataset_file_names)
    dataset_file_name = dataset_file_names(j).name;
    
    [filepath, name, ext] = fileparts(dataset_file_name);
            
    [gene_expression, standardized_gene_expression, time_points, gene_ID_type, smooth_gene_trajectories] = step_2_from_csv(dataset_file_name, false);
            
    [gene_expression_sorted_by_F_value, standardized_gene_expression_sorted_by_F_value, number_of_statistically_significant_DRGs, smooth_gene_expression, fd_smooth_coefficients] = step_3(gene_expression, time_points, smooth_gene_trajectories, 3000, false);
            
    writetable(gene_expression_sorted_by_F_value(1:number_of_statistically_significant_DRGs,1), [name, '.drgs_pipeline'], 'WriteVariableNames', false, 'FileType', 'text');
end

cd(Dynamics4GenomicBigData_HOME);
