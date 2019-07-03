set_paths_and_imports;

folder_with_simulated_datasets = [Dynamics4GenomicBigData_HOME filesep 'Input' filesep 'Output'];

cd(folder_with_simulated_datasets);

subfolders = dir('*');
for j = 1:numel(subfolders)
    subfolder = subfolders(j).name;
    
    if strcmp(subfolder, '.')==0 & strcmp(subfolder, '..')==0
    
        cd(subfolder);
        subsubfolders = dir('*.csv');
        
        for k = 1:numel(subsubfolders)
            dataset_file_name = subsubfolders(k).name;
            
            [filepath, name, ext] = fileparts(dataset_file_name);
            
            [gene_expression, standardized_gene_expression, time_points, gene_ID_type, smooth_gene_trajectories] = step_2_from_csv(dataset_file_name, false);
            
            [gene_expression_sorted_by_F_value, standardized_gene_expression_sorted_by_F_value, number_of_statistically_significant_DRGs, smooth_gene_expression, fd_smooth_coefficients] = step_3(gene_expression, time_points, smooth_gene_trajectories, 3000, false);
            
            writetable(gene_expression_sorted_by_F_value(1:number_of_statistically_significant_DRGs,1), [name, '_Pipeline_DRGs.txt'], 'WriteVariableNames', false);
            
        end
        cd('..');
    end
end

cd(Dynamics4GenomicBigData_HOME);
