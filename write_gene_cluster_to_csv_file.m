function write_gene_cluster_to_csv_file(ids_of_genes_in_current_cluster,output_file_name)

  fid = fopen(output_file_name,'w');
  for row = 1:size(ids_of_genes_in_current_cluster,1)
    ids_of_genes_in_current_cluster = ids_of_genes_in_current_cluster;
    fprintf(fid, repmat('%s\t',1,size(ids_of_genes_in_current_cluster,2)-1), ids_of_genes_in_current_cluster{row,1:end-1});
    fprintf(fid, '%s\n', ids_of_genes_in_current_cluster{row,end});
  end
  fclose(fid);