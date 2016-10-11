function index_of_gpl_column_with_gene_ids = capture_index_of_gpl_column_with_gene_ids(platform_struct)

    prompt = ['Check platform ' platform_struct.Accession '''s record on GEO and enter the index of the column where the gene IDs/names are located: '];
   
    index_of_gpl_column_with_gene_ids = input(['\n' prompt]);
    if(~isnumeric(index_of_gpl_column_with_gene_ids))
      index_of_gpl_column_with_gene_ids = 1;
    else if(index_of_gpl_column_with_gene_ids < 1 | index_of_gpl_column_with_gene_ids > size(platform_struct.Data,2))
      index_of_gpl_column_with_gene_ids = 1;
    end
end