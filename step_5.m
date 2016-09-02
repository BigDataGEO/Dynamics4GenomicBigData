function step_5(list_of_genes, list_of_gene_clusters, indices_of_DRGs, gene_ID_type)

  global Dynamics4GenomicBigData_HOME;

  currentFolder = pwd;
  cd(Dynamics4GenomicBigData_HOME);
  py.importlib.import_module('DAVIDWS');
  cd(currentFolder);

  gene_annotation(gene_ID_type, list_of_genes, indices_of_DRGs, list_of_gene_clusters, 'Step_5', Dynamics4GenomicBigData_HOME, true, true);

end