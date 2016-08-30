function step_5(fidxcluster, gid, IND_DRG, gene_ID_type)

  global Dynamics4GenomicBigData_HOME;

  currentFolder = pwd;
  cd(Dynamics4GenomicBigData_HOME);
  py.importlib.import_module('DAVIDWS');
  cd(currentFolder);

%    [uselessVariable, cluster_indexes_by_size] = sort(cellfun('size', fidxcluster, 1), 'descend');
%    clusters_sorted_by_size = fidxcluster(cluster_indexes_by_size);

  gene_annotation(gene_ID_type, gid, IND_DRG, fidxcluster, 'Step_5', Dynamics4GenomicBigData_HOME, true, true);

end