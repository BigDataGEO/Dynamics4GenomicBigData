function [cluster_indexes_by_size, clusters_sorted_by_size] = step_5(i, fidxcluster, gid, IND_DRG, path)

  [uselessVariable, cluster_indexes_by_size] = sort(cellfun('size', fidxcluster{i}, 1), 'descend');
  clusters_sorted_by_size = fidxcluster{i}(cluster_indexes_by_size);

  gene_annotation(gid, IND_DRG{i}, clusters_sorted_by_size, 'Annotation', path, true, true);