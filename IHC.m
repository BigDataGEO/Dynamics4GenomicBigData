function [fidxcluster, rmclusters, c, mean_clusters_mat, clusters] = IHC(data, alpha, dist, link)

tic;

n = size(data,1);

if nargin < 3, dist='Spearman'; end
if nargin < 4, link='average'; end

idxclustertmp = clusterdata(data,'criterion','distance','cutoff',1-alpha,'distance',dist,'linkage',link);

for j = 1:max(idxclustertmp)
idxcluster{1}{j}=find(idxclustertmp==j);
end

orig_idxcluster{1} = idxcluster{1};

clusters      = cell(length(idxcluster{1}),1);
rmclusters    = cell(length(idxcluster{1}),1);
for j = 1:length(idxcluster{1})
    clusters{j} = data(orig_idxcluster{1}{j},:);
end

mean_clusters_mat = cellfun(@mod_mean,clusters,'UniformOutput', false); 
mean_clusters_mat = cell2mat(mean_clusters_mat);

k=2;
con =1;
crit = 0;

while(con)
    
     clear('idxclustertmp');
     
     idxclustertmp = clusterdata(mean_clusters_mat,'criterion','distance','cutoff',1-alpha,'distance',dist,'linkage',link);

     for j = 1:max(idxclustertmp)
     idxcluster{k}{j}=find(idxclustertmp==j);
     end
    
     clusters      = cell(length(idxcluster{k}),1);
     for j = 1:length(idxcluster{k})
         orig_idxcluster{k}{j}  = vertcat(orig_idxcluster{k-1}{idxcluster{k}{j}});
         clusters{j} = data(orig_idxcluster{k}{j},:);
     end
     
    %% Prune
    for j = 1:length(clusters)
        PCorR = 1-pdist2(clusters{j},mod_mean(clusters{j}),dist);
        if(any(PCorR<(alpha)))
            indrm = find(PCorR<alpha);
            rmclusters{k}{j}= clusters{j}(indrm,:);
            clusters{j}(indrm,:) = [];
            rmclustersidx{k}{j} = orig_idxcluster{k}{j}(indrm);
            orig_idxcluster{k}{j}(indrm) = [];
        end
    end
    
    if(~isempty(rmclusters{k}))      
          orig_idxcluster{k}  = horzcat(orig_idxcluster{k},num2cell(vertcat(rmclustersidx{k}{:}))');
          clusters = vertcat(clusters,num2cell(vertcat(rmclusters{k}{:}),2));
    end

    mean_clusters_mat = cellfun(@mod_mean,clusters,'UniformOutput', false);
    mean_clusters_mat = cell2mat(mean_clusters_mat);
    
    %% Stopping Conditions
    
    %1. The indices of the clusters are identical
    %2. The number of iterations is >100
    %3. The number of jumping genes is less than 5%.
    
    for l = 1:length(orig_idxcluster{k})
    markov_ind(orig_idxcluster{k}{l}) = l;
    end

    for p = 1:length(orig_idxcluster{k-1})
    mcl_ind(orig_idxcluster{k-1}{p}) = p;
    end
    
    mem = zeros(n,1);
    mem_a = zeros(n,1);
    D = zeros(n,1);
    for i= 1:n
    ind_g_i = find(markov_ind(i)==markov_ind);
    ind_g_i_2 = find(mcl_ind(i)==mcl_ind);
    mem(ind_g_i) = 1;
    mem_a(ind_g_i_2) = 1;
    D(i) = pdist2(mem',mem_a','hamming');
    W(i) = sum(mem)+sum(mem_a);
    end
    
    c(k) = sum(1 - (D'*n)./W)./n;
         
    if(k>3)
    
    m = mean(abs(diff(c((k-3):k))));
         
    if((k>100) || m<0.02)
         
         con1 = 1;
         while(con1)
             
             k=k+1;
             
             idxclustertmp = clusterdata(mean_clusters_mat,'criterion','distance','cutoff',1-alpha,'distance',dist,'linkage',link);

             for j = 1:max(idxclustertmp)
             idxcluster{k}{j}=find(idxclustertmp==j);
             end
             
             clusters      = cell(length(idxcluster{k}),1);
             for j = 1:length(idxcluster{k})
                 orig_idxcluster{k}{j}  = vertcat(orig_idxcluster{k-1}{idxcluster{k}{j}});
                 clusters{j} = data(orig_idxcluster{k}{j},:);
             end
             mean_clusters_mat = cellfun(@mod_mean,clusters,'UniformOutput', false);
             mean_clusters_mat = cell2mat(mean_clusters_mat);
             
             C = 1-pdist(mean_clusters_mat,dist);
             
             con1 = max(max(C))>(alpha);
      
         end
         break;
    end
   
    end
     
     %display(['No of cycles: ',num2str(k)]);

     k=k+1;
    
end

fidxcluster = orig_idxcluster{k};

end


function[y] = mod_mean(x)

if(size(x,1)>1) 
    y = mean(x);
else
    y = x ;
end

end
