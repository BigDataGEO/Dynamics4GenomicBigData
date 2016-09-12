load('paper.mat');
%% Introduction

%% Materials and methods

%% Experimental design

% We use a pipeline for high-dimensional time course gene expression data. This pipeline
% utilizes functional data analysis (FDA) techniques in order to study the dynamics of gene
% networks.

% The pipeline includes:

%% Preprocessing
%
% Affymetrix Genechip arrays are currently among the most widely used high-throughput technologies
% for the genome-wide measurement of expression profiles. To minimize mis- and cross-hybridization
% problems, this technology includes both perfect match (PM) and mismatch (MM) probe pairs as well
% as multiple probes per gene (Lipshutz et al., 1999). As a result, significant preprocessing is
% required before an absolute expression level for a specific gene may be accurately assessed. In
% general, preprocessing probe-level expression data consists of three steps: background adjustment
% (remove local artifacts and ``noise"), normalization (remove array effects), and summarization at
% the probe set level (combine probe intensities across arrays to obtain a measure of the
% expression level of corresponding mRNA). The Pipeline allows the user to select from the following
% popular preprocessing techniques: Microarray suite 5 (MAS5), Robust Multi-array Average (RMA) and
% Guanine Cytosine Robust Multi-Array Analysis (GCRMA).
%%

[gene_expression, time_points, smooth_gene_trajectories] = step_2(raw_gene_expression, raw_time_points, true);

%%
% Lets assume that the centred expression profile of the $i^{th}$ gene, belonging to subject $j$,
% $X_{i,j}(t)$, is a smooth function over time and that the time course gene expression
% measurements are discrete observations from this smooth function, which have been distorted by
% noise, i.e., $Y_{i,j}(t_{k}) -\mu_{i,j} =X_{i,j}(t_{k})+\epsilon_{i,j}(t_{k})$, for
% $i=1,\ldots,n,$ $j=1,\ldots,N$ and $k=1,\ldots,K_{i,j},$ where $n$ is the number of genes, $N$ is
% the number of subjects, and $K_{i,j}$ is the number of time points observed for the $i^{th}$
% gene, belonging to subject $j$. The noise denoted by $\epsilon_{i,j}(t_{k})$ is assumed to be an
% independently identically distributed (i.i.d.) random variable with mean $0$ and variance
% $\sigma^{2}$. Figure (1) provides an illustration the time course gene expression measurements
% over time for each subject.
%%

%%
% The pipeline incorporates the option to obtain the functional entity $X_{i,j}(t)$ by spline
% smoothing (\cite{green1993nonparametric,silverman2005functional}). Spline smoothing involves minimizing the following
% penalized residual sum of squares $$(Y_{i,j}(t_{k}) -\mu_{i,j} - X_{i,j}(t))^{2} + \lambda \int \left[D^{2} X_{i,j}(t) \textrm{d}t\right]^{2}.$$
%
% The first term defines the discrepancy between the observed centred gene expression measurements
% $Y_{i,j}(t_{k}) -\mu_{i,j}$ and the functional entity $X_{i,j}(t),$ while the second term
% requires $X_{i,j}(t)$ to be sufficiently smooth. The parameter $\lambda$ controls the trade-off
% between these two competing effects and hence ensures that $X_{i,j}(t)$ has an appropriate amount
% of regularity or smoothness. As smoothing each individual gene is computationally expensive and
% is likely to obtain biased estimates as our data does not contain replicates and is expected to
% contain a low signal-to-noise ratio. We follow an approach similar to
% (Yao et al. 2005, Wu and Wu, 2013) and pool the genes together in order to obtain an estimate of
% the amount of regularity required.
%%

%%
% As the majority of the genes have a relatively flat expression profile over time (see Figure 1) estimating the smoothing parameter by pooling all
% the genes together is not ideal as it will tend to over smooth the data. Therefore, we choose a subset of the genes that exhibit time course
% patters that have relatively smooth trajectories that do not fluctuate widely. Then we rank these genes by their interquartile range and select
% the top genes for our estimation subset. Figure 2 illustartes the
% genes expresion curves in the estimation subset over time for each subject.
% The majority of the genes have a relatively flat expression levels over time thus estimating the regularity parameter $\lambda$ using all of the genes together is not ideal as it will tend select a $\lambda$ that is too large in order to minimize the prediction error for the majority of the unresponsive genes. As we are interested in obtaining an appropriate amount of regularity for the responsive gene we apply an approach similar to \cite{yao2005functional} and \cite{wu2013more} we choose a subset of the genes that exhibit time course patterns that have relatively smooth trajectories that do not fluctuate widely. Then we rank these genes by their interquartile range and select the top genes for our subset. The regularity parameter is estimated by minimizing the generalized cross validation (the prediction error) of the responsive genes in our estimation subset. 
%%

%%
% The smoothing parameter is estimated by the conventional generalized cross validation approach using the genes in our estimation subset
% (see Figure 2). Once the smoothing parameter is selected one can proceed to the estimate the functional entity $X_{i,j}(t)$ using all gene
% expression profiles. Table 1 provides the optimal smoothing parameters, the corresponding degrees of freedom,
% the generalized cross validation and the sum of squared errors for the fitted curves $X_{i,j}(t)$ for each subject in the study respectively.
%%

%% Hypothesis Testing to Identify the Dynamic Response Genes (DRGs)
% In a typical gene expression experiment, tens of thousands of genes are measured simultaneously, but only a fraction
% of them are associated with the biological process of interest or a particular stimulus, such as a
% therapeutic treatment or virus infection. Since it is reasonable to
% include only these responsive
% genes in our analysis, here we identify dynamic response genes (DRGs), i.e., genes with expressions that have changed significantly over time.
%
% The pipeline uses the following hypothesis to identify the DRGs,
%
% $H_{0}$: $X_{i,j}(t)$ equal to $0$
%
%
% $H_{a}$: $X_{i,j}(t)$ not equal to $0$.
%
%
% The test statistic is the conventional F-statistic, which compares the goodness-of-fit of the null model to the alternative model:
%
% $F_{i,j}=\frac{RSS_{i,j}^0-RSS_{i,j}^1}{RSS_{i,j}^1},$
%
% where $RSS_{i,j}^0= (Y_{i,j}(t_{k}) -\mu_{i,j} )^{2}$ and $RSS_{i,j}^1=(Y_{i,j}(t_{k}) -\mu_{i,j} - X_{i,j}(t))^{2}$ are the
% residual sum of squares under the null and the alternative models for the $i$-th gene, belonging to subject $j$. The genes with large F-ratios
% can be considered as exhibiting notable changes with respect to time. As we wish to have an equal amount of DRGs for each subject we rank the F-ratios and select 3000 of the top ranking genes (TRGs).
%%


[list_of_DRGs, indices_of_DRGs, indices_of_genes_sorted_by_F_value, smooth_gene_expression, fd_smooth_coefficients] = step_3(list_of_genes, gene_expression, time_points, number_of_top_DRGs_considered, smooth_gene_trajectories, true);



%% Cluster the DRG's into temporal gene response modules
% As many of the DRG's exhibit similar expression patterns over time we
% wish to cluster them into temporal gene response modules
% (groups of genes which have similar gene expression profiles over time). In order to achieve this we adopt the Iterative Hierarchal
% Clustering (IHC) method \cite{careycluster}. The Iterated Hierarchical clustering (IHC) algorithm can identify inhomogeneous
% clusters, capture both the large and very small clusters and provides an automated selection of the optimal number of clusters. The IHC
% algorithm is only reliant on a single parameter $\alpha,$ which controls the trade-off for the between- and within-cluster correlations and
% is outlined below:

%%
%
% #  Initialization: Cluster the gene expression curves using hierarchical clustering. Let the distance metric be the Spearman correlation with a threshold of $\alpha$.
% #  Merge: Treat each of the cluster centres (exemplars) as `new genes', use the same rule as in Step 1 to merge the exemplars into new clusters.
% #  Prune: Let $c_{i}$ be the centre of cluster $i$. If the correlation between the cluster centre and gene $j$, which will be denoted by $\rho_{i,j},$ is less than $\alpha$ remove $gene_{j}$ from the cluster $i$. Let $M$ be the number of genes removed from the existing $N$ clusters. Assign all $M$ genes into single-clusters. Hence there is now $(N+M)$ clusters in total.
% #  Repeat Steps 2 and 3 until either the index of clusters fully or partially converges or the number of iterations reaches the maximum number of iterations.
% #  Repeat Step 2 until the between-cluster correlations are less than $\alpha$.
%

[list_of_gene_clusters, gene_expression_by_cluster, list_of_cluster_means] = step_4(gene_expression, time_points, list_of_DRGs, indices_of_DRGs, indices_of_genes_sorted_by_F_value, smooth_gene_expression, number_of_top_DRGs_considered, true);

%%
% Moreover, these temporal gene response modules can be classified into single-gene modules (SGM) with only one gene in each cluster,
% small-size modules (SSM) that contain between 2-10 genes in each cluster, medium-size modules (MSM) that consist of 11-99 genes in each of
% the clusters and large-size modules (LSM) which contain over 100 genes in each cluster.

%% Functional Annotation of the genes in the temporal gene response modules
% Using the DAVID Functional Annotation Tool \cite{huang2009systematic,huang2009bioinformatics} we produce the following annotaion reports.
%
% # Table Report
% # Chart Report
% # Cluster Report
%
% The cluster report provides the gene-enrichment analysis, pathway mapping, gene/term similarity search, homologue match, ID translation, etc.
% To reduce redundancy, the Functional Annotation Clustering report groups/displays similar
% annotations together which makes the biology clearer and more focused to be read vs.
% the traditional chart report. The grouping algorithm is based on the hypothesis that similar
% annotations should have similar gene members. The Functional Annotation Clustering integrates the
% same techniques of  Kappa statistics to measure the degree of the common genes between two
% annotations, and fuzzy heuristic clustering (used in Gene Functional Classification Tool) to
% classify the groups of similar annotations according kappa values. In this sense, the more common
% genes annotations share, the higher chance they will be grouped together.
%
% The table report provides the ID, Gene Name, Species, BBID, BIOCARTA, COG_ONTOLOGY, GOTERM_BP_FAT, GOTERM_CC_FAT, GOTERM_MF_FAT, INTERPRO, KEGG_PATHWAY, OMIM_DISEASE, PIR_SUPERFAMILY, SMART and SP_PIR_KEYWORDS for all the genes passed as input.
%
% The chart report is an annotation-term-focused view which lists annotation terms and their associated genes under study.
% The threshold of EASE Score, a modified Fisher Exact P-Value, for gene-enrichment analysis. It ranges from 0 to 1.
% Fisher Exact P-Value = 0 represents perfect enrichment. Usually P-Value is equal or smaller than 0.05 to be considered strongly enriched in the annotation categories. Default is 0.1.

step_5(list_of_genes, list_of_gene_clusters, indices_of_DRGs, gene_ID_type);

%% Construct high-dimensional gene regulation networks (GRNs) using differential equation models
%
% There has been an abundance of literature regarding the use of ordinary differential equation (ODE) modelling for constructing high-dimensional gene regulation networks (GRN) e.g., \cite{hecker2009gene}, \cite{lu2011high} and \cite{wu2013high}. A gene regulatory network, attempts to map how different genes control the expression of other genes. The gene regulations can be modelled by rate equations,
% $$
% D M_{q,j} = \alpha_{0,q,j} + \sum_{p=1}^{Q} \alpha_{p,q,j}M_{p,j}, \qquad \textrm{for} \quad q =1,\ldots,Q,
% $$
% where $\alpha_{0,q,j}$ is the intercept for the $q^{th}$ gene response module, belonging to subject $j$ and the coefficients $\{\alpha_{p,q,j}\}_{p=1}^{Q}$ quantify the regulation effects of the $p^{th}$ gene response module on the instantaneous rate of change in $q^{th}$ gene response module. This model can appropriately capture both up and down regulations as well as up and down self-regulations. Typically, only a few gene response modules will effect the instantaneous rate of change in $q^{th}$ gene response module thus only a few of the $\{\alpha_{p,q,j}\}_{p=1}^{Q}$ will be non-zero. We first perform a model selection which determines which $\{\alpha_{p,q,j}\}_{p=1}^{Q}$ are non-zero and then we estimate their coefficients to determine the regulation effects. 
% Here we use the ordinary differential equation (ODE) modelling approach in order to reconstruct the high-dimensional gene regulation networks (GRN).
% In ODE network models, gene regulations are modeled by rate equations, which quantify the rate of change (derivative) of the expression level
% of one temporal gene response module in the system as a function of expression levels of all related temporal gene response modules.
% Both up and down regulations as well as self-regulations can be appropriately captured by the ODE model.
% The general form of the ODE model can be written as
%
% $\frac{d(M_{q,j}}{dt} = \alpha_{0,q,j} + \sum_{p=1}^{P} \alpha_{p,q,j}M_{p,j},$
%
% for i =q,..,Q where $M_{q,j}$ is the $q^{th}$ gene response module, belonging to subject $j$, $\alpha_{0,q,j}$ is the intercept for the $q^{th}$
% gene response module, belonging to subject $j$ and the coefficients $\{\alpha_{p,q,j}\}_{i=1}^{P}$ quantify the regulation effects of other
% gene response modules on the rate of expression change of the q-th gene response module, belonging to subject $j$.
%
% Although the model size Q is comparably smaller than that of the original model N,
% simultaneous model selection and parameter optimization of ODE parameters $\{\alpha{p,q,j}\}_{q,p=1}^{Q,P}$ are still computationally very expensive,
% because it involves costly numerical integration and complicated parameter regularization.
% The two-stage smoothing-based estimation method (\cite{voit2004decoupling,liang2008parameter}), decouples the system of
% differential equations into a set of pseudo-regression models.  This method avoids numerically solving the differential equations directly and 
% does not require the initial or boundary conditions of the state variables. More importantly, it allows us to perform model selection
% and parameter estimation for one equation at a time, which significantly reduces the computational cost. 


%%
% _Step one_
%
% # Obtain the estimates of the temporal gene response modules $M_{q,j}$ and their derivatives $\frac{d(M_{q,j}}{dt}$ using the smoothing splines estimates obtained in step one.
%
% _Step two_
%
% We plug the estimated mean expression curves $\hat{M}_{q,j}$ and their derivatives $\frac{d(\hat{M}_{q,j}}{dt}$ into the ODE model to obtain the following set of pseudo linear regression models.
%
% $\frac{d(\hat{M}_{q,j}}{dt} = \beta_{0,q,j} + \sum_{p=1}^{P} \beta_{p,q,j}\hat{M}_{p,j} + \epsilon_{p,j}$
%
% for q=1,...,Q. For linear regression models, many penalized methods have been proposed in the regularization framework to conduct variable selection and estimation, such as the
% least absolute shrinkage and selection operator (LASSO) (\cite{tibshirani1996regression}), smoothly clipped absolute deviation (SCAD) (Fan and Li, 2001) and
% so on. Here we use LASSO regularization to conduct variable selection and estimation using the the follwing following penalized objective function
%
% $\epsilon_{p,j} + \lambda \sum_{p=0}^{P} \| \beta_{p,q,j} \|$
%
% where $\sum_{p=0}^{P} \| \beta_{p,q,j} \|$ is the LASSO penalty which
% shrinks the parameters $\beta_{p,q,j}$ to zero and $\lambda$ is the sparisty parameter which controls the trade-off between minimizing the error of the
% pseudo linear regression model and requiring the parameters $\beta_{p,q,j}$ to be zero. A large $\lambda$ sets all the parameters $\beta_{p,q,j}$ to zero.
% A low $\lambda$ minmizes the error of the pseudo linear regression model.
%
% _Refinement Step_
%
% The two-stage smoothing-based estimation method is employed to simplify the computation
% of the ODE models and also facilitate the variable selection procedure. However, the
% parameter estimates from the two-stage method are not efficient in terms of estimation accuracy,
% as there can be considerable approximation error in the estimates of the
% modules expression curves $\hat{M}_{q,j}$ and their derivatives
% $\frac{d(\hat{M}_{q,j}}{dt}$ and decoupling the system of
% differential equations into a set of pseudo-regression models only accounts for the direct effect, that is, how the $q^{th}$ gene response module
% is regulated by the $p^{th}$ gene response module (q $\rightarrow$ p) and does not account for the indirect effect, that is, how the $r^{th}$ gene response module
% is regulated by the $q^{th}$ gene response module which in turn regululates the $p^{th}$ gene response module (r $\rightarrow$ q $\rightarrow$ p) which 
% is a fundamental property of any system of differential equations.
% 
% To overcome the estimation deficiency and the decoupling effect of the two-stage method, we propose to refine the parameter estimates for the
% selected ODE model using the nonlinear least squares (NLS) method. The parameter estimates
% from the two-stage method are used to estimate the sturcture of the
% estimates only, that is, to inform us which parameters are non-zero.

% Non-linear least squares attributable to \cite{hemker1972numerical,bard1974nonlinear} estimates the parameters $\{\hat{\alpha}_{p,q,j}\}_{p=0}^{Q}$ by minimizing the dependency between the observed measurements of the temporal gene response modules and the numerical approximation to the solution of (\ref{ODE}).
% The initial parameter estimates are given by the non-zero $\{\hat{\beta}_{p,q,j}\}_{p=0}^{Q}$ parameters from the two-stage method the remaining parameters are set to zero. The initial states of the temporal gene response modules are $M_{q,j}(0)$ and $DM_{q,j}(0)$ for $q =1,\ldots,Q.$  

adjacency_matrix_of_gene_regulatory_network = step_7(list_of_gene_clusters, time_points, indices_of_DRGs, fd_smooth_coefficients, true);

%% Obtain Network Analysis of the gene regulation networks (GRNs).
% Graph theorists and network analysts have developed a number of metrics to characterise biological networks for an overview see \cite{huber2007graphs} and \cite{lee2004coexpression}. These metrics facilitate drug target identification and insight on potential strategies for treating various diseases. The pipeline uses the SBE Toolbox description of the metrics which are produced by the pipeline are listed in Table (\ref{Graph Metrics}).
%   Graph Statistics
%   ----------------
%
%       graph_clustercoeff       - Calculate overall clustering coefficient of the graph.
%       graph_diameter           - Calculate diameter of the graph.
%       graph_meandist           - Calculate mean distance of the graph.
%       graph_density            - Calculate density of the graph.

%%
%   Node Statistics
%   ---------------
%
%       bridging_centrality      - Calculate Bridging Centrality.
%       closeness_centrality     - Calculate Closeness Centrality.
%       eccentricity_centrality  - Calculate Eccentricity Centrality.
%       delta_centrality         - Calculate Delta Centrality.
%       current_info_flow        - Calculate Current Information Flow.
%       assortativitycoeff       - Calculate Assortativity Coefficient.
%       bridgingcoeff            - Calculate Bridging Coefficient.
%       brokeringcoeff           - Calculate Brokering Coefficient.
%       clusteringcoeff          - Calculate Clustering Coefficient.
%       order2cc                 - Calculate Clustering Coefficient considering 2nd level of connections.
%       hierarchy                - Calculate Degree of Hierarchy.
%       kclique                  - Caldulate K-clique.
%       kcore                    - Calculate K-core of network.
%       locavgcon                - Calculate Local Average Connectivity.
%       neighborhood_conn        - Calculate Neighborhood Connectivity.     
%       participationcoeff       - Calculate Participation Coefficient.
%       richclubcoeff            - Calculate Rich Club Coefficient.
%       smallworldindex          - Calculate Small World Index.
%       within_module_deg        - Calculate Within Module Degree.
%

step_8(adjacency_matrix_of_gene_regulatory_network);

%% Results

%% Discussion

%% Conclusion

%% Supporting Information

%%
% \bibliography{bibliography}
%
