
%% Author Summary

%% Introduction
%
% Read the papers opened in google scholar to help with introduction.

%% Materials and Methods

%% Expermental Design
if ~isempty(geoStruct.Header.Samples.description)
for i = 1:length(cond)
    wraptext(geoStruct.Header.Samples.description{pr_ind(cond(i))},60);
    fprintf('%\n');
end
end

%  if ~isempty(geoStruct.Header.Samples.label_protocol_ch1)
if(isfield(geoStruct.Header.Samples,'label_protocol_ch1'))
wraptext(geoStruct.Header.Samples.label_protocol_ch1{pr_ind(cond)},60)
end

%% Step 3. Detect the Dynamic Response Genes (DRGs)
%
% Lets assume that the centred expression levels of the $i^{th}$ gene, belonging to subject $j$, $X_{i,j}(t)$, is a smooth function over time $t$ 
% and that each measurement represents a discrete observations from this smooth function, which have been distorted by noise, i.e., 
% $$Y_{i,j}(t_{k}) -\mu_{i,j} = X_{i,j}(t_{k})+\epsilon_{i,j}(t_{k}),$$ 
% for $i=1,\ldots,n,$ $j=1,\ldots,N$ and $k=1,\ldots,K_{i,j},$ where $n$ is the number of genes, $N$ is the number of subjects, 
% $K_{i,j}$ is the number of time points observed for the $i^{th}$ gene, belonging to subject $j$ and $\mu_{i,j}$ is the average 
% gene expression measurement over time for the $i^{th}$ gene, belonging to subject $j$. The noise $\epsilon_{i,j}(t_{k})$ is assumed to be 
% an independently identically distributed (i.i.d.) random variable with mean $0$ and variance $\sigma^{2}$. 
%
% The estimated centred expression levels $\hat{X}_{i,j}(t)$ are produced by spline smoothing attributable to \cite{green1993nonparametric,silverman2005functional}. This approach approximates $X_{i,j}(t)$ by a linear combination of $L$ independent basis functions, $X_{i,j}(t) \approx \sum_{l=1}^{L} c_{i,j,l}B_{i,j,l}(t)$, where $B_{i,j,l}(t)$ are the basis functions and $c_{i,j,l}$ are the corresponding coefficients. The estimated coefficients $\hat{c}_{i,j,l}$ are obtained by minimising 
% $$
% (Y_{i,j}(t_{k}) -\mu_{i,j} - X_{i,j}(t))^{2} + \lambda \int \left[D^{2} X_{i,j}(t) \textrm{d}t\right]^{2}.
% $$
% The first term defines the squared discrepancy between the observed centred gene expression measurement $Y_{i,j}(t_{k}) -\mu_{i,j},$ and the estimated function $X_{i,j}(t),$ while the second term is the integral of squared second derivative of the estimated function which requires $X_{i,j}(t)$ to be sufficiently smooth. The parameter $\lambda$ controls the trade-off between the fit to the data and the smoothness requirement and hence ensures that $X_{i,j}(t)$ has an appropriate amount of regularity. 
%
% The majority of the genes have a relatively flat expression levels over time thus estimating the regularity parameter $\lambda$ using all of the genes together is not ideal as it will tend select a $\lambda$ that is too large in order to minimize the prediction error for the majority of the unresponsive genes. As we are interested in obtaining an appropriate amount of regularity for the responsive gene we apply an approach similar to \cite{yao2005functional} and \cite{wu2013more} we choose a subset of the genes that exhibit time course patterns that have relatively smooth trajectories that do not fluctuate widely. Then we rank these genes by their interquartile range and select the top genes for our subset. The regularity parameter is estimated by minimizing the generalized cross validation (the prediction error) of the responsive genes in our estimation subset. 
%
% The majority of the genes have a relatively flat expression levels over time thus estimating the regularity parameter $\lambda$ using all of the genes together is not ideal as it will tend select a $\lambda$ that is too large in order to minimize the prediction error for the majority of the unresponsive genes. As we are interested in obtaining an appropriate amount of regularity for the responsive gene we apply an approach similar to \cite{yao2005functional} and \cite{wu2013more} we choose a subset of the genes that exhibit time course patterns that have relatively smooth trajectories that do not fluctuate widely. Then we rank these genes by their interquartile range and select the top genes for our subset. The regularity parameter is estimated by minimizing the generalized cross validation (the prediction error) of the responsive genes in our estimation subset. 
%%

tmp2 = zeros(length(cond),4);
myVars = {'dfgenens','gcvgenens','lambdagenes','STDERR'};
for i = 1:length(cond)
    load(strcat(flder,'/',strcat(GEO_number,conditions_analyzed{cond(i)},date),'.mat'),myVars{:});
    tmp2(i,:) = round2([cell2mat(dfgenens), mean(cell2mat(gcvgenens),2),log10(cell2mat(lambdagenes)),sum(cell2mat(STDERR),2)]);
end

colLab = {'Df','GCV','log10($\lambda$)','Std Error'};
rowLab = strcat(repmat({'Subject '},N,1),cellstr(arrayfun(@num2str, cond, 'UniformOutput', false))');
Caption = 'The degrees of freedom (Df), the generalized cross validation (GCV), the smoothing paramter $\lambda$ and the standard error of the fitted curves produced by spline smoothing';
Label = 'fit1';
disp(char(Generate_Latex_Tabels(tmp2,rowLab,colLab,Caption,Label)));

ind= 0 ;
myVars = {'yhat','Time'};
for i = 1:length(cond)
    load(strcat(flder,'/',strcat(GEO_number,conditions_analyzed{cond(i)},date),'.mat'),myVars{:});
    for sub = 1:N
        surf(yhat{sub}','FaceColor','interp','EdgeColor','none');
        
        axisLabelFontSize = 30;
        
        set(gcf, 'units', 'centimeters');
        set(gcf, 'position', [0, 0, 30, 24]);
        
        set(gcf, 'PaperPositionMode', 'manual');
	set(gcf, 'PaperUnits', 'centimeters');
	set(gcf, 'PaperPosition', [0 0 30 24]);
	set(gcf, 'PaperUnits', 'centimeters');
	set(gcf, 'PaperSize', [30 24]);        
        
        set(gca,'XTick',1:length(Time{sub}),'Xticklabel',Time{sub})
        set(gca,'FontSize',11);
        
        xlim([Time{sub}(1),length(Time{sub})])
        ylim([1,size(yhat{sub},2)])
        zlim([min(min(yhat{sub})),max(max(yhat{sub}))])
        
        xlabel('Time', 'FontSize', axisLabelFontSize)
        ylabel('Probe Set', 'FontSize', axisLabelFontSize)
        zlabel('Probe Set Expression Level', 'FontSize', axisLabelFontSize)
        
        title([char(Subject_name),' smooth gene expresion curves'], 'FontSize', axisLabelFontSize)
    end
end

 
%%
% DRGs can be defined as genes with expressions that change significantly with respect to time. We use the following hypothesis to detect the DRGs,
% $$H_{0}: \hat{X}_{i,j}(t) = 0$$
% $$H_{a}: \hat{X}_{i,j}(t) \neq 0.$$
%
% The test statistic is the conventional F-statistic, which compares the goodness-of-fit of the null model to the alternative model:
% $$
% F_{i,j}=\frac{\textrm{RSS}_{i,j}^0-\textrm{RSS}_{i,j}^1}{\textrm{RSS}_{i,j}^1},
% $$
% where $\textrm{RSS}_{i,j}^0= (Y_{i,j}(t_{k}) -\mu_{i,j} )^{2}$ and $\textrm{RSS}_{i,j}^1=(Y_{i,j}(t_{k}) -\mu_{i,j} - \hat{X}_{i,j}(t))^{2}$ are the residual sum of squares under the null and the alternative models for the $i$-th gene, belonging to subject $j$. The genes with large F-ratios can be considered as exhibiting notable changes with respect to time. As we wish to have an equal amount of DRGs for each subject, we rank the F-ratios and select 3000 of the top ranking dynamic response genes.
%
% Look at the annotation.xlsx file and report the most enriched pathways,
% GO BP, MF and CC temrs. 
%

%% Cluster these DRGs into temporal gene response modules (GRMs)
%
% As many of the DRGs exhibit similar expression patterns over time, we wish to cluster them into co-expressed modules (groups of genes which have similar gene expression levels over time). It is widely recognized that many co-expressed genes may follow similar temporal patterns, but at the same time, some genes may have very few or even no co-expressed genes and thus may exhibit unique temporal response patterns. Consequently, the temporal gene response modules or clusters can be inhomogeneous, i.e., some clusters are very large and contain many genes while others are small or even only contain a single gene. In order to obtain inhomogeneous clusters we adopt the Iterative Hierarchical Clustering (IHC) method introduced in \cite{careycluster}. This approach is only reliant on a single parameter $\alpha,$ which controls the trade-off of the between- and within-cluster correlations, in particular, the average within-cluster correlation will be approximately $\alpha$ and the between- cluster correlation will be below $\alpha$. The IHC algorithm is outlined below:
% 
% Initialization: Cluster the DRG using hierarchical clustering. Let the distance metric be the Spearman rank correlation with a threshold of $\alpha$ and the linkage method be the average of the genes in the clusters.
% Merge: Treat each of the cluster centres (exemplars) as `new genes', use the same rule as in Step 1 to merge the exemplars into new clusters. The cluster centres provide the average time-course pattern of the cluster members.
% Prune: Let $c_{i}$ be the centre of cluster $i$. If the correlation between the cluster centre and gene $j$, which will be denoted by $\rho_{i,j},$ is less than $\alpha,$ then remove $gene_{j}$ from the cluster $i$. Let $M$ be the number of genes removed from the existing $N$ clusters. Assign all $M$ genes into single-element clusters. Hence, there is now $(N+M)$ clusters in total.
% Repeat Steps 2-3 until the index of clusters converges.
% Repeat Step 2 until the between-cluster correlations are less than $\alpha$.
% 
% The IHC method typically identifies four types of temporal gene response modules: single-gene modules (SGM) with only one gene in each cluster, small-size modules (SSM) that contain between 2-10 genes in each cluster, medium-size modules (MSM) that consist of 11-99 genes in each of the clusters and large-size modules (LSM) which contain over 100 genes in each cluster.  Let $M_{q,j}$ denote the centre of the $q^{th}$ temporal gene response module, belonging to subject $j$, for $q=1,\ldots,Q$. An estimate of $M_{q,j}$ is obtained by averaging the estimated gene expression levels provided in Step 3, $\hat{X}_{i,j}(t),$ for all the genes contained in the $q^{th}$ gene response module. 



tmp2 = zeros(length(cond),4);
myVars ={'dfgenens','gcvgenens','lambdagenes','STDERR'};
for i = 1:length(cond)
    load(strcat(flder,'/',strcat(GEO_number,conditions_analyzed{cond(i)},date),'.mat'),myVars{:});
    tmp2(i,:) = round2([cell2mat(dfgenens), mean(cell2mat(gcvgenens),2),log10(cell2mat(lambdagenes)),sum(cell2mat(STDERR),2)]);
end

colLab = {'No. of Modules','No. of LSM','No. of MSM','No. of SSM','No. of SGM'};
rowLab = strcat(repmat({'Subject '},N,1),cellstr(arrayfun(@num2str, cond, 'UniformOutput', false))');
Caption = 'The number of gene response modules, number of LSM, MSM, SMS and SGM';
Label = 'fit1';
disp(char(Generate_Latex_Tabels(tmp2,rowLab,colLab,Caption,Label)));



%% Obtain the functional enrichment analysis of the GRMs
%
% The pipeline uses the DAVID (the database for annotation, visualization and integrated discovery) functional annotation tool developed by 
% \cite{huang2009systematic,huang2009bioinformatics} to produce its annotation report. The functional annotation tool systematically extracts 
% biological meaning from large genomic data. Thus given a list of probe set ids the Pipeline can utilise the DAVID database to obtain the 
% associated biological annotation (e.g., gene ontology terms, associated pathways, functional analysis, e.c.t). 

%% Construct the high-dimensional gene regulation network (GRN) that determines the interactions between the GRMs
%
% There has been an abundance of literature regarding the use of ordinary differential equation (ODE) modelling for constructing high-dimensional gene regulation networks (GRN) e.g., \cite{hecker2009gene}, \cite{lu2011high} and \cite{wu2013high}. A gene regulatory network, attempts to map how different genes control the expression of other genes. The gene regulations can be modelled by rate equations,
% $$
% D M_{q,j} = \alpha_{0,q,j} + \sum_{p=1}^{Q} \alpha_{p,q,j}M_{p,j}, \qquad \textrm{for} \quad q =1,\ldots,Q,
% $$
% where $\alpha_{0,q,j}$ is the intercept for the $q^{th}$ gene response module, belonging to subject $j$ and the coefficients $\{\alpha_{p,q,j}\}_{p=1}^{Q}$ quantify the regulation effects of the $p^{th}$ gene response module on the instantaneous rate of change in $q^{th}$ gene response module. This model can appropriately capture both up and down regulations as well as up and down self-regulations. Typically, only a few gene response modules will effect the instantaneous rate of change in $q^{th}$ gene response module thus only a few of the $\{\alpha_{p,q,j}\}_{p=1}^{Q}$ will be non-zero. We first perform a model selection which determines which $\{\alpha_{p,q,j}\}_{p=1}^{Q}$ are non-zero and then we estimate their coefficients to determine the regulation effects. 
%
% (a) Model Selection
% 
% The two-stage smoothing-based estimation method (\cite{voit2004decoupling,liang2008parameter}), decouples the system of
% differential equations into a set of pseudo-regression models.  This method significantly reduces the computational cost, avoids numerically solving the differential equations directly and does not require initial conditions $M_{q,j}$. The two-stage method is outlined below,
% Let $D \hat{M}_{q,j}=\sum_{l=1}^{L} \hat{c}_{i,j,l}DB_{i,j,l}(t)$, where $DB_{i,j,l}(t)$ are the derivatives of the basis functions and $\hat{c}_{i,j,l}$ are the corresponding coefficients.
% The corresponding set of pseudo linear regression models are
%	$$D\hat{M}_{q,j} = \beta_{0,q,j} + \sum_{p=1}^{Q} \beta_{p,q,j}\hat{M}_{p,j} + \epsilon_{p,j}$$
%	for $q=1,\ldots,Q$. In order to conduct variable selection we use the least absolute shrinkage and selection operator (LASSO) \citep{tibshirani1996regression} penalised regression procedure. LASSO uses the the following penalized objective function 
%	$$\epsilon_{p,j} + \omega \sum_{p=0}^{Q} \| \beta_{p,q,j} \|,$$ 
%	where $\sum_{p=0}^{Q} \| \beta_{p,q,j} \|$ is the LASSO penalty which shrinks the parameters $\beta_{p,q,j}$ to zero and $\omega$ is the sparsity parameter which controls the trade-off between minimising the error of the pseudo linear regression model and requiring the parameters $\beta_{p,q,j}$ to be zero. 
%
% As the parameter estimates from the two-stage method are not efficient in terms of estimation accuracy, as there can be considerable approximation error in the estimates of $M_{q,j}$ and its derivatives. The $\{\hat{\beta}_{p,q,j}\}_{p=0}^{Q}$ estimates are only used to inform us which parameters are non-zero. 
%
% (b) Estimation of the regulation effects $\{\hat{\alpha}_{p,q,j}\}_{p=1}^{Q}$\\}
%
% Non-linear least squares attributable to \cite{hemker1972numerical,bard1974nonlinear} estimates the parameters $\{\hat{\alpha}_{p,q,j}\}_{p=0}^{Q}$ by minimizing the dependency between the observed measurements of the temporal gene response modules and the numerical approximation to the solution of (\ref{ODE}).
% The initial parameter estimates are given by the non-zero $\{\hat{\beta}_{p,q,j}\}_{p=0}^{Q}$ parameters from the two-stage method the remaining parameters are set to zero. The initial states of the temporal gene response modules are $M_{q,j}(0)$ and $DM_{q,j}(0)$ for $q =1,\ldots,Q.$  

%% Perform a network feature and dynamic property analysis on the GRN

%Graph theorists and network analysts have developed a number of metrics to characterise biological networks for an overview see \cite{huber2007graphs} and \cite{lee2004coexpression}. These metrics facilitate drug target identification and insight on potential strategies for treating various diseases. The pipeline uses the SBE Toolbox description of the metrics which are produced by the pipeline are listed in Table (\ref{Graph Metrics}). 

%%
% \begin{table}[!h]
% 	\caption{\textbf{Network metrics}\label{Graph Metrics}}
% 		\resizebox{\columnwidth}{!}{%
% 	\begin{tabular}{|l|l|}
% 		\hline
% 		Assortativity             & Assortativity coefficient of a graph\\
% 		Response2TNremoval    & Response to targeted node removal\\
% 		Betweenness            & Betweenness centrality of all network nodes\\
% 		Centralization            & Compute connectivity centralization\\
% 		Charpath                  & Characteristic path length\\
% 		Closeness                 &Compute closeness centrality\\
% 		Clusteringcoefs           & Clustering coefficients of all network nodes\\
% 		Degreecentrality          & Compute degree centrality of a node\\
% 		Degreedistribution        & Return node degree distribution vector\\
% 		Diameter                 & Return network diameter\\
% 		Shortest Path             & Return all shortest paths\\
% 		Distance                  &Global network distances\\
% 		Eccentricity              &Eccentricity of each vertex\\
% 		Efficiency                &Local node efficiency\\
% 		Eigencentrality           & Eigenvector centrality\\
% 		Response2RNremoval      & Response to random node removal\\
% 		Gethubs                   & Find network hub nodes\\
% 		Heterogeneity             & Network heterogeneity\\
% 		Maxspantree               &Extract a maximum spanning tree\\
% 		Minspantree               & Extract a minimum spanning tree\\
% 		Modularity                &Clustering modularity\\
% 		Nodedegree                & Return node degree of each vertex\\
% 		Nodeneighborstats         & Node neighbour topological properties\\
% 		Participationcoef         & Compute participation coefficient\\
% 		Radiality                 & Compute vertex radiality centrality\\
% 		\hline
% 	\end{tabular}
% }
% \end{table}


%% Results

%% Discussion

%% Conclusion

%% Supporting Information

%%
% \bibliography{bibliography}
%