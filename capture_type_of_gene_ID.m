function gene_ID_type = capture_type_of_gene_ID(GEO_number)
    idTypes = {
	      'AFFYMETRIX_3PRIME_IVT_ID'
	      'AFFYMETRIX_EXON_GENE_ID'
	      'AFFYMETRIX_SNP_ID'
	      'AGILENT_CHIP_ID'
	      'AGILENT_ID'
	      'AGILENT_OLIGO_ID'
	      'ENSEMBL_GENE_ID'
	      'ENSEMBL_TRANSCRIPT_ID'
	      'ENTREZ_GENE_ID'
	      'FLYBASE_GENE_ID'
	      'FLYBASE_TRANSCRIPT_ID'
	      'GENBANK_ACCESSION'
	      'GENOMIC_GI_ACCESSION'
	      'GENPEPT_ACCESSION'
	      'ILLUMINA_ID'
	      'IPI_ID'
	      'MGI_ID'
	      'PFAM_ID'
	      'PIR_ID'
	      'PROTEIN_GI_ACCESSION'
	      'REFSEQ_GENOMIC'
	      'REFSEQ_MRNA'
	      'REFSEQ_PROTEIN'
	      'REFSEQ_RNA'
	      'RGD_ID'
	      'SGD_ID'
	      'TAIR_ID'
	      'UCSC_GENE_ID'
	      'UNIGENE'
	      'UNIPROT_ACCESSION'
	      'UNIPROT_ID'
	      'UNIREF100_ID'
	      'WORMBASE_GENE_ID'
	      'WORMPEP_ID'
	      'ZFIN_ID'};
    

    fprintf('\n\n');
    
    for indx = 1:size(idTypes,1)
      display([num2str(indx) ': ' idTypes{indx} '']);
    end
    
    prompt = ['Enter the type of gene ID used in the study associated to ' GEO_number ' (e.g., enter 9 for ENSEMBL_GENE_ID): '];
    index_of_gene_ID_type = input(['\n\n' prompt]);
    
    gene_ID_type = idTypes{1};
    if(isnumeric(index_of_gene_ID_type) & index_of_gene_ID_type > 0 & index_of_gene_ID_type <= size(idTypes,1))
      gene_ID_type = idTypes{index_of_gene_ID_type};
    end
end
