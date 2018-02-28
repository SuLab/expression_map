library(WGCNA)
library(topGO)

setwd('/gpfs/group/su/lhgioia/map/')

wgcna.results <- readRDS('results/allen/wgcna/consensus_phenotype_clustering.RDS')
gene.names <- readRDS('results/allen/variable_genes_log_4k.RDS')
gene.names <- gsub('X','',gene.names)

all.genes <- read.table('data/allen/tmp/gene_list.txt')
all.genes <- all.genes$V1

geneId.to.go <- readMappings('data/go/entrez_to_go.tsv')

all.results <- list()
var.results <- list()

for(col in unique(wgcna.results$merged.colors)){

    gene.slice <- gene.names[wgcna.results$merged.colors == col]
    
    gene.list <- factor(as.integer(gene.names %in% gene.slice))
    names(gene.list) <- gene.names
    
    go.obj <- new('topGOdata',
                  ontology = 'BP',
                  allGenes = gene.list,
                  annot = annFUN.gene2GO,
                  gene2GO = geneId.to.go)

    
    fisher.test <- new('elimCount',testStatistic = GOFisherTest)
    result.fisher <- getSigGroups(go.obj,fisher.test)
    var.results[[col]] <- result.fisher
    

    gene.list <- factor(as.integer(all.genes %in% gene.slice))
    names(gene.list) <- sapply(all.genes,toString)

    go.obj <- new('topGOdata',
                  ontology = 'BP',
                  allGenes = gene.list,
                  annot = annFUN.gene2GO,
                  gene2GO = geneId.to.go)

    
    fisher.test <- new('elimCount',testStatistic = GOFisherTest)
    result.fisher <- getSigGroups(go.obj,fisher.test)
    all.results[[col]] <- result.fisher


}
    


saveRDS(all.results,'results/allen/go/cons_pheno_all_bg_fisher_elim.RDS')
saveRDS(var.results,'results/allen/go/cons_pheno_var_bg_fisher_elim.RDS')
