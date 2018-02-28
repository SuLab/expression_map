library(topGO)

setwd('/gpfs/group/su/lhgioia/map/')

marker.results <- readRDS('results/allen/markers/filtered_markers_pairwise_wilvain.RDS')

var.genes <- c()

for(marker in 1:length(marker.results)){

    marker.genes <- names(marker.results[[marker]])
    marker.genes <- gsub('X','',marker.genes)

    marker.results[[marker]] <- marker.genes

    var.genes <- c(var.genes,marker.genes)

}

gene.names <- read.table('data/allen/gene_list.txt')$V1

## wgcna.results <- readRDS('results/allen/wgcna/consensus_phenotype_clustering.RDS')
## gene.names <- readRDS('results/allen/variable_genes_log_4k.RDS')
## gene.names <- gsub('X','',gene.names)

## all.genes <- read.table('data/allen/tmp/gene_list.txt')
## all.genes <- all.genes$V1

geneId.to.go <- readMappings('data/go/entrez_to_go.tsv')

all.results <- list()
var.results <- list()

for(marker in 1:length(marker.results)){

    gene.slice <- marker.results[[marker]]
    
    gene.list <- factor(as.integer(var.genes %in% gene.slice))
    names(gene.list) <- var.genes
    
    go.obj <- new('topGOdata',
                  ontology = 'BP',
                  allGenes = gene.list,
                  annot = annFUN.gene2GO,
                  gene2GO = geneId.to.go)

    
    fisher.test <- new('elimCount',testStatistic = GOFisherTest)
    result.fisher <- getSigGroups(go.obj,fisher.test)
    var.results[[marker]] <- result.fisher
    

    gene.list <- factor(as.integer(gene.names %in% gene.slice))
    names(gene.list) <- sapply(gene.names,toString)

    go.obj <- new('topGOdata',
                  ontology = 'BP',
                  allGenes = gene.list,
                  annot = annFUN.gene2GO,
                  gene2GO = geneId.to.go)

    
    fisher.test <- new('elimCount',testStatistic = GOFisherTest)
    result.fisher <- getSigGroups(go.obj,fisher.test)
    all.results[[marker]] <- result.fisher


}
    


saveRDS(all.results,'results/allen/go/wilvain_all_bg_fisher_elim.RDS')
saveRDS(var.results,'results/allen/go/wilvain_var_bg_fisher_elim.RDS')
