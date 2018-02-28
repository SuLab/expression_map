library(WGCNA)
library(topGO)

setwd('/gpfs/group/su/lhgioia/map/')

sep.results <- readRDS('results/allen/wgcna/merged_results_separate_phenotype.RDS')
gene.names <- readRDS('results/allen/variable_genes_log_4k.RDS')
gene.names <- gsub('X','',gene.names)

all.genes <- read.table('data/allen/tmp/gene_list.txt')
all.genes <- all.genes$V1

geneId.to.go <- readMappings('data/go/entrez_to_go.tsv')

all.results <- list()

var.results <- list()

for(reg in names(sep.results)){

    all.results[[reg]] <- list()
    var.results[[reg]] <- list()

    for(col in unique(sep.results[[reg]]$colors)){

        gene.slice <- gene.names[sep.results[[reg]]$colors == col]
        
        gene.list <- factor(as.integer(gene.names %in% gene.slice))
        names(gene.list) <- gene.names
        
        go.obj <- new('topGOdata',
                      ontology = 'BP',
                      allGenes = gene.list,
                      annot = annFUN.gene2GO,
                      gene2GO = geneId.to.go)

        
        fisher.test <- new('elimCount',testStatistic = GOFisherTest)
        result.fisher <- getSigGroups(go.obj,fisher.test)
        var.results[[reg]][[col]] <- result.fisher
        

        gene.list <- factor(as.integer(all.genes %in% gene.slice))
        names(gene.list) <- sapply(all.genes,toString)

        go.obj <- new('topGOdata',
                      ontology = 'BP',
                      allGenes = gene.list,
                      annot = annFUN.gene2GO,
                      gene2GO = geneId.to.go)

        
        fisher.test <- new('elimCount',testStatistic = GOFisherTest)
        result.fisher <- getSigGroups(go.obj,fisher.test)
        all.results[[reg]][[col]] <- result.fisher


    }
    
}

saveRDS(all.results,'results/allen/go/sep_pheno_all_bg_fisher_elim.RDS')
saveRDS(var.results,'results/allen/go/sep_pheno_var_bg_fisher_elim.RDS')
