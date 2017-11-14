#### Merge phenotype, sex, and age info ####
## read phenotype data
gtex_pheno <- read.table("/gpfs/group/su/lhgioia/map/data/recount/gtex/pheno_formatted_withKey.tsv", stringsAsFactors = F, sep = "\t", quote = "", header=T)
names(gtex_pheno)[21] <- "sample_key"

## read sex & age data
gtex_sexAge <- read.table("/gpfs/group/su/lhgioia/map/data/recount/gtex_sex_age.tsv", stringsAsFactors = F, header=T)
names(gtex_sexAge) <- c("sample_key", "sex", "age")

## merge by common id
gtex <- merge(gtex_pheno, gtex_sexAge, by="sample_key")

## remove old rows
gtex$sample_key <- NULL
gtex$age.x <- NULL
gtex$sex.x <- NULL

## reset column order
names(gtex)[19:20] <- c("sex", "age")
gtex <- gtex[, match(names(gtex_pheno[,1:20]), names(gtex))]

## write output
write.table(gtex, "/gpfs/group/su/lhgioia/map/data/recount/gtex/pheno_final.tsv", row.names = F, quote = F, sep = "\t")
