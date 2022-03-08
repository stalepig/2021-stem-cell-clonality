setwd("/Volumes/SHAQ/2019/01_15_19_DSS_single_cell_post/monocle3/")

library(monocle3)
library(ggplot2)
library(magrittr)
library(reshape2)
library(ggthemes)
library(gplots)
library(RColorBrewer)
source("camMonocleHelper.R")

# load(file = "cds_full.rds")
# d <- read.csv(file = "table_stem_repro_dynamics.csv",header = T)

## computes abundance of FPC signature in each sample in each cell type
# cds.full.subset <- cds.full[map.gene.short.names(cds = cds.full,gene.short.names = d$gene_short_name),]
# colData(cds.full.subset)$condition <- with(colData(cds.full.subset),paste(annotation,sample.name,sep = "-"))
# agg <- aggregate.expression.by.factor(cds = cds.full.subset,grouping = "condition",do.print = T)

## divides signature into genes that changed and those that didn't
# summ <- data.frame(agg)
# summ$id <- row.names(agg)
# summ10 <- melt(data = summ,id.vars = "id",variable.name = "group",value.name = "abundance")
# summ10$sample <- sapply(X = as.character(summ10$group),FUN = function(txt) {
#   retlist <- strsplit(x = txt,split = "\\.")[[1]]
#   return(retlist[length(retlist)])
# })
# summ10$annotation <- sapply(X = as.character(summ10$group),FUN = function(txt) {
#   retlist <- strsplit(x = txt,split = "\\.")[[1]]
#   return(paste(retlist[1:(length(retlist)-1)],sep = " ",collapse = " "))
# })
# summa <- split(x = summ10,f = list(summ10$id,summ10$annotation),drop = T) %>%
#   lapply(FUN = function(df) {
#     id <- df[1,"id"]
#     annotation <- df[1,"annotation"]
#     minval <- df$abundance[order(df$abundance)][2]
#     maxval <- max(df$abundance)
#     return(data.frame(id = id, annotation = annotation, minval = minval, maxval = maxval))
#   }) %>% do.call(what = "rbind",args = .)
# summa10 <- subset(summa,summa$maxval > 0.05)
# summa10$log2ratio <- log2(summa10$maxval / summa10$minval)
# summa10$upper <- ifelse(test = summa10$log2ratio > 1.18,yes = "yes",no = "no")
# summa20 <- merge(x = rowData(cds.full.subset),y = summa10,by = "id",all = F) %>% data.frame
# summa30 <- summa20[order(summa20$upper,summa20$annotation,summa20$log2ratio,decreasing = T),]

## writes the data table
# write.table(x = summa30,file = "table_convergence_single_genes.csv",sep = ",",row.names = F,quote = F)

## sees what is unique to epithelium
# summa40 <- summa30[!duplicated(summa30$gene_short_name),]

## loads the enrichr data and makes the plot
# file.list <- list.files(path = "enrichr/",pattern = "txt",all.files = F)
# da <- lapply(X = file.list,FUN = function(theFile) {
#   result <- strsplit(x = theFile,split = "\\_")[[1]]
#   annotation <- paste(result[1:(length(result)-2)],sep = " ",collapse = " ")
#   database <- tail(x = result,n = 1)
#   thePath <- paste(getwd(),"/",theFile,sep="")
#   fileContents <- read.table(file = paste(getwd(),"/enrichr/",theFile,sep=""),sep = "\t",header = T)
#   fileContents$annotation <- rep(x = annotation,times = dim(fileContents)[[1]])
#   fileContents$database <- rep(x = database,times = dim(fileContents)[[1]])
#   return(fileContents)
# }) %>% do.call(what = "rbind",args = .)

## NA handling
## if odds ratio or Combined Score, then each NA = 0
## if p or adjusted p value, then each NA = 1
# da10 <- dcast(data = da,formula = annotation~Term,value.var = "Combined.Score",fun.aggregate = median,fill = 0)
# m <- data.matrix(da10[,-1])
# row.names(m) <- da10$annotation

## filtering to make the heatmap more readable
m05 <- t(apply(X = m,MARGIN = 1,FUN = function(x) { x/sum(x) }))
m10 <- m05[,which(colSums(m05) > 0.1)]
colnames(m10)[grepl(pattern = "cytokine and cytokine receptor",x = colnames(m10))] <- "Cytokine and cytokine receptor"
my.palette <- colorRampPalette(colors = c("black","magenta"))(n = 300)
heatmap.2(x = m10,trace = "none",margins = c(15,10),scale = "row",col = my.palette,tracecol = NA)