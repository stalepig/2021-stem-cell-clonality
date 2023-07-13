library(monocle3)
library(magrittr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(enrichR)
library(pheatmap)

assign.clusters <- function(cds,assignments,column = 2,new.col.name="annotation") {
  assigner <- assignments[,column]
  names(assigner) <- as.character(assignments[,1])
  colData(cds)[,new.col.name] <- sapply(X = clusters(cds),FUN = function(x) {
    assigner[as.character(x)]
  })
  return(cds)
}

recluster.cds <- function(cds,cluster.method="leiden") {
  cds <- preprocess_cds(cds = cds,num_dim = 100,verbose = T)
  cds <- reduce_dimension(cds = cds,verbose = T)
  cds <- cluster_cells(cds = cds,verbose = T,cluster_method = cluster.method)
  cds <- learn_graph(cds = cds,verbose = T)
  return(cds)
}

filter.marker.results <- function(marker_test_res,num.to.keep=3,criterion="pseudo_R2") {
  spl <- split(x = marker_test_res,f = list(marker_test_res$cell_group))
  top_specific_markers <- do.call("rbind",lapply(X = spl,FUN = function(df) {
    df10 <- df[order(df[,criterion],decreasing = T),]
    df10[1:num.to.keep,]
  }))
  top_specific_markers_ids <- unique(top_specific_markers$gene_id)
  return(top_specific_markers_ids)
}

add.gene.sum.column <- function(cds,genes.of.interest,col.name="gene.sum",use.raw=F,use.norm=T) {
  feats <- rowData(x = cds) %>% data.frame %>% filter(gene_short_name %in% genes.of.interest)
  if (use.norm) {
    exprs <- normalized_counts(cds = cds[feats$id],norm_method = "log",pseudocount = 1)  
  } else {
    exprs <- exprs(cds[feats$id,])
  }
  
  if (use.raw) {
    exprs.norm <- exprs
  } else {
    exprs.norm <- t(apply(X = exprs,MARGIN = 1,FUN = function(x) {
      x / max(x)
    }))
  }
  exprs.sum <- apply(X = exprs.norm,MARGIN = 2,FUN = sum)
  colData(cds)[,col.name] <- exprs.sum
  return(cds)
}

find.closest.row <- function(input.matrix,center.point=c(0,0)) {
  ss <- apply(X = input.matrix,MARGIN = 1,FUN = function(x) {
    return((center.point[1]-x[1])^2 + (center.point[2]-x[2])^2)
  })
  return(which.min(ss))
}

aggregate.expression.by.factor <- function(cds,grouping,do.print=F,drop.zeros=F,pct = 1,use.norm = T) {
  group.call <- as.factor(colData(cds)[,grouping])
  names(group.call) <- row.names(colData(cds))
  if (use.norm) {
    exprs.subset <- normalized_counts(cds = cds,norm_method = "log",pseudocount = 1)
  } else {
    exprs.subset <- exprs(cds)  
  }
  aggm <- matrix(nrow = dim(exprs.subset)[1],ncol = length(levels(group.call)))
  for (i in c(1:dim(exprs.subset)[1])) {
    if (do.print) {
      print(i)
    }
    theVec <- as.vector(exprs.subset[i,])
    names(theVec) <- colnames(exprs.subset)

    if (drop.zeros) {
      aggvec <- tapply(X = theVec,INDEX = group.call,FUN = function(x) {
        x.trunc <- subset(x = x,subset = x>0)
        if (length(x.trunc) > pct * 0.01 * length(x)) {
          return(mean(x.trunc))
        } else {
          return(0)
        }
      })
    } else {
      aggvec <- tapply(X = theVec,INDEX = group.call,FUN = mean)
    }
    aggm[i,] <- aggvec
  }
  row.names(aggm) <- row.names(exprs.subset)
  colnames(aggm) <- levels(group.call)
  return(aggm)
}

map.gene.short.names <- function(cds,gene.short.names) {
  feats <- rowData(cds) %>% data.frame
  feats.subset <- filter(feats,gene_short_name %in% gene.short.names)
  return(feats.subset$id)
} 

plot.expression.with.dots <- function(cds.subset,grouping,yjitter = 0.1,logticks=F,pseudocount=0.0,use.norm=T,subsample=F,dotsize=0.5) {
  if (use.norm) {
    expr.df <- cbind(colData(cds.subset),data.frame(t(normalized_counts(cds.subset,"log",1) %>% as.matrix))) %>% data.frame
  } else {
    expr.df <- cbind(colData(cds.subset),data.frame(t(exprs(cds.subset) %>% as.matrix))) %>% data.frame
  }
  startcol <- dim(colData(cds.subset))[[2]]+1
  endcol <- dim(expr.df)[[2]]
  colnames(expr.df)[startcol:endcol] <- rowData(cds.subset)$gene_short_name
  expr.df10 <- melt(data = expr.df,measure.vars = c(startcol:endcol),variable.name = "gene",value.name = "expression")
  expr.df10$expression <- expr.df10$expression + pseudocount
  if (subsample) {
    max.size <- table(expr.df10[,grouping]) %>% min
    expr.df10 <- expr.df10 %>% group_by(gene) %>% group_by_(grouping) %>% sample_n(max.size)
  }
  if (logticks == T) {
    dplot <- ggplot(data = expr.df10,mapping = aes_string(x = grouping,y = "expression",colour = grouping)) +
      facet_wrap(facets = .~gene,scales = "free") +
      geom_point(position = position_jitter(height = yjitter),size = dotsize,alpha=0.7) +
      scale_y_log10() +
      theme_classic(base_size = 16) +
      theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))    
  } else {
    dplot <- ggplot(data = expr.df10,mapping = aes_string(x = grouping,y = "expression",colour = grouping)) +
      facet_wrap(facets = .~gene,scales = "free") +
      geom_point(position = position_jitter(height = yjitter),size = dotsize,alpha=0.7) +
      theme_classic(base_size = 16) +
      theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
  }
  return(dplot)
}

query.enrichr <- function(genes,grouping,dbs,limit.result = 100,q.value = 0.05) {
  setEnrichrSite(site = "Enrichr")
  indf <- data.frame(gene_short_name = genes,type = grouping)
  
  retdf <- data.frame()
  for (gset in unique(indf$type)) {
    goi <- subset(indf,indf$type==gset)$gene_short_name
    enr <- enrichr(genes = goi,databases = dbs)
    for (i in 1:length(enr)) {
      db.name <- names(enr)[i]
      df <- enr[[i]]
      if (dim(df)[[1]] > 0) {
        df$database <- db.name
        df$type <- gset
        holder <- subset(df,df$Adjusted.P.value < q.value)
        if (dim(holder)[[1]] > limit.result) {
          holder <- holder[1:limit.result,]
        }
        if (dim(holder)[[1]] > 0) {
          retdf <- rbind(retdf,holder)
        }
      }
    }
  }
  
  retdf10 <- dcast(data = retdf,formula = Term~type,value.var = "Adjusted.P.value",fun.aggregate = min,fill = 1)
  mat <- data.matrix(retdf10[,-1])
  row.names(mat) <- retdf10$Term
  mat10 <- -log10(mat)
  
  ph <- pheatmap(mat = mat10,cluster_cols = F)
  return(retdf)
}

subset.cds.features <- function(cds,gene.short.names) {
  goi <- map.gene.short.names(cds = cds,gene.short.names = gene.short.names)
  return(cds[goi,])
}

subset.cds.cells <- function(cds,column,elements) {
  barcodes <- colData(cds) %>% data.frame
  cds.sub <- cds[,subset(barcodes,barcodes[,column] %in% elements) %>% row.names]
  return(cds.sub)
}