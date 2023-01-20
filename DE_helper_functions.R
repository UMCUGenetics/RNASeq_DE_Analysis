#  install.packages("tidyverse")
#  install.packages("purrr")
#  install.packages("generics")
#  install.packages("vctrs")
#  install.packages("clusterProfiler")
# install.packages("plotly")
# install.packages("DT")
# BiocManager::install("clusterProfiler")

library("biomaRt")#, lib=custom_libs)
library("purrr")#, lib=custom_libs)
library("dplyr")#, lib=custom_libs)
library("DESeq2")#, lib=custom_libs)
library("pheatmap")#, lib=custom_libs)
library("RColorBrewer")#, lib=custom_libs)
library("ggplot2")#, lib=custom_libs)
library("ggrepel")#, lib=custom_libs)

library(plotly)
library(DT)

#pathway analysis libs (gage, etc)
library(AnnotationDbi)
library(BiocManager)
library(pathview)
library(gage)
library(gageData)

#goterm analysis (goseq)
library(goseq)
library(geneLenDataBase)
library(GO.db)

library(clusterProfiler)

setCores <- function(nrcores=2){
  #set appropriate multicore using path as indicator for unix/windows
  library("BiocParallel")
  if(.Platform$OS.type == "unix") {
    register(MulticoreParam(nrcores))  
  } else{
    register(SnowParam(nrcores))
  }
}

save_pheatmap_pdf <- function(x, filename, width=10, height=10) {
       stopifnot(!missing(x))
       stopifnot(!missing(filename))
       pdf(filename, width=width, height=height)
       grid::grid.newpage()
       grid::grid.draw(x$gtable)
       dev.off()
}

save_pdf <- function(x, filename, width=10, height=10) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x)
  dev.off()
}


save_svg <- function(x, filename, width=10, height=10) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  svg(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x)
  dev.off()
}

#initializes for use of specific annotation database. Returns a list object with anno_db, gene_id_column and gene_name_column objects.
init_annotation_db <- function (use_biomart=T, species=NULL, annotation_data=NULL){
  anno_db=NULL
  kegg.genesets=NULL
  go.seq.genome=NULL
  gene_id_column=NULL
  gene_name_column=NULL
  
  #load annotation data if specified
  if( !is.null(annotation_data) &&annotation_data!=""){
    load(annotation_data)
    anno_db = anno_db_data$anno_db
    gene_id_column = anno_db_data$gene_id_column
    gene_name_column = anno_db_data$gene_name_column
  }
  
  if(is.null(anno_db)){
    #general init
    if(use_biomart){
      gene_id_column = "ensembl_gene_id"
      gene_name_column = "external_gene_name"
      library("biomaRt")
      
    } else {
      gene_id_column = "ENSEMBL"
      gene_name_column= "SYMOL"
    }
  }
  
  #genome specific inits
  if(species=="human_grch37"){
    go.seq.genome="hg19"
    data(kegg.sets.hs)
    kegg.genesets=kegg.sets.hs
    
    if(is.null(anno_db)){
      if(use_biomart){
        grch37 = useEnsembl(biomart="ensembl",GRCh=37)
        anno_db = useDataset("hsapiens_gene_ensembl", grch37)
      }
      else{
        library("org.Hs.eg.db")
        anno_db <-org.Hs.eg.db
      }
    }
  } else if(species=="human_grch38"){
    data(kegg.sets.hs)
    kegg.genesets=kegg.sets.hs
    go.seq.genome="hg38"
    
    if(is.null(anno_db)){
      if(use_biomart){
        anno_db = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
      }
    }
  } else if(species=="mouse_grcm38") {
    #mm10
    data(kegg.sets.mm)
    kegg.genesets=kegg.sets.mm
    go.seq.genome="mm10"
    
    if(is.null(anno_db)){
      if(use_biomart){
        mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl",host="http://aug2020.archive.ensembl.org")
        anno_db = useDataset("mmusculus_gene_ensembl", mart)
      } else {
        library("org.Mm.eg.db")
        anno_db <- org.Mm.eg.db 
      }
    }
    
  } else if(species=="mouse_grcm39"){
    #mm11
    data(kegg.sets.mm)
    kegg.genesets=kegg.sets.mm
    go.seq.genome="mm39"
    
    if(is.null(anno_db)){
      if(use_biomart){
        anno_db = useDataset("mmusculus_gene_ensembl", useMart("ensembl"))  
      }
    }
  } else if(species=="dog_canfam3.1"){
    
    if(is.null(anno_db)){
      if(use_biomart){
        anno_db = useDataset("clfamiliaris_gene_ensembl", useMart("ensembl"))
      }
    }
    
  } else if(species=="cat_gFelis_catus_9.0"){
    
    if(is.null(anno_db)){
      if(use_biomart){
        anno_db = useDataset("fcatus_gene_ensembl", useMart("ensembl"))
      }
    }
  }
  
  result <- list("gene_id_column"=gene_id_column, "gene_name_column"=gene_name_column, "anno_db"=anno_db, go.seq.genome=go.seq.genome, kegg.genesets=kegg.genesets)
  
  return(result)  
}


# read in metadata file
readMetadataFile <- function(metafile= NULL){
  md <-
    data.frame(read.table(
      metafile,
      sep = '\t',
      header = T,
      comment.char = ""
    ))
  
  return (md)
}

# read counts file (tab delimited). exclude is a Vector containing the column names to exclude from the data.frame
read_in_feature_counts_custom <- function(file, exclude = c("Chr", "Start", "End", "Strand", "Length", "gene_name")){
  cnt<- read.table(file, sep="\t", header=T, comment.char="#")
  
  if(!is.null(exclude)){
    for (excl_col in exclude){
      if(excl_col %in% colnames(cnt)){
        cnt<- cnt %>% dplyr::select(-all_of(excl_col))
      }
    }
  }
  return(cnt)
}

removeColumnsFromCountTable <- function(count.table, remove.columns= c("Chr", "Start", "End", "Strand", "Length", "gene_name")){
  if(!is.null(remove.columns)){
    for (excl_col in remove.columns){
      if(excl_col %in% colnames(count.table)){
        count.table<- count.table %>% dplyr::select(-all_of(excl_col))
      }
    }
  }
  return(count.table)
}

readCountsFile <- function(countfile = NULL, sep='\t'){
  rawcounts <-
    data.frame(read.table(
      countfile,
      sep = sep,
      header = T,
      row.names = 1
    ))
  return(rawcounts)
}

#Correct sampleids for proper processing in R, samples starting with a digit get an 'X' prepended. Also '-' are replaced with a '.'
correctIds <- function(x){
  x = as.character(x)
  if(is.list(x) || is.vector(x)){
    for(i in 1:length(x)){
      c = substring(x[i], 1, 1)
      if(c>="0" && c<="9"){
        x[i] = paste0("X", x[i])
      }
    }
  }else{
    c = substring(x, 1, 1)
    if(c>="0" && c<="9"){
      x = paste0("X", x)
    }
  }
  gsub("\\-", ".", x)
}

generatePCARepel <- function(pcadata,
                             color = "condition",
                             shape = NULL,
                             title = NULL
) {
  percentPCA = round(100*attr(pcadata, "percentVar"))
  
  if(is.null(shape)){
    print(
      ggplot(pcadata, aes(x=PC1, y=PC2, label = name)) + 
        geom_point(aes(color = !!sym(color)), size = 6, alpha = .8) + 
        geom_text_repel(size=5) +
        ggtitle(paste0("PCA ", title)) + 
        xlab(paste0("PC1: ", percentPCA[1],"%variance")) +
        ylab(paste0("PC2: ", percentPCA[2],"%variance")) +
        theme(plot.title = element_text(hjust = 0.5, face ="bold"), legend.title = element_text(size=14), legend.text=element_text(size=14))
    )
  }
  else{
    print(
      ggplot(pcadata, aes(x=PC1, y=PC2, label = name)) + 
        geom_point(aes(shape = !!sym(shape), color = !!sym(color)), size = 6, alpha = .8) + 
        geom_text_repel(size=5) +
        ggtitle(paste0("PCA ", title)) + 
        xlab(paste0("PC1: ", percentPCA[1],"%variance")) +
        ylab(paste0("PC2: ", percentPCA[2],"%variance")) +
        theme(plot.title = element_text(hjust = 0.5, face ="bold"), legend.title = element_text(size=14), legend.text=element_text(size=14))
    )
  }
}

generatePCARepelDeprecated <- function(pcadata,
                                       color = "condition",
                                       shape = NULL,
                                       title = NULL
) {
  if(is.null(shape)){
    print(
      ggplot(pcadata, aes(x=PC1, y=PC2, label = name)) + 
        geom_point(aes_string(color = color), size = 10, alpha = .8) + 
        #geom_text(aes(label = name), vjust = 2, size = 2) + 
        geom_text_repel(size=3) +
        ggtitle(paste0("DE-PCA ", title)) + 
        theme(plot.title = element_text(hjust = 0.5, face ="bold"))
    )
  }
  else{
    print(
      ggplot(pcadata, aes(x=PC1, y=PC2, label = name)) + 
        geom_point(aes_string(shape = shape, color = color), size = 10, alpha = .8) + 
        #geom_text(aes(label = name), vjust = 2, size = 2) + 
        geom_text_repel(size=3) +
        ggtitle(paste0("DE-PCA ", title)) + 
        theme(plot.title = element_text(hjust = 0.5, face ="bold"))
    )
  }
}

generateHeatmap <- function(dds, title=NULL, rowlabeling="condition", filename=NA) {
  #Heatmap
  rld <- vst(dds, blind=FALSE)
  sampleDists <- dist(t(assay(rld)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(rownames(sampleDistMatrix), rld[[rowlabeling]], sep = " | ")
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors,
           main = title,
           fontsize=14,
           filename=filename)
}

getDataForContrast <- function(dds=NULL, contrast=NULL, fcval=1.5, fdrval=0.05, anno_db=NULL) {
  #Generate comparision of interest
  title <- paste0("DESeq2_Volcano_",contrast[1],"_",contrast[2],"_versus_",contrast[3])
  
  #Filter hits according
  #results_contrast are all results for contrast
  results_contrast <- results(dds, contrast = contrast)
  
  #Filter adjusted p-val and L2FC threshold
  results_contrast$filter <- 0
  results_contrast$filter[results_contrast$padj < fdrval & results_contrast$log2FoldChange > log2(fcval)] <- 1
  results_contrast$filter[results_contrast$padj < fdrval & results_contrast$log2FoldChange < -log2(fcval)] <- -1
  results_contrast$filter <- as.factor(results_contrast$filter)
  #filtered only contains the results matching tresholds
  filtered <- results_contrast
  
  #Annotate gene names
  if (!is.null(anno_db)) {
    select <- AnnotationDbi::select
    
    genesymbols <- select(anno_db,  rownames(rawcounts), c(gene_id_column, gene_name_column), gene_id_column)
    symbols <- genesymbols[!duplicated(genesymbols[[gene_name_column]]),]
    symbols.subset <- symbols[symbols[[gene_id_column]] %in% rownames(filtered),]
    rownames(symbols.subset) <- symbols.subset[[gene_id_column]]
    symbols.subset.2 <- subset(symbols.subset,select=c(gene_name_column))
    final <- merge(symbols.subset.2 ,as(filtered, "data.frame"), by='row.names', all=T)
    rownames(final) <- final[,1]
    # Write Filtered genes to output file
  } else {
    final <- filtered
  } 
  final <- final[order(final$padj),]
  
  #write only genes passing filter
  filtered = final[final$filter!=0,]
  filtered <- subset(filtered, select=-c( `filter`))
  write.csv(data.frame("GENEID"=rownames(filtered),filtered), file = paste0(contrast[1],"_",contrast[2],"_vs_",contrast[3],"_padj_",fdrval,".csv"), quote = F, row.names = F)
  
  #write full set
  final <- subset(final, select=-c(`filter`))
  write.csv(data.frame("GENEID"=rownames(final),final), file = paste0(contrast[1],"_",contrast[2],"_vs_",contrast[3],".csv"), quote = F, row.names = F)
  
  return(results_contrast) 
}

plotVolcano <- function(results_contrast=NULL, contrast=NULL, fcval=1.5, fdrval=0.05, fc.max=6){
  # -------- QC figures --------------#
  par(mfrow=c(1,1))
  with(results_contrast, plot(log2FoldChange, -log10(pvalue), pch=20, main=paste0("Volcano-Plot. ",contrast[1],": ",contrast[2],"_versus_",contrast[3]),xlim = c(-1*fc.max,fc.max)))
  with(subset(results_contrast, padj<fdrval & log2FoldChange < -log2(fcval)), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
  with(subset(results_contrast, padj<fdrval & log2FoldChange > log2(fcval)), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
  abline(v = -log2(fcval), col="blue", lwd=3, lty=2)
  abline(v = log2(fcval), col="red", lwd=3, lty=2)
  legend("bottomright", xjust=1, yjust=1, legend=c("Up","Down"), pch=20, col=c("red","blue"))
}

plotMA <- function(results_contrast=NULL, contrast=NULL,fdrval=0.05, fc.max=2){
  DESeq2::plotMA(results_contrast, ylim = c(-1*fc.max,fc.max), alpha = fdrval, main=paste0("MA-plot. ",contrast[1],": ",contrast[2],"_versus_",contrast[3]))
}

#This function translates ensemble ids to entrez ids which are needed for kegg pathway analysis. This can be done either using a biomart object or using annodb...
ensemblToEntrezIds <- function(ensembl.ids=NULL, biomart.db=NULL, anno.db=NULL){
  entrezIds = NULL
  
  if(!is.null(biomart.db)){
    #retrieve entrezIds for all ensembl gene ids using biomart
    entrezIds =  getBM(filters="ensembl_gene_id", attributes=c("ensembl_gene_id","entrezgene_id"), values=ensembl.ids, mart=biomart.db)       
    
    #concatenate multiple entrezIds for an ensembl geneid using |
    entrezIdsDF = entrezIds %>% 
      group_by(ensembl_gene_id) %>%
      summarize(across(entrezgene_id, ~ifelse(all(is.na(.x)), NA, paste0(na.omit(.x), collapse = "|"))), .groups = "drop") %>%
      arrange(factor(ensembl_gene_id,levels=row.names(results_contrast)))
    
    entrezIdsDF = entrezIdsDF[!duplicated(entrezIdsDF$ensembl_gene_id),]
    
    
  } else {
    # Alternative method using annodb, retrieves less kegg ids? OFten used?
    #library(org.Hs.eg.db)
    #entrezIdsDF = mapIds(org.Hs.eg.db,
    entrezIdsDF = mapIds(anno.db,
                         keys=ensembl.ids, 
                         column="ENTREZID",
                         keytype="ENSEMBL",
                         multiVals="first")
  }
  return(entrezIdsDF)
}

#instead of printing a pathway plot to just a png file, this will also plot it in the Rmd environment...
plot_pathview <- function(genelist, pathwayid, species="hsa", out.suffix = "genesinfo", save_image=F){
  
  grid::grid.newpage()
  
  viewPath <- pathview(gene.data  = genelist,
                       pathway.id = pathwayid,
                       species    = species,
                       out.suffix = out.suffix, kegg.native = T, same.layer = F)
  filename=paste(pathwayid, ".", out.suffix, ".png", sep = "")
  img <- png::readPNG(filename)
  if(!save_image) invisible(file.remove(filename))
  
  grid::grid.raster(img)
}

plot_pathviews <- function(genelist, pathwayids, species="hsa", out.suffix = "genesinfo", save_image=F){
  for(pathwayid in pathwayids){
    plot_pathview(genelist, pathwayids, species, out.suffix, save_image)
  }
}

filter.gage.pathway.results <- function(fc.kegg.p, cutoff.field="q.val", cutoff.val=0.05){
  sel <- fc.kegg.p$greater[, cutoff.field] < cutoff.val & !is.na(fc.kegg.p$greater[, cutoff.field])
  path.ids.g <- rownames(fc.kegg.p$greater)[sel]
  #only used when same.dir is set to TRUE
  sel.l <- fc.kegg.p$less[, cutoff.field] < cutoff.val & !is.na(fc.kegg.p$less[,cutoff.field])
  path.ids.l <- rownames(fc.kegg.p$less)[sel.l]
  #return combined results
  path.ids <- substr(c(path.ids.g, path.ids.l), 1, 8)
  
  return(path.ids)
  
}

#example: https://github.com/plagnollab/RNASeq_pipeline/blob/master/goseq_script.R
#info: https://support.bioconductor.org/p/68076/
go.analysis.goseq.init <- function(){
  genome=go.seq.genome
  geneid.source="ensGene"
  bias.data=NULL
  test.categories=NULL
  enriched.file=NULL
  plot.fit=TRUE
}


process.go.results <- function(go.data=NULL, p.value=0.05, display.top.n=10, enriched.file=NULL){
  if(is.null(go.data)){  return()  }
  
  #get only statistically significant ones
  enriched.go = go.data$category[go.data$over_represented_pvalue<0.05]
  
  #print enriched results details if file is specified
  if(!is.null(enriched.file) && length(enriched.go)>0){
    #save relevant info for top x 
    #library(GO.db)
    print.count=min(250, length(enriched.go))
    capture.output(
      for(go in enriched.go[1:print.count]){
        print(GOTERM[[go]])
        cat("-------------------------\n")
      }
      , file = enriched.file)
  }
  #top_n(10, wt=-over_represented_pvalue) %>% 
  top.hits = go.data %>% 
    top_n(n=display.top.n, wt=-over_represented_pvalue) %>% 
    mutate(hitsPerc=numDEInCat*100/numInCat)
  
  if(length(top.hits) > display.top.n){
    cat("  **No overexpressed GO terms found..**\n\n")
  } else{
    print(ggplot(top.hits, aes(x=hitsPerc, 
                               y=term, 
                               colour=over_represented_pvalue, 
                               size=numDEInCat)) +
            geom_point() +
            expand_limits(x=0) +
            labs(x="Hits (%)", y="GO term", colour="p value", size="Count"))
  }
  
}

go.analysis.goseq <- function(de.genes=NULL, all.genes=NULL, genome=NULL, geneid.source=NULL, test.categories=NULL,  bias.data=NULL, plot.fit=TRUE){
  
  #convert to vectors
  de.genes.vector = c(t(de.genes))
  all.genes.vector = c(t(all.genes))
  
  #construct genevector with de genes marked
  gene.vector = as.integer(all.genes.vector%in%de.genes.vector)
  names(gene.vector)=all.genes.vector
  
  #weigh gene vector by length of our genes
  pwf=nullp(gene.vector, genome, geneid.source, bias.data, plot.fit)
  
  #Find enriched go terms
  #test.cast = combination of "GO:CC", "GO:BP", "GO:MF" & "KEGG". Default: "GO:CC", "GO:BP", "GO:MF"
  if(is.null(test.categories)){
    test.categories=c("GO:CC", "GO:BP", "GO:MF")
  }
  go.wall=goseq(pwf=pwf, genome=genome, id=geneid.source, test.cats=test.categories)
  
  return(go.wall)
}

datatable_custom <- function(data, row.names=T, round.cols=NULL, signif.cols=NULL, digits=3){
  x=datatable(data, rownames = row.names, 
              extensions = 'Buttons', options = list (
                dom = 'Bftrip',
                autoWidth = TRUE,
                buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                fnDrawCallback = JS('function(){ HTMLWidgets.staticRender(); }'),
                deferRender=TRUE
              ),
              filter='top') %>% 
    formatRound(columns=round.cols, digits=digits) %>%  
    formatSignif(columns=signif.cols, digits = digits)
  
  return(x)
  
}
