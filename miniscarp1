#!/usr/bin/Rscript env
source("/hwfs2/RD/Pipes/Project/ST/pipeline/scRNA/SCRAP/v3.0/bin/base.R")
data.path <- "/hwfs2/RD/Pipes/Project/ST/pipeline/scRNA/SCRAP/v2.1/refData/cellCycle"
# parameter
suppressMessages(library(getopt))
spec <- matrix(c(
  'scoutdir'    , 'i', 2, 'character', 'SCRAP output directory',
  'od'          , 'o', 1, 'character', 'output dir name',
  'method'      , 'm', 2, 'character'  ,'Other remove batch method expect RPCA',
  'data'        , 'r', 2, 'character', 'scRNA_upload.Rds',
  'json'        , 'j', 2, 'character', 'json file',
  'help'        , 'h', 0, 'logical'  , 'print usage'
), byrow = TRUE, ncol = 5)
opt <- getopt(spec)

PrintUsage <- function(){
  cat("Usage:nobatch.r -i /output_of_SCRAP -m cca -o /output \n
  or nobatch.r -r /scRNA_upload.Rds -j /parameter.json -m cca-o /output \n")
  cat(paste(getopt(spec, usage = TRUE), sep = "\n"))
  q(status = 1)		
}

if(!is.null(opt$help)) {PrintUsage()}
if(is.null(opt$od)) {opt$od <- "./output"}
if(!dir.exists(opt$od)) {dir.create(opt$od)}

#判断并修改输入文件格式
if (!is.null(opt$scoutdir) && !is.null(opt$json) && !is.null(opt$data)) {
  stop("Error: Please provide either --data and --json parameters or --od parameter.")
}
if ((!is.null(opt$data) && is.null(opt$json)) || (is.null(opt$data) && !is.null(opt$json))) {
  stop("Error: Please provide both --data and --json parameters.")
}

if (is.null(opt$scoutdir)) {
    opt$data <- paste0(opt$data,'/result/Backup/04_singleSampleData/scRNA_upload.Rds')
    opt$json <- paste0(opt$data,'/JSON/parameter.json')
}
if (!is.null(opt$method)){othermethod<-opt$method}



time1 <- proc.time()
# load packages
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(RJSONIO))
suppressMessages(library(ggplot2))
suppressMessages(library(parallel))
suppressMessages(library(RColorBrewer))

##增加包用来批次矫正
suppressMessages(library(harmony))

kanchorval <- 5
# add 20220213
.Kparameter <- function(object){
  kVal = min(unlist(sapply(1:length(object), function(i) dim(object[[i]])[2])))
  kVal = plyr::round_any(kVal, 10, f = floor)
  kfilterVal = ifelse(kVal > 200, 200, kVal); cat("Using FindIntegrationAnchors(k.filter = ", kfilterVal, ") for downstream analysis.
\n")
  kweightVal = ifelse(kVal > 100, 100, ifelse(kVal/2 > 30, kVal, 30)); cat("Using IntegrateData(k.weight =", kweightVal, ") for downstream analysis.\n")
  p <- c(kfilterVal, kweightVal)
  return(p)
}

DealWithNormalized <- function(object, method = "SCT", removeCycle = "no"){
  if(method == "Normalize"){
    object <- NormalizeData(object, normalization.method =  "LogNormalize", scale.factor = 10000)
    object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = as.numeric(parameter_json$fea_n))
    if(removeCycle == "yes"){
      object <- ScaleData(object, vars.to.regress = c("CellCycle.diff"))
    }
    if(removeCycle == "no"){
      object <- ScaleData(object)
    }
  }
  if(method == "SCT"){
    if(removeCycle == "yes"){
      object <- SCTransform(object = object, do.scale = TRUE, verbose = FALSE, vars.to.regress = c("CellCycle.diff"))
    }
    if(removeCycle == "no"){
      object <- SCTransform(object = object, do.scale = TRUE, verbose = FALSE)
    }
  }
  return(object)
}

MergeByBatch <- function(object, batchs){
  listData <- list()
  uniq.batch <- unique(batchs)
  cat("use batch is: ", batchs, "\n")
  cat("unique batch is: ", uniq.batch, "\n")
  k <- 1
  for(i in 1:length(uniq.batch)){
    #listData[[i]] <- list()
    temp <- list()
    index <- sapply(1:length(object), function(j){
      #cat("current sample is: ", unique(object[[j]]$sample),"\n")
      #cat("sample batch is: ", batchs[unique(object[[j]]$sample)], "\n")
      #cat("current batch is", uniq.batch[i],"\n")
      if(batchs[as.character(unique(object[[j]]$sample))] == uniq.batch[i]){
        return(j)
      }else{
        return(NULL)
      }
    })
    if(is.null(unlist(index))){
      next
    }
    cat("merge index is: ", unlist(index), "\n")
    temp <- object[unlist(index)]
    if(length(temp) > 1){
      listData[[k]] <- merge(x = temp[[1]], y = temp[2:length(temp)])
      listData[[k]]$batch <- uniq.batch[i]
    }else{
      listData[[k]] <- temp[[1]]
      listData[[k]]$batch <- uniq.batch[i]
    }
    k <- k + 1
  }
  return(listData)
}

NotRemoveBatch <- function(singleList){
  #removeBatch == "no"
  merge.y <- function(object.list) {
    y = c(object.list[[2]])
    if (length(object.list) > 2L) {
      for (i in 3L:length(object.list)) {
        y <- c(y, object.list[[i]])
      }
    }
    return(y)
  }
  print("start merge multi-sample...")
  #sceObjMerge <- merge(x = singleList[[1]], y = merge.y(singleList), add.cell.ids = mySamples)
  if(length(singleList) == 1){
    sceObjMerge <- singleList[[1]]
  }else{
    sceObjMerge <- merge(x = singleList[[1]], y = merge.y(singleList))
  }
  print("start normalized...")
  sceObjMerge <- DealWithNormalized(object = sceObjMerge, method = normalize_method)
  sceObjMerge$batch <- NA
  for(i in 1:length(unique(sceObjMerge@meta.data$sample))){
    s <- unique(sceObjMerge@meta.data$sample)[i]
    sceObjMerge@meta.data$batch[which(sceObjMerge@meta.data$sample == s)] <- my.batchs[s]
  }
  print("not removebatch, just merge sample is done...")
  return(sceObjMerge)
}

#removeBatch == "yes"
RemoveBatchAndIntegrate <- function(singleList){
  print("start integrated multi-sample...")
  print("merge the same batch of data")
  singleList2 <- MergeByBatch(singleList, batchs = my.batchs)
  print(singleList2)
  print("merge the same batch is done...")
  myK <- .Kparameter(object = singleList)
  cat("k parameter is: ", myK, "\n")
  kfilterVal <- myK[1]
  kweightVal <- myK[2]
  kanchorval <- 5
  cat("data normalized, method is: ", normalize_method, "\n")
  cat("kfilterVal is: ", kfilterVal, "\n")
  cat("kanchorval is: ", kanchorval, "\n")
  cat("kweightVal is: ", kweightVal, "\n")
  singleList2 <- lapply(X = singleList2, FUN = DealWithNormalized, method = normalize_method)
  
  print("integrated data")
  if(normalize_method == "SCT"){
    pancreas.features <- SelectIntegrationFeatures(object.list = singleList2, nfeatures = as.numeric(parameter_json$fea_n))
    save(pancreas.features, file = file.path(opt$od, "pancreas_features.Rdata"))
    pancreas.list <- PrepSCTIntegration(object.list = singleList2, anchor.features = pancreas.features)
    if(rpca == "yes"){
      pancreas.list <- lapply(X = pancreas.list, FUN = RunPCA, features = pancreas.features)
      save(pancreas.list, file = file.path(opt$od, "pancreas_list.Rdata"))
      single.anchors <- FindIntegrationAnchors(object.list = pancreas.list, anchor.features = pancreas.features, dims = 1:30, normalization.method = "SCT", reduction = "rpca", k.filter = kfilterVal, k.anchor = kanchorval)
    }else{
      single.anchors <- FindIntegrationAnchors(object.list = pancreas.list, anchor.features = pancreas.features, dims = 1:30, normalization.method = "SCT", k.filter = kfilterVal, k.anchor = kanchorval)
    }
    sceObjInt <- IntegrateData(anchorset = single.anchors, normalization.method = "SCT", dims = 1:30, k.weight = kweightVal)
  }
  if(normalize_method == "Normalize"){
    pancreas.features <- SelectIntegrationFeatures(object.list = singleList2, nfeatures = as.numeric(parameter_json$fea_n))
    if(rpca == "yes"){
      singleList2 <- lapply(X = singleList2, FUN = RunPCA, features = pancreas.features)
      single.anchors <- FindIntegrationAnchors(object.list = singleList2, anchor.features = pancreas.features, dims = 1:30, normalization.method = "LogNormalize", reduction = "rpca", k.filter = kfilterVal, k.anchor = kanchorval)
    }else{
      single.anchors <- FindIntegrationAnchors(object.list = singleList2, anchor.features = pancreas.features, dims = 1:30, normalization.method = "LogNormalize", k.filter = kfilterVal, k.anchor = kanchorval)
    }
    sceObjInt <- IntegrateData(anchorset = single.anchors, dims = 1:30, normalization.method = "LogNormalize", k.weight = kweightVal)
  }
  print("data integrated done...")
  #save(sceObjInt, file = file.path(opt$od, "sceObjInt.Rdata"))
  return(sceObjInt)
}

CalculateCellCycle <- function(object, s.genes, g2m.genes){
  # cell cycle
  object <- CellCycleScoring(object, s.features = s.genes, g2m.features = g2m.genes, assay = "RNA")
  object$CellCycle.diff <- object$S.Score - object$G2M.Score
  
  #object2 <- object %>% ScaleData() %>% RunPCA()
  #p.cycle.pca <- DimPlot(object2, reduction = "pca", group.by = "Phase")
  #SavePlot(od = opt$od, filename = "cell_cycle_noremove_pca", data = p.cycle.pca)
  
  print("Start remove cell cycle...")
  object <- DealWithNormalized(object = object,  method = normalize_method, removeCycle = parameter_json$removeCellcycle)
  # object <- ScaleData(object, vars.to.regress = c("CellCycle.diff"))
  print("End...")
  #object2 <- object %>% RunPCA()
  #p.cycle.pca <- DimPlot(object2, reduction = "pca", group.by = "Phase")
  #SavePlot(od = opt$od, filename = "cell_cycle_remove_pca", data = p.cycle.pca)
  
  return(object)
}

## 增加去批次的函数
#removeBatch == "yes"
RemoveBatchAndIntegrate_cca <- function(singleList){
  print("start integrated multi-sample...")
  print("merge the same batch of data")
  singleList2 <- MergeByBatch(singleList, batchs = my.batchs)
  print(singleList2)
  print("merge the same batch is done...")
  myK <- .Kparameter(object = singleList)
  cat("k parameter is: ", myK, "\n")
  kfilterVal <- myK[1]
  kweightVal <- myK[2]
  kanchorval <- 5
  cat("data normalized, method is: ", normalize_method, "\n")
  cat("kfilterVal is: ", kfilterVal, "\n")
  cat("kanchorval is: ", kanchorval, "\n")
  cat("kweightVal is: ", kweightVal, "\n")
  singleList2 <- lapply(X = singleList2, FUN = DealWithNormalized, method = normalize_method)
  
  print("integrated data")
  if(normalize_method == "SCT"){
    pancreas.features <- SelectIntegrationFeatures(object.list = singleList2, nfeatures = as.numeric(parameter_json$fea_n))
    save(pancreas.features, file = file.path(opt$od, "pancreas_features.Rdata"))
    pancreas.list <- PrepSCTIntegration(object.list = singleList2, anchor.features = pancreas.features)
    if(cca == "yes"){
      pancreas.list <- lapply(X = pancreas.list, FUN = RunPCA, features = pancreas.features)
      save(pancreas.list, file = file.path(opt$od, "pancreas_list.Rdata"))
      single.anchors <- FindIntegrationAnchors(object.list = pancreas.list, anchor.features = pancreas.features, dims = 1:30, normalization.method = "SCT", reduction = "cca", k.filter = kfilterVal, k.anchor = kanchorval)
    }else{
      single.anchors <- FindIntegrationAnchors(object.list = pancreas.list, anchor.features = pancreas.features, dims = 1:30, normalization.method = "SCT", k.filter = kfilterVal, k.anchor = kanchorval)
    }
    sceObjInt <- IntegrateData(anchorset = single.anchors, normalization.method = "SCT", dims = 1:30, k.weight = kweightVal)
  }
  if(normalize_method == "Normalize"){
    pancreas.features <- SelectIntegrationFeatures(object.list = singleList2, nfeatures = as.numeric(parameter_json$fea_n))
    if(cca == "yes"){
      singleList2 <- lapply(X = singleList2, FUN = RunPCA, features = pancreas.features)
      single.anchors <- FindIntegrationAnchors(object.list = singleList2, anchor.features = pancreas.features, dims = 1:30, normalization.method = "LogNormalize", reduction = "cca", k.filter = kfilterVal, k.anchor = kanchorval)
    }else{
      single.anchors <- FindIntegrationAnchors(object.list = singleList2, anchor.features = pancreas.features, dims = 1:30, normalization.method = "LogNormalize", k.filter = kfilterVal, k.anchor = kanchorval)
    }
    sceObjInt <- IntegrateData(anchorset = single.anchors, dims = 1:30, normalization.method = "LogNormalize", k.weight = kweightVal)
  }
  print("data integrated done...")
  #save(sceObjInt, file = file.path(opt$od, "sceObjInt.Rdata"))
  return(sceObjInt)
}
RemoveBatchAndIntegrate_harmony <- function(singleList){
  print("start integrated multi-sample...")
  print("merge the same batch of data")
  singleList2 <- MergeByBatch(singleList, batchs = my.batchs)
  print(singleList2)
  cat("data normalized, method is: ", normalize_method, "\n")
  cat("kfilterVal is: ", kfilterVal, "\n")
  cat("kanchorval is: ", kanchorval, "\n")
  cat("kweightVal is: ", kweightVal, "\n")
  singleList2 <- lapply(X = singleList2, FUN = DealWithNormalized, method = normalize_method)
  print("integrated data")
  RunHarmony
  singleList3 <- lapply(X = singleList2, FUN = RunHarmony, ,reduction = "pca",group.by.vars = "group",reduction.save = "harmony")
  print("data integrated done...")
  #save(sceObjInt, file = file.path(opt$od, "sceObjInt.Rdata"))
  return(singleList3)
}


set.seed(1)
parameter_json <- RJSONIO::fromJSON(opt$json)
if(class(parameter_json) != "list") { parameter_json <- as.list(parameter_json) }
Species <- parameter_json$Species
normalize_method <- parameter_json$normalize_method
removeBatch <- tolower(parameter_json$removeBatch)
rpca <- tolower(parameter_json$run_rpca)
my.batchs <- strsplit(parameter_json$Batch, split = ",") %>% unlist 
my.fsample <- strsplit(parameter_json$SampleName, split = ",") %>% unlist
if(length(my.batchs) != length(my.fsample)){
  cat("Error:: There is no one-to-one correspondence between batch and sample,please check json!!!")
  q()
}
names(my.batchs) <- my.fsample
print("initial batch")
print(my.batchs)
my.batchs <- my.batchs[str_order(my.batchs, numeric = T)]
print("sorted batch")
print(my.batchs)
if(length(unique(my.batchs)) == 1){
  removeBatch <- "no"
  parameter_json$removeBatch <- "no"
  parameter_json$removeCellcycle <- "no"
}
removeCellcycle <- parameter_json$removeCellcycle
cat("remove batch:", removeBatch, "\n")
cat("remove cell cycle:", removeCellcycle, "\n")
outjson <- toJSON(parameter_json, pretty = TRUE)
cat(outjson, file = opt$json)

if(!is.null(opt$data)){
  #load(opt$data)
  singleList <- readRDS(opt$data)
  
  if(length(singleList) == 1 && !is.list(singleList)){
    mySamples <- unique(singleList$sample)
    singleList <- list(singleList)
    names(singleList) <- mySamples
  }else{
    mySamples <- names(singleList)
    if(length(my.fsample) <= length(mySamples)){
      singleList <- singleList[my.fsample]
    }
  }
  
  cat("start deal with data...")
  if(removeBatch == "no"){
    sceObj <- NotRemoveBatch(singleList)
  }
  ##修改
  
  if(is.null(opt$method)){removeBatch == "no"}

  if(removeBatch == "yes"){
    sceObj <- RemoveBatchAndIntegrate(singleList)
  }
  if(othermethod == "cca"){
    cca = "yes"
    sceObj <- RemoveBatchAndIntegrate_cca(singleList)
  }
  if(othermethod == "harmony"){
    sceObj <- RemoveBatchAndIntegrate_harmony(singleList)
  }
  


  ##
  print(sceObj) 
  saveRDS(sceObj, file = file.path(opt$od, "seurat_object.Rds"))
  
  if(length(intersect(c("nGene", "nUMI", "percent.mt", "percent.ribosome"), colnames(sceObj@meta.data))) == 4){
  # plot
    P.vlnplot.pre <- VlnPlot(sceObj, group.by = "sample", features = c("nGene", "nUMI", "percent.mt", "percent.ribosome"), ncol = 1, pt.size = 0)  * theme(axis.title.x = element_blank(), axis.text = element_text(size = 10))
    SavePlot(od = opt$od, filename = "All_sample_filter.vln", data = P.vlnplot.pre, width = 6, height = 8)
    P1.scatter.pre <- FeatureScatter(sceObj, group.by = "sample", feature1 = "nUMI", feature2 = "percent.mt", pt.size = 0.5) + theme(legend.position = "none", axis.text = element_text(size = 10))
    P2.scatter.pre <- FeatureScatter(sceObj, group.by = "sample", feature1 = "nUMI", feature2 = "nGene", pt.size = 0.5) + theme(axis.text = element_text(size = 10))
    P.scatter.pre <- P1.scatter.pre + P2.scatter.pre
    SavePlot(od = opt$od, filename = "All_sample_filter.scatter", data = P.scatter.pre, width = 4, height = 2, scale = 1.3)  
  }
  # remove cell cycle effect
  if(removeCellcycle == "yes"){
    print("Calculate the cell cycle score...")
    if(length(grep("^human|^Homo_sapiens|^mouse|^Mus_musculus", Species, ignore.case = TRUE, perl = TRUE)) != 0){
      if(length(grep("^mouse|^Mus_musculus", Species, ignore.case = TRUE, perl = TRUE)) != 0) { mykey <- "Mouse" }
      if(length(grep("^human|^Homo_sapiens", Species, ignore.case = TRUE, perl = TRUE)) != 0) { mykey <- "Human" }
      cellcycle_file <- paste(mykey, "cell_cycle_gene.Rda", sep = "_")
      load(file.path(data.path, cellcycle_file))
      sceObj <- CalculateCellCycle(object = sceObj, s.genes = s_genes, g2m.genes = g2m_genes)
    }else{
      cat("there is no cell cycle genes in species:", Species, "\n")
    }
    file1 <- file.path(opt$od, "seurat_object.Rds")
    file2 <- file.path(opt$od, "seurat_object_noremoveCellcycle.Rds")
    cmd <- paste("mv", file1, file2, sep = " ")
    print(cmd)
    system(cmd)
    saveRDS(sceObj, file = file.path(opt$od, "seurat_object.Rds"))
  }
}

time2 <- proc.time()
run.time <- time2 - time1
print(paste0('执行时间：',run.time[3][[1]],'秒'))

