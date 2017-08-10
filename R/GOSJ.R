# #These functions were (loosely) based on similar functions created for the DEXSeq package.
# #
# # Note that DEXSeq is licensed under the GPL v3. Therefore this
# #   code packaged together is licensed under the GPL v3, as noted in the LICENSE file.
# # Some code snippets are "united states government work" and thus cannot be
# #   copyrighted. See the LICENSE file for more information.
# #
# # The current versions of the original functions upon which these were based can be found
# #    here: http://github.com/Bioconductor-mirror/DEXSeq
# #
# # Updated Authorship and license information can be found here:
# #   here: http://github.com/Bioconductor-mirror/DEXSeq/blob/master/DESCRIPTION
#
#
#
# setClass( "JunctionSeqCountSet",
#           contains = "eSet",
#           representation = representation(
#             designColumns = "character",
#             dispFitCoefs = "numeric",
#             fittedMu = "matrix",
#             dispFunctionType = "list",
#             dispFunction = "function",
#             dispFunctionJct  = "function",
#             dispFunctionExon = "function",
#             formulas = "list",
#             annotationFile = "character",
#             geneCountData = "matrix",
#             countVectors = "matrix",
#             altSizeFactors = "data.frame",
#             plottingEstimates = "list",
#             plottingEstimatesVST = "list", #DEPRECATED! VST-xform is fast enough that it's better to calculate them as needed.
#             geneLevelPlottingEstimates = "list",
#             modelFitForHypothesisTest = "list", #USUALLY unused.
#             modelFitForEffectSize = "list", #USUALLY unused.
#             flatGffData = "data.frame",
#             flatGffGeneData = "list",
#             analysisType = "character",
#             DESeqDataSet = "DESeqDataSet",
#             modelCoefficientsSample    = "list", #USUALLY unused.
#             modelCoefficientsGene      = "list"  #USUALLY unused.
#           ),
#           prototype = prototype( new( "VersionedBiobase",
#                                       versions = c( classVersion("eSet"), JunctionSeqCountSet = "0.0.5" ) ) )
# )
#
#
# makeDESeqDataSetFromJSCS <- function(jscs, test.formula1){
#   countData <- jscs@countVectors
#   colData <- rbind.data.frame(
#     cbind.data.frame(data.frame(sample = rownames(pData(jscs))) , pData(jscs) ),
#     cbind.data.frame(data.frame(sample = rownames(pData(jscs))) , pData(jscs) )
#   )
#   colData$countbin <- factor(c(  rep("this",ncol(countData)/2)  ,   rep("others",ncol(countData)/2)   ), levels = c("this","others"))
#   for(i in 1:ncol(colData)){
#     colData[[i]] <- factor(colData[[i]])
#   }
#   colData <- DataFrame(colData)
#   rownames(colData) <- colnames(countData)
#
#   jscs@DESeqDataSet <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = test.formula1, ignoreRank = TRUE)
#   return(jscs)
# }
#
# getModelFit <- function(jscs, featureID, geneID, countbinID, modelFitType = c("fitH0","fitH1","fitEffect")){
#   if(is.null(featureID)){
#     if(is.null(geneID) | is.null(countbinID)){
#       stop("ERROR: getModelFit(): either featureID or geneID AND countbinID must be set!")
#     } else {
#       featureID <- paste0(geneID, ":", countbinID)
#     }
#   }
#
#   out <- list()
#   if(any(modelFitType == "fitH0")){
#     if(is.null(jscs@modelFitForHypothesisTest)){
#       stop("ERROR: getModelFit(): no modelFitForHypothesisTest found. Run testJunctionsForDJU().")
#     }
#     out[["fitH0"]] <- jscs@modelFitForHypothesisTest[[featureID]][["fitH0"]]
#   }
#   if(any(modelFitType == "fitH1")){
#     if(is.null(jscs@modelFitForHypothesisTest)){
#       stop("ERROR: getModelFit(): no modelFitForHypothesisTest found. Run testJunctionsForDJU().")
#     }
#     out[["fitH1"]] <- jscs@modelFitForHypothesisTest[[featureID]][["fitH1"]]
#   }
#   if(any(modelFitType == "fitEffect")){
#     if(is.null(jscs@modelFitForEffectSize)){
#       stop("ERROR: getModelFit(): no modelFitForHypothesisTest found. Run testJunctionsForDJU().")
#     }
#     out[["fitEffect"]] <- jscs@modelFitForEffectSize[[featureID]]
#   }
#   return(out)
# }
#
#
# newJunctionSeqCountSet <- function( countData, geneCountData, design, geneIDs, countbinIDs, featureIntervals=NULL, transcripts=NULL){
#
#   countData <- as.matrix( countData )
#   if( any( round( countData ) != countData ) ){
#     stop( "The countData is not integer." )}
#   mode( countData ) <- "integer"
#
#   geneCountData <- as.matrix( geneCountData )
#   if( any( round( geneCountData ) != geneCountData ) ){
#     stop( "The geneCountData is not integer." )}
#   mode( geneCountData ) <- "integer"
#
#   if( is( design, "matrix" ) ){
#     design <- as.data.frame( design )}
#
#   rownames(countData) <- paste(geneIDs, countbinIDs, sep=":")
#   if( any( duplicated(rownames(countData) ) ) ) {
#     stop("The geneIDs or countbinIDs are not unique")
#   }
#   if(any(duplicated(rownames(geneCountData)))){
#     stop("The geneIDs in the geneCountData are not unique")
#   }
#   if(any(duplicated(colnames(geneCountData)))){
#     stop("The sample ID's in the geneCountData are not unique")
#   }
#
#   phenoData <- annotatedDataFrameFrom( countData, byrow=FALSE )
#   featureData <- annotatedDataFrameFrom( countData, byrow=TRUE )
#
#   phenoData$sizeFactor <- rep( NA_real_, ncol(countData) )
#   varMetadata( phenoData )[ "sizeFactor", "labelDescription" ] <- "size factor (relative estimate of sequencing depth)"
#
#   geneIDs <- as.factor( geneIDs )
#   if( length(geneIDs) != nrow(countData) )
#     stop( "geneIDs must be of the same length as the number of columns in countData")
#
#   featureData$geneID <- geneIDs
#   varMetadata( featureData )[ "geneID", "labelDescription" ] <- "ID of gene to which the exon belongs"
#
#   countbinIDs <- as.character( countbinIDs )
#   if( length(countbinIDs) != nrow(countData) )
#     stop( "countbinIDs must be of the same length as the number of columns in countData")
#
#   featureData$countbinID <- countbinIDs
#   varMetadata( featureData )[ "countbinID", "labelDescription" ] <- "feature ID (unique only within a gene)"
#
#   if( is.null(featureIntervals) ){
#     featureIntervals <- data.frame(
#       chr    = rep( NA_character_, nrow( countData ) ),
#       start  = rep( NA_integer_,   nrow( countData ) ),
#       end    = rep( NA_integer_,   nrow( countData ) ),
#       strand = rep( NA_character_, nrow( countData ) ) ) }
#
#   featureData$testable <- rep( NA_real_, nrow( countData ) )
#   varMetadata( featureData )[ "testable", "labelDescription" ] <- "slot indicating if an feature should be considered in the test."
#   featureData$status <- rep( "TBD", nrow( countData ) )
#
#   featureData$allZero <- rep( NA_real_, nrow( countData ) )
#   varMetadata( featureData )[ "allZero", "labelDescription" ] <- "slot indicating if the feature count is zero across all samples."
#
#
#   featureData$status <- rep( "TBD", nrow( countData ) )
#   varMetadata( featureData )[ "status", "labelDescription" ] <- "Feature status (either 'OK' or the reason that the feature is untestable)."
#
#   featureData$baseMean <- rep( NA_real_, nrow( countData ) )
#   varMetadata( featureData )[ "baseMean", "labelDescription" ] <- "Mean normalized counts across all samples."
#   featureData$baseVar <- rep( NA_real_, nrow( countData ) )
#   varMetadata( featureData )[ "baseVar", "labelDescription" ] <- "Simple variance of normalized counts across all samples."
#
#   featureData$dispBeforeSharing <- rep( NA_real_, nrow( countData ) )
#   varMetadata( featureData )[ "dispBeforeSharing", "labelDescription" ] <- "feature dispersion (Cox-Reid estimate)"
#
#   featureData$dispFitted <- rep( NA_real_, nrow( countData ) )
#   varMetadata( featureData )[ "dispFitted", "labelDescription" ] <- "Fitted mean-variance estimate."
#
#   featureData$dispersion <- rep( NA_real_, nrow( countData ) )
#   varMetadata( featureData )[ "dispersion", "labelDescription" ] <- "Final dispersion estimate."
#
#   featureData$pvalue <- rep( NA_real_, nrow( countData ) )
#   varMetadata( featureData )[ "pvalue", "labelDescription" ] <- "p-value from testForDEU"
#
#   featureData$padjust <- rep( NA_real_, nrow( countData ) )
#   varMetadata( featureData )[ "padjust", "labelDescription" ] <- "BH adjusted p-value"
#
#
#   featureIntervals <- as.data.frame( featureIntervals )
#
#   # in case it was a GRanges object before, change the colname:
#   if( "seqnames" %in% colnames(featureIntervals) ){
#     colnames(featureIntervals)[ colnames(featureIntervals) == "seqnames" ] <- "chr"   }
#
#   if( !all( c( "chr", "start", "end", "strand" ) %in% colnames(featureIntervals) ) ){
#     stop( "featureIntervals must be a data frame with columns 'chr', 'start', 'end', and 'strand'." )}
#
#   if(is.null(transcripts)){
#     transcripts <- rep(NA_character_, nrow( countData ) )}
#
#   if(!is(transcripts, "character")){
#     stop("transcript information must be a character vector")}
#
#   featureData$chr    <- factor( featureIntervals$chr )
#   featureData$start  <- featureIntervals$start
#   featureData$end    <- featureIntervals$end
#   featureData$strand <- factor( featureIntervals$strand )
#   featureData$transcripts <- transcripts
#   varMetadata( featureData )[ "chr",    "labelDescription" ] <- "chromosome of feature"
#   varMetadata( featureData )[ "start",  "labelDescription" ] <- "start of feature"
#   varMetadata( featureData )[ "end",    "labelDescription" ] <- "end of feature"
#   varMetadata( featureData )[ "strand", "labelDescription" ] <- "strand of feature"
#   varMetadata( featureData )[ "transcripts", "labelDescription" ] <- "transcripts in which this feature is contained"
#
#   featureData$baseMean <- rep( NA_real_, nrow( countData ) )
#   varMetadata( featureData )[ "baseMean", "labelDescription" ] <- "The mean normalized counts across all samples."
#
#
#   if( is( design, "data.frame" ) || is( design, "AnnotatedDataFrame" ) ) {
#     stopifnot( nrow( design ) == ncol( countData ) )
#     stopifnot( all( unlist( lapply(design, class) ) == "factor" ) )
#     design <- as( design, "AnnotatedDataFrame" )
#     dimLabels(design) <- dimLabels(phenoData)
#     rownames( pData(design) ) <- rownames( pData(phenoData) )
#     phenoData <- combine( phenoData, design )
#     rvft <- c( `_all` = NA_character_ )
#     designColumns <- varLabels(design)
#   } else {
#     design <- factor( design, levels=unique(design))
#     stopifnot( length( design ) == ncol( countData ) )
#     phenoData$`condition` <- factor( design )
#     varMetadata( phenoData )[ "condition", "labelDescription" ] <- "experimental condition, treatment or phenotype"
#     designColumns <- "condition"
#   }
#   jscs <- new( "JunctionSeqCountSet",
#                assayData = assayDataNew( "environment", counts=countData ),
#                phenoData = phenoData,
#                featureData = featureData,
#                designColumns = designColumns,
#                dispFitCoefs = c( NA_real_, NA_real_ ),
#                geneCountData = geneCountData
#   )
#   jscs
# }
#
# setValidity( "JunctionSeqCountSet", function( object ) {
#
#   if( !all( object@designColumns %in% names(pData(object)) ) )
#     return( "Not all designColumns appear in phenoData." )
#
#   if( ! "sizeFactor" %in% names(pData(object)) )
#     return( "phenoData does not contain a 'sizeFactor' column.")
#   if( ! is( pData(object)$`sizeFactor`, "numeric" ) )
#     return( "The 'sizeFactor' column in phenoData is not numeric." )
#
#   if( ! "geneID" %in% names(fData(object)) )
#     return( "featureData does not contain a 'geneID' column.")
#   if( ! is( fData(object)$geneID, "factor" ) )
#     return( "The 'geneID' column in fData is not a factor." )
#
#   if( ! "countbinID" %in% names(fData(object)) )
#     return( "featureData does not contain an 'countbinID' column.")
#   if( ! is( fData(object)$countbinID, "character" ) )
#     return( "The 'countbinID' column in fData is not a character vector." )
#
#   if( ! "chr"  %in% names(fData(object)) )
#     return( "featureData does not contain a 'chr' column.")
#   if( ! is( fData(object)$chr, "factor" ) )
#     return( "The 'chr' column in fData is not a factor." )
#
#   if( ! "start"  %in% names(fData(object)) )
#     return( "featureData does not contain a 'start' column.")
#   if( ! is( fData(object)$start, "integer" ) )
#     return( "The 'start' column in fData is not integer." )
#
#   if( ! all(featureNames(object) == paste(geneIDs(object), countbinIDs(object), sep=":") ) )
#     return( "The featureNames do not match with the geneIDs:countbinIDs" )
#
#   if( ! all(rownames( counts(object) ) == featureNames(object) ) )
#     return( "The rownames of the count matrix do not coincide with the featureNames" )
#
#   if( ! all(rownames( fData( object ) ) == featureNames( object ) ) )
#     return( "The rownames of the featureData do not coincide with the featureNames" )
#
#
#
#   if( ! "end"  %in% names(fData(object)) )
#     return( "featureData does not contain a 'end' column.")
#   if( ! is( fData(object)$end, "integer" ) )
#     return( "The 'end' column in fData is not integer." )
#
#   if( ! "strand"  %in% names(fData(object)) )
#     return( "featureData does not contain a 'strand' column.")
#   if( ! is( fData(object)$strand, "factor" ) )
#     return( "The 'strand' column in fData is not a factor." )
#   if( !is(fData(object)$dispersion, "numeric")){
#     return( "The 'dispersion' is not numeric")}
#   if( !is(fData(object)$dispFitted, "numeric")){
#     return( "The 'dispFitted' is not numeric")}
#   if( !is(fData(object)$dispBeforeSharing, "numeric")){
#     return( "The 'dispBeforeSharing' column is not numeric")}
#   if( !is(fData(object)$pvalue, "numeric")){
#     return( "The 'pvalue' values are not numeric")}
#   if( !is(fData(object)$padjust, "numeric")){
#     return( "The 'padjust' values are not numeric")}
#   if( !is.integer( assayData(object)[["counts"]] ) )
#     return( "The count data is not in integer mode." )
#
#   if( any( assayData(object)[["counts"]] < 0 ) )
#     return( "The count data contains negative values." )
#
#   if( length( object@dispFitCoefs ) != 2 )
#     return( "dispFitCoefs is not a vector of length 2." )
#
#   TRUE
# } )
#
#
# setMethod("counts", signature(object="JunctionSeqCountSet"),
#           function( object, normalized=FALSE) {
#             cds <- object
#             if(!normalized){
#               assayData(cds)[["counts"]]
#             } else {
#               if(any(is.na( sizeFactors(cds)))) {
#                 stop( "Please first calculate size factors or set normalized=FALSE")
#               } else {
#                 t(t( assayData( cds )[["counts"]] ) / sizeFactors(cds) )
#               }
#             }
#           })
# setReplaceMethod("counts", signature(object="JunctionSeqCountSet", value="matrix"),
#                  function( object, value ) {
#                    cds <- object
#                    assayData(cds)[[ "counts" ]] <- value
#                    validObject(cds)
#                    cds
#                  })
#
# setMethod("sizeFactors",  signature(object="JunctionSeqCountSet"),
#           function(object) {
#             cds <- object
#             sf <- pData(cds)$sizeFactor
#             names( sf ) <- colnames( counts(cds) )
#             sf
#           })
#
# setReplaceMethod("sizeFactors",  signature(object="JunctionSeqCountSet", value="numeric"),
#                  function(object, value ) {
#                    cds <- object
#                    pData(cds)$sizeFactor <- value
#                    validObject( cds )
#                    cds
#                  })
#
# setMethod("design", signature(object="JunctionSeqCountSet"),
#           function( object, drop=TRUE, asAnnotatedDataFrame=FALSE ) {
#             cds <- object
#             if( asAnnotatedDataFrame )
#               return( phenoData(cds)[, cds@designColumns ] )
#             ans <- pData(cds)[, cds@designColumns, drop=FALSE ]
#             if( ncol(ans) == 1 && drop ) {
#               ans <- ans[,1]
#               names(ans) <- colnames( counts(cds) ) }
#             else
#               rownames( ans ) <- colnames( counts(cds) )
#             ans
#           })
# setReplaceMethod("design", signature(object="JunctionSeqCountSet"),
#                  function( object, value ) {
#                    cds <- object
#                    ## Is it multivariate or just a vector?
#                    if( ncol(cbind(value)) > 1 )
#                      value <- as( value, "AnnotatedDataFrame" )
#                    else {
#                      value <- new( "AnnotatedDataFrame",
#                                    data = data.frame( condition = value ) )
#                      varMetadata( value )[ "condition", "labelDescription" ] <-
#                        "experimental condition, treatment or phenotype" }
#
#                    rownames( pData(value) ) <- rownames( pData(cds) )
#                    dimLabels( value ) <- dimLabels( phenoData(cds) )
#                    phenoData(cds) <- combine(
#                      phenoData(cds)[ , !( colnames(pData(cds)) %in% cds@designColumns ), drop=FALSE ],
#                      value )
#                    cds@designColumns <- colnames( pData(value) )
#                    validObject(cds)
#                    cds
#                  })
#
# geneIDs <- function( ecs ) {
#   stopifnot( is( ecs, "JunctionSeqCountSet" ) )
#   g <- fData(ecs)$geneID
#   names(g) <- rownames( counts(ecs) )
#   g
# }
#
# `geneIDs<-` <- function( ecs, value ) {
#   stopifnot( is( ecs, "JunctionSeqCountSet" ) )
#   fData(ecs)$geneID <- value
#   validObject(ecs)
#   ecs
# }
#
# countbinIDs <- function( ecs ) {
#   stopifnot( is( ecs, "JunctionSeqCountSet" ) )
#   g <- fData(ecs)$countbinID
#   names(g) <- rownames( counts(ecs) )
#   g
# }
#
# `countbinIDs<-` <- function( ecs, value ) {
#   stopifnot( is( ecs, "JunctionSeqCountSet" ) )
#   fData(ecs)$countbinID <- value
#   validObject(ecs)
#   ecs
# }
#
#
# subsetByGenes <- function( ecs, genes ) {
#   stopifnot( is( ecs, "JunctionSeqCountSet" ) )
#   stopifnot( all( genes %in% levels(geneIDs(ecs)) ) )
#   ecs2 <- ecs[ as.character(geneIDs(ecs)) %in% genes, ]
#   ecs2
# }
#
# geneCountTable <- function( ecs ) {
#   stopifnot( is( ecs, "JunctionSeqCountSet" ) )
#   do.call( rbind,
#            tapply( seq_len(nrow(ecs)), geneIDs(ecs), function(rows)
#              colSums( counts(ecs)[rows,,drop=FALSE] ) ) )
# }
#
# DEUresultTable <- function(ecs)
# {
#   result <- data.frame(
#     geneID=geneIDs(ecs),
#     countbinID=countbinIDs(ecs),
#     dispersion=featureData(ecs)$dispersion,
#     pvalue=fData(ecs)$pvalue,
#     padjust=fData(ecs)$padjust,
#     baseMean=rowMeans(counts(ecs, normalized=TRUE)))
#
#   extracol <- regexpr("log2fold", colnames(fData(ecs)))==1
#   if(any(extracol)){
#     w <- which(extracol)
#     result <- data.frame(result, fData(ecs)[,w])
#     colnames(result)[7:(6+length(w))] <- colnames(fData(ecs))[w]
#   }
#   result
# }
#
#' codes for simulation and plot PWF
#' simulation code from xiaobei
#' PWF plot code from goseq
#' Get data set for simulation studies
#'
#' @param counts
#' @param drop.extreme.dispersion
#' @param drop.low.lambda
#'
#' @return
#' @export
#'
#' @examples
getDataset <- function(counts, drop.extreme.dispersion = 0.1, drop.low.lambda = TRUE) {
  ## this function generates NB parameters from real dataset ##
  ## it is low-level function of NBsim ##
  d <- DGEList(counts)
  d <- calcNormFactors(d)
  cp <- round(cpm(d,normalized.lib.sizes=TRUE),1)
  if(drop.low.lambda) d <- d[rowSums(cp>1) >= 2, ]
  d$AveLogCPM <- log2(rowMeans(cpm(d, prior.count = 1e-5)))
  d <- estimateGLMCommonDisp(d)
  d <- estimateGLMTrendedDisp(d)
  d <- estimateGLMTagwiseDisp(d)
  dispersion <- d$tagwise.dispersion
  AveLogCPM <- d$AveLogCPM
  if(is.numeric(drop.extreme.dispersion))
  {
    bad <- quantile(dispersion, 1-drop.extreme.dispersion, names = FALSE)
    ids <- dispersion <= bad
    AveLogCPM <- AveLogCPM[ids]
    dispersion <- dispersion[ids]
  }
  dataset.AveLogCPM <- AveLogCPM
  dataset.dispersion <- dispersion
  dataset.lib.size <- d$samples$lib.size
  dataset.nTags <- nrow(d)
  list(dataset.AveLogCPM = dataset.AveLogCPM, dataset.dispersion = dataset.dispersion, dataset.lib.size = dataset.lib.size, dataset.nTags = dataset.nTags)
}
#' Combine3Re
#'
#' This function combines feature based differential usage results, the derived gene based results from feature based results and
#' the resutls from rMAT to identify the common results between 3 methods
#'
#' @param Re.PJ: Output from JunctionSeq
#' @param re.PJ.gene.based: Gene based results from JunctionSeq
#' @param re.rMAT: Results from rMAT
#' @param outputfile_DGE_2_FC_P_geneWise
#' @param gene.model
#' @param Output_GO_BP_feature.gene.rMAT.2.FC
#' @param Re.PJ.selected.feature.FC.p
#' @param Output_GO_BP_gene
#' @param Output_GO_BP_feature
#' @param Output_GO_BP_feature_gene
#' @param Output_GO_BP_feature_gene_redefined
#'
#' @return
#' @export
#'
#' @examples
#' Re.combine.3.methods<-Combine3Re(Re.PJ, re.PJ.gene.based, re.rMAT,gene.model,output.dir.name.new)
#'
#' Re.combine.3.methods.rMAT.SE<-Combine3Re(Re.PJ, re.PJ.gene.based, re.rMAT,gene.model,output.dir.name.rMAT_SE_based)
#'
#'
Combine3Re <- function(Re.PJ, re.PJ.gene.based, re.rMAT,gene.model, Output_files_dir) {

  Re.PJ.selected.feature.2.FC.p<-Select_DE_gene_basd_on_Feature(Re.PJ,re.PJ.gene.based,re.rMAT,"SE",1,0.05,gene.model,Output_files_dir)

  DE_type<-c("Feature","GeneWise","rMAT","FeatureGeneWise","FeaturerMAT","GeneWiserMAT","FeatureGeneWiseRMAT")

  Re.combine.1<-OutputGOBasedDEfromFeatures3(Re.PJ.selected.feature.2.FC.p[[1]],"Feature",gene.model,Output_files_dir)
  Re.combine.2<-OutputGOBasedDEfromFeatures3(Re.PJ.selected.feature.2.FC.p[[1]],"GeneWise",gene.model,Output_files_dir)
  Re.combine.3<-OutputGOBasedDEfromFeatures3(Re.PJ.selected.feature.2.FC.p[[1]],"rMAT",gene.model,Output_files_dir)
  Re.combine.4<-OutputGOBasedDEfromFeatures3(Re.PJ.selected.feature.2.FC.p[[1]],"FeatureGeneWise",gene.model,Output_files_dir)
  Re.combine.5<-OutputGOBasedDEfromFeatures3(Re.PJ.selected.feature.2.FC.p[[1]],"FeaturerMAT",gene.model,Output_files_dir)
  Re.combine.6<-OutputGOBasedDEfromFeatures3(Re.PJ.selected.feature.2.FC.p[[1]],"GeneWiserMAT",gene.model,Output_files_dir)
  Re.combine.7<-OutputGOBasedDEfromFeatures3(Re.PJ.selected.feature.2.FC.p[[1]],"FeatureGeneWiseRMAT",gene.model,Output_files_dir)

  #Re.combine<-list(rMAT=Re.combine.3)

   Re.combine<-list(Feature=Re.combine.1,
                    GeneWise=Re.combine.2,
                    rMAT=Re.combine.3,
                    FeatureGeneWise=Re.combine.4,
                    FeaturerMAT=Re.combine.5,
                    GeneWiserMAT=Re.combine.6,
                    FeatureGeneWiseRMAT=Re.combine.7)

  return(Re.combine)

}
#' CombineCpTftGoGeneID
#'
#' @return
#' @export
#'
#' @examples
#'
#' gene.id.all<-CombineCpTftGoGeneID("/media/H_driver/Annotation/MsigDB/c2.cp.Mouse.v5.1.symbols.gmt",
#' "/media/H_driver/Annotation/MsigDB/c3.tft.Mouse.v5.1.symbols.gmt","/media/H_driver/Annotation/GO_gene_mouse/gene_GO.mgi",
#' "/media/H_driver/Annotation/mm10/genes_table_02052016.csv")
#'
#'
CombineCpTftGoGeneID <- function(cp_input,tft_input,Go_input,gene_anno){

  gene.2.cat.cp.mouse<-Gmt2GeneCat(cp_input,gene_anno)

  gene.2.cat.tft.mouse<-Gmt2GeneCat(tft_input,gene_anno)

  gene.id.from.go<-ExtractGeneIDFromGo(Go_input,gene_anno)

  gene.id.all<-unique(c(names(gene.2.cat.cp.mouse),names(gene.2.cat.tft.mouse),gene.id.from.go))

  return(gene.id.all)

}

#' CompareDEfromGeneWithDEfromFeatures
#'
#' This function compares gene-based DE with feature-based DE
#'
#'
#' @param Re.PJ.gene: gene-based DE results
#' @param re.PJ.gene.based: feature-based DE results
#' @param Output_venn_file: file name for venn
#'
#' @return
#' @export
#'
#' @examples
#'
#' Gene.based.DE.feature.based.DE<-CompareDEfromGeneWithDEfromFeatures(Re.PJ.gene,re.PJ.gene.based,"/media/H_driver/PJ/gene_feature_DE_2.tiff")
#'
CompareDEfromGeneWithDEfromFeatures<-function(Re.PJ.gene,re.PJ.gene.based,Output_venn_file){

re<-merge(Re.PJ.gene,pData(re.PJ.gene.based),by=0)

no.re.testable.index<-which(as.character(re$mostSigID)=="character(0)")

re2<-re[-no.re.testable.index,]

gene.based.DE<-re2[which(re2$padj<0.05),]$geneID

feature.based.DE<-re2[which(re2$geneWisePadj<0.05),]$geneID


all.gene.index<-rep(0,length(re2$geneID))

names(all.gene.index)<-re2$geneID

all.gene.index.gene.based<-all.gene.index
all.gene.index.gene.based[which(names(all.gene.index.gene.based) %in% gene.based.DE)]=1

exonplussj=re2$numExons+re2$numKnown

pwf.DE_interest.gene.based.using.gene.length=nullp(all.gene.index.gene.based,"mm10","ensGene",plot.fit = FALSE)

pwf.DE_interest.gene.based.using.exonplussj=nullp(all.gene.index.gene.based,"mm10","ensGene",bias.data = exonplussj,
                                                  plot.fit = FALSE)


all.gene.index.feature.based<-all.gene.index
all.gene.index.feature.based[which(names(all.gene.index.feature.based) %in% feature.based.DE)]=1

pwf.DE_interest.feature.based.using.gene.length=nullp(all.gene.index.feature.based,"mm10","ensGene",plot.fit = FALSE)

pwf.DE_interest.feature.based.using.exonplussj=nullp(all.gene.index.feature.based,"mm10","ensGene",bias.data = exonplussj,
                                                  plot.fit = FALSE)


Re3<-list(GeneBased=gene.based.DE,FeatureBased=feature.based.DE,GenePlusFeature=re2,pwfGeneGL=pwf.DE_interest.gene.based.using.gene.length,
          pwfGeneFeature=pwf.DE_interest.gene.based.using.exonplussj,pwfFeatureGL=pwf.DE_interest.feature.based.using.gene.length,
          pwfFeatureFeature=pwf.DE_interest.feature.based.using.exonplussj)

venn.plot <- venn.diagram(
  x = Re3[c(1,2)],
  filename = Output_venn_file,
  height = 3000,
  width = 3500,
  resolution = 1000,
  col = "black",
  lty = "dotted",
  lwd = 1,
  fill = c("red","blue"),
  alpha = 0.50,
  label.col = c(rep("white",3)),
  cex = 0.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("red","blue"),
  cat.cex = 0.5,
  cat.pos = 0.5,
  cat.dist = 0.05,
  cat.fontfamily = "serif"
)

return(Re3)

}
#' CompareGOResults
#'
#' @param Gene.based.DE.feature.based.DE
#' @param gene_model
#'
#' @return
#' @export
#'
#' @examples
#'
#'
CompareGOResults <- function(Gene.based.DE.feature.based.DE, gene_model)
  {
  GO.wall.DE_interest.geneGL=goseq2(Gene.based.DE.feature.based.DE$pwfGeneGL,"mm10","ensGene",gene.model=gene_model)

  GO.wall.DE_interes.geneFT=goseq2(Gene.based.DE.feature.based.DE$pwfGeneFeature,"mm10","ensGene",gene.model=gene_model)

  GO.wall.DE_interest.FtFT=goseq2(Gene.based.DE.feature.based.DE$pwfFeatureFeature,"mm10","ensGene",gene.model=gene_model)

  GO.wall.DE_interest.FtFT=goseq2(Gene.based.DE.feature.based.DE$pwfFeatureFeature,"mm10","ensGene",gene.model=gene_model)
}
#!/usr/bin/env Rscript
#Usage: Rscript Run_DESeq_3.R inputfile_count inputfile_sample_information outputfile_prefix

#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")

#.libPaths("~/R_local_libs")

#install.packages("installr",repos="http://cran.rstudio.com/")
#library(installr)
#updateR()

#.libPaths()

#library("DESeq2")
#library("base")
#library("heatmap3")
#library("lattice")
#library("reshape")
#library("ggplot2")
#library("grid")
#library(gplots)
#library(RColorBrewer)
#library(survival)
#library(limma)
#library(edgeR)
#library(multtest)
#library(biomaRt)
#library(goseq)

#args = commandArgs(trailingOnly=T)

DE_based_on_gene_count <- function(input.file,sample.info.file,output.file) {
  # input.file = args[1]
  # sample.info.file=args[2]
  # output.file = args[3]

  cat(input.file,"\n")
  cat(sample.info.file,"\n")

  data.count<-read.table(input.file,sep="\t",header = T,row.names=1)
  data.info<-read.table(sample.info.file,sep="\t",header = F)

  print(head(data.count))
  print(data.info)

  #data.sample<-unique(trimws(as.character(data.info[,2])))
  #data.number.sample<-length(data.sample)

  #print(data.sample)
  #print(data.number.sample)


  #Function.Check.condition<-function(data.sample,data.info){


  Re<-table(data.info[,2:3])
  print(Re)

  Re1<-as.data.frame.matrix(Re)
  print(Re1)

  print(rownames(Re1))
  print(colnames(Re1))

  data.sample<-rownames(Re1)
  data.condition<-colnames(Re1)

  data.number.sample<-length(data.sample)

  sample.index<-as.list(c(seq(1,data.number.sample)))

  Function.Get.DE.FC<-function(i){

    count.4.sample.control<-data.info[which(data.info[,2] %in% data.sample[i]&trimws(data.info[,3]) %in% c(data.condition[1])),1]
    count.4.sample.treatment<-data.info[which(data.info[,2] %in% data.sample[i]&trimws(data.info[,3]) %in% c(data.condition[2])),1]

    index.control=which(as.character(gsub("X","",colnames(data.count))) %in% c(as.character(count.4.sample.control)))
    index.treatment=which(as.character(gsub("X","",colnames(data.count))) %in% c(as.character(count.4.sample.treatment)))

    data.count.reformated<-cbind(data.count[,index.control],data.count[,index.treatment])

    numControl<-length(count.4.sample.control)
    numTreatment<-length(count.4.sample.treatment)

    condition <- factor(c(rep(data.condition[1],numControl),rep(data.condition[2],numTreatment)))

    dds <- DESeqDataSetFromMatrix(data.count.reformated, DataFrame(condition), ~ condition)
    re<-results(DESeq(dds))
    re.FC<-cbind(as.data.frame(re),2^re[,2])
    colnames(re.FC)[7]="FoldChange"
    write.csv(cbind(re.FC,as.data.frame(counts(dds))),file=paste(output.file,"_",data.sample[i],"_DE.csv",sep=""))
    write.csv(as.data.frame(colData(dds)),file=paste(output.file,"_",data.sample[i],"_info.csv",sep=""))

  }

  lapply(sample.index, Function.Get.DE.FC)
}
#' Draw4Cat
#'
#' @param RE.exon.sj.all.gene
#'
#' @return
#' @export
#'
#' @examples
#' Draw4Cat(RE.exon.sj.all.gene)
#'
Draw4Cat <- function(RE.exon.sj.all.gene.enriched) {
  x=RE.exon.sj.all.gene.enriched[1:10,1:5]
  x$category <- factor(x$category,levels = x$category[order(x$over_represented_pvalue_adjusted,decreasing = TRUE)])
  Function<-x$category
  negative_log10p=-log10(x$over_represented_pvalue_adjusted)
  ggplot(x, aes(x=Function, y=negative_log10p,fill=factor(x$category)))+geom_bar(stat="identity")+geom_hline(yintercept = -log10(0.05))+coord_flip()
}
#' Draw4GO
#'
#' @param RE.exon.sj.all.gene
#'
#' @return
#' @export
#'
#' @examples
#' Draw4GO(RE.exon.sj.all.gene)
#'
Draw4GO <- function(RE.exon.sj.all.gene.enriched) {
  x=RE.exon.sj.all.gene.enriched[1:10,1:7]
  x$term <- factor(x$term,levels = x$term[order(x$over_represented_pvalue_adjusted,decreasing = TRUE)])
  Function<-x$term
  negative_log10p=-log10(x$over_represented_pvalue_adjusted)
  ggplot(x, aes(x=Function, y=negative_log10p,fill=factor(x$category)))+geom_bar(stat="identity")+geom_hline(yintercept = -log10(0.05))+coord_flip()
}

#' DrawFeature
#'
#' @param Re.PJ
#' @param out.dir
#' @param gene_id
#'
#' @return
#' @export
#'
#' @examples
#'
#' DrawFeature(Re.PJ,output.dir.name,"ENSMUSG00000026563+ENSMUSG00000040596+ENSMUSG00000089853+ENSMUSG00000103400")
#'
#'
DrawFeature <- function(Re.PJ,out.dir,gene_id) {
buildAllPlots(Re.PJ,outfile.prefix=paste0(out.dir,gene_id,"/"),gene.list=gene_id,use.plotting.device="png",plot.gene.level.expression=TRUE)
}
#' ExtractGeneIDFromGo
#'
#' @return
#' @export
#'
#' @examples
#' gene.id.from.go<-ExtractGeneIDFromGo("/media/H_driver/Annotation/GO_gene_mouse/gene_GO.mgi","/media/H_driver/Annotation/mm10/genes_table_02052016.csv")
#'
#'
ExtractGeneIDFromGo<- function(gene_go_file,gene_anno_file) {

  GO.gene.mouse <- read.table(gene_go_file,header=F)

  colnames(GO.gene.mouse)=c("gene_id","GO")

  gene.ID.conversion<-read.csv(gene_anno_file)

  names.gene.gmt.2<-match(GO.gene.mouse$gene_id,gene.ID.conversion$gene_id)

  gene.ID.conversion.2<-gene.ID.conversion[names.gene.gmt.2,]

  #gene.2.cat.gmt.2<-gene.2.cat.gmt

  ensembl_gene_id.from.go<-unique(unlist(split(gene.ID.conversion.2[,3],";")))

  ensembl_gene_id.from.go.2<-ensembl_gene_id.from.go[-which(is.na(ensembl_gene_id.from.go))]

  ensembl_gene_id.from.go.2

}
#' Title
#'
#' @return
#' @export
#'
#' @examples
FDR_genewise <- function() {
  cat("tab[notZero]\n")
  print(tab[notZero])
  cat(tab[notZero],file="mm.txt")

  mm=read.table("mm.txt")
  mm=as.numeric(mm)

  cat("which(notZero)\n")
  print(which(notZero))
  cat(which(notZero),file="nn.txt")
  cat(sort(which(notZero)),file="nn2.txt")

  nn=read.table("nn.txt")
  nn=as.numeric(nn)

  cat("theta\n")
  print(theta)
  cat(theta,file="theta.txt")

  theta=read.table("theta.txt")
  theta=as.numeric(theta)

  exon.all=read.table("exon_all.txt")
  exon.all=as.numeric(exon.all)


  cat(pGene,file="pGene.txt",sep="\n",append=TRUE)

  stopifnot(all(is.finite(pGene)))

  ## Determine the thetas to be used
  theta = unique(sort(pGene))

  cat(theta,file="theta.txt",sep="\n",append=TRUE)

  pGene=read.table("pGene.txt")
  pGene=as.numeric(pGene[,1])
  theta = unique(sort(pGene))

  bins   = cut(pGene, breaks=c(-Inf, as.vector(theta)), right = TRUE, include.lowest = TRUE)

  head(bins)

  counts = tabulate(bins, nbins = nlevels(bins))

  denom  = cumsum(counts)

  #theta=read.table("theta.txt")
  #exon.all=as.numeric(ex)

  cut(c(2,3,5,1,6), breaks=c(-Inf, as.vector(c(1,3,2))), right = TRUE, include.lowest = TRUE)



  length(which(exon.all==6))

  cat(length(mm),"\t",length(nn),"\t",length(theta),"\n")



  mm2=mm[1:2]
  nn2=nn[1:2]

  mm3=c(2,3,4,5)
  nn3=c(4,5,2,3)
  theta=c(0.1,0.2,0.3,0.4,0.005,0.3)

  # 2*(1-(1-theta)^4)
  #
  # 3*(1-(1-theta)^5)
  #
  # 4*(1-(1-theta)^2)
  #
  # 5*(1-(1-theta)^3)

  numerator.test    = mapply(function(m, n) {

    #cat("m\n")
    #print(m)

    #cat("n\n")
    #print(n)

    #cat("theta","\n")
    #cat(length(theta),"\n")
    #cat(length(unique(theta)),"\n")
    #print(theta)

    re<-m * (1 - (1-theta)^n)

    #cat("re\n")
    #print(re)
    #cat(re,file="re.txt")

    re

  },
  m = mm3,
  n = nn3)
}
#' GSA.read.gmt.2
#' Descrption: In addition to GO terms, the GOSJ package can analyze other types of gene sets defined by users, for example, the gene sets in MSigDB database.The GSA.read.gmt.2 function takes input the gene sets in .gmt format, i.e. the first column is the pathway or gene set name, and the rest columns include gene names of the genes within the pathway. The gene names are specified in gene symbols.
#' @param gmt filename
#'
#' @return
#' @export
#'
#' @examples
#'
#' re.gsa<-GSA.read.gmt.2("/media/H_driver/Annotation/MsigDB/c2.cp.Mouse.v5.1.symbols.gmt")
#'
#' re.ggsa.2<-GSA.read.gmt("/media/H_driver/Annotation/MsigDB/c2.cp.Mouse.v5.1.symbols.gmt")

GSA.read.gmt.2<-function (filename)
{
  a = scan(filename, what = list("", ""), sep = "\t", quote = NULL,
           fill = T, flush = T, multi.line = F)
  geneset.names = a[1][[1]]
  geneset.descriptions = a[2][[1]]
  dd = scan(filename, what = "", sep = "\t", quote = NULL)
  nn = length(geneset.names)
  n = length(dd)
  ox = rep(NA, nn)
  ii = 1
  for (i in 1:nn) {
    #cat(i)
    while ((dd[ii] != geneset.names[i]) | (dd[ii + 1] !=
                                           geneset.descriptions[i])) {
      ii = ii + 1
    }
    ox[i] = ii
    ii = ii + 1
  }
  genesets = vector("list", nn)
  for (i in 1:(nn - 1)) {
    #cat(i, fill = T)
    i1 = ox[i] + 2
    i2 = ox[i + 1] - 1
    geneset.descriptions[i] = dd[ox[i] + 1]
    genesets[[i]] = dd[i1:i2]
  }
  geneset.descriptions[nn] = dd[ox[nn] + 1]
  genesets[[nn]] = dd[(ox[nn] + 2):n]
  out = list(genesets = genesets, geneset.names = geneset.names,
             geneset.descriptions = geneset.descriptions)
  class(out) = "GSA.genesets"
  return(out)
}
#' GeneBasedAnalysis
#'
#' This function use gene-based DE to get GO
#'
#'
#' @param Re.PJ.gene: gene-based DE results
#' @param gene_model: gene model
#' @param output_file: file name for venn
#'
#' @return
#' @export
#'
#' @examples
#'
#' Gene.non.coding.GO<-GeneBasedAnalysis(Re.PJ.gene.rm.non.coding,gene.model,"/media/H_driver/PJ/geneGL_rm_non_coding.xls")
#'
#'
GeneBasedAnalysis<-function(Re.PJ.gene,gene_model,output_file){

  re2<-Re.PJ.gene

  gene.based.DE<-rownames(re2[which(re2$padj<0.05),])

  all.gene.index<-rep(0,length(rownames(re2)))

  names(all.gene.index)<-rownames(re2)

  all.gene.index.gene.based<-all.gene.index
  all.gene.index.gene.based[which(names(all.gene.index.gene.based) %in% gene.based.DE)]=1

  Re3=nullp(all.gene.index.gene.based,"mm10","ensGene",plot.fit = FALSE)

  GO.wall.DE_interest=goseq2(Re3,"mm10","ensGene",gene.model=gene_model,use_genes_without_cat=TRUE)

  OuputGO(GO.wall.DE_interest,output_file)

  re4<-list(pwf=Re3,GO=GO.wall.DE_interest)

  return(re4)

}
#' @title Compare gene based p value with exon based p value
#'
#' @param input_data_set data set for all genes
#' @param gene_wise_p_c genewise p value column
#' @param feature_wise_p_c feature(exon or SJ) p value column
#'
#' @return
#' @export
#'
#' @examples
#'
#' load("/media/H_driver/2015/Nimer_Cheng/Data_set_two_methods.RData")
#' GeneWisePvsExonWiseP(data.table.gene.based.all.jscs4,7,9)
#'
#'
#' Re<-GeneWisePvsExonWiseP("/media/H_driver/2015/Nimer_Cheng/GeneWise_jscs3_all_with_anno_2_24_2016.csv",
#'"/media/H_driver/2015/Nimer_Cheng/DE_cheng_output_sample1_DE.csv")
#'
#'
#'
#'
#'
#'
GeneWisePvsExonWiseP<-function(subfeature_based_input_file,gene_based_input_file){
#load(input_file)

  Re.Jun<-read.csv(subfeature_based_input_file)
  Re.DE<-read.csv(gene_based_input_file)

  print(head(Re.Jun))
  print(head(Re.DE))


  Re.Jun.2<-Re.Jun[,-1]

  colnames(Re.Jun.2)[1]="ensembl_gene_id"

  colnames(Re.DE)[1]="geneID"

  mart = useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", host="uswest.ensembl.org", dataset="mmusculus_gene_ensembl")
  geneMaps<-getBM(attributes=c("ensembl_gene_id","mgi_symbol"), filters="mgi_symbol",values=Re.DE[,1],mart=mart)

  colnames(geneMaps)=c("ensembl_gene_id","geneID")

  Re.DE.2<-merge(Re.DE,geneMaps,by="geneID")

  print(head(Re.DE.2))

  Re.DE.Jun<-merge(Re.DE.2,Re.Jun.2,by="ensembl_gene_id")


  print(dim(Re.DE.Jun))

  Re.DE.Jun.2<-Re.DE.Jun[-which(is.na(Re.DE.Jun[,8])),]

  n<-dim(Re.DE.Jun.2)[1]
  ran.p<-runif(n)


  Re<-list(Re.DE.Jun.with.NA=Re.DE.Jun,Re.DE.Jun.without.NA=cbind(Re.DE.Jun.2,ran.p))

  #plot()

  pairs(~padj+geneWisePadj+mostSigPadjust+numKnown+ran.p,data=Re[[2]][,c(8,22,24,26,33)],main="4 p values vs SJ")
  #boxplot(cbind(Re[[2]][,c(8,22,24)],ran.p))


  Re


}
#' @ Generate gene annotations based on the Ensembl ID of genes
#'
#' @param TobeAnno: gene list
#'
#' @return
#' @export
#'
#' @examples
#'
#' Cheng.gene.all.anno.3<-do.call(rbind,lapply(re.PJ.gene.based.testable[,1],GenerateGeneAnno,gene.model))
#'
#'
GenerateGeneAnno<-function(TobeAnno,data.set){
  #First element
  gene.split<-unlist(strsplit(TobeAnno,"\\+"))
  #print(gene.split)

  #use bioMart data base
  #gene.anno<-GetMgiSymbolDescription(gene.split[1])

  #use local file
  gene.anno<-GetMgiSymbolUsingLocalDataBase(gene.split[1],data.set)

  if(length(gene.split)>=2){
    for(j in 2:length(gene.split)){
      #gene.anno.temp=GetMgiSymbolDescription(gene.split[j])
      gene.anno.temp<-GetMgiSymbolUsingLocalDataBase(gene.split[j],data.set)
      gene.anno<-paste(gene.anno,gene.anno.temp,sep="+")
    }
  }

  Annotatedgene<-as.data.frame(cbind(TobeAnno,gene.anno))
  colnames(Annotatedgene)=c("geneID","geneAnno")

  return(Annotatedgene)

}

#' Title
#'
#' @param Re
#'
#' @return
#' @export
#'
#' @examples
#'
#'
GetCCBetweenMostSigPadjustSJ <- function(Re) {
  #take some times to get this calculate doine
  re.gene.based<-makeGeneWiseTable(Re,gene.list=unique(as.character(fData(Re)$geneID)))

  no.testable.index<-which(as.character(pData(re.gene.based)$mostSigID)=="character(0)")
  re.gene.based.testable<-pData(re.gene.based)[-no.testable.index,]

  cc<-cor(as.numeric(re.gene.based.testable$numKnown),
      as.numeric(re.gene.based.testable$mostSigPadjust))
  return(cc)
}



#' Title
#'
#' @param sample.name
#' @param condition
#' @param use.covars
#'
#' @return
#' @export
#'
#' @examples
#'
#' dir.name="/Volumes/Bioinformatics$/2015/Nimer_Cheng/"
#' file.sample="decoder.bySample.txt"
#' path.file.sample<-paste0(dir.name,file.sample)
#' decoder.bySample<-read.table(path.file.sample,header=T,stringsAsFactors = F)
#' print(decoder.bySample)
#' Re<-GetDesign(decoder.bySample$sample.ID,decoder.bySample$group.ID)
#'
GetDesign <- function(sample.names,condition,use.covars=NULL) {

  if(! is.factor(condition)){
    condition <- factor(condition, levels = sort(unique(condition)))
  }

  design <- data.frame(condition = condition)
  if(! is.null(use.covars)){
    message(paste0("> rJSA: using covars:"," ",date()))
    if(class(use.covars) != "data.frame"){
      stop(paste0("FATAL ERROR: use.covars must be a data.frame! Instead it appears to be: ",class(use.covars)))
    }
    for(i in 1:ncol(use.covars)){
      if(! is.factor(use.covars[[i]]) ){
        use.covars[[i]] <- factor(use.covars[[i]], levels = sort(unique(use.covars[[i]]))  )
      }
    }

    design <- data.frame(cbind(design,use.covars))
    for(i in 1:length(names(use.covars))){
      message(paste0("      covar: ",names(use.covars)[i]))
      message(paste0(c("      ",paste0(use.covars[,i],collapse=", "))))
    }
    names(design) <- c("condition",names(use.covars))
  }
  row.names(design) <- sample.names

  print(design)

}

#' Title
#'
#' @param data.gene
#'
#' @return
#' @export
#'
#' @examples
#'
#'
GetMgiSymbolDescription<-function(data.gene){

  gene.anno<-as.character(paste(getBM(attributes=c("mgi_symbol","description"), filters="ensembl_gene_id",
                                      values=data.gene,
                                      mart=mart)$mgi_symbol,
                                getBM(attributes=c("mgi_symbol","description"), filters="ensembl_gene_id",
                                      values=data.gene,
                                      mart=mart)$description,sep="->"))
  if(length(gene.anno)==0)
  {gene.anno="unmatched"}

  return(gene.anno)
}

#' Title
#'
#' @param data.gene
#'
#' @return
#' @export
#'
#' @examples
#'
#'
GetMgiSymbolUsingLocalDataBase<-function(data.gene,data.set){

  gene.model<-data.set

  gene.anno<-gene.model[which(gene.model[,3]==data.gene),1]

  if(length(gene.anno)==0)
  {gene.anno="unmatched"}

  return(gene.anno)
}
#' GetPwfUseReformatedData
#'
#' Use reformated data to calculate probability weight function
#'
#' @param re.gene.based
#' @param ad
#' @param sub_feature
#' @param threshold
#'
#' @return
#' @export
#'
#' @examples
#' Re.pwf.exon.sj<-GetPwfUseReformatedData(ReformatData(re.PJ.gene.based),ad="exon_SJ",sub_feature=NULL,0.05)
#'
GetPwfUseReformatedData<-function(re.gene.based,ad="GL",sub_feature=NULL,threshold){

  Data4Goterm<-re.gene.based

  if(is.null(sub_feature)){Data4Goterm.sub_feature<-Data4Goterm}
  else{Data4Goterm.sub_feature<-Data4Goterm[grep(sub_feature,Data4Goterm[,8]),]}

  Data4Goterm.sub_feature.geneID.NumOfJunctions<-Data4Goterm.sub_feature[,c(1,11)]

  Data4Goterm.sub_feature.Sig<-Data4Goterm.sub_feature[which(Data4Goterm.sub_feature[,7]<threshold),]

  All.gene.id.based.on.sub_feature<-unique(Data4Goterm.sub_feature[,1])
  All.gene.id.index<-rep(0,length(All.gene.id.based.on.sub_feature))
  names(All.gene.id.index)=All.gene.id.based.on.sub_feature

  All.genes.based.on.Sig.sub_feature<-unique(Data4Goterm.sub_feature.Sig[,1])
  gene.DE_interest<-as.integer(which( All.gene.id.based.on.sub_feature %in% All.genes.based.on.Sig.sub_feature ))

  All.gene.id.index[gene.DE_interest]<-1
  #print(length(All.gene.id.index))

  gene.with.matched.junction<-which(Data4Goterm.sub_feature.geneID.NumOfJunctions[,1] %in% c(names(All.gene.id.index)))
  num.junction.4.matched.gene<-as.numeric(Data4Goterm.sub_feature.geneID.NumOfJunctions[gene.with.matched.junction,2])

  All.gene.id.index.2<-All.gene.id.index

  #print(All.gene.id.index.2)

  if(ad=="GL"){
    pwf.DE_interest=nullp(All.gene.id.index.2,"mm10","ensGene",plot.fit = FALSE)
  }
  else
  {
    pwf.DE_interest=nullp(All.gene.id.index.2,"mm10","ensGene",bias.data = num.junction.4.matched.gene,plot.fit = FALSE)
  }

  re<-pwf.DE_interest

  return(re)
}
#' GetResults4GeneLevel
#'
#' This function applies DESeq2 to gene level count to identify differentially expressed genes
#'
#' @param dir.name: path that count file and sample information file is in
#' @param file.sample: sample information file
#' @param file.count: count files
#'
#'
#' @return results from DESeq2
#' @export
#'
#' @examples
#' # For PJ project
#' dir.name.PJ.gene="/media/H_driver/PJ/"
#' file.sample.PJ.gene="decoder.bySample.rtf"
#' file.count.PJ.gene="/QC.geneCounts.formatted.for.DESeq.txt"
#'
#' Re.PJ.gene<-GetResults4GeneLevel(dir.name.PJ.gene,file.sample.PJ.gene,file.count.PJ.gene)

GetResults4GeneLevel<-function(dir.name,file.sample,file.count){

  suppressPackageStartupMessages(library(DESeq2))

  path.file.sample<-paste0(dir.name,file.sample)
  decoder.bySample<-read.table(path.file.sample,header=T,stringsAsFactors = F)
  #print(decoder.bySample)

  sampleCondition <- decoder.bySample$group.ID;
  sampleName <- decoder.bySample$sample.ID;


  directory <- dir.name

  sampleFiles <- paste0(decoder.bySample$sample.ID,file.count);

  sampleTable <- data.frame(sampleName = sampleName,fileName = sampleFiles,condition = sampleCondition);

  #print(sampleTable)

  dds <-  DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,directory = directory,design = ~ condition)

  dds2 <- DESeq(dds);
  res <- results(dds2);

  return(res)

}
#' calculate the p-value for GO terms
#'
#' @param num_de_incat
#' @param num_incat
#' @param num_genes
#' @param num_de
#' @param weight
#'
#' @return
#' @export
#'
#' @examples
#'
#'28 	 186 	 1852 	 125 	 1.03043
#'
#'use the gene pooled from all GO terms
#'Get_Over_under_represented_pvalue(28,186,1852,125,1.03043)
#'
#'use the genes that have matched GO terms
#'
#'Get_Over_under_represented_pvalue(28,186,1840,125,1.026767)
#'
Get_Over_under_represented_pvalue<-function(num_de_incat,num_incat,num_genes,num_de,weight){
c(dWNCHypergeo(num_de_incat,num_incat,num_genes-num_incat,num_de,weight)
  +pWNCHypergeo(num_de_incat,num_incat,num_genes-num_incat,num_de,weight,lower.tail=FALSE),
  pWNCHypergeo(num_de_incat,num_incat,num_genes-num_incat,num_de,weight))
}
#' Gmt2GeneCat
#'
#' Read a gmt file, and return a list with the name of element being a gene id based on gene_anno_file, and each element
#' being the pathways that this gene corresponds to
#'
#' @param gmt_input_file
#' @param gene_anno_file
#' @param based_by
#'
#' @return
#' @export
#'
#' @examples
#'
#' gene.2.cat.cp.mouse<-Gmt2GeneCat("/media/H_driver/Annotation/MsigDB/c2.cp.Mouse.v5.1.symbols.gmt",
#' "/media/H_driver/Annotation/mm10/genes_table_02052016.csv")
#'
#' gene.2.cat.tft.mouse<-Gmt2GeneCat("/media/H_driver/Annotation/MsigDB/c3.tft.Mouse.v5.1.symbols.gmt",
#' "/media/H_driver/Annotation/mm10/genes_table_02052016.csv")
#'

Gmt2GeneCat <- function(gmt_input_file,gene_anno_file) {

  gene.2.cat.gmt<-gene2cat2(gmt_input_file)

  names.gene.gmt<-as.data.frame(names(gene.2.cat.gmt))
  colnames(names.gene.gmt)<-"gene_id"

  gene.ID.conversion<-read.csv(gene_anno_file)
  names.gene.gmt.2<-match(names.gene.gmt$gene_id,gene.ID.conversion$gene_id)
  gene.ID.conversion.2<-gene.ID.conversion[names.gene.gmt.2,]
  gene.2.cat.gmt.2<-gene.2.cat.gmt

  names(gene.2.cat.gmt.2)<-gene.ID.conversion.2[,3]

  gene.2.cat.gmt.2
}
#'This function is to subset data set based on feature
#'
#'
#' @param Data4Goterm
#' @param ad
#' @param sub_feature
#' @param threshold
#'
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#' Re.Go.adjusted.by.number.junction.2<-GotermAnalysis2GetAllGOTerms_AdjustedByNumOfJunctionWithinOneGene(
#' re.PJ.gene.based,ad="J","J",0.05,"Splice_junction_based")
#'
#' data.pwf2.SJs<-plotPWF2(Re.Go.adjusted.by.number.junction.2[[2]],binsize=30,xlab = "Number of SJs(<binsize> gene bins)")
#'
#' Re.Go.adjusted.by.exon.SJ<-GotermAnalysis2GetAllGOTerms_AdjustedByNumOfJunctionWithinOneGene(
#' re.PJ.gene.based,ad="exon_SJ",sub_feature=NULL,0.05,file_prefix="Exon_Splice_junction_based.xls",gene_model=gene.model)
#'
#'
#'
GotermAnalysis2GetAllGOTerms_AdjustedByNumOfJunctionWithinOneGene<-function(re.gene.based,ad="GL",sub_feature=NULL,threshold,gene_model){

  #Data4Goterm<-pData(re.gene.based)

  Data4Goterm<-re.gene.based

  if(is.null(sub_feature)){Data4Goterm.sub_feature<-Data4Goterm}
  else{Data4Goterm.sub_feature<-Data4Goterm[grep(sub_feature,Data4Goterm[,8]),]}


  print(dim(Data4Goterm.sub_feature))

  Data4Goterm.sub_feature.geneID.NumOfJunctions<-Data4Goterm.sub_feature[,c(1,11)]

  print(dim(Data4Goterm.sub_feature.geneID.NumOfJunctions))

  Data4Goterm.sub_feature.Sig<-Data4Goterm.sub_feature[which(Data4Goterm.sub_feature[,7]<threshold),]

  #GO term analysis using GOSeq
  All.gene.id.based.on.sub_feature<-unique(unlist(strsplit(Data4Goterm.sub_feature[,1],"\\+")))
  #length(All.gene.id.based.on.sub_feature)
  All.gene.id.index<-rep(0,length(All.gene.id.based.on.sub_feature))
  names(All.gene.id.index)=All.gene.id.based.on.sub_feature

  All.genes.based.on.Sig.sub_feature<-unique(unlist(strsplit(Data4Goterm.sub_feature.Sig[,1],"\\+")))
  gene.DE_interest<-as.integer(which( All.gene.id.based.on.sub_feature %in% All.genes.based.on.Sig.sub_feature ))

  All.gene.id.index[gene.DE_interest]<-1
  print(length(All.gene.id.index))

  gene.with.matched.junction<-which(Data4Goterm.sub_feature.geneID.NumOfJunctions[,1] %in% c(names(All.gene.id.index)))
  num.junction.4.matched.gene<-Data4Goterm.sub_feature.geneID.NumOfJunctions[gene.with.matched.junction,2]

  names.4.matched.gene<-Data4Goterm.sub_feature.geneID.NumOfJunctions[gene.with.matched.junction,1]

  All.gene.id.index.2<-All.gene.id.index[which(names(All.gene.id.index) %in% c(names.4.matched.gene))]

  print(length(All.gene.id.index.2))

  print(All.gene.id.index.2)

if(ad=="GL"){
pwf.DE_interest=nullp(All.gene.id.index.2,"mm10","ensGene",plot.fit = FALSE)
}
else
{
pwf.DE_interest=nullp(All.gene.id.index.2,"mm10","ensGene",bias.data = num.junction.4.matched.gene,plot.fit = FALSE)
}

  #GO.wall.DE_interest=goseq2(pwf.DE_interest,"mm10","ensGene",gene.model=gene_model,use_genes_without_cat=TRUE)

  GO.wall.DE_interest=goseq2(pwf.DE_interest,"mm10","ensGene",gene.model=gene_model)
  enriched.GO.DE_interest=GO.wall.DE_interest[p.adjust(GO.wall.DE_interest$over_represented_pvalue,method="BH")<threshold,]
  re<-list()

  re[[1]]<-GO.wall.DE_interest
  re[[2]]<-pwf.DE_interest
  return(re)

  }
#' GotermAnalysisUseFeatureDefineDE
#' Descrption: This function computes statistical significance of GO terms, given a list of significant genes. These significant genes can be directly from differntial splicing analysis software such as rMat, DEXSeq or JunctionSeq. 

#' Based on genewise results to perform Go term analysis using DGEs defined by different criterion
#' @param Data4Goterm
#' @param ad
#' @param sub_feature
#' @param gene_model
#'
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#' Re.Go.adjusted.by.number.junction.2<-GotermAnalysis2GetAllGOTerms_AdjustedByNumOfJunctionWithinOneGene(
#' re.PJ.gene.based,ad="J","J",0.05,"Splice_junction_based")
#'
#' data.pwf2.SJs<-plotPWF2(Re.Go.adjusted.by.number.junction.2[[2]],binsize=30,xlab = "Number of SJs(<binsize> gene bins)")
#'
#' Re.Go.adjusted.by.exon.SJ<-GotermAnalysis2GetAllGOTerms_AdjustedByNumOfJunctionWithinOneGene(
#' re.PJ.gene.based,ad="exon_SJ",sub_feature=NULL,0.05,file_prefix="Exon_Splice_junction_based.xls",gene_model=gene.model)
#'
GotermAnalysisUseFeatureDefineDE<-function(re.gene.based,ad="GL",sub_feature=NULL,DE_define=c("Feature","GeneWise","rMAT","FeatureGeneWise","FeaturerMAT","GeneWiserMAT","FeatureGeneWiseRMAT"),gene_model,Output_file_dir){

  #Data4Goterm<-pData(re.gene.based)

  Data4Goterm<-re.gene.based

  if(is.null(sub_feature)){Data4Goterm.sub_feature<-Data4Goterm}
  else{Data4Goterm.sub_feature<-Data4Goterm[grep(sub_feature,Data4Goterm[,8]),]}

  #print(dim(Data4Goterm.sub_feature))

  Data4Goterm.sub_feature.geneID.NumOfJunctions<-Data4Goterm.sub_feature[,c(1,11)]

  #print(dim(Data4Goterm.sub_feature.geneID.NumOfJunctions))

  #Data4Goterm.sub_feature.Sig<-Data4Goterm.sub_feature[which(Data4Goterm.sub_feature[,7]<threshold),]

  if(DE_define=="Feature"){
  Data4Goterm.sub_feature.Sig<-Data4Goterm.sub_feature[which(Data4Goterm.sub_feature$DE_or_not_feature==1),]
  }else if(DE_define=="GeneWise"){
  Data4Goterm.sub_feature.Sig<-Data4Goterm.sub_feature[which(Data4Goterm.sub_feature$DE_or_not_geneWise==1),]
  }else if(DE_define=="rMAT"){
    Data4Goterm.sub_feature.Sig<-Data4Goterm.sub_feature[which(Data4Goterm.sub_feature$DE_or_not_rMAT_based==1),]
  }else if(DE_define=="FeatureGeneWise"){
  Data4Goterm.sub_feature.Sig<-Data4Goterm.sub_feature[which(Data4Goterm.sub_feature$DE_or_not_feature==1&Data4Goterm.sub_feature$DE_or_not_geneWise==1),]
  }else if(DE_define=="FeaturerMAT"){
    Data4Goterm.sub_feature.Sig<-Data4Goterm.sub_feature[which(Data4Goterm.sub_feature$DE_or_not_feature==1&
                                                                 Data4Goterm.sub_feature$DE_or_not_rMAT_based==1),]
  }else if(DE_define=="GeneWiserMAT"){
    Data4Goterm.sub_feature.Sig<-Data4Goterm.sub_feature[which(Data4Goterm.sub_feature$DE_or_not_geneWise==1&
                                                               Data4Goterm.sub_feature$DE_or_not_rMAT_based==1),]
  }else if(DE_define=="FeatureGeneWiseRMAT"){
  Data4Goterm.sub_feature.Sig<-Data4Goterm.sub_feature[which(Data4Goterm.sub_feature$DE_or_not_feature==1&
                                                               Data4Goterm.sub_feature$DE_or_not_geneWise==1&
                                                               Data4Goterm.sub_feature$DE_or_not_rMAT_based==1),]}

  #GO term analysis using GOSeq

  All.gene.id.based.on.sub_feature<-unique(Data4Goterm.sub_feature[,1])

  cat("How many DE genes?","\n")
  cat(dim(Data4Goterm.sub_feature.Sig)[1],"\n")
  cat(length(unique(as.character(Data4Goterm.sub_feature.Sig$geneID))),"\n")

  All.gene.id.index<-rep(0,length(All.gene.id.based.on.sub_feature))
  names(All.gene.id.index)=All.gene.id.based.on.sub_feature

  All.genes.based.on.Sig.sub_feature<-unique(Data4Goterm.sub_feature.Sig[,1])
  gene.DE_interest<-as.integer(which( All.gene.id.based.on.sub_feature %in% All.genes.based.on.Sig.sub_feature ))

  All.gene.id.index[gene.DE_interest]<-1
  #print(length(All.gene.id.index))

  gene.with.matched.junction<-which(Data4Goterm.sub_feature.geneID.NumOfJunctions[,1] %in% c(names(All.gene.id.index)))
  num.junction.4.matched.gene<-as.numeric(Data4Goterm.sub_feature.geneID.NumOfJunctions[gene.with.matched.junction,2])

  #names.4.matched.gene<-Data4Goterm.sub_feature.geneID.NumOfJunctions[gene.with.matched.junction,1]

  #All.gene.id.index.2<-All.gene.id.index[which(names(All.gene.id.index) %in% c(names.4.matched.gene))]

  #print(length(All.gene.id.index.2))

  All.gene.id.index.2<-All.gene.id.index

  #print(All.gene.id.index.2)

  if(ad=="GL"){
    pwf.DE_interest=nullp(All.gene.id.index.2,"mm10","ensGene",plot.fit = FALSE)
  }
  else
  {
    pwf.DE_interest=nullp(All.gene.id.index.2,"mm10","ensGene",bias.data = num.junction.4.matched.gene,plot.fit = FALSE)
  }

  GO.wall.DE_interest=goseq2(pwf.DE_interest,"mm10","ensGene",gene_model,test.cats=c("GO:BP"),use_genes_without_cat=TRUE)

  #GO.wall.DE_interest=goseq(pwf.DE_interest,"mm10","ensGene",test.cats=c("GO:BP"),use_genes_without_cat=TRUE)

  #GO.wall.DE_interest=goseq2(pwf.DE_interest,"mm10","ensGene",gene.model=gene_model)
  #enriched.GO.DE_interest=GO.wall.DE_interest[p.adjust(GO.wall.DE_interest$over_represented_pvalue,method="BH")<threshold,]
  re<-list()

  re[[1]]<-GO.wall.DE_interest[[1]]
  re[[2]]<-pwf.DE_interest
  re[[3]]<-GO.wall.DE_interest[[2]]


  DE.gene.symbol<-gene.model[which(as.character(gene.model[,3]) %in% unique(as.character(Data4Goterm.sub_feature.Sig$geneID))),1]

  Re4_temp<-list(A=DE.gene.symbol,B=re[[3]])

  venn.plot <- venn.diagram(
    x = Re4_temp[c(1,2)],
    filename = paste0(Output_file_dir,DE_define,"_overlap_with_DE_from_GO.tiff"),
    col = "black",
    lty = "dotted",
    lwd = 2,
    fill = c("red","blue"),
    alpha = 0.50,
    label.col = c(rep("white",3)),
    cex = 1,
    fontfamily = "serif",
    fontface = "bold",
    cat.col = c("red","blue"),
    cat.cex = 0.8,
    cat.fontfamily = "serif"
  )

  return(re)
}
#' GotermAnalysisUseReformatedData
#'
#'
#' @param Data4Goterm
#' @param ad
#' @param sub_feature
#' @param threshold
#'
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#' Re.Go.adjusted.by.number.junction.2<-GotermAnalysis2GetAllGOTerms_AdjustedByNumOfJunctionWithinOneGene(
#' re.PJ.gene.based,ad="J","J",0.05,"Splice_junction_based")
#'
#' data.pwf2.SJs<-plotPWF2(Re.Go.adjusted.by.number.junction.2[[2]],binsize=30,xlab = "Number of SJs(<binsize> gene bins)")
#'
#' Re.Go.adjusted.by.exon.SJ<-GotermAnalysis2GetAllGOTerms_AdjustedByNumOfJunctionWithinOneGene(
#' re.PJ.gene.based,ad="exon_SJ",sub_feature=NULL,0.05,file_prefix="Exon_Splice_junction_based.xls",gene_model=gene.model)
#'
GotermAnalysisUseReformatedData<-function(re.gene.based,ad="GL",sub_feature=NULL,threshold,genomeID,geneID,gene_model){

  #Data4Goterm<-pData(re.gene.based)

  Data4Goterm<-re.gene.based

  if(is.null(sub_feature)){Data4Goterm.sub_feature<-Data4Goterm}
  else{Data4Goterm.sub_feature<-Data4Goterm[grep(sub_feature,Data4Goterm[,8]),]}


  #print(dim(Data4Goterm.sub_feature))

  Data4Goterm.sub_feature.geneID.NumOfJunctions<-Data4Goterm.sub_feature[,c(1,11)]

  #print(dim(Data4Goterm.sub_feature.geneID.NumOfJunctions))

  Data4Goterm.sub_feature.Sig<-Data4Goterm.sub_feature[which(Data4Goterm.sub_feature[,7]<threshold),]

  #GO term analysis using GOSeq
  All.gene.id.based.on.sub_feature<-unique(Data4Goterm.sub_feature[,1])
  #length(All.gene.id.based.on.sub_feature)
  All.gene.id.index<-rep(0,length(All.gene.id.based.on.sub_feature))
  names(All.gene.id.index)=All.gene.id.based.on.sub_feature

  All.genes.based.on.Sig.sub_feature<-unique(Data4Goterm.sub_feature.Sig[,1])
  gene.DE_interest<-as.integer(which( All.gene.id.based.on.sub_feature %in% All.genes.based.on.Sig.sub_feature ))

  All.gene.id.index[gene.DE_interest]<-1
  print(length(All.gene.id.index))

  gene.with.matched.junction<-which(Data4Goterm.sub_feature.geneID.NumOfJunctions[,1] %in% c(names(All.gene.id.index)))
  num.junction.4.matched.gene<-as.numeric(Data4Goterm.sub_feature.geneID.NumOfJunctions[gene.with.matched.junction,2])

  #names.4.matched.gene<-Data4Goterm.sub_feature.geneID.NumOfJunctions[gene.with.matched.junction,1]

  #All.gene.id.index.2<-All.gene.id.index[which(names(All.gene.id.index) %in% c(names.4.matched.gene))]

  #print(length(All.gene.id.index.2))

  All.gene.id.index.2<-All.gene.id.index

  print(All.gene.id.index.2)

  if(ad=="GL"){
    pwf.DE_interest=nullp(All.gene.id.index.2,genomeID,geneID,plot.fit = FALSE)
  }
  else
  {
    pwf.DE_interest=nullp(All.gene.id.index.2,genomeID,geneID,bias.data = num.junction.4.matched.gene,plot.fit = FALSE)
  }

  GO.wall.DE_interest=goseq2(pwf.DE_interest,genomeID,geneID,gene.model=gene_model,use_genes_without_cat=TRUE)

  #GO.wall.DE_interest=goseq2(pwf.DE_interest,"mm10","ensGene",gene.model=gene_model)
  enriched.GO.DE_interest=GO.wall.DE_interest[p.adjust(GO.wall.DE_interest$over_represented_pvalue,method="BH")<threshold,]
  re<-list()

  re[[1]]<-GO.wall.DE_interest
  re[[2]]<-pwf.DE_interest

  return(re)
}
# installed.packages()
# install.packages("RCurl")
# source("http://bioconductor.org/biocLite.R")
# source("https://bioconductor.org/biocLite.R")
# biocLite("biomaRt")
# library(biomaRt)
# install.packages("roxygen2")
# library(roxygen2)
# biocLite("JunctionSeq")
# library(JunctionSeq)
# biocLite("goseq")
# library(goseq)
# biocLite("org.Mm.eg.db")
# library(org.Mm.eg.db)
# biocLite("org.Hs.eg.db")
# library(org.Hs.eg.db)

# require(GO.db)
# biocLite("geneLenDataBase")
# install.packages("geneLenDataBase")
# require(geneLenDataBase)
# install.packages("GSA")
# library(GSA)
# install.packages("VennDiagram")
# library(VennDiagram)
# install.packages("openssl")
# biocLite("cummeRbund")
# library(cummeRbund)
# install.packages("popbio")
# library(popbio)
# source("https://bioconductor.org/biocLite.R")
# biocLite("affycoretools")
# library(affycoretools)
# source("https://bioconductor.org/biocLite.R")
# biocLite("ReportingTools")
# library(ReportingTools)
# biocLite("hgu95av2.db")
# library(hgu95av2.db")
# install.packages("devtools")
# library(devtools)
# install.packages("rgl")
# library(rgl)
# install.packages("rgl")
# install.packages("qpcR")
# library(qpcR)
# biocLite("ggbio")

# gene.model<-read.table("/media/H_driver/Annotation/mm10/genes_table_02052016.csv",header = TRUE, sep = ",", as.is=TRUE)
# save(gene.model,file="~/GOSJ/data/gene.model.RData")

# gene.model.hg38<-read.table("/media/H_driver/Annotation/hg38/genes_table_02092016.csv",header = TRUE, sep = ",", as.is=TRUE)
# save(gene.model.hg38,file="~/GOSJ/data/gene.model.hg38.RData")

#' JS.perGeneQValue
#'
#' @param pvals
#' @param wTest
#' @param geneID
#' @param method
#'
#' @return
#' @export
#'
#' @examples
#'

#'data.perGeneQ.use.large.FC.featutre<-JS.perGeneQValue(Re.PJ.selected$pvalue,
#'Re.PJ.selected$testable,Re.PJ.selected$geneID,method = JS.perGeneQValueExact)
#'
#'length(unique(as.character(Re.PJ.selected$geneID)))
#'
#'head(Re.PJ.selected)
#'
#'
#'data.perGeneQ.use.large.FC.featutre.2<-list_to_df(data.perGeneQ.use.large.FC.featutre)
#'colnames(data.perGeneQ.use.large.FC.featutre.2)<-c("geneID","geneWisePadj_FC")
#'
#'Re.PJ.selected.2<-merge(Re.PJ.selected,data.perGeneQ.use.large.FC.featutre.2,by="geneID")
#'
#'
#'head(Re.PJ.selected.2)
#'
#'re.PJ.gene.based.selected<-merge(pData(re.PJ.gene.based),data.perGeneQ.use.large.FC.featutre.2,by="geneID")

JS.perGeneQValue = function(pvals, wTest, geneID, method = JS.perGeneQValueExact) {

  ## use only those exons that were testable
  pvals     = pvals[wTest]
  ## 'factor' removes ununsed levels
  cat(geneID,file="geneID_0.txt",sep="\n",append=TRUE)

  cat(wTest,file="wTest.txt",sep="\n",append=TRUE)

  cat(geneID[wTest],file="geneID_wTest.txt",sep="\n",append=TRUE)

  geneID    = factor(geneID[wTest])

  cat(geneID,file="geneID.txt",sep="\n",append=TRUE)

  geneSplit = split(seq(along=geneID), geneID)

  cat(unlist(geneSplit),file="geneSplit.txt",sep="\n",append=TRUE)

  ## summarise p-values of exons for one gene: take the minimum
  pGene = sapply(geneSplit, function(i) min(pvals[i]))

  #print(pGene)
  cat(pGene,file="pGene.txt",sep="\n",append=TRUE)

  stopifnot(all(is.finite(pGene)))

  ## Determine the thetas to be used
  theta = unique(sort(pGene))

  cat(theta,file="theta.txt",sep="\n",append=TRUE)

  ## compute q-values associated with each theta
  q = method(pGene, theta, geneSplit)

  ## return a named vector of q-values per gene
  res        = rep(NA_real_, length(pGene))
  res        = q[match(pGene, theta)]
  res = pmin(1, res)
  names(res) = names(geneSplit)
  #stopifnot(!any(is.na(res)))
  return(res)
}
#' JS.perGeneQValueExact
#'
#' Exact computation - see methods part of the paper
#'
#'
#' @param pGene:p value
#' @param theta:cutoff
#' @param geneSplit:split
#'
#' @return
#' @export
#'
#' @examples
#'
#'
JS.perGeneQValueExact = function(pGene, theta, geneSplit) {
  stopifnot(length(pGene)==length(geneSplit))

  ## Compute the numerator \sum_{i=1}^M 1-(1-theta)^{n_i}
  ## Below we first identify the summands which are the same
  ## (because they have the same n_i), then do the sum via the
  ## mapply
  numExons     = listLen(geneSplit)
  tab          = tabulate(numExons)
  notZero      = (tab>0)
  numerator    = mapply(function(m, n) m * (1 - (1-theta)^n),
                        m = tab[notZero],
                        n = which(notZero))
  numerator    = rowSums(numerator)

  ## Compute the denominator: for each value of theta, the number
  ## of genes with pGene <= theta[i].
  ## Note that in cut(..., right=TRUE), the intervals are
  ## right-closed (left open) intervals.
  bins   = cut(pGene, breaks=c(-Inf, as.vector(theta)), right = TRUE, include.lowest = TRUE)
  counts = tabulate(bins, nbins = nlevels(bins))
  denom  = cumsum(counts)
  stopifnot(denom[length(denom)]==length(pGene))

  return(numerator/denom)
}
##--------------------------------------------------
## Exact computation - see methods part of the paper
##---------------------------------------------------
#' Title
#'
#' @param pGene
#' @param theta
#' @param geneSplit
#'
#' @return
#' @export
#'
#' @examples
#' source("https://bioconductor.org/biocLite.R")
#' biocLite("Biobase")
#' library(Biobase)
#'
JS.perGeneQValueExact.test = function(pGene, theta, geneSplit) {
  stopifnot(length(pGene)==length(geneSplit))

  ## Compute the numerator \sum_{i=1}^M 1-(1-theta)^{n_i}
  ## Below we first identify the summands which are the same
  ## (because they have the same n_i), then do the sum via the
  ## mapply
  numExons     = listLen(geneSplit)



  cat("numExons\n")

  cat("How many unique\n")

  cat(numExons,file="exon_all.txt","\n")

  cat(length(unique(numExons)),"\n")

  cat(unique(numExons),file="exon.txt","\n")
  cat(sort(unique(numExons)),file="exon2.txt")

  print(numExons)

  tab          = tabulate(numExons)

  cat("tab\n")
  print(tab)

  notZero      = (tab>0)

  cat("tab[notZero]\n")
  print(tab[notZero])
  cat(tab[notZero],file="mm.txt")

  cat("which(notZero)\n")
  print(which(notZero))
  cat(which(notZero),file="nn.txt")
  cat(sort(which(notZero)),file="nn2.txt")


  cat("theta\n")
  print(theta)
  cat(theta,file="theta.txt")

  numerator    = mapply(function(m, n) {

    cat("m\n")
    print(m)

    cat("n\n")
    print(n)

    cat("theta","\n")
    cat(length(theta),"\n")
    cat(length(unique(theta)),"\n")
    print(theta)

    re<-m * (1 - (1-theta)^n)

    cat("re\n")
    print(re)
    cat(re,file="re.txt")

    re

    },
                        m = tab[notZero],
                        n = which(notZero))

  print(numerator)

  cat(dim(numerator))

  numerator    = rowSums(numerator)

  ## Compute the denominator: for each value of theta, the number
  ## of genes with pGene <= theta[i].
  ## Note that in cut(..., right=TRUE), the intervals are
  ## right-closed (left open) intervals.

  bins   = cut(pGene, breaks=c(-Inf, as.vector(theta)), right = TRUE, include.lowest = TRUE)

  counts = tabulate(bins, nbins = nlevels(bins))
  denom  = cumsum(counts)
  stopifnot(denom[length(denom)]==length(pGene))

  return(numerator/denom)

}
#' LabelGeneBasedFeature
#'
#' @param Re.PJ
#'
#' @return
#' @export
#'
#' @examples
#' Re.PJ.selected<-SelectFeature(Re.PJ)
#'
#'dim(Re.PJ.selected)
#'
#'LabelGeneBasedFeature(Re.PJ,cutoff_FC,cutoff_p,Output_file)
#'
#'
LabelGeneBasedFeature<-function(Re.PJ,cutoff_FC,cutoff_p,Output_file){

  feature.based.Re.PJ<-fData(Re.PJ)

  re2<-feature.based.Re.PJ[which(feature.based.Re.PJ[,20]>cutoff_FC&feature.based.Re.PJ[,11]<cutoff_p),]

  DE.gene<-unique(re2$geneID)

  DE_or_not<-rep(0,dim(feature.based.Re.PJ)[1])

  re3<-cbind(feature.based.Re.PJ,DE_or_not)

  re3[which(re3$geneID %in% DE.gene),]$DE_or_not<-1

  dataset2<- re3

  dataset2[sapply(dataset2, is.list)] <-
    sapply(dataset2[sapply(dataset2, is.list)],
           function(x)sapply(x, function(y) paste(unlist(y),collapse=", ") ) )

  write.table(dataset2,file=Output_file,row.names = FALSE,quote=FALSE,sep="\t")


  return(re3)

}
#' This function is to generate the gene-based results
#'
#' @param jscs
#' @param gene.list
#' @param FDR.threshold
#' @param verbose
#' @param debug.mode
#'
#' @return
#' @export
#'
#' @examples
#'
#' # Generate genewise table for PJ data
#' re.PJ.gene.based<-makeGeneWiseTable(Re.PJ,gene.list=unique(as.character(fData(Re.PJ)$geneID)))
#'
makeGeneWiseTable <- function(jscs, gene.list, FDR.threshold = 0.05, verbose = TRUE, debug.mode = FALSE){
  if(verbose) message("   Compiling data table. ",date())

  mainTable <- data.frame(geneID = as.character(gene.list), stringsAsFactors=FALSE)
  row.names(mainTable) <- gene.list
  mainTable <- AnnotatedDataFrame(mainTable)
  varMetadata(mainTable)["geneID", "labelDescription"] <- "Gene Unique Identifier"

  noGenes <- (length(gene.list) == 0)

  if(! noGenes){
    geneAnno <- as.data.frame(t(sapply(gene.list, function(g){
      geneRows <- which(fData(jscs)$geneID == g)
      c(as.character(fData(jscs)$chr[geneRows[1]]),
        as.numeric(min(fData(jscs)$start[geneRows])),
        as.numeric(max(fData(jscs)$end[geneRows])),
        as.character(fData(jscs)$strand[geneRows[1]])
      )
    })))
    colnames(geneAnno) <- c("chr","start","end","strand")
  } else {
    geneAnno <- data.frame(chr = character(), start = numeric(), end = numeric(), strand = character())
  }

  mainTable$chr <- as.character(geneAnno$chr)
  mainTable$start <- geneAnno$start
  mainTable$end <- geneAnno$end
  mainTable$strand <- geneAnno$strand
  varMetadata(mainTable)[c("chr","start","end","strand"), "labelDescription"] <-
    c("Gene chromosome",
      "Gene start",
      "Gene end",
      "Gene strand")

  #message("2")
  geneBaseMeans <- if(noGenes){numeric()} else { rowMeans(jscs@geneCountData[match(gene.list,rownames(jscs@geneCountData)),, drop=FALSE] / sizeFactors(jscs))}
  mainTable$baseMean <- if(noGenes){character()} else {sprintf("%.1f",geneBaseMeans)}
  varMetadata(mainTable)["baseMean", "labelDescription"] <- "Gene BaseMean (simple normalized mean read or read-pair count per sample)"

  if(! is.null(fData(jscs)$geneWisePadj)){
    mainTable$geneWisePadj <- if(noGenes){ numeric()} else {sapply(gene.list, function(g){
      geneRows <- which(fData(jscs)$geneID == g)
      min( fData(jscs)$geneWisePadj[geneRows] , na.rm = TRUE)
    })}
    varMetadata(mainTable)["geneWisePadj", "labelDescription"] <- "Gene-level adjusted p-value. P-value for the hypothesis that one or more features are DU."
  }

  mainTable$mostSigID <- if(noGenes){character()} else {sapply(gene.list, function(g){
    geneRows <- which(fData(jscs)$geneID == g)
    fData(jscs)$countbinID[ geneRows[which.min( fData(jscs)$padjust[geneRows])] ]
  })}
  varMetadata(mainTable)["mostSigID", "labelDescription"] <- "Feature ID of the most singificant feature."

  mainTable$mostSigPadjust <- if(noGenes){numeric()} else {sapply(gene.list, function(g){
    geneRows <- which(fData(jscs)$geneID == g)
    fData(jscs)$padjust[ geneRows[which.min( fData(jscs)$padjust[geneRows])] ]
  })}
  mainTable$mostSigPadjust <- sprintf("%.3g",mainTable$mostSigPadjust)
  varMetadata(mainTable)["mostSigPadjust", "labelDescription"] <- "Adjusted p-value of the most singificant feature."

  gene.row.list <- if(noGenes){ integer() } else {lapply(gene.list, function(g){  which(fData(jscs)$geneID == g) })}

  mainTable$numExons <- if(noGenes){ integer() } else { sapply(gene.row.list, function(geneRows){
    sum(fData(jscs)$featureType[geneRows] == "exonic_part", na.rm = TRUE)
  })}
  varMetadata(mainTable)["numExons", "labelDescription"] <- "Number of distinct exonic regions belonging to the gene."

  mainTable$numKnown <- if(noGenes){ integer() } else {sapply(gene.row.list, function(geneRows){
    sum(fData(jscs)$featureType[geneRows] == "splice_site", na.rm = TRUE)
  })}
  varMetadata(mainTable)["numKnown", "labelDescription"] <- "Number of distinct known splice sites belonging to the gene."
  mainTable$numNovel <- if(noGenes){ integer() } else {sapply(gene.row.list, function(geneRows){
    sum(fData(jscs)$featureType[geneRows] == "novel_splice_site", na.rm = TRUE)
  })}
  varMetadata(mainTable)["numNovel", "labelDescription"] <- "Number of distinct novel splice sites belonging to the gene."

  mainTable$exonsSig <- if(noGenes){ integer() } else {sapply(gene.row.list, function(geneRows){
    sum(fData(jscs)$padjust[geneRows] < FDR.threshold & fData(jscs)$featureType[geneRows] == "exonic_part", na.rm = TRUE)
  })}
  varMetadata(mainTable)["exonsSig", "labelDescription"] <- paste0("Number of signficant exonic regions at p-adjust < ", FDR.threshold)
  mainTable$knownSig <- if(noGenes){ integer() } else {sapply(gene.row.list, function(geneRows){
    sum(fData(jscs)$padjust[geneRows] < FDR.threshold & fData(jscs)$featureType[geneRows] == "splice_site", na.rm = TRUE)
  })}
  varMetadata(mainTable)["knownSig", "labelDescription"] <- paste0("Number of signficant known splice junctions at p-adjust < ", FDR.threshold)
  mainTable$novelSig <- if(noGenes){ integer() } else {sapply(gene.row.list, function(geneRows){
    sum(fData(jscs)$padjust[geneRows] < FDR.threshold & fData(jscs)$featureType[geneRows] == "novel_splice_site", na.rm = TRUE)
  })}
  varMetadata(mainTable)["novelSig", "labelDescription"] <- paste0("Number of signficant novel splice junctions at p-adjust < ", FDR.threshold)

  mainTable$numFeatures = if(noGenes){ character()} else {paste0(mainTable$numExons,"/",mainTable$numKnown,"/",mainTable$numNovel)}
  varMetadata(mainTable)["numFeatures", "labelDescription"] <- "Number exonic regions / num known SJ / num novel SJ"
  mainTable$numSig = if(noGenes){ character()} else {paste0(mainTable$exonsSig,"/",mainTable$knownSig,"/",mainTable$novelSig)}
  varMetadata(mainTable)["numSig", "labelDescription"] <- "Number sig exonic regions / num sig known SJ / num sig novel SJ"

  return(mainTable)
}
#' Generate simulated count based on NB distribution
#'
#' @param dataset
#' @param group
#' @param nTags
#' @param nlibs
#' @param fix.dispersion
#' @param lib.size
#' @param drop.low.lambda
#' @param drop.extreme.dispersion
#' @param add.outlier
#' @param outlierMech
#' @param pOutlier
#' @param min.factor
#' @param max.factor
#' @param pDiff
#' @param pUp
#' @param foldDiff
#' @param name
#' @param save.file
#' @param file
#' @param only.add.outlier
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
NBsim <-
  function(dataset, group, nTags = 10000, nlibs = length(group), fix.dispersion = NA, lib.size = NULL, drop.low.lambda = TRUE, drop.extreme.dispersion = 0.1,  add.outlier = FALSE, outlierMech = c("S", "R", "M"), pOutlier = 0.1, min.factor = 1.5, max.factor = 10, pDiff=.1, pUp=.5, foldDiff=3, name = NULL, save.file = FALSE, file = NULL, only.add.outlier = FALSE, verbose=TRUE)

  {
    ## NBsim generate simulated count from the real dataset followed by the NB model ##
    require(edgeR)
    group = as.factor(group)

    sample.fun <- function(object)
    {
      ## it is low-level function of NBsim ##
      ## it samples from the real dataset ##
      nlibs <- object$nlibs
      nTags <- object$nTags
      AveLogCPM <-object$dataset$dataset.AveLogCPM
      dispersion <- object$dataset$dataset.dispersion

      id_r <- sample(length(AveLogCPM), nTags, replace = TRUE)
      Lambda <- 2^(AveLogCPM[id_r])
      Lambda <- Lambda/sum(Lambda)
      Dispersion <- dispersion[id_r]
      id_0<- Lambda == 0
      Lambda <- Lambda[!id_0]
      Dispersion <- Dispersion[!id_0]
      Lambda <- expandAsMatrix(Lambda, dim = c(nTags, nlibs))
      object$Lambda <- Lambda
      if(!is.na(fix.dispersion))
        Dispersion <- expandAsMatrix(fix.dispersion, dim = c(nTags, nlibs))
      else Dispersion <- expandAsMatrix(Dispersion, dim = c(nTags, nlibs))
      object$Dispersion <- Dispersion
      object

    }
    diff.fun <- function(object)
    {

      ## it is low-level function of NBsim ##
      ## it creates diff genes according to foldDiff ##
      group <- object$group
      pDiff <- object$pDiff
      pUp <-  object$pUp
      foldDiff <- object$foldDiff
      Lambda <- object$Lambda
      nTags <- object$nTags
      g <- group == levels(group)[1]
      ind <- sample(nTags, floor(pDiff*nTags))
      if(length(ind)>0 & !foldDiff == 1 ) {
        fcDir <- sample(c(-1,1), length(ind), prob=c(1-pUp,pUp), replace=TRUE)
        Lambda[ind,g] <- Lambda[ind,g]*exp(log(foldDiff)/2*fcDir)
        Lambda[ind,!g] <- Lambda[ind,!g]*exp(log(foldDiff)/2*(-fcDir))
        #Lambda <- t(t(Lambda)/colSums(Lambda))
        object$Lambda <- Lambda
        object$indDE <- ind
        object$indNonDE <- (1:nTags)[-ind]
        object$mask_DEup <- object$mask_DEdown <- object$mask_nonDE <- expandAsMatrix(FALSE, dim = dim(Lambda))
        object$mask_DEup[ind[fcDir == 1], g] <- TRUE
        object$mask_DEup[ind[fcDir == -1], !g] <- TRUE
        object$mask_DEdown[ind[fcDir == -1], g] <- TRUE
        object$mask_DEdown[ind[fcDir == 1], !g] <- TRUE
        object$mask_nonDE[-ind,] <- TRUE
        object$mask_DE <- object$mask_DEup | object$mask_DEdown}
      if(foldDiff == 1| pDiff == 0)
        object$indDE <- NA
      object
    }
    sim.fun <- function(object)
    {
      ## it is low-level function of NBsim ##
      ## it simulate counts using rnbinom ##
      Lambda <- object$Lambda
      Dispersion <- object$Dispersion
      nTags <- object$nTags
      nlibs <- object$nlibs
      lib.size <- object$lib.size
      counts <- matrix(rnbinom(nTags*nlibs, mu = t(t(Lambda)*lib.size), size = 1/Dispersion), nrow = nTags, ncol = nlibs)
      rownames(counts) <- paste("ids", 1:nTags, sep = "")
      object$counts <- counts
      object
    }

    outlier.fun <- function(object, outlierMech, pOutlier, min.factor = 2, max.factor = 5)
    {
      ## it is low-level function of NBsim ##
      ## it makes outlier ##
      outlierMech <- match.arg(outlierMech, c("S", "M", "R"))
      dim <- dim(object$counts)
      outlier.factor <- function() runif(1, min.factor, max.factor)
      countAddOut <- object$counts
      LambdaAddOut <- object$Lambda
      DispersionAddOut <- object$Dispersion
      switch(outlierMech,
             S = {
               mask_outlier <- expandAsMatrix(FALSE, dim = dim)
               id_r <- which(runif(dim[1]) < pOutlier)
               id_c <- sample(dim[2], length(id_r), replace = TRUE)
               for(i in seq(id_r))
                 mask_outlier[id_r[i], id_c[i]] <- TRUE
               countAddOut[mask_outlier] <- sapply(countAddOut[mask_outlier], function(z) round(z*outlier.factor()))
             },
             R = {
               mask_outlier <- matrix(runif(dim[1]*dim[2]) < pOutlier, dim[1], dim[2])
               countAddOut[mask_outlier] <- sapply(countAddOut[mask_outlier], function(z) round(z*outlier.factor()))
             },

             M = {
               mask_outlier <- matrix(runif(dim[1]*dim[2]) < pOutlier, dim[1], dim[2])
               LambdaAddOut[mask_outlier] <- sapply(LambdaAddOut[mask_outlier], function(z) z*outlier.factor())
               countAddOut[mask_outlier] <- rnbinom(sum(mask_outlier), mu = t(t(LambdaAddOut)*object$lib.size)[mask_outlier], size = 1/DispersionAddOut[mask_outlier])
             }
      )
      if(!object$foldDiff == 1 & !pDiff == 0)
      {
        indDEupOutlier <- which(apply(object$mask_DEup & mask_outlier, 1, any))
        indDEdownOutlier <- which(apply(object$mask_DEdown & mask_outlier, 1, any))
        indDEnoOutlier <- which(apply((object$mask_DE & !mask_outlier) , 1, all))
        indNonDEOutlier <- which(apply(object$mask_nonDE & mask_outlier, 1, any))
        indNonDEnoOutlier <- which(apply((object$mask_nonDE & !mask_outlier) , 1, all))
        indDEbothOutlier <- NA
        o <- indDEupOutlier %in% indDEdownOutlier
        q <-  indDEdownOutlier %in% indDEupOutlier
        if(any(o))
        {
          indDEupOutlier <- indDEupOutlier[!o]
          indDEbothOutlier <- indDEupOutlier[o]
          indDEdownOutlier <- indDEdownOutlier[!q]
        }
      }
      else
      {
        indDEupOutlier <- indDEdownOutlier <- indDEnoOutlier <- indNonDEOutlier <- indNonDEnoOutlier <- indDEbothOutlier <- NA
      }
      out <- list(countAddOut = countAddOut, outlierMech = outlierMech, pOutlier = pOutlier, mask_outlier = mask_outlier, indDEupOutlier = indDEupOutlier,
                  indDEdownOutlier = indDEdownOutlier, indDEbothOutlier = indDEbothOutlier, indDEnoOutlier = indDEnoOutlier, indNonDEOutlier = indNonDEOutlier,
                  indNonDEnoOutlier = indNonDEnoOutlier, LambdaAddOut = LambdaAddOut, DispersionAddOut = DispersionAddOut)

    }

    calProb <- function(x, l) round(1 -(1 - x)^(1/l), 2) ## calculate probability to make sure all the outlierMech produce the same amount of outliers ##


    if(verbose) message("Preparing dataset.\n")
    if(class(dataset) == "DGEList")
    {
      dat <- dataset
      dat[["R"]] <- dat[["S"]] <- dat[["M"]] <- dat[["pOutlier"]] <- dat[["outlierMech"]]<- NULL
    }
    else if(is.character(dataset))
    {
      load(dataset)
      dat <- get(gsub("(\\.)(\\w+)", "", basename(dataset)))
      dat[["R"]] <- dat[["S"]] <- dat[["M"]] <- dat[["pOutlier"]] <- dat[["outlierMech"]]<- NULL
    }
    else if(is.matrix(dataset))
    {
      if(is.null(name)) name <- deparse(substitute(dataset))
      dataset <- getDataset(counts =dataset, drop.extreme.dispersion = drop.extreme.dispersion, drop.low.lambda = drop.low.lambda)
      dat <- new("DGEList", list(dataset = dataset, nTags = nTags, lib.size = lib.size, nlibs = nlibs, group = group, design = model.matrix(~group), pDiff= pDiff, pUp = pUp, foldDiff = foldDiff, outlierMech = outlierMech, min.factor = min.factor, max.factor = max.factor, name = name))
    }
    else
      dat <- new("DGEList", list(dataset = dataset, nTags = nTags, lib.size = lib.size, nlibs = nlibs, group = group, design = model.matrix(~group), pDiff= pDiff, pUp = pUp, foldDiff = foldDiff, outlierMech = outlierMech, min.factor = min.factor, max.factor = max.factor, name = name))

    if(!only.add.outlier)
    {
      if(is.null(lib.size))
        dat$lib.size <- runif(nlibs, min = 0.7*median(dat$dataset$dataset.lib.size), max = 1.3*median(dat$dataset$dataset.lib.size))

      if(is.null(nTags))
        dat$nTags <- dat$dataset$dataset.nTags
      if(verbose) message("Sampling.\n")
      dat <- sample.fun(dat)
      if(verbose) message("Calculating differential expression.\n")
      dat <- diff.fun(dat)
      if(verbose) message("Simulating data.\n")
      dat <- sim.fun(dat)
    }
    if(add.outlier){
      outlierMech <- match.arg(outlierMech,  c("S", "R", "M"), several.ok = TRUE)
      if(length(pOutlier)== 1 & length(outlierMech) > 1 & any(outlierMech == "S"))
      {
        prob <- calProb(pOutlier, length(group))
        pOutlier <- rep(pOutlier, length = length(outlierMech))
        pOutlier[!outlierMech == "S"] <- prob
      }
      else if(!length(pOutlier) == length(outlierMech))
        stop("pOutlier is not equal to outlierMech")
      if(verbose) message("Adding outliers.\n")
      dat[outlierMech] <- mapply(function(x, y) outlier.fun(dat, outlierMech = x, pOutlier = y, min.factor = min.factor, max.factor = max.factor), x = outlierMech, y = pOutlier, SIMPLIFY = FALSE)
      dat$pOutlier <- pOutlier
    }

    if(save.file)
    {

      ## save file for shiny app ##
      if(verbose) message("Saving file.\n")
      if(is.null(file))
      { g <- paste0("g", sum(levels(group)[1] == group), "v", sum(levels(group)[2] == group))
      f <- paste0("f", foldDiff)
      if(add.outlier) o <- paste0("o", sprintf( "%02d",100*pOutlier[1L]))
      else o <- paste0("o", sprintf( "%02d", 0 ))
      file <- paste0(dat$name, "/", g, f, o, ".Rdata")
      dir.create(dat$name, showWarnings = FALSE)
      }
      filenm <- eval(gsub("(\\.)(\\w+)", "", basename(file)))
      assign(filenm, dat)
      save(list = filenm, file = file)
    }
    dat
  }
#' OutputCatBasedPwf
#'
#' Reformat genewised results from JunctionSeq, and use this one to get enriched GO terms
#'
#' @param Re.pwf.exon.sj
#' @param gene_model
#' @param gene_2_cat
#' @param Output_file
#'
#' @return
#' @export
#'
#' @examples
#'
#' RE.cp<-OutputCatBasedPwf(Re.pwf.exon.sj,gene_model=gene.model,gene_2_cat=gene.2.cat.cp.mouse,Output_file="Cp.xls")
#'
#' RE.cft<-OutputCatBasedPwf(Re.pwf.exon.sj,gene_model=gene.model,gene_2_cat=gene.2.cat.tft.mouse,Output_file="C_tft.xls")
#'
#'
#' RE.cp.gene<-OutputCatBasedPwf(Gene.non.coding.GO$pwf,gene_model=gene.model,gene_2_cat=gene.2.cat.cp.mouse,
#' Output_file="/media/H_driver/PJ/geneGL_rm_non_coding_Cp.xls")
#'
#' RE.cft.gene<-OutputCatBasedPwf(Gene.non.coding.GO$pwf,gene_model=gene.model,gene_2_cat=gene.2.cat.tft.mouse,
#' Output_file="/media/H_driver/PJ/geneGL_rm_non_coding_tft.xls")
#'
#'

OutputCatBasedPwf<-function(Re.pwf.exon.sj,gene_model,gene_2_cat,Output_file){

  Re<-goseq2(Re.pwf.exon.sj,"mm10","ensGene",gene.model=gene_model,gene2cat=gene_2_cat,use_genes_without_cat=TRUE)

  #select GO term(10<=numInCat<=300 and BP only)

  index.select<-which(Re$numInCat>=10&Re$numInCat<=300)

  Re.select<-Re[index.select,]

  over_represented_pvalue_adjusted<-p.adjust(Re.select$over_represented_pvalue,method="BH")

  Re.select.with.adjP<-cbind(Re.select[,c(1,2)],over_represented_pvalue_adjusted,Re.select[,-c(1,2,3)])


  dataset2<- Re.select.with.adjP

  Re2<-list(EnrichedGO=Re,SelectedEnrichedGO=dataset2)

  dataset2[sapply(dataset2, is.list)] <-
    sapply(dataset2[sapply(dataset2, is.list)],
           function(x)sapply(x, function(y) paste(unlist(y),collapse=", ") ) )

  write.table(dataset2,file=Output_file,row.names = FALSE,quote=FALSE,sep="\t")

  return(Re2)

}
#' Title
#'
#' @param jscs.object
#' @param gene.name
#' @param outfile_prefix
#'
#' @return
#' @export
#'
#' @examples
#' OutputFigure(jscs.object=Re.PJ,gene.name="ENSMUSG00000000028+ENSMUSG00000005262",
#' outfile_prefix="/media/H_driver/PJ/plot_ENSMUSG00000000028_ENSMUSG00000005262/")
#'
OutputFigure <- function(jscs.object,gene.name,outfile_prefix) {
  buildAllPlots(jscs=Re.PJ,outfile.prefix=outfile_prefix,gene.list=gene.name,
                use.plotting.device="png",plot.gene.level.expression=TRUE);
}
#' OutputGOBasedDEfromFeatures3
#'
#' Reformat genewised results from JunctionSeq, and use this one to get enriched GO terms
#'
#' @param re.PJ.gene.based
#' @param gene.model
#' @param Output_file
#'
#' @return
#' @export
#'
#' @examples
#'
#' RE.exon.sj.all.gene<-OutputGOBasedDEfromFeatures3(re.PJ.gene.based,gene.model,"GO_exon_sj_use_all_gene.xls")
#'
#'
#'
OutputGOBasedDEfromFeatures3<-function(re.PJ.gene.based,DE_type,gene.model,Output_file_dir){

  #re2<-ReformatDataUseTable(re.PJ.gene.based)

  re2<-re.PJ.gene.based

  Re.Go.adjusted.by.exon.SJ<-GotermAnalysisUseFeatureDefineDE(re2,ad="exon_SJ",sub_feature=NULL,DE_define=DE_type,gene_model=gene.model,Output_file_dir)

  #head(Re.Go.adjusted.by.exon.SJ[[1]])

  #select GO term(10<=numInCat<=300 and BP only)

  index.select<-which(Re.Go.adjusted.by.exon.SJ[[1]]$numInCat>=10&Re.Go.adjusted.by.exon.SJ[[1]]$numInCat<=300&Re.Go.adjusted.by.exon.SJ[[1]]$ontology=="BP")

  Re.Go.adjusted.by.exon.SJ.select<-Re.Go.adjusted.by.exon.SJ[[1]][index.select,]

  over_represented_pvalue_adjusted<-p.adjust(Re.Go.adjusted.by.exon.SJ.select$over_represented_pvalue,method="BH")

  Re.Go.adjusted.by.exon.SJ.select.with.adjP<-cbind(Re.Go.adjusted.by.exon.SJ.select[,c(1,2)],over_represented_pvalue_adjusted,Re.Go.adjusted.by.exon.SJ.select[,-c(1,2,3)])


  dataset2<- Re.Go.adjusted.by.exon.SJ.select.with.adjP

  dataset2[sapply(dataset2, is.list)] <-
    sapply(dataset2[sapply(dataset2, is.list)],
           function(x)sapply(x, function(y) paste(unlist(y),collapse=", ") ) )

  write.table(dataset2,file=paste0(Output_file_dir,"/",DE_type,".xls"),row.names = FALSE,quote=FALSE,sep="\t")

  write.table(Re.Go.adjusted.by.exon.SJ[[2]],file=paste0(Output_file_dir,"/",DE_type,"pwf.xls"),row.names = T,quote=FALSE,sep="\t")

  return(Re.Go.adjusted.by.exon.SJ)

}
#' OutputGOBasedDEfromFeatures
#'
#' Use genewised results from JunctionSeq to get enriched GO terms
#'
#' @param re.PJ.gene.based
#' @param gene.model
#' @param Output_file
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#' RE.exon.sj<-OutputGOBasedDEfromFeatures(re.PJ.gene.based,gene.model,"GO_exon_sj_3.xls")
#'
#' sink("test.txt")
#' RE.exon.sj<-OutputGOBasedDEfromFeatures(re.PJ.gene.based,gene.model,"GO_exon_sj_4.xls")
#' sink()
#'
OutputGOBasedDEfromFeatures<-function(re.PJ.gene.based,gene.model,Output_file){

  re<-pData(re.PJ.gene.based)

  no.re.testable.index<-which(as.character(re$mostSigID)=="character(0)")

  re2<-re[-no.re.testable.index,]

  print(names(re2))

  Re.Go.adjusted.by.exon.SJ<-GotermAnalysis2GetAllGOTerms_AdjustedByNumOfJunctionWithinOneGene(re2,ad="exon_SJ",sub_feature=NULL,0.05,gene_model=gene.model)

  head(Re.Go.adjusted.by.exon.SJ[[1]])

  #select GO term(10<=numInCat<=300 and BP only)

  index.select<-which(Re.Go.adjusted.by.exon.SJ[[1]]$numInCat>=10&Re.Go.adjusted.by.exon.SJ[[1]]$numInCat<=300&Re.Go.adjusted.by.exon.SJ[[1]]$ontology=="BP")

  Re.Go.adjusted.by.exon.SJ.select<-Re.Go.adjusted.by.exon.SJ[[1]][index.select,]

  over_represented_pvalue_adjusted<-p.adjust(Re.Go.adjusted.by.exon.SJ.select$over_represented_pvalue,method="BH")

  Re.Go.adjusted.by.exon.SJ.select.with.adjP<-cbind(Re.Go.adjusted.by.exon.SJ.select[,c(1,2)],over_represented_pvalue_adjusted,Re.Go.adjusted.by.exon.SJ.select[,-c(1,2,3)])


  dataset2<- Re.Go.adjusted.by.exon.SJ.select.with.adjP

  dataset2[sapply(dataset2, is.list)] <-
    sapply(dataset2[sapply(dataset2, is.list)],
           function(x)sapply(x, function(y) paste(unlist(y),collapse=", ") ) )

  write.table(dataset2,file=Output_file,row.names = FALSE,quote=FALSE,sep="\t")

  return(dataset2)

}
#' OutputGOBasedDEfromFeatures2
#'
#' Reformat genewised results from JunctionSeq, and use this one to get enriched GO terms
#'
#' @param re.PJ.gene.based
#' @param gene.model
#' @param Output_file
#'
#' @return
#' @export
#'
#' @examples
#'
#' RE.exon.sj.all.gene<-OutputGOBasedDEfromFeatures2(re.PJ.gene.based,gene.model,"GO_exon_sj_use_all_gene.xls")
#'
#'
#'
OutputGOBasedDEfromFeatures2<-function(re.PJ.gene.based,gene.model,Output_file){

  re2<-ReformatData(re.PJ.gene.based)

  Re.Go.adjusted.by.exon.SJ<-GotermAnalysisUseReformatedData(re2,ad="exon_SJ",sub_feature=NULL,0.05,"hg19","ensGene",gene_model=gene.model)

  Re<-Re.Go.adjusted.by.exon.SJ

  head(Re.Go.adjusted.by.exon.SJ[[1]])

  #select GO term(10<=numInCat<=300 and BP only)

  index.select<-which(Re.Go.adjusted.by.exon.SJ[[1]]$numInCat>=10&Re.Go.adjusted.by.exon.SJ[[1]]$numInCat<=300&Re.Go.adjusted.by.exon.SJ[[1]]$ontology=="BP")

  Re.Go.adjusted.by.exon.SJ.select<-Re.Go.adjusted.by.exon.SJ[[1]][index.select,]

  over_represented_pvalue_adjusted<-p.adjust(Re.Go.adjusted.by.exon.SJ.select$over_represented_pvalue,method="BH")

  Re.Go.adjusted.by.exon.SJ.select.with.adjP<-cbind(Re.Go.adjusted.by.exon.SJ.select[,c(1,2)],over_represented_pvalue_adjusted,Re.Go.adjusted.by.exon.SJ.select[,-c(1,2,3)])


  dataset2<- Re.Go.adjusted.by.exon.SJ.select.with.adjP

  dataset2[sapply(dataset2, is.list)] <-
    sapply(dataset2[sapply(dataset2, is.list)],
           function(x)sapply(x, function(y) paste(unlist(y),collapse=", ") ) )

  Re2<-list(EnrichedGO=Re[[1]],Pwf=Re[[2]],SelectedEnrichedGO=dataset2)

  write.table(dataset2,file=Output_file,row.names = FALSE,quote=FALSE,sep="\t")

  return(Re2)

}
#' Title
#'
#' @param re.PJ.gene.based
#' @param output_file
#'
#' @return
#' @export
#'
#' @examples
#'
#' data(gene.model)
#'
#' OutputGeneWiseTable(re.PJ.gene.based,gene.model,output_file="/media/H_driver/PJ/GeneWise_Re_annotated_using_new_annotation_3.csv")
#'
#' OutputGeneWiseTable(re.PJ.gene.based,gene.model,output_file="/media/H_driver/PJ/GeneWise_Re_annotated_using_new_annotation_4.csv")
#'
#' OutputGeneWiseTable(re.PJ.gene.based,gene.model,output_file="/media/H_driver/PJ/GeneWise_41540.csv")
#'
#' OutputGeneWiseTable(Re.PJ,output_file="/media/H_driver/PJ/AllFeature_Based.csv")
#' write.table(fData(Re.PJ),row.names = FALSE,file="/media/H_driver/PJ/AllFeature_Based.csv", quote=FALSE, sep="\t")
#'
#'
OutputGeneWiseTable <- function(re.PJ.gene.based,gene.model,output_file) {

   no.PJ.testable.index<-which(as.character(pData(re.PJ.gene.based)$mostSigID)=="character(0)")
   re.PJ.gene.based.testable<-pData(re.PJ.gene.based)[-no.PJ.testable.index,]
   dim(re.PJ.gene.based.testable)

   cor(as.numeric(re.PJ.gene.based.testable$numKnown),
   as.numeric(re.PJ.gene.based.testable$mostSigPadjust))

  table.gene.based.all.3<-re.PJ.gene.based.testable
  x <- vapply(table.gene.based.all.3$mostSigID, length, 1L) ## How many items per list element
  table.gene.based.all.3<- table.gene.based.all.3[rep(rownames(table.gene.based.all.3), x), ] ## Expand the data frame
  table.gene.based.all.3$mostSigID <- unlist(table.gene.based.all.3$mostSigID, use.names = FALSE)  ## Replace w

  gene.annotation.4.geneID<-do.call(rbind,lapply(table.gene.based.all.3[,1],GenerateGeneAnno,gene.model))

  table.gene.based.all.4<-merge(gene.annotation.4.geneID,table.gene.based.all.3,by="geneID")

  write.csv(table.gene.based.all.3,file = output_file)

}


# test.goseq2.2.gene.symbol<-sapply(test.goseq2[,8],function(u,gene.model){
#    xx<-gene.model[match(u,as.character(gene.model[,3])),1]
#    names(xx)<-
#
#   },gene.model)
#
# head(test.goseq2.2.gene.symbol)
# dataset2<-test.goseq33
#
# dataset2[sapply(dataset2, is.list)] <-
#   sapply(dataset2[sapply(dataset2, is.list)],
#          function(x)sapply(x, function(y) paste(unlist(y),collapse=", ") ) )
#
# head(dataset2)
# write.table(dataset2,row.names = FALSE,"myfile.xls", quote=FALSE, sep="\t")
#' Title
#'
#' @param dir.name
#' @param file.sample
#' @param file.count
#' @param file.gff
#'
#' @return
#' @export
#'
#' @examples
#'
#' dir.name="/Volumes/Bioinformatics$/2015/Nimer_Cheng/"
#' file.sample="decoder.bySample.txt"
#' file.count="_junction_seq_new_gtf_7/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt"
#' file.gff="Mus_musculus.GRCm38.83.JunctionSeq.flat.gff"
#' Re<-PermutattionOfRawData(dir.name,file.sample,file.count,file.gff)
#'
PermutattionOfRawData<-function(dir.name,file.sample,file.count,file.gff){
  #Get sample file
  path.file.sample<-paste0(dir.name,file.sample)
  decoder.bySample<-read.table(path.file.sample,header=T,stringsAsFactors = F)

  print(decoder.bySample)

  #Get count file
  path.file.count<-paste0(dir.name,decoder.bySample$sample.ID,file.count)
  countFiles<-paste0(path.file.count)

  print(countFiles)

  #Get annotation file
  path.file.gff<-paste0(dir.name,file.gff)

  print(path.file.gff)

  #Analysis
#  jscs.2<-runJunctionSeqAnalyses(sample.files= countFiles,sample.names= decoder.bySample$sample.ID,condition= decoder.bySample$group.ID,
#                                 flat.gff.file=path.file.gff,nCores=1,verbose=TRUE,debug.mode=TRUE,use.multigene.aggregates = TRUE)

#  return(jscs.2)
}
#' Plot the probability Weighting Function
#'
#' @param pwf: The first step of the pathway analysis is to assess if there is any systematic bias of the significant genes. The plotPWF function takes input a vector of genes (with label 1 as significant and 0 otherwise), and the bias factor to be adjusted (e.g. number of splice junctions in the genes). The output of the function is a figure that illustrates the relationship between bias factor and proportion of significant genes. The size of the gene set size (i.e. number of genes in each gene bin) on the x-axis can be adjusted by user. 
#' @param binsize
#' @param pwf_col
#' @param pwf_lwd
#' @param xlab
#' @param ylab
#' @param ...
#'
#' @return
#' @export
#'
#' @examples

plotPWF2<-function (pwf, binsize = "auto", pwf_col = 3, pwf_lwd = 2, xlab = "Biased Data in <binsize> gene bins.",
          ylab = "Proportion DE", ...)
{
  w = !is.na(pwf$bias.data)
  print(w)
  o = order(pwf$bias.data[w])
  print(o)

  rang = max(pwf$pwf, na.rm = TRUE) - min(pwf$pwf, na.rm = TRUE)
  if (rang == 0 & binsize == "auto")
    binsize = 1000
  if (binsize == "auto") {
    binsize = max(1, min(100, floor(sum(w) * 0.08)))
    resid = rang
    oldwarn = options()$warn
    options(warn = -1)
    while (binsize <= floor(sum(w) * 0.1) & resid/rang >
           0.001) {
      binsize = binsize + 100
      splitter = ceiling(1:length(pwf$DEgenes[w][o])/binsize)
      de = sapply(split(pwf$DEgenes[w][o], splitter), mean)
      binlen = sapply(split(as.numeric(pwf$bias.data[w][o]),
                            splitter), mean)
      resid = sum((de - approx(pwf$bias.data[w][o], pwf$pwf[w][o],
                               binlen)$y)^2)/length(binlen)
    }
    options(warn = oldwarn)
  }
  else {
    splitter = ceiling(1:length(pwf$DEgenes[w][o])/binsize)
    print(splitter)
    de = sapply(split(pwf$DEgenes[w][o], splitter), mean)
    print(de)
    binlen = sapply(split(as.numeric(pwf$bias.data[w][o]),
                          splitter), median)
    print(binlen)
  }
  xlab = gsub("<binsize>", as.character(binsize), xlab)
  if ("xlab" %in% names(list(...))) {
    if ("ylab" %in% names(list(...))) {
      plot(binlen, de, ...)
    }
    else {
      plot(binlen, de, ylab = ylab, ...)
    }
  }
  else if ("ylab" %in% names(list(...))) {
    plot(binlen, de, xlab = xlab, ...)
  }
  else {
    plot(binlen, de, xlab = xlab, ylab = ylab, ...)
  }
  lines(pwf$bias.data[w][o], pwf$pwf[w][o], col = pwf_col,
        lwd = pwf_lwd)

 return(de)

}
#' ProcessCuffLinkResults
#'
#' @return
#' @export
#'
#' @examples
ProcessCuffLinkResults <- function() {
  getwd()

  cuff_data<-readCufflinks()

  diffGeneIDs <- getSig(cuff_data,level="genes",alpha=0.05)
  diffGenes<-getGenes(cuff_data,diffGeneIDs)

  names<-featureNames(diffGenes)
  row.names(names)=names$tracking_id
  diffGenesNames<-as.matrix(names)
  diffGenesNames<-diffGenesNames[,-1]

  names<-featureNames(diffGenes)
  row.names(names)=names$tracking_id
  diffGenesNames<-as.matrix(names)
  diffGenesNames<-diffGenesNames[,-1]


  diffGenesData<-diffData(diffGenes)
  row.names(diffGenesData)=diffGenesData$gene_id
  diffGenesData<-diffGenesData[,-1]

  # merge the two matrices by row names
  diffGenesOutput<-merge(diffGenesNames,diffGenesData,by="row.names")

  head(diffGenesOutput)

  featureNames(diffGenes)

  csDensity(genes(cuff_data))
  csScatter(genes(cuff_data),'Mut','HC')
  csVolcano(genes(cuff_data),'Mut','HC')
  mygene<-getGene(cuff_data, 'XLOC_011342')
  expressionBarplot(mygene)

  gene_diff_data<-diffData(genes(cuff_data))
  sig_gene_data<-subset(gene_diff_data,(significant=='yes'))
  nrow(sig_gene_data)

  isoform_diff_data<-diffData(isoforms(cuff_data), 'Mut', 'HC')
  sig_isoform_data<-subset(isoform_diff_data, (significant=='yes'))
  nrow(sig_isoform_data)

  tss_diff_data<-diffData(TSS(cuff_data), 'Mut', 'HC')
  sig_tss_data<-subset(tss_diff_data, (significant=='yes'))
  nrow(sig_tss_data)

  cds_diff_data<-diffData(CDS(cuff_data), 'Mut', 'HC')
  sig_cds_data<-subset(cds_diff_data, (significant=='yes'))
  nrow(sig_cds_data)

  promoter_diff_data<-distValues(promoters(cuff_data))
  sig_promoter_data<-subset(promoter_diff_data, (significant=='yes'))
  nrow(sig_promoter_data)

  splicing_diff_data<-distValues(splicing(cuff_data))
  sig_splicing_data<-subset(splicing_diff_data, (significant=='yes'))
  nrow(sig_splicing_data)

  relCDS_diff_data<-distValues(relCDS(cuff_data))
  sig_relCDS_data<-subset(relCDS_diff_data, (significant=='yes'))
  nrow(sig_relCDS_data)

  gene_diff_data<-diffData(genes(cuff_data))
  sig_gene_data<-subset(gene_diff_data, (significant=='yes'))
  write.table(sig_gene_data, "diff_genes.txt", sep = "\t", row.names = F,col.names = T, quote = F)

  dim(gene_diff_data)

  data.isoform<-read.table("isoform_exp.diff",header=T)
}
#' @ ProcessOutputFilesFrom_rMATS
#'
#' @ read output from rMAT
#'
#' @param input_file
#'
#' @return
#' @export
#'
#' @examples
#'
#' dir.name="/media/H_driver/PJ/rMAT/pHamard_MATS-selected-re-run-exon-centric/pHamard_MATS-selected-re-run-exon-centric/Neg_WTvKO_c0.01/"
#'
#' input.file.pattern="*Neg_WTvKO_c0.01.txt"
#'
#' sink(paste0(output.dir.name,"rMAT_output_from_PJ.txt"))
#' re.rMAT<-ProcessOutputFilesFrom_rMATS(dir.name,input.file.pattern)
#' sink()
#'
#'
ProcessOutputFilesFrom_rMATS<-function(dir.name,input.file.pattern,De_defined_by_what){

  ProcessOutputFilesFrom_rMATS_read<-function(input_file){

    re=read.table(input_file,header=T)

    return(re)

  }

  file.name=paste0(dir.name,dir(dir.name,recursive = TRUE,pattern=input.file.pattern))
  file.name.2<-as.list(file.name)

  names(file.name.2)=sapply(strsplit(file.name,split="\\/"),"[[",9)

  re.out<-lapply(file.name.2,ProcessOutputFilesFrom_rMATS_read)

  #Define DE by FDR
  if(De_defined_by_what=="FDR"){
  re.out.2<-do.call(c,lapply(re.out, function(u){
    x<-as.character(u[which(u$FDR<0.05),]$GeneID)
    x
    #x<-as.data.frame(t(x))
    #colnames(x)<-colnames(Data4Goterm)
    #x
  }))
  }else if(De_defined_by_what=="P_and_inclusion"){
    re.out.2<-do.call(c,lapply(re.out, function(u){
      x<-as.character(u[which(u$PValue<0.05&abs(u$IncLevelDifference)>0.20),]$GeneID)
      x
      #x<-as.data.frame(t(x))
      #colnames(x)<-colnames(Data4Goterm)
      #x
    }))
    }

# print(names((re.out.2)))

  n1=length(re.out.2[grep("A3SS.MATS.ReadsOnTargetAndJunctionCounts",names(re.out.2))])
  n2=length(re.out.2[grep("A5SS.MATS.ReadsOnTargetAndJunctionCounts",names(re.out.2))])
  n3=length(re.out.2[grep("MXE.MATS.ReadsOnTargetAndJunctionCounts",names(re.out.2))])
  n4=length(re.out.2[grep("RI.MATS.ReadsOnTargetAndJunctionCounts",names(re.out.2))])
  n5=length(re.out.2[grep("SE.MATS.ReadsOnTargetAndJunctionCounts",names(re.out.2))])

  cat(n1,"\t",n2,"\t",n3,"\t",n4,"\t",n5,"\n")

  re.out.3<-list(JunctionCountOnly=unique(re.out.2[grep("JunctionCountOnly",names(re.out.2))]),
                 ReadsOnTargetAndJunctionCounts=unique(re.out.2[grep("ReadsOnTargetAndJunctionCounts",names(re.out.2))]),
                 SEMATSJunctionCountOnly=unique(re.out.2[grep("SE.MATS.JunctionCountOnly",names(re.out.2))]),
                 SEReadsOnTargetAndJunctionCounts=unique(re.out.2[grep("SE.MATS.ReadsOnTargetAndJunctionCounts",names(re.out.2))]),
                 SEs=unique(re.out.2[grep("SE.MATS",names(re.out.2))])
                 )

  return(re.out.3)

}
#' Title
#'
#' @param re.PJ.gene.based
#'
#' @return
#' @export
#'
#' @examples
#' re.PJ.gene.based.testable.reformat<-ReformatData(re.PJ.gene.based)
#'
ReformatData <- function(re.PJ.gene.based) {

  re<-pData(re.PJ.gene.based)

  no.re.testable.index<-which(as.character(re$mostSigID)=="character(0)")
  re2<-re[-no.re.testable.index,]

  All.gene.id.based.on.sub_feature<-unique(unlist(strsplit(re2[,1],"\\+")))
  All.gene.id.index<-rep(0,length(All.gene.id.based.on.sub_feature))
  names(All.gene.id.index)=All.gene.id.based.on.sub_feature


  reformat.gene.p<-do.call(rbind,sapply(All.gene.id.based.on.sub_feature, function(u,re2){
    x<-re2[grep(u,re2[,1]),-1]
    x<-as.data.frame(t(x))
    #colnames(x)<-colnames(Data4Goterm)
    #x
  },re2))

  re3<-as.data.frame(reformat.gene.p)
  re3<-cbind(All.gene.id.based.on.sub_feature,re3)
  colnames(re3)[1]="geneID"

  return(re3)

}
#' ReformatDataUseTable
#'
#' @param re.PJ.gene.based
#'
#' @return
#' @export
#'
#' @examples
#' re.PJ.gene.based.testable.reformat.2<-ReformatDataUseTable(Re.PJ.selected.feature.FC.p)
#'
#' dim(re.PJ.gene.based.testable.reformat.2)
#'
#'head(re.PJ.gene.based.testable.reformat.2)
#'
ReformatDataUseTable <- function(re.PJ.gene.based) {

  re<-re.PJ.gene.based

  no.re.testable.index<-which(as.character(re$mostSigID)=="character(0)")
  re2<-re[-no.re.testable.index,]

  All.gene.id.based.on.sub_feature<-unique(unlist(strsplit(re2[,1],"\\+")))
  All.gene.id.index<-rep(0,length(All.gene.id.based.on.sub_feature))
  names(All.gene.id.index)=All.gene.id.based.on.sub_feature


  reformat.gene.p<-do.call(rbind,sapply(All.gene.id.based.on.sub_feature, function(u,re2){
    x<-re2[grep(u,re2[,1]),-1]
    x<-as.data.frame(t(x))
    #colnames(x)<-colnames(Data4Goterm)
    #x
  },re2))

  re3<-as.data.frame(reformat.gene.p)
  re3<-cbind(All.gene.id.based.on.sub_feature,re3)
  colnames(re3)[1]="geneID"

  return(re3)

}
#  GetResultsFromJunctionSeq
#' Get analysis results using JunctionSeq
#'
#' @param dir.name
#' @param file.sample
#' @param file.count
#' @param file.gff
#'
#' @return
#' @export
#'
#' @examples
#'
#' # For Guoyan Project
#' load("/Volumes/Bioinformatics\$/2015/Nimer_Cheng/1_29_2016.RData")
#' dir.name="/Volumes/Bioinformatics$/2015/Nimer_Cheng/"
#' file.sample="decoder.bySample.txt"
#' file.count="_junction_seq_new_gtf_7/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt"
#' file.gff="Mus_musculus.GRCm38.83.JunctionSeq.flat.gff"
#' Re<-GetResultsFromJunctionSeq(dir.name,file.sample,file.count,file.gff)
#' head(fData(Re))
#' save(Re,file="Re_Run_test_GOSJ.RData")
#'

#' # For PJ project
#' dir.name.PJ="/media/H_driver/PJ/"
#' file.sample.PJ="decoder.bySample.rtf"
#' file.count.PJ="/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt"
#'
#' file.gff=paste0("/media/H_driver/2015/Nimer_Cheng/",file.gff)
#' Re.PJ<-GetResultsFromJunctionSeq(dir.name.PJ,file.sample.PJ,file.count.PJ,file.gff)
#'
#' save(Re.PJ,file="/media/H_driver/PJ/PJ_jscs.RData")
#' buildAllPlots(jscs=jscs,outfile.prefix="./plots_based_on_DE_splice_site_gene1/",
#' gene.list=gene.based.de.splice.site[1],use.plotting.device="png",plot.gene.level.expression=TRUE,sequencing.type="single-end");
#'
#'
#
#'
GetResultsFromJunctionSeq<-function(dir.name,file.sample,file.count,file.gff){
#Get sample file
path.file.sample<-paste0(dir.name,file.sample)
decoder.bySample<-read.table(path.file.sample,header=T,stringsAsFactors = F)

print(decoder.bySample)

#Get count file
path.file.count<-paste0(dir.name,decoder.bySample$sample.ID,file.count)
countFiles<-paste0(path.file.count)

print(countFiles)

#Get annotation file
#path.file.gff<-paste0(dir.name,file.gff)
path.file.gff<-file.gff

print(path.file.gff)

#Analysis
jscs.2<-runJunctionSeqAnalyses(sample.files= countFiles,sample.names= decoder.bySample$sample.ID,condition= decoder.bySample$group.ID,
                             flat.gff.file=path.file.gff,nCores=1,verbose=TRUE,debug.mode=TRUE,use.multigene.aggregates = TRUE)

return(jscs.2)
}
#' SelectFeature
#'
#' @param Re.PJ
#'
#' @return
#' @export
#'
#' @examples
#' Re.PJ.selected<-SelectFeature(Re.PJ)
#'
#'dim(Re.PJ.selected)

SelectFeature<-function(Re.PJ){

  table.based.RE.PJ<-fData(Re.PJ)

  re2<-table.based.RE.PJ[table.based.RE.PJ$testable,]

  re2<-table.based.RE.PJ[which(table.based.RE.PJ[,20]>=1),]

  return(re2)

}
#' Select_DE_gene_based_on_Feature
#'
#' @param Re.PJ
#' @param re.PJ.gene.based
#' @param re.rMAT
#' @param splicing_type: there are 5 splicing types from rMAT, need to choose one
#'
#' @return
#' @export
#'
#' @examples
#'
#' Re.PJ.selected.feature.FC.p<-Select_DE_gene_basd_on_Feature(Re.PJ,re.PJ.gene.based,re.rMAT,0.58,0.05,outputfile_DGE_FC_P_geneWise)
#' dim(Re.PJ.selected.feature.FC.p)
#' head(Re.PJ.selected.feature.FC.p)
#' length(which(Re.PJ.selected.feature.FC.p$DE_or_not==1))
#'
#' Re.PJ.selected.feature.FC.p.check<-Select_DE_gene_basd_on_Feature(Re.PJ,re.PJ.gene.based,re.rMAT,"SE",2,0.05,outputfile_DGE_FC_P_geneWise)

Select_DE_gene_basd_on_Feature<-function(Re.PJ,re.PJ.gene.based,re.rMAT,splicing_type,cutoff_FC,cutoff_P_value,gene_model,venn_output){

  feature.based.RE.PJ<-fData(Re.PJ)

  gene.based.RE.PJ<-ReformatDataUseTable(pData(re.PJ.gene.based))

  cat(dim(feature.based.RE.PJ),"\t",dim(gene.based.RE.PJ),"\n")

  re2<-feature.based.RE.PJ[which(feature.based.RE.PJ[,20]>cutoff_FC&feature.based.RE.PJ[,11]<cutoff_P_value),]

  DE.gene.based.on.FC.p.of.feature<-unique(unlist(strsplit(as.character(re2$geneID),"\\+")))

  re2.geneWise.p.based<-feature.based.RE.PJ[which(feature.based.RE.PJ$geneWisePadj<cutoff_P_value),]

  DE.gene.based.on.geneWise.p.only<-unique(unlist(strsplit(as.character(re2.geneWise.p.based$geneID),"\\+")))

  #DE.gene.based.on.FC.p.of.feature.and.geneWise.p<-intersect(DE.gene.based.on.FC.p.of.feature,DE.gene.based.on.geneWise.p.only)



  DE.gene.feature<-DE.gene.based.on.FC.p.of.feature
  DE.gene.geneWisePadj<-DE.gene.based.on.geneWise.p.only

  if(splicing_type=="All_5_Types"){
    cat(length(DE.gene.based.on.FC.p.of.feature),"\t",length(DE.gene.based.on.geneWise.p.only),"\t",
        length(re.rMAT$ReadsOnTargetAndJunctionCounts),"\n")
        DE.gene.rMAT<-re.rMAT$ReadsOnTargetAndJunctionCounts
  }else if(splicing_type=="SE"){
    cat(length(DE.gene.based.on.FC.p.of.feature),"\t",length(DE.gene.based.on.geneWise.p.only),"\t",
        length(re.rMAT$SEReadsOnTargetAndJunctionCounts),"\n")
      DE.gene.rMAT<-re.rMAT$SEReadsOnTargetAndJunctionCounts
    }

  DE_or_not_feature<-rep(0,dim(gene.based.RE.PJ)[1])
  DE_or_not_geneWise<-rep(0,dim(gene.based.RE.PJ)[1])
  DE_or_not_rMAT_based<-rep(0,dim(gene.based.RE.PJ)[1])

  re3<-cbind(gene.based.RE.PJ,DE_or_not_feature,DE_or_not_geneWise,DE_or_not_rMAT_based)


  re3[which(re3$geneID %in% DE.gene.feature),]$DE_or_not_feature<-1
  re3[which(re3$geneID %in% DE.gene.geneWisePadj),]$DE_or_not_geneWise<-1
  re3[which(re3$geneID %in% DE.gene.rMAT),]$DE_or_not_rMAT_based<-1

  #re4<-unique(re3[which(re3$DE_or_not==1),]$geneID)

  Re4<-list(GeneAll=re3,DGE.geneWise=DE.gene.geneWisePadj,
            DGE.featureWise=DE.gene.feature,
            DGE.rMAT=DE.gene.rMAT)

  # venn.plot <- venn.diagram(
  #   x = Re4[c(1,2)],
  #   filename = venn_output,
  #   height = 3000,
  #   width = 3500,
  #   resolution = 1000,
  #   col = "black",
  #   lty = "dotted",
  #   lwd = 1,
  #   fill = c("red","blue"),
  #   alpha = 0.50,
  #   label.col = c(rep("white",3)),
  #   cex = 0.5,
  #   fontfamily = "serif",
  #   fontface = "bold",
  #   cat.col = c("red","blue"),
  #   cat.cex = 0.5,
  #   cat.pos = 0.5,
  #   cat.dist = 0.05,
  #   cat.fontfamily = "serif"
  # )

  venn.plot <- venn.diagram(
    x = Re4[c(2,3,4)],
    filename = paste0(venn_output,"DE_gene_overlap_rMAT_SE.tiff"),
    col = "black",
    lty = "dotted",
    lwd = 2,
    fill = c("red", "orange", "blue"),
    alpha = 0.50,
    label.col = c(rep("white",7)),
    cex = 1,
    fontfamily = "serif",
    fontface = "bold",
    cat.col = c("red", "orange", "blue"),
    cat.cex = 0.8,
    cat.fontfamily = "serif"
  )

  write.table(as.data.frame(intersect(Re4$DGE.featureWise,Re4$DGE.rMAT)),file=paste0(venn_output,"DGE_overlap_rMAT_SE_with_feature.txt"),row.names = FALSE,quote=FALSE,sep="\t")

  FindElementsBetweenSets<-function(xx,gene_model,venn_output){

    a1<-xx[[1]]
    a2<-xx[[2]]
    a3<-xx[[3]]

    a12<-intersect(a1,a2)
    a13<-intersect(a1,a3)
    a23<-intersect(a2,a3)

    a123<-intersect(a12,a3)
    a12not3<-setdiff(a12,a123)
    a13not2<-setdiff(a13,a123)
    a23not1<-setdiff(a23,a123)

    a1not23<-setdiff(a1,union(a12,a13))
    a2not13<-setdiff(a2,union(a12,a23))
    a3not12<-setdiff(a3,union(a13,a23))

    DGE.geneWise<-c(a1not23,a12not3,a13not2,a123)
    DGE.featureWise<-c(a12not3,a2not13,a123,a23not1)
    DGE.rMAT<- c(a13not2,a123,a23not1,a3not12)

    re<-qpcR:::cbind.na(a123,a12not3,a13not2,a23not1,a1not23,a2not13,a3not12,DGE.geneWise,DGE.featureWise,DGE.rMAT)
    re[is.na(re)] <- ""
    re<-data.frame(apply(re,2,sort,decreasing=T))

    colnames(re)<-c(paste0(length(a123),"_genes"),paste0(length(a12not3),"_genes"),paste0(length(a13not2),"_genes"),
                    paste0(length(a23not1),"_genes"),paste0(length(a1not23),"_genes"),paste0(length(a2not13),"_genes"),
                    paste0(length(a3not12),"_genes"),paste0(length(DGE.geneWise),"_DGE.geneWise"),
                    paste0(length(DGE.featureWise),"_DGE.featureWise"),paste0(length(DGE.rMAT),"_DGE.rMAT"))

     a123.gene.symbol<-sapply(a123,function(u,gene.model){
      y<-gene.model[which(as.character(gene.model[,3]) %in% u),1]
      if(length(as.character(y))==0){
        y<-"NotMatched"
      }
      y
    },gene.model)

     a12not3.gene.symbol<-sapply(a12not3,function(u,gene.model){
       y<-gene.model[which(as.character(gene.model[,3]) %in% u),1]
       if(length(as.character(y))==0){
         y<-"NotMatched"
       }
       y
     },gene.model)

     a13not2.gene.symbol<-sapply(a13not2,function(u,gene.model){
       y<-gene.model[which(as.character(gene.model[,3]) %in% u),1]
       if(length(as.character(y))==0){
         y<-"NotMatched"
       }
       y
     },gene.model)

     a23not1.gene.symbol<-sapply(a23not1,function(u,gene.model){
       y<-gene.model[which(as.character(gene.model[,3]) %in% u),1]
       if(length(as.character(y))==0){
         y<-"NotMatched"
       }
       y
     },gene.model)

     a1not23.gene.symbol<-sapply(a1not23,function(u,gene.model){
       y<-gene.model[which(as.character(gene.model[,3]) %in% u),1]
       if(length(as.character(y))==0){
         y<-"NotMatched"
       }
       y
     },gene.model)

     a2not13.gene.symbol<-sapply(a2not13,function(u,gene.model){
       y<-gene.model[which(as.character(gene.model[,3]) %in% u),1]
       if(length(as.character(y))==0){
         y<-"NotMatched"
       }
       y
     },gene.model)

     a3not12.gene.symbol<-sapply(a3not12,function(u,gene.model){
       y<-gene.model[which(as.character(gene.model[,3]) %in% u),1]
       if(length(as.character(y))==0){
         y<-"NotMatched"
       }
       y
     },gene.model)


     DGE.geneWise.gene.symbol<-sapply(DGE.geneWise,function(u,gene.model){
       y<-gene.model[which(as.character(gene.model[,3]) %in% u),1]
       if(length(as.character(y))==0){
         y<-"NotMatched"
       }
       y
     },gene.model)

     DGE.featureWise.gene.symbol<-sapply(DGE.featureWise,function(u,gene.model){
       y<-gene.model[which(as.character(gene.model[,3]) %in% u),1]
       if(length(as.character(y))==0){
         y<-"NotMatched"
       }
       y
     },gene.model)

     DGE.rMAT.gene.symbol<-sapply(DGE.rMAT,function(u,gene.model){
       y<-gene.model[which(as.character(gene.model[,3]) %in% u),1]
       if(length(as.character(y))==0){
         y<-"NotMatched"
       }
       y
     },gene.model)


     re2<-qpcR:::cbind.na(a123.gene.symbol,a12not3.gene.symbol,a13not2.gene.symbol,
                a23not1.gene.symbol,a1not23.gene.symbol,a2not13.gene.symbol,a3not12.gene.symbol,DGE.geneWise.gene.symbol,
                DGE.featureWise.gene.symbol,DGE.rMAT.gene.symbol)
     re2[is.na(re2)] <- ""
     re2<-data.frame(apply(re2,2,sort,decreasing=T))

     colnames(re2)<-c(paste0(length(a123.gene.symbol),"_genes"),paste0(length(a12not3.gene.symbol),"_genes"),paste0(length(a13not2.gene.symbol),"_genes"),
                     paste0(length(a23not1.gene.symbol),"_genes"),paste0(length(a1not23.gene.symbol),"_genes"),paste0(length(a2not13.gene.symbol),"_genes"),
                     paste0(length(a3not12.gene.symbol),"_genes"),
                     paste0(length(DGE.geneWise),"_DGE.geneWise"),
                     paste0(length(DGE.featureWise),"_DGE.featureWise"),paste0(length(DGE.rMAT),"_DGE.rMAT"))


    write.table(re,file=paste0(venn_output,"DGE_list_Ensemble.xls"),row.names = FALSE,quote=FALSE,sep="\t")

    write.table(re2,file=paste0(venn_output,"DGE_list_gene_symbol.xls"),row.names = FALSE,quote=FALSE,sep="\t")

   re3<-list(re=re,re2=re2)

   return(re3)

  }

  Ree4<-FindElementsBetweenSets(Re4[c(2,3,4)],gene_model,venn_output)


  return(Re4)

}
#' @title Simulate gene sets based on the number of splicing junctions
#
#' @description
#'We sample the mean of the number of splicing junction, and use this mean to simulate an array A including
#'splicing junction numbers with certain amounts. we convert this number into an array AAA between 0 and 1. Use the median of
#'AAA array as the probability for 1, we sample from [0,1] to get a list with certain length. This list is used as
#'a gene set. For example, if this gene set has 30 genes, each gene has its splicing junction number, we will
#'have 30 splicing junction number. After differential gene expression analysis, some genes are differentiallly expressed,
#'and are labeled as 1,otherwise labeled as 0.
#'Here we simulate two scenarios:
#'Scenario1: the probability of being 1 is dependent on the median of AAA
#'Scenario2: the probability of being 1 is not dependent on the median of AAA
#'
#' @param min_num_splicing minium splicing junction number
#' @param max_num_splicing max splicing junction number
#' @param num_gene number of gene
#'
#' @return a list with two set of genes
#' @export
#' @examples
#' re.random.DE2SJ<-SimulationSJ2DE(20,200,30)

SimulationSJ2DE<- function(min_num_splicing,max_num_splicing,num_gene) {
  crr.random<-array()
  for(j in 1:1:100){
    ans<-list()
    for(i in 1:50){
      mean_p<-sample(seq(min_num_splicing,max_num_splicing,10),1)
      temp<-data.frame()
      a=mean_p
      A=rpois(num_gene,a)
      AA=max(A)
      AAA=A/AA
      #AAA

      #p=median(AAA)
      p=runif(1)
      q=1-p
      SJNum<-A
      De.list<-sample(c(0,1),num_gene,replace = TRUE,prob=c(q,p))
      temp.data<-cbind(SJNum,De.list)

      De.prop<-length(which(De.list==1))/num_gene
      SJ.De.prop<-cbind(median(A),De.prop)

      ans[[i]]<-list(SJ.De.prop=SJ.De.prop,SJ.De.prop.data=temp.data)
    }

    ans.2<-do.call(rbind,ans)
    SJ.DE.prop<-matrix(unlist(ans.2[,1]),50,2,byrow=T)
    crr.random[j]<-cor(SJ.DE.prop)[1,2]

  }

  crr.force<-array()

  for(j in 1:100){
    ans<-list()
    for(i in 1:50){
      mean_p<-sample(seq(min_num_splicing,max_num_splicing,10),1)
      temp<-data.frame()
      a=mean_p
      A=rpois(num_gene,a)
      AA=max(A)
      AAA=A/AA
      #AAA

      p=median(AAA)
      #p=runif(1)
      q=1-p
      SJNum<-A
      De.list<-sample(c(0,1),num_gene,replace = TRUE,prob=c(q,p))
      temp.data<-cbind(SJNum,De.list)

      De.prop<-length(which(De.list==1))/30
      SJ.De.prop<-cbind(median(A),De.prop)

      ans[[i]]<-list(SJ.De.prop=SJ.De.prop,SJ.De.prop.data=temp.data)
    }

    ans.2<-do.call(rbind,ans)
    SJ.DE.prop<-matrix(unlist(ans.2[,1]),50,2,byrow=T)

    crr.force[j]<-cor(SJ.DE.prop)[1,2]

  }

  re<-list(Unrelated_SJ_DE=crr.random,related_SJ_DE=crr.force)
  re

}
#' Title
#'
#' @param jscs_object
#'
#' @return
#' @export
#'
#' @examples
#' re.PJ.gene.based.testable.reformat<-ReformatData(re.PJ.gene.based)
#' UseLogistic2CKBias(re.PJ.gene.based.testable.reformat)
#'
UseLogistic2CKBias <- function(jscs_object) {
  #mydata <- read.csv(input_file)

  mydata<-jscs_object
  #mydata<-mydata[-which(as.numeric(mydata$numKnown)==548),]

  ## view the first few rows of the data
  head(mydata)

  print(colnames(mydata))

  n.gene<-dim(mydata)[1]

  DE.out<-rep(0,n.gene)

  de.index<-which(mydata$geneWisePadj<0.05)

  DE.out[de.index]<-1

  mydata.2<-cbind(mydata,DE.out)

  print(head(mydata.2))

  hist(as.numeric(mydata.2$geneWisePadj))

  print(colnames(mydata.2))

  mylogit.2 <- glm(DE.out ~ as.numeric(numKnown), data = mydata.2, family = "binomial")

  logi.hist.plot(as.numeric(mydata.2$numKnown),mydata.2$DE.out,boxp=TRUE,type="hist",col="gray",xlabel="Number of splicing junctions"
                 ,counts=T)

  print(summary(mylogit.2))

}
#' Title
#'
#' @param input_file
#' @return
#' @export
#'
#' @examples
#'
#' UseLogistic2_adDE("/Volumes/Bioinformatics$/2015/Nimer_Cheng/
#' GeneWise_jscs3_all_with_anno_2_24_2016.csv")
#'
#'
UseLogistic2_adDE <- function(input_file) {
  mydata <- read.csv(input_file)
  ## view the first few rows of the data
  #head(mydata)

  #print(colnames(mydata))

  n.gene<-dim(mydata)[1]

  DE.out<-rep(0,n.gene)

  de.index<-which(mydata[,11]<0.05)

  DE.out[de.index]<-1

  mydata.2<-cbind(mydata,DE.out)

  #print(head(mydata.2))

  #hist(mydata.2$geneWisePadj)

  print(colnames(mydata.2))

  mylogit.sj.exon <- glm(DE.out ~ numKnown+numExons, data = mydata.2, family = "binomial")

  mylogit.sj <- glm(DE.out ~ numKnown, data = mydata.2, family = "binomial")

  mylogit.exon <- glm(DE.out ~ numExons, data = mydata.2, family = "binomial")

  #print(summary(mylogit.sj.exon))

  print(summary(mylogit.sj))

  print(summary(mylogit.exon))

}
#' Title
#'
#' @param input_file
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#'
#'
#' Re<-UseOtherDataSet("/media/H_driver/DataSet_SJ/GSE66793_cmp1-3-geo-juncs.tsv")
#'
#'
UseOtherDataSet<-function(input_file){
  #load(input_file)

  Re.Jun<-read.table(input_file,sep="\t",header = TRUE)



  # Re.DE<-read.csv(gene_based_input_file)
  #
  print(colnames(Re.Jun))

  cat(dim(Re.Jun)[1],"\n")

  cat(length(unique(Re.Jun[,13])),"\n")


  # print(head(Re.DE))
  #
  #
  # Re.Jun.2<-Re.Jun[,-1]
  #
  # colnames(Re.Jun.2)[1]="ensembl_gene_id"
  #
  # colnames(Re.DE)[1]="geneID"
  #
  # mart = useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", host="uswest.ensembl.org", dataset="mmusculus_gene_ensembl")
  # geneMaps<-getBM(attributes=c("ensembl_gene_id","mgi_symbol"), filters="mgi_symbol",values=Re.DE[,1],mart=mart)
  #
  # colnames(geneMaps)=c("ensembl_gene_id","geneID")
  #
  # Re.DE.2<-merge(Re.DE,geneMaps,by="geneID")
  #
  # print(head(Re.DE.2))
  #
  # Re.DE.Jun<-merge(Re.DE.2,Re.Jun.2,by="ensembl_gene_id")
  #
  #
  # print(dim(Re.DE.Jun))
  #
  # Re.DE.Jun.2<-Re.DE.Jun[-which(is.na(Re.DE.Jun[,8])),]
  #
  # n<-dim(Re.DE.Jun.2)[1]
  # ran.p<-runif(n)
  #
  #
  # Re<-list(Re.DE.Jun.with.NA=Re.DE.Jun,Re.DE.Jun.without.NA=cbind(Re.DE.Jun.2,ran.p))
  #
  # #plot()
  #
  # pairs(~padj+geneWisePadj+mostSigPadjust+numKnown+ran.p,data=Re[[2]][,c(8,22,24,26,33)],main="4 p values vs SJ")
  # #boxplot(cbind(Re[[2]][,c(8,22,24)],ran.p))
  #
  #
  # Re

}
#' Title
#'
#' @param libname
#' @param pkgname
#'
#' @return
#' @export
#'
#' @examples

.onAttach <- function(libname, pkgname)
{
  if(.Platform$OS.type=="windows" && .Platform$GUI=="Rgui" ) {
    winMenuAddItem("Vignettes","goseq","shell.exec(system.file(\"doc\",\"goseq.pdf\",package=\"goseq\"))")
  }
}

#These two variables are required for automatic fetching of categories to function.  Their purpose is to take the UCSC genome and gene ID values given when looking up length data and convert them to the names used for the same organism and gene identifier in the organism packages.

#Mappings that are primarily required by getgo, the purpose of this is to convert the UCSC genome IDs, to the bioconductor organism names, e.g. "mm"->"org.Mm."
.ORG_PACKAGES=paste("org.",c("Ag.eg","At.tair","Bt.eg","Ce.eg","Cf.eg","Dm.eg","Dr.eg","EcK12.eg","EcSakai.eg","Gg.eg","Hs.eg","Mm.eg","Mmu.eg","Pf.plasmo","Pt.eg","Rn.eg","Sc.sgd","Ss.eg","Xl.eg"),sep='')
names(.ORG_PACKAGES)=c("anoGam","Arabidopsis","bosTau","ce","canFam","dm","danRer","E. coli K12","E. coli Sakai","galGal","hg","mm","rheMac","Malaria","panTro","rn","sacCer","susScr","xenTro")

#These are the only formats supported by getgo at the moment, the purpose is to convert the USCC gene
#ID formats, to the shorthand used by the bioconductor organism packages, .e.g. "refGene"->"ENSEMBL"
.ID_MAP=c("eg","eg","ENSEMBL","SYMBOL","sgd","plasmo","tair")
names(.ID_MAP)=c("knownGene","refGene","ensGene","geneSymbol","sgd","plasmo","tair")

#Below are the exceptions to the function name for gene to go term mappings
.ORG_GOMAP_FUNCTION=c("GO2ALLEGS","GO2ALLTAIRS","GO2ALLORFS","GO2ALLORFS")
names(.ORG_GOMAP_FUNCTION)=c("default","org.At.tair","org.Pf.plasmo","org.Sc.sgd")

#TxDb Length databases
.TXDB_ORGS=c("ce6","dm3","hg18","hg19","hg38","mm10","mm9","rn4","rn5","sacCer2","sacCer3")
#' @ Use QoRTs to get the count of subfeatures in each gene
#'
#' @ cmd1
#' @ inputfile:bam file, gene annotation file,
#' @ outfile
#'
#' @return
#' @export
#'
#' @examples
#'
#' dir.name="/media/H_driver/Aimin_project/GOSJ_example_data/"
#' dir.name="/media/H_driver/Aimin_project/GOSJ_STAR_Bam/"
#'
#' file.name=dir(dir.name,recursive = TRUE,pattern="sorted.bam")
#' file.name.whole<-paste0(dir.name,file.name)
#' file.name.selected<-file.name.whole
#' file.name.selected.2<-as.list(file.name.selected)
#' names(file.name.selected.2)=sapply(strsplit(file.name.selected,split="\\/"),"[[",6)
#'
#' file.name.selected.3<-file.name.selected.2[-c(1,3)]
#'
#' cmd1="java -Xmx4000M -jar /home/aiminyan/QoRTs/QoRTsFullExampleData/QoRTsRelease/QoRTs.jar QC --noGzipOutput --keepMultiMapped"
#'
#' cmd2="java -Xmx5000M -jar /home/aiminyan/QoRTs/QoRTsFullExampleData/QoRTsRelease/QoRTs.jar QC --noGzipOutput"
#'
#' gtf1="/media/H_driver/Aimin_project/GTF_Files/Homo_sapiens.GRCh38.84.processed.sorted.2.gtf"
#'
#' re.out<-lapply(file.name.selected.3,callQoRT,gtf_file=gtf1,runing_cmd=cmd2)
#'

callQoRT<-function(input_file,runing_cmd,gtf_file){

  inputfile=paste(input_file,gtf_file,sep=" ")

  outfile=paste("",sapply(strsplit(input_file,split="\\/"),"[[",2),sapply(strsplit(input_file,split="\\/"),"[[",3),
                sapply(strsplit(input_file,split="\\/"),"[[",4),sapply(strsplit(input_file,split="\\/"),"[[",6),sep="/")


  cmd2=paste(runing_cmd,inputfile,outfile,sep=" ")

  print(cmd2)


  system(cmd2, intern = TRUE, ignore.stderr = TRUE)

  return(cmd2)

}
#' gene2cat
#'
#' Given a gene name and the object built from reading Gmt file, and find the pathways that this gene corresponds to
#'
#' @param gene_name
#' @param re
#'
#' @return
#' @export
#'
#' @examples
#'
#'
gene2cat <- function(gene_name,re) {
  z<-re$genesets
  res <- lapply(z, function(ch) grep(gene_name, ch))
  res2<-sapply(res, function(x) length(x) > 0)
  gene2cat<-list(re$geneset.names[res2])
  gene2cat
}
#' gene2cat2
#'
#' gene2cat2 is called in Gmt2GeneCat to convert a gmt file to a list with the name of each element being gene name
#' and each element being the pathways that this gene corresponds to
#'
#' @param gmt_input_file
#'
#' @return
#' @export
#'
#' @examples
#'
#' gene.2.cat.hallmark<-gene2cat2("/media/H_driver/2015/Nimer_Cheng/h.all.v5.1.symbols.gmt")
#'
gene2cat2 <- function(gmt_input_file) {

  re<-GSA.read.gmt(gmt_input_file)
  gene.name<-unique(do.call(c,re$genesets))

  gene.2.cat<-sapply(gene.name,gene2cat,re)
  names(gene.2.cat)<-gene.name
  gene.2.cat

}
#'#############################################################################
#'#Description: Attempts to fetch the categories specified for the given genes, genome and gene ID format
#'#Notes: Relies on the bioconductor organism packages being installed for whatever genome you specify
#'#Author: Matthew Young
#'#Date Modified: 17/2/2015
#' Title
#'
#' @param genes
#' @param genome
#' @param id
#' @param fetch.cats
#'
#' @return
#' @export
#'
#' @examples
getgo=function(genes,genome,id,fetch.cats=c("GO:CC","GO:BP","GO:MF")){
	#Check for valid input
	if(any(!fetch.cats%in%c("GO:CC","GO:BP","GO:MF","KEGG"))){
		stop("Invaled category specified.  Categories can only be GO:CC, GO:BP, GO:MF or KEGG")
	}
	#Convert from genome ID to org.__.__.db format

  print(gsub("[0-9]+",'',genome))
  cat(".ORG_PACKAGES")
  print(names(.ORG_PACKAGES))
	orgstring=as.character(.ORG_PACKAGES[match(gsub("[0-9]+",'',genome),names(.ORG_PACKAGES))])
	#Multimatch or no match
	if(length(orgstring)!=1){
		stop("Couldn't grab GO categories automatically.  Please manually specify.")
	}
	#Load the library
	library(paste(orgstring,"db",sep='.'),character.only=TRUE)
	#What is the default ID that the organism package uses?
	coreid=strsplit(orgstring,"\\.")[[1]][3]

	#Now we need to convert it into the naming convention used by the organism packages
	userid=as.character(.ID_MAP[match(id,names(.ID_MAP))])
	#Multimatch or no match
	if(is.na(userid) | (length(userid)!=1)){
		stop("Couldn't grab GO categories automatically.  Please manually specify.")
	}
	#The (now loaded) organism package contains a mapping between the internal ID and whatever
        #the default is (usually eg), the rest of this function is about changing that mapping to
	#point from categories to the ID specified
	#Fetch the mapping in its current format
	#Because GO is a directed graph, we need to get not just the genes associated with each ID,
	#but also those associated with its children.  GO2ALLEGS does this.
	core2cat=NULL
	if(length(grep("^GO",fetch.cats))!=0){
	        #Get the name of the function which maps gene ids to go terms
		#usually this will be "GO2ALLEG"
	        gomapFunction=.ORG_GOMAP_FUNCTION[orgstring]
		if(is.na(gomapFunction)) gomapFunction=.ORG_GOMAP_FUNCTION["default"]
		x=toTable(get(paste(orgstring,gomapFunction,sep='')))
		#Keep only those ones that we specified and keep only the names
#		core2cat=x[x$Ontology%in%gsub("^GO:",'',fetch.cats),1:2]
		x[!x$Ontology%in%gsub("^GO:",'',fetch.cats),2]<-"Other"
		core2cat=x[,1:2]
		colnames(core2cat)=c("gene_id","category")
	}
	if(length(grep("^KEGG",fetch.cats))!=0){
		x=toTable(get(paste(orgstring,"PATH",sep='')))
		#Either add it to existing table or create a new one
		colnames(x)=c("gene_id","category")
		if(!is.null(core2cat)){
			core2cat=rbind(core2cat,x)
		}else{
			core2cat=x
		}
	}

	#Now we MAY have to convert the "gene_id" column to the format we are using
	if(coreid!=userid){
		#The mapping between user id and core id, don't use the <USER_ID>2<CORE_ID> object as the naming is not always consistent
		user2core=toTable(get(paste(orgstring,userid,sep='')))
		#Throw away any user ID that doesn't appear in core2cat
		user2core=user2core[user2core[,1]%in%core2cat[,1],]
		#Make a list version of core2cat, we'll need it
		list_core2cat=split(core2cat[,2],core2cat[,1])
		#Now we need to replicate the core IDs that need replicating
		list_core2cat=list_core2cat[match(user2core[,1],names(list_core2cat))]
		#Now we can replace the labels on this list with the user ones from user2core,
		#but there will be duplicates, so we have to unlist, label, then relist
		user2cat=split(unlist(list_core2cat,FALSE,FALSE),rep(user2core[,2],sapply(list_core2cat,length)))
		#Now we only want each category listed once for each entry...
		user2cat=sapply(user2cat,unique)
		###In case you don't believe that this works as it should, here is the slow as all hell way for comparison...
		###Make first list
		##list_user2core=split(user2core[,1],user2core[,2])
		###Make the second
		##list_core2cat=split(core2cat[,2],core2cat[,1])
		###Go through each entry in first list and expand using second...
		##user2cat=sapply(list_user2core,function(u){unique(unlist(list_core2cat[u],FALSE,FALSE))})

	}else{
		#We don't need to convert anything (WOO!), so just make it into a list
		user2cat=split(core2cat[,2],core2cat[,1])
		user2cat=sapply(user2cat,unique)
	}
	#remove any empty strings
	user2cat=lapply(user2cat,function(x){
	        if(length(x)>1) x=x[x!="Other"]
		x })

	## we don't like case sensitivity
	names(user2cat)<-toupper(names(user2cat))

	#Now look them up
	return(user2cat[toupper(genes)])
}
#' Title
#'
#' @param genes
#' @param genome
#' @param id
#' @param fetch.cats
#'
#' @return
#' @export
#'
#' @examples
#'
#' getgo3(rownames(Re.Go.adjusted.by.number.junction.2[[2]]),"mm10","ensGene",fetch.cats=c("GO:BP"))
#'
getgo3=function(genes,genome,id,fetch.cats=c("GO:CC","GO:BP","GO:MF")){
  #Check for valid input
  if(any(!fetch.cats%in%c("GO:CC","GO:BP","GO:MF","KEGG"))){
    stop("Invaled category specified.  Categories can only be GO:CC, GO:BP, GO:MF or KEGG")
  }
  #Convert from genome ID to org.__.__.db format
  orgstring=as.character(.ORG_PACKAGES[match(gsub("[0-9]+",'',genome),names(.ORG_PACKAGES))])
  #Multimatch or no match
  if(length(orgstring)!=1){
    stop("Couldn't grab GO categories automatically.  Please manually specify.")
  }
  #Load the library
  library(paste(orgstring,"db",sep='.'),character.only=TRUE)
  #What is the default ID that the organism package uses?
  coreid=strsplit(orgstring,"\\.")[[1]][3]

  #Now we need to convert it into the naming convention used by the organism packages
  userid=as.character(.ID_MAP[match(id,names(.ID_MAP))])
  #Multimatch or no match
  if(is.na(userid) | (length(userid)!=1)){
    stop("Couldn't grab GO categories automatically.  Please manually specify.")
  }
  #The (now loaded) organism package contains a mapping between the internal ID and whatever
  #the default is (usually eg), the rest of this function is about changing that mapping to
  #point from categories to the ID specified
  #Fetch the mapping in its current format
  #Because GO is a directed graph, we need to get not just the genes associated with each ID,
  #but also those associated with its children.  GO2ALLEGS does this.
  core2cat=NULL
  if(length(grep("^GO",fetch.cats))!=0){
    #Get the name of the function which maps gene ids to go terms
    #usually this will be "GO2ALLEG"
    gomapFunction=.ORG_GOMAP_FUNCTION[orgstring]
    if(is.na(gomapFunction)) gomapFunction=.ORG_GOMAP_FUNCTION["default"]
    x=toTable(get(paste(orgstring,gomapFunction,sep='')))
    #Keep only those ones that we specified and keep only the names
    #		core2cat=x[x$Ontology%in%gsub("^GO:",'',fetch.cats),1:2]
    x[!x$Ontology%in%gsub("^GO:",'',fetch.cats),2]<-"Other"
    core2cat=x[,1:2]
    colnames(core2cat)=c("gene_id","category")
  }
  if(length(grep("^KEGG",fetch.cats))!=0){
    x=toTable(get(paste(orgstring,"PATH",sep='')))
    #Either add it to existing table or create a new one
    colnames(x)=c("gene_id","category")
    if(!is.null(core2cat)){
      core2cat=rbind(core2cat,x)
    }else{
      core2cat=x
    }
  }

  #Now we MAY have to convert the "gene_id" column to the format we are using
  if(coreid!=userid){
    #The mapping between user id and core id, don't use the <USER_ID>2<CORE_ID> object as the naming is not always consistent
    user2core=toTable(get(paste(orgstring,userid,sep='')))
    #Throw away any user ID that doesn't appear in core2cat
    user2core=user2core[user2core[,1]%in%core2cat[,1],]
    #Make a list version of core2cat, we'll need it
    list_core2cat=split(core2cat[,2],core2cat[,1])
    #Now we need to replicate the core IDs that need replicating
    list_core2cat=list_core2cat[match(user2core[,1],names(list_core2cat))]
    #Now we can replace the labels on this list with the user ones from user2core,
    #but there will be duplicates, so we have to unlist, label, then relist
    user2cat=split(unlist(list_core2cat,FALSE,FALSE),rep(user2core[,2],sapply(list_core2cat,length)))
    #Now we only want each category listed once for each entry...
    user2cat=sapply(user2cat,unique)
    ###In case you don't believe that this works as it should, here is the slow as all hell way for comparison...
    ###Make first list
    ##list_user2core=split(user2core[,1],user2core[,2])
    ###Make the second
    ##list_core2cat=split(core2cat[,2],core2cat[,1])
    ###Go through each entry in first list and expand using second...
    ##user2cat=sapply(list_user2core,function(u){unique(unlist(list_core2cat[u],FALSE,FALSE))})

  }else{
    #We don't need to convert anything (WOO!), so just make it into a list
    user2cat=split(core2cat[,2],core2cat[,1])
    user2cat=sapply(user2cat,unique)
  }
  #remove any empty strings
  user2cat=lapply(user2cat,function(x){
    if(length(x)>1) x=x[x!="Other"]
    x })

  ## we don't like case sensitivity
  names(user2cat)<-toupper(names(user2cat))

  #Now look them up
  return(user2cat[toupper(genes)])
}
#' goseq2
#'Description: goseq2 function computes p-values for each gene set and output the list of differentially expressed genes within the gene set. The user can select serveral statistical methods for computing p-values - xxx. 

#' Modifying goseq to generate GO term with the listed DEgene
#'
#' @param pwf: probability weight function
#' @param genome: genome you use
#' @param id: which gene id
#' @param gene2cat: a list with gene as name and category as value
#' @param test.cats: the category you choose
#' @param method: which method for calculting overpresentation p-value
#' @param repcnt: the number of replications
#' @param use_genes_without_cat: whether using genes without category or not
#'
#' @return
#' @export
#'
#' @examples
#' data(gene.model)
#' gene_model<-gene.model
#'
#' GO.wall.DE_interest.geneGL=goseq2(Gene.based.DE.feature.based.DE$pwfGeneGL,"mm10","ensGene",gene.model=gene_model)
#'
#' GO.wall.DE_interes.geneFT=goseq2(Gene.based.DE.feature.based.DE$pwfGeneFeature,"mm10","ensGene",gene.model=gene_model)
#'
#' GO.wall.DE_interest.FtFT=goseq2(Gene.based.DE.feature.based.DE$pwfFeatureFeature,"mm10","ensGene",gene.model=gene_model)
#'
#'
#'
#'
goseq2=function(pwf,genome,id,gene.model,gene2cat=NULL,test.cats=c("GO:CC","GO:BP","GO:MF"),method="Wallenius",repcnt=2000,use_genes_without_cat=FALSE){
  ################# Input pre-processing and validation ###################
  #Do some validation of input variables
  if(any(!test.cats%in%c("GO:CC","GO:BP","GO:MF","KEGG"))){
    stop("Invalid category specified.  Valid categories are GO:CC, GO:BP, GO:MF or KEGG")
  }
  if((missing(genome) | missing(id))){
    if(is.null(gene2cat)){
      stop("You must specify the genome and gene ID format when automatically fetching gene to GO category mappings.")
    }
    #If we're using user specified mappings, this obviously isn't a problem
    genome='dummy'
    id='dummy'
  }
  if(!any(method%in%c("Wallenius","Sampling","Hypergeometric"))){
    stop("Invalid calculation method selected.  Valid options are Wallenius, Sampling & Hypergeometric.")
  }
  if(!is.null(gene2cat) && (!is.data.frame(gene2cat) & !is.list(gene2cat))){
    stop("Was expecting a dataframe or a list mapping categories to genes.  Check gene2cat input and try again.")
  }

  #Factors are evil
  pwf=unfactor(pwf)
  gene2cat=unfactor(gene2cat)

  ###################### Data fetching and processing ########################
  if(is.null(gene2cat)){
    #When we fetch the data using getgo it will be in the list format
    message("Fetching GO annotations...")
    gene2cat=getgo3(rownames(pwf),genome,id,fetch.cats=test.cats)
    names(gene2cat)=rownames(pwf)

    #cat("OK")
    #Do the two rebuilds to remove any nulls
    cat2gene=reversemapping(gene2cat)
    gene2cat=reversemapping(cat2gene)

    #print(cat2gene)
    #print(gene2cat)

  }else{
    #The gene2cat input accepts a number of formats, we need to check each of them in term
    message("Using manually entered categories.")
    #The options are a flat mapping (that is a data frame or matrix) or a list, where the list can be either gene->categories or category->genes
    if(class(gene2cat)!="list"){
      #it's not a list so it must be a data.frame, work out which column contains the genes
      genecol_sum=as.numeric(apply(gene2cat,2,function(u){sum(u%in%rownames(pwf))}))
      genecol=which(genecol_sum!=0)
      if(length(genecol)>1){
        genecol=genecol[order(-genecol_sum)[1]]
        warning(paste("More than one possible gene column found in gene2cat, using the one headed",colnames(gene2cat)[genecol]))
      }
      if(length(genecol)==0){
        genecol=1
        warning(paste("Gene column could not be identified in gene2cat conclusively, using the one headed",colnames(gene2cat)[genecol]))
      }
      othercol=1
      if(genecol==1){othercol=2}
      #Now put it into our delicious listy format
      gene2cat=split(gene2cat[,othercol],gene2cat[,genecol])
      #Do the appropriate builds
      cat2gene=reversemapping(gene2cat)
      gene2cat=reversemapping(cat2gene)
    }
    #!!!!
    #The following conditional has been flagged as a potential issue when using certain
    #types of input where the category names are the same as gene names (which seems like
    #something you should avoid anyway...).  Leave it for now
    #!!!!
    #We're now garunteed to have a list (unless the user screwed up the input) but it could
    #be category->genes rather than the gene->categories that we want.
    if(sum(unique(unlist(gene2cat,use.names=FALSE))%in%rownames(pwf))>sum(unique(names(gene2cat))%in%rownames(pwf))){
      gene2cat=reversemapping(gene2cat)
    }
    #Alright, we're garunteed a list going in the direction we want now.  Throw out genes which we will not use
    gene2cat=gene2cat[names(gene2cat)%in%rownames(pwf)]

    #Rebuild because it's a fun thing to do
    cat2gene=reversemapping(gene2cat)
    gene2cat=reversemapping(cat2gene)

    ## make sure we remove duplicate entries .. e.g. see
    ## http://permalink.gmane.org/gmane.science.biology.informatics.conductor/46876
    cat2gene=lapply(cat2gene,function(x){unique(x)})
    gene2cat=lapply(gene2cat,function(x){unique(x)})
  }

  nafrac=(sum(is.na(pwf$pwf))/nrow(pwf))*100
  if(nafrac>50){
    warning(paste("Missing length data for ",round(nafrac),"% of genes.  Accuarcy of GO test will be reduced.",sep=''))
  }
  #Give the genes with unknown length the weight used by the median gene (not the median weighting!)
  pwf$pwf[is.na(pwf$pwf)]=pwf$pwf[match(sort(pwf$bias.data[!is.na(pwf$bias.data)])[ceiling(sum(!is.na(pwf$bias.data))/2)],pwf$bias.data)]

  ###################### Calculating the p-values ########################
  # Remove all the genes with unknown GOterms
  unknown_go_terms=nrow(pwf)-length(gene2cat)
  if((!use_genes_without_cat) && unknown_go_terms>0 ){
    message(paste("For",unknown_go_terms,"genes, we could not find any categories. These genes will be excluded."))
    message("To force their use, please run with use_genes_without_cat=TRUE (see documentation).")
    message("This was the default behavior for version 1.15.1 and earlier.")
    pwf=pwf[rownames(pwf) %in% names(gene2cat),]
  }
  #A few variables are always useful so calculate them
  cats=names(cat2gene)
  DE=rownames(pwf)[pwf$DEgenes==1]
  num_de=length(DE)
  num_genes=nrow(pwf)
  pvals=data.frame(category=cats,over_represented_pvalue=NA,under_represented_pvalue=NA,stringsAsFactors=FALSE,numDEInCat=NA,numInCat=NA)
  if(method=="Sampling"){
    #We need to know the number of DE genes in each category, make this as a mask that we can use later...
    num_DE_mask=rep(0,length(cats))
    a=table(unlist(gene2cat[DE],FALSE,FALSE))

    num_DE_mask[match(names(a),cats)]=as.numeric(a)
    num_DE_mask=as.integer(num_DE_mask)
    #We have to ensure that genes not associated with a category are included in the simulation, to do this they need an empty entry in the gene2cat list
    gene2cat=gene2cat[rownames(pwf)]
    names(gene2cat)=rownames(pwf)
    message("Running the simulation...")
    #Now do the actual simulating
    lookup=matrix(0,nrow=repcnt,ncol=length(cats))
    for(i in 1:repcnt){
      #A more efficient way of doing weighted random sampling without replacment than the built in function
      #The order(runif...)[1:n] bit picks n genes at random, weighting them by the PWF
      #The table(as.character(unlist(...))) bit then counts the number of times this random set occured in each category
      a=table(as.character(unlist(gene2cat[order(runif(num_genes)^(1/pwf$pwf),decreasing=TRUE)[1:num_de]],FALSE,FALSE)))
      lookup[i,match(names(a),cats)]=a
      pp(repcnt)
    }
    message("Calculating the p-values...")
    #The only advantage of the loop is it uses less memory...
    #for(i in 1:length(cats)){
    #	pvals[i,2:3]=c((sum(lookup[,i]>=num_DE_mask[i])+1)/(repcnt+1),(sum(lookup[,i]<=num_DE_mask[i])+1)/(repcnt+1))
    #	pp(length(cats))
    #}
    pvals[,2]=(colSums(lookup>=outer(rep(1,repcnt),num_DE_mask))+1)/(repcnt+1)
    pvals[,3]=(colSums(lookup<=outer(rep(1,repcnt),num_DE_mask))+1)/(repcnt+1)
  }
  if(method=="Wallenius"){
    message("Calculating the p-values...")
    #All these things are just to make stuff run faster, mostly because comparison of integers is faster than string comparison
    degenesnum=which(pwf$DEgenes==1)
    #Turn all genes into a reference to the pwf object

    cat2genenum=relist(match(unlist(cat2gene),rownames(pwf)),cat2gene)
    #This value is used in every calculation, by storing it we need only calculate it once
    alpha=sum(pwf$pwf)

    #Each category will have a different weighting so needs its own test
    pvals[,2:3]=t(sapply(cat2genenum,function(u){

     #The number of DE genes in this category
      num_de_incat=sum(degenesnum%in%u)

      #The total number of genes in this category
      num_incat=length(u)

      #This is just a quick way of calculating weight=avg(PWF within category)/avg(PWF outside of category)
      avg_weight=mean(pwf$pwf[u])
      weight=(avg_weight*(num_genes-num_incat))/(alpha-num_incat*avg_weight)
      if(num_incat==num_genes){ weight=1 } #case for the root GO terms
      #Now calculate the sum of the tails of the Wallenius distribution (the p-values)

    #cat(num_de_incat,"\t",num_incat,"\t",num_genes,"\t",num_de,"\t",weight,"\n")

      c(dWNCHypergeo(num_de_incat,num_incat,num_genes-num_incat,num_de,weight)
        +pWNCHypergeo(num_de_incat,num_incat,num_genes-num_incat,num_de,weight,lower.tail=FALSE),
        pWNCHypergeo(num_de_incat,num_incat,num_genes-num_incat,num_de,weight))
    }))
    }
  if(method=="Hypergeometric"){
    message("Calculating the p-values...")
    #All these things are just to make stuff run faster, mostly because comparison of integers is faster than string comparison
    degenesnum=which(pwf$DEgenes==1)
    #Turn all genes into a reference to the pwf object
    cat2genenum=relist(match(unlist(cat2gene),rownames(pwf)),cat2gene)
    #Simple hypergeometric test, one category at a time
    pvals[,2:3]=t(sapply(cat2genenum,function(u){
      #The number of DE genes in this category
      num_de_incat=sum(degenesnum%in%u)
      #The total number of genes in this category
      num_incat=length(u)
      #Calculate the sum of the tails of the hypergeometric distribution (the p-values)
      c(dhyper(num_de_incat,num_incat,num_genes-num_incat,num_de)+phyper(num_de_incat,num_incat,num_genes-num_incat,num_de,lower.tail=FALSE),phyper(num_de_incat,num_incat,num_genes-num_incat,num_de))
    }))
  }

  #Populate the count columns...
  degenesnum=which(pwf$DEgenes==1)
  cat2genenum=relist(match(unlist(cat2gene),rownames(pwf)),cat2gene)
  pvals[,4:5]=t(sapply(cat2genenum,function(u){
    c(sum(degenesnum%in%u),length(u))
  }))

  #Got the name of DE gene in each GO
  #print(head(cat2gene))

  DE_pwf=rownames(pwf[degenesnum,])

  #cat2degenenum=relist(match(unlist(cat2gene),rownames(pwf[degenesnum,])),cat2gene)
  #print(cat2degenenum)

  pvals.6<-sapply(cat2gene,function(u,DE_pwf){
    #c(sum(degenesnum%in%u),length(u))
    #c(rownames(pwf)[u[-which(is.na(u))]])
    x<-u[which(u %in% DE_pwf)]
    x
  },DE_pwf)

  #cat("After matching\n")

  #cat("DE_ensemble\n")

  #print(unique(as.character(unlist2(pvals.6))))
  #cat(length(unique(as.character(unlist2(pvals.6)))),"\n")

  pvals.6.gene.symbol<-sapply(pvals.6,function(u,gene.model){
    y<-gene.model[which(as.character(gene.model[,3]) %in% u),1]
    y
    },gene.model)

  #cat("DE_symbol\n")
  #print(unique(as.character(unlist2(pvals.6.gene.symbol))))
  #cat(length(unique(as.character(unlist2(pvals.6.gene.symbol)))),"\n")

  #Convert list to data frame
  pvals.6.df<-list_to_df(pvals.6)
  #cat("pvals_6_df","\n")
  #cat(dim(pvals.6.df)[1],"\n")

  pvals.6.gene.symbol.df<-list_to_df(pvals.6.gene.symbol)
  #cat("pvals_6_gene_symbol_df","\n")
  #cat(dim(pvals.6.gene.symbol.df)[1],"\n")

  #cat(length(unique(as.character(pvals.6.gene.symbol.df[,2]))),"\n")
  #print(unique(as.character(pvals.6.gene.symbol.df[,2])))

  dataset2<- pvals.6.gene.symbol.df
  dataset2[sapply(dataset2, is.list)] <-
    sapply(dataset2[sapply(dataset2, is.list)],
           function(x)sapply(x, function(y) paste(unlist(y),collapse=", ") ) )

  #print(dataset2)
  #df_args <- c(pvals.6.gene.symbol.df[,2], sep=",")
  #df_args.data<-do.call(paste, df_args)

  #print(class(dataset2))
  #cat(dim(dataset2),"\n")

  temp.gene.name=unique(apply(dataset2[,2],1,c))
  temp.gene.name.2=unique(trim(unlist(strsplit(temp.gene.name,split=","))))
  #print(class(temp.gene.name.2))

  #cat(length(temp.gene.name.2),"\n")
  #print(temp.gene.name.2)
  DE_from_GO<-temp.gene.name.2

  colnames(pvals.6.df)=c("category","DEgene_ID")
  colnames(pvals.6.gene.symbol.df)=c("category","DEgene_symbol")

  #Finally, sort by p-value
  pvals=pvals[order(pvals$over_represented_pvalue),]

  # Supplement the table with the GO term name and ontology group
  # but only if the enrichment categories are actually GO terms
  if(any(grep("^GO:",pvals$category))){
    GOnames=select(GO.db,keys=pvals$category,columns=c("TERM","ONTOLOGY"))[,2:3]
    colnames(GOnames)<-tolower(colnames(GOnames))
    pvals=cbind(pvals,GOnames)
  }

  # And return
  pvals.2<-merge(pvals,pvals.6.df,by="category",sort = FALSE)

  pvals.3<-merge(pvals.2,pvals.6.gene.symbol.df,by="category",sort = FALSE)

  pvals.4<-list(GO=pvals.3,DE_GO=DE_from_GO)

  return(pvals.4)

}
#' Title
#'
#' @param list_for_df
#'
#' @return
#' @export
#'
#' @examples
#' test.goseq3.df<-list_to_df(test.goseq3)
#'
list_to_df <- function(list_for_df) {
  list_for_df <- as.list(list_for_df)

  nm <- names(list_for_df)
  if (is.null(nm))
    nm <- seq_along(list_for_df)

  df <- data.frame(name = nm, stringsAsFactors = FALSE)
  df$value <- unname(list_for_df)
  df
}
#' @title Plot the probability Weighting Function
#' @description
#' @usage
#' plotPWF2(pwf, binsize = "auto", pwf_col = 3, pwf_lwd = 2,
#' xlab = "Biased Data in <binsize> gene bins.", ylab = "Proportion DE", ...)
#' @param pwf probability weigth function
#' @param binsize the number of gene in each bin(gene set)
#' @param pwf_col the color for the fitted line of probability weigth function
#' @param pwf_lwd the font for the fitted line of probability weigth function
#' @param xlab the label for x-axis
#' @param ylab the label for y-axis
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#'
#' png("~/GOSJ/Figure/pwfGeneGL.tiff")
#' plotPWF2(Gene.based.DE.feature.based.DE$pwfGeneGL)
#' dev.off()
#'
#' png("~/GOSJ/Figure/pwfGeneFeature.tiff")
#' plotPWF2(Gene.based.DE.feature.based.DE$pwfGeneFeature)
#' dev.off()
#'
#' png("~/GOSJ/Figure/pwfFeatureGL.tiff")
#' plotPWF2(Gene.based.DE.feature.based.DE$pwfFeatureGL)
#' dev.off()
#'
#' png("~/GOSJ/Figure/pwFeatureFeature.tiff")
#' plotPWF2(Gene.based.DE.feature.based.DE$pwfFeatureFeature)
#' dev.off()
#'
#'
plotPWF2<-function (pwf, binsize = "auto", pwf_col = 3, pwf_lwd = 2, xlab = "Biased Data in <binsize> gene bins.",
                    ylab = "Proportion DE", ...)
{
  w = !is.na(pwf$bias.data)
  print(w)
  o = order(pwf$bias.data[w])
  print(o)

  rang = max(pwf$pwf, na.rm = TRUE) - min(pwf$pwf, na.rm = TRUE)
  if (rang == 0 & binsize == "auto")
    binsize = 1000
  if (binsize == "auto") {
    binsize = max(1, min(100, floor(sum(w) * 0.08)))
    resid = rang
    oldwarn = options()$warn
    options(warn = -1)
    while (binsize <= floor(sum(w) * 0.1) & resid/rang >
           0.001) {
      binsize = binsize + 100
      splitter = ceiling(1:length(pwf$DEgenes[w][o])/binsize)
      de = sapply(split(pwf$DEgenes[w][o], splitter), mean)
      binlen = sapply(split(as.numeric(pwf$bias.data[w][o]),
                            splitter), mean)
      resid = sum((de - approx(pwf$bias.data[w][o], pwf$pwf[w][o],
                               binlen)$y)^2)/length(binlen)
    }
    options(warn = oldwarn)
  }
  else {
    splitter = ceiling(1:length(pwf$DEgenes[w][o])/binsize)
    print(splitter)
    de = sapply(split(pwf$DEgenes[w][o], splitter), mean)
    print(de)
    binlen = sapply(split(as.numeric(pwf$bias.data[w][o]),
                          splitter), median)
    print(binlen)
  }
  xlab = gsub("<binsize>", as.character(binsize), xlab)
  if ("xlab" %in% names(list(...))) {
    if ("ylab" %in% names(list(...))) {
      plot(binlen, de, ...)
    }
    else {
      plot(binlen, de, ylab = ylab, ...)
    }
  }
  else if ("ylab" %in% names(list(...))) {
    plot(binlen, de, xlab = xlab, ...)
  }
  else {
    plot(binlen, de, xlab = xlab, ylab = ylab, ...)
  }
  lines(pwf$bias.data[w][o], pwf$pwf[w][o], col = pwf_col,
        lwd = pwf_lwd)

  return(de)

}
#' Title
#'
#' @param countfiles
#' @param countdata
#' @param samplenames
#' @param design
#' @param flat.gff.file
#' @param test.formula1
#' @param analysis.type
#' @param nCores
#' @param use.exons
#' @param use.junctions
#' @param use.known.junctions
#' @param use.novel.junctions
#' @param use.multigene.aggregates
#' @param gene.names
#' @param verbose
#' @param method.countVectors
#' @param noDESeqMatrix
#'
#' @return
#' @export
#'
#' @examples
#' dir.name="/Volumes/Bioinformatics$/2015/Nimer_Cheng/"
#' file.sample="decoder.bySample.txt"
#' file.count="_junction_seq_new_gtf_7/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt"
#' file.gff="Mus_musculus.GRCm38.83.JunctionSeq.flat.gff"

#' path.file.sample<-paste0(dir.name,file.sample)
#' decoder.bySample<-read.table(path.file.sample,header=T,stringsAsFactors = F)
#' print(decoder.bySample)

#' #Get count file
#' path.file.count<-paste0(dir.name,decoder.bySample$sample.ID,file.count)
#' countFiles<-paste0(path.file.count)
#' print(countFiles)

#' #Get annotation file
#' flat.file.gff<-paste0(dir.name,file.gff)
#' print(flat.file.gff)
#'
#' method.countVectors = c("geneLevelCounts","sumOfAllBinsForGene","sumOfAllBinsOfSameTypeForGene")
#' method.countVectors <- match.arg(method.countVectors)
#'
#' jscs = readJunctionSeqCounts(countfiles = as.character(sample.files),
#' samplenames = sample.names,
#' design = Re,
#' flat.gff.file = flat.gff.file,
#'  verbose = verbose,
#'  use.junctions = use.junctions,
#'  use.novel.junctions = use.novel.junctions,
#'  use.known.junctions = use.known.junctions,
#'  use.exons = use.exons,
#'  use.multigene.aggregates = use.multigene.aggregates,
#'  nCores = nCores,
#'  method.countVectors = method.countVectors,
#'  test.formula1 = test.formula1,
#'  gene.names = gene.names
#'  )
#'
#'
#'
readJunctionSeqCounts <- function(countfiles = NULL, countdata = NULL,
                                  samplenames,  design,
                                  flat.gff.file=NULL,
                                  test.formula1 = formula(~ sample + countbin + condition : countbin),
                                  analysis.type = c("junctionsAndExons","junctionsOnly","exonsOnly"),
                                  nCores = 1,
                                  use.exons = NULL, use.junctions = NULL,
                                  use.known.junctions = TRUE,
                                  use.novel.junctions = TRUE,
                                  use.multigene.aggregates = FALSE,
                                  gene.names = NULL,
                                  verbose = TRUE,
                                  method.countVectors = c("geneLevelCounts","sumOfAllBinsForGene","sumOfAllBinsOfSameTypeForGene"),
                                  noDESeqMatrix = FALSE)
{
  method.countVectors <- match.arg(method.countVectors)

  if(isTRUE(verbose)) {
    message("-> STARTING readJunctionSeqCounts (",date(),")")
  }


  analysis.type <- match.arg(analysis.type)
  if(is.null(use.junctions) && is.null(use.exons)){
    if(analysis.type == "junctionsAndExons"){
      use.junctions <- TRUE
      use.exons <- TRUE
    } else if(analysis.type == "junctionsOnly"){
      use.junctions <- TRUE
      use.exons <- FALSE
    } else if(analysis.type == "exonsOnly"){
      use.junctions <- FALSE
      use.exons <- TRUE
    }
  } else {
    if(is.null(use.junctions) || is.null(use.exons)){
      stop(paste0("Illegal syntax! If parameter use.junctions or use.exons are used, then BOTH must be set!\n use.junctions = '",use.junctions,"', use.exons = '",use.exons,"'"))
    }

    if(use.junctions && use.exons){
      analysis.type <- "junctionsAndExons"
    } else if(use.junctions && (! use.exons)){
      analysis.type <- "junctionsOnly"
    } else  if((! use.junctions) && use.exons){
      analysis.type <- "exonsOnly"
    } else {
      stop("Illegal syntax! Parameters use.exons and use.junctions cannot both be false!")
    }
  }
  if(isTRUE(verbose)){
    message("---> RJSC; (v",packageVersion("JunctionSeq"),")")
    message("---> RJSC: samplenames: ",paste0(samplenames, collapse=","))
    message("---> RJSC: flat.gff.file: ",flat.gff.file)
    message("---> RJSC: use.exons:",use.exons)
    message("---> RJSC: use.junctions:",use.junctions)
    message("---> RJSC: use.novel.junctions:",use.novel.junctions)
  }


  if((is.null(countfiles) && is.null(countdata))){
    stop("Fatal error: Either countfiles OR countdata must be set! Both are null!")
  }
  if(  (!is.null(countfiles)) && (!is.null(countdata))   ){
    stop("Fatal error: Either countfiles OR countdata must be set! Both are non-null!")
  }

  stopifnot( class(design) == "data.frame" )

  for(i in 1:ncol(design)){
    if( ! is.factor(design[[i]])){
      stop("ERROR: design must be a data.frame composed entirely of factors!")
    }
  }

  if(! is.null(countfiles)){
    lf <- lapply( countfiles, function(x)
      read.table( x, header=FALSE,stringsAsFactors=FALSE ) )
  } else {
    lf <- countdata
  }

  if( !all( sapply( lf[-1], function(x) all( x$V1 == lf[1]$V1 ) ) ) )
    stop( "Count files have differing gene ID column." )
  if(isTRUE(verbose)) message("---> File read complete.");

  dcounts <- sapply( lf, `[[`, "V2" )
  rownames(dcounts) <- lf[[1]][,1]
  dcounts <- dcounts[ substr(rownames(dcounts),1,1)!="_", ]

  bin.type <- sapply( rownames(dcounts),
                      function(x){
                        substr(strsplit(x, ":",fixed=TRUE)[[1]][2],1,1)
                      })
  raw.geneID <- sapply( rownames(dcounts),
                        function(x){
                          strsplit(x, ":",fixed=TRUE)[[1]][1]
                        })

  if(isTRUE(verbose)) message(paste0("---> Extracted counts. Found ",dim(dcounts)[1]," features so far."));

  geneCountTable <- dcounts[bin.type == "A",, drop=FALSE]
  rownames(geneCountTable) <- sapply(strsplit(rownames(geneCountTable), ":"),"[[",1)
  colnames(geneCountTable) <- as.character(samplenames)
  use.bins <- bin.type != "A"

  if(isTRUE(verbose)) message(paste0("---> Extracted gene-level counts. Found: ",dim(geneCountTable)[1], " genes and aggregate-genes."))
  if(isTRUE(verbose)) message(paste0("---> Removed gene features. Found: ",sum(use.bins), " features to be included so far."))

  if(isTRUE(! use.junctions) && isTRUE(! use.exons)){
    stop("FATAL ERROR: At least one of: use.junctions or use.exons must be set to TRUE. Otherwise you've got no data to test!")
  }

  if(isTRUE(! use.junctions)){
    use.bins <- use.bins & bin.type != "J" & bin.type != "N"
    if(isTRUE(verbose)) message(paste0("---> Removed splice junction features. Found: ",sum(use.bins), " features to be included so far."))
  }
  if(isTRUE(! use.novel.junctions)){
    use.bins <- use.bins & bin.type != "N"
    if(isTRUE(verbose)) message(paste0("---> Removed novel splice junction features. Found: ",sum(use.bins), " features to be included so far."))
  }
  if(isTRUE(! use.known.junctions)){
    use.bins <- use.bins & bin.type != "J"
    if(isTRUE(verbose)) message(paste0("---> Removed known splice junction features. Found: ",sum(use.bins), " features to be included so far."))
  }
  if(isTRUE(! use.exons)){
    use.bins <- use.bins & bin.type != "E"
    if(isTRUE(verbose)) message(paste0("---> Removed exon features. Found: ",sum(use.bins), " features to be included so far."))
  }

  is.multiGene.aggregate <- grepl("+", raw.geneID, fixed=TRUE)
  multiGene.aggregate.IDs <- unique(raw.geneID[is.multiGene.aggregate & use.bins])
  multiGene.aggregate.ct <- length(multiGene.aggregate.IDs)
  multiGene.aggregate.geneCt <- length(unlist( strsplit(multiGene.aggregate.IDs, "+",fixed=TRUE)))

  if(isTRUE(verbose)) message("---> Note: ",sum(is.multiGene.aggregate[use.bins])," counting bins from overlapping genes")
  if(isTRUE(verbose)) message("--->          There are ",multiGene.aggregate.ct,    " multigene aggregates.")
  if(isTRUE(verbose)) message("--->          There are ",multiGene.aggregate.geneCt," genes that are part of an aggregate.")
  if(isTRUE(! use.multigene.aggregates)){
    use.bins <- use.bins & (! is.multiGene.aggregate)
    if(isTRUE(verbose)) message(paste0("---> Removed multigene-aggregate features. Found: ",sum(use.bins), " features to be included so far."))
  }


  if(isTRUE(verbose)) message(paste0("---> Final feature count: ",sum(use.bins), " features to be included in the analysis."))
  dcounts <- dcounts[use.bins,, drop=FALSE]
  bin.type <- bin.type[use.bins]

  if(isTRUE(verbose)) message("---> Extracted feature counts.");

  colnames(dcounts) <- as.character(samplenames)
  splitted <- strsplit(rownames(dcounts), ":")
  exons <- sapply(splitted, "[[", 2)
  genesrle <- sapply( splitted, "[[", 1)

  if(isTRUE(verbose)) message("---> counts complete.");



  if(! is.null(flat.gff.file)){
    if(isTRUE(verbose)) message("-----> reading annotation...");
    anno.data <- readAnnotationData(flat.gff.file)
    if(isTRUE(verbose)) message("-----> formatting annotation...");
    #featureName     featureType     chrom   start   end     strand  gene_id part_number     transcripts

    exoninfo <- data.frame(chr = anno.data$chrom, start = anno.data$start, end = anno.data$end, strand = anno.data$strand)
    if(isTRUE(verbose)) message("-----> initial generation...");
    rownames(exoninfo) <- anno.data$featureName

    transcripts <- anno.data$transcripts
    transcripts <- gsub("\\+", ";", transcripts)
    names(transcripts) <- rownames(exoninfo)

    matching <- match(rownames(dcounts), rownames(exoninfo))
    if(any(is.na(matching))){
      stop("FATAL ERROR! Annotation file appears to be missing information! Are you sure you are using the correct flattened annotation file, created by prepare_annotation_with_splices.py?")
    }
    if(isTRUE(verbose)) message("-----> creating jscs...");
    jscs <- newJunctionSeqCountSet(countData=dcounts, geneCountData = geneCountTable, design=design, geneIDs=genesrle, countbinIDs=exons, featureIntervals=exoninfo[matching,], transcripts=transcripts[matching])
    jscs@annotationFile <- flat.gff.file
    jscs@flatGffData <- anno.data

    jscs@flatGffGeneData <- readGeneInfo(flat.gff.file)
    jscs <- mapGeneNames(jscs, gene.names)
  } else {
    if(isTRUE(verbose)) message("-> FINISHED readJunctionSeqCounts (",date(),")");
    message("Warning: flat gff annotation not set (via parameter flat.gff.file)! While technically optional, running without the annotation data may make interpretation of the data difficult. Much of the plotting functionality will not work!")
    warning("Warning: flat gff annotation not set (via parameter flat.gff.file)! While technically optional, running without the annotation data may make interpretation of the data difficult. Much of the plotting functionality will not work!")
    jscs <- newJunctionSeqCountSet(countData=dcounts, geneCountData = geneCountTable, design=design, geneIDs=genesrle, countbinIDs=exons);

  }
  attr(jscs,"AltMethods") <- c(attr(jscs,"AltMethods"), method.countVectors = method.countVectors)
  attr(jscs,"CallStack") <- list(deparse(match.call()))

  jscs@analysisType <- analysis.type
  featureChar <- substr(fData(jscs)$countbinID,1,1)
  fData(jscs)[["featureType"]] <- ifelse(featureChar == "E","exonic_part",ifelse(featureChar == "J", "splice_site", "novel_splice_site"))
  varMetadata( featureData(jscs) )[ "featureType", "labelDescription" ] <- "The type of feature (exonic_part,, splice_site, or novel_splice_site)."
  pData(jscs)$countfiles <- countfiles
  if(isTRUE(verbose)) message("-----> generating count vectors... (",date(),")")
  jscs@countVectors <- getAllJunctionSeqCountVectors(jscs, nCores = nCores, method.countVectors); #use.alternate.method = use.alternate.method)
  if(isTRUE(verbose)) message("-----> count vectors generated (",date(),")");

  if(! noDESeqMatrix){
    if(isTRUE(verbose)) message("-----> generating DESeqDataSet... (",date(),")")
    jscs <- makeDESeqDataSetFromJSCS(jscs, test.formula1 = test.formula1)
    if(isTRUE(verbose)) message("-----> DESeqDataSet generated (",date(),")")
    fData(jscs)[["allZero"]] <- (rowSums(counts(jscs)) == 0) |
      (rowSums(counts(jscs@DESeqDataSet)[, colData(jscs@DESeqDataSet)$countbin == "others"]) ==0)
    mcols(jscs@DESeqDataSet)$allZero <- fData(jscs)[["allZero"]]
  }

  return(jscs)
  if(isTRUE(verbose)) message("-> FINISHED readJunctionSeqCounts (",date(),")");
}
#' Title
#'
#' @param map
#'
#' @return
#' @export
#'
#' @examples
reversemapping=function(map){
  tmp=unlist(map,use.names=FALSE)
  names(tmp)=rep(names(map),times=as.numeric(summary(map)[,1]))
  return(split(names(tmp),as.vector(tmp)))
}
