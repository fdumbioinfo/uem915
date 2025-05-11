#' @title Omic bioanalysis workflow
#' @description
#'  Biostatistics analysis: QCs, ANOVA, threshold filtering, venn diagramm, cluster analysis for significant rows, pattern search
#'  
#'  Graphics: volcanoplots, heatmaps, lineplots, boxplots (with Kruskal-Wallis test)
#'  
#'  Functional analysis: MSigDB enrichment analysis, stringDB protein interaction network, topGO analysis (gene ontology)
#'   
#'  See example section to test workflow with internal data.
#'   
#' @param dat data.frame normalize data table
#' @param sif data.frame sample information file including model factors
#' @param annot data.frame annotation with Symbol column for functional analysis
#' @param species character available species: bt ce dr dm gg hs mm pt rn ss xt
#' @param qcs logical quality controls
#' @param model character anova model factors (see details)
#' @param paired character factor for paired design
#' @param nested character factor for nested design
#' @param batch character factor for batch effect design
#' @param addfactor character additionnal factors
#' @param threshold numeric vector from 1 to 160 (see details)
#' @param pattern logical search relevant pattern accross comparisons (see details)
#' @param venn logical venn diagram
#' @param cluster logical row hierarchical clustering using pearson correlation
#' @param nc numeric number of clusters to cut in dendrogramm
#' @param padj character fdr by defaut for Benjamini-Hochberg false discovery correction
#' @param heatmap logical do heatmaps for all lists
#' @param maxheatmap numeric max rows for heatmap
#' @param volcanoplot logical make volcanoplot for each threshold
#' @param lineplot logical do lineplot for significant features
#' @param boxplotrow logical do boxplot for significant features with Kruskal
#' @param ena logical msigdb enrichement analysis (over-representation analysis)
#' @param bg numeric background used for over-representation test
#' @param filtergeneset character regular expression to filter collection geneset (e.g. "reactome|tft")
#' @param enaScoremin numeric for minimum ratio ena
#' @param enamin numeric min list size for functional analysis
#' @param sample numeric subset analysis
#' @param dopar numeric core number
#' @param path character results directory path
#' @param dirname character results directory name
#' @param zip logical compress results directory if TRUE
#' @param remove logical remove uncompress results directory if TRUE
#' @return omic results directory
#' @details
#' Use uem915::env() to load required libraries before uem915::omic() (see example)
#' 
#' Accepted values for threshold param (1 to 150): see list -> uem915:::thresholdlist %>% lapply("[",c(1,2)) %>% unlist %>% matrix(ncol=2,byrow = T) %>% data.frame %>% setNames(c("pval","fc"))
#' 
#'  Pattern: search relevant profiles among up and down comparison combinations
#' 
#' filtergeneset param: see geneset collections -> moalannotgene::genesetdb %>% lapply(names)
#' 
#' Experimental design examples:
#' 
#'  - model = "TREAMTENT" for 1-way anova
#'
#'  - model = "TREATMENT+TIME+TREATMENT*TIME" for 3-ways anova with interaction
#'  
#'  - model = "TREATMENT", paired = "CASE" for 2-ways paired anova
#'  
#'  - model = "TREATMENT", batch = "BATCH" for 2-ways anova with remove batch effect (<=> paired anova) 
#'  
#'  - model = "TREATMENT+PHENOTYPE" for 2-ways anova
#'  
#'  - model = "TREATMENT", addfactor = "PHENOTYPE" for 2-ways anova but venn, cluster and pattern are not applied to addfactor
#'    
#'  - model =  "TREATMENT", nested = "TREATMENTinCASE" for 2-ways nested design
#'  
#'  - model =  "TREATMENT", paired = "CASE", batch = "BATCH" for 3-ways paired anova with remove batch effect
#'  
#'  
#' Limitations: 
#' 
#'  - Complete block designs only
#'   
#'  - Use dopar = 2 to decrease computing resources
#'  
#'  - sample param will subsets random rows in dat and decrease analysis time.
#' 
#' Annotation updates: 05112023 for gene and ensembl, MSigDB 7.5.1, StringDB 12.0
#' 
#' Input format:
#'  
#'  - Use uem915::input() to load external data from tsv files
#' 
#'  - norm data table must contains IDs in first column
#'  
#'  - norm data (in columns) and sample information (in rows) must have same order
#'  
#'  - norm data and annotation rows must have same order
#'  
#'  - use uem915::annot() function to annotate IDs (Symbols, ensembl and gene ids accepted)
#'  
#' @examples
#' # # test omic() with internal dataset GSE65055:
#' # # loading libraries
#' # library(uem915)
#' # uem915::env()
#' # # loading norm data
#' # moal:::GSE65055normdata -> normdata
#' # normdata %>% head
#' # # loading sample information file
#' # moal:::GSE65055sampledata -> sampledata
#' # sampledata %>% head
#' # # ordering factor levels
#' # sampledata$ANEUPLOIDY %>% ordered( c("Control","T13","T18","T21") ) -> sampledata$ANEUPLOIDY
#' # sampledata$TISSUE %>% as.factor -> sampledata$TISSUE
#' # # annotation
#' # normdata$rowID %>% moal::annot( species = "hs", idtype = "GENE" ) -> annotdata
#' # # omic analysis
#' # moal::omic(
#' #  dat = normdata, sif = sampledata, annot = annotdata, species = "hs",
#' #  model = "ANEUPLOIDY", batch = "TISSUE", threshold = c(6,126),
#' #  heatmap = T, lineplot = T, boxplotrow = T,
#' #  venn = F, cluster = F, pattern = F,
#' #  ena = T, network = F, topgo = F,
#' #  sample = NULL, dopar = NULL, zip = F,
#' #  dirname = "test", path = "." )
#' 
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom dplyr select inner_join left_join
#' @importFrom tidyselect all_of
#' @importFrom rlang .data
#' @importFrom ggplot2 ggsave
#' @importFrom graphics hist
#' @importFrom grDevices pdf graphics.off
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom limma removeBatchEffect
#' @importFrom stats model.matrix setNames p.adjust
#' @importFrom gridExtra marrangeGrob
#' @importFrom utils packageVersion
#' @importFrom tidyselect all_of
#' @export
omic <- function(
  dat, sif, annot = NULL, species = "hs",
  model = NULL, paired = NULL, nested = NULL, batch = NULL , addfactor = NULL,
  qcs = TRUE, threshold = c(c(2,6,11),c(2,6,11)+120), padj = "none",
  pattern = TRUE, venn = TRUE, cluster = TRUE, nc = c(2,3,6,12),
  heatmap = TRUE, maxheatmap = NULL, volcanoplot = TRUE, lineplot = TRUE, boxplotrow = TRUE,
  ena = TRUE, enamin = 2, filtergeneset = "all", enaScoremin = 1.1, bg = 25000,
  sample = NULL , dopar = NULL, path = ".", dirname = NULL, zip = FALSE, remove = FALSE )
{
  Sys.time() -> start
  list() -> log
  i=j=k=1
  #
  # data preprocessing
  #
  paste( "\n#######################\n", sep="" ) %>% cat
  paste( "###  uem915",packageVersion("uem915") %>% gsub( "\\.","",.),"  omic  ###\n", sep="" ) %>% cat
  paste( "#######################\n", sep="" ) %>% cat
  paste("#\n" , sep = "" ) %>% cat
  paste("data preprocessing...\n" , sep = "" ) %>% cat
  # normdata
  colnames(dat)[1] <- "rowID" ; mode(dat[,1]) <- "character"
  # annotdata
  if( is.null(annot) ){ dat[,1] %>% as.character %>% data.frame( "rowID"=. , "Symbol"=. , stringsAsFactors=F ) -> Annot0 ; ena <- FALSE ; topgo <- FALSE ; network <- FALSE }
  if( !is.null(annot) & all(colnames(annot)!="Symbol") ){ annot %>% cbind( "Symbol"=dat[,1] %>% as.character ) -> Annot0 ; ena <- FALSE ; topgo <- FALSE ; network <- FALSE }
  if( !is.null(annot) & any(colnames(annot)=="Symbol") )
  {
    colnames(annot)[1] <- "rowID" ; mode(annot[,1]) <- "character"
    annot -> Annot0 ; Annot0$Symbol -> Symbol
    # NCBI annotation
    if(any(colnames(annot)=="GeneID")){annot %>% dplyr::select(-.data$GeneID) -> annot}
    Symbol -> symbollist
    Symbol %>% annot(species=species) -> AnnotEnaSpecies
    AnnotEnaSpecies %>% setNames( c("Symbol","GeneID","descriptionNCBI","chromoseNCBI","GeneTypeNCBI","SpeciesNCBI","SynonymsNCBI") ) -> AnnotEnaSpecies
    data.frame( rowID = Annot0$rowID %>% as.character , AnnotEnaSpecies , Annot0 %>% dplyr::select( -.data$rowID , -.data$Symbol ) ) -> Annot0
  }
  #
  # metadata
  #
  if( !is.null(sample) )
  {
    sample(x=rownames(dat),size=sample) %>% as.numeric -> sel
    dat %>% dplyr::slice(sel) -> dat
    Annot0 %>% dplyr::slice(sel) -> Annot0
    if(!is.null(annot) & any(colnames(annot)=="Symbol")){ AnnotEnaSpecies %>% data.frame %>% dplyr::slice(sel) %>% data.frame -> AnnotEnaSpecies }
  }
  # output
  paste("output done.\n",sep="") %>% cat
  ifelse(is.null(dirname),"omic" -> DirName,paste("omic_",dirname,sep="" ) -> DirName)
  paste( DirName,"_",ncol(dat)-1,"_",nrow(dat),sep="" ) -> DirName
  path %>% file.path( DirName ) -> Path ; Path %>% dir.create
  # normdata
  Path %>% file.path( "inputdata" ) %>% dir.create
  paste("datanorm_",ncol(dat)-1,"_",nrow(dat),".tsv",sep="") %>% file.path( Path , "inputdata" , . ) -> FileName
  dat %>% output( FileName )
  # sample metadata
  paste("metadata_",ncol(sif),"_",nrow(sif),".tsv",sep="") %>% file.path( Path, "inputdata", . ) -> FileName
  sif %>% output( FileName )
  # annotdata
  paste( "annotations_",ncol(Annot0),"_",nrow(Annot0),".tsv",sep="" ) %>% file.path( Path , "inputdata" , . ) -> FileName
  Annot0 %>% output( FileName )
  paste("Preprocessing done.\n" , sep = "" ) %>% cat
  #
  # QCs all
  #
  if( qcs )
  {
    paste( "QC_all_",ncol(dat)-1,"_",nrow(dat),sep="") -> DirNameQC
    dat %>% qc( sif = sif, dirname = DirNameQC, path = Path, inputdata = F )
    paste("QC all done.\n") %>% cat
  }
  #
  # anova
  #
  # anova preprocessing
  if( !is.null(model) )
  {
    paste( "anova preprocessing...\n" ) %>% cat
    # model factors
    model %>% gsub( " " , "" , . ) -> model
    model %>% gsub( " " , "" , . ) -> RbeModel
    model %>% strsplit( "\\+" ) %>% unlist %>% unique -> ModelFactors
    ModelFactors -> RbeFactors
    ModelFactors -> AnovaFactors
    ModelFactors -> CompFactors
    # interaction factor
    INTERACTION <- F
    if( any( grepl( "\\*" , ModelFactors ) ) )
    {
      # remove interaction in AnovaFactors
      AnovaFactors %>% grep("\\*" , . , invert = T ) %>% AnovaFactors[.] -> AnovaFactors
      RbeFactors %>% grep("\\*" , . , invert = T ) %>% RbeFactors[.] -> RbeFactors
      # create interaction factor
      ModelFactors %>% grep("\\*" , . , invert = F ) %>% ModelFactors[.] -> IntFactor0
      IntFactor0 %>% strsplit("\\*") %>% unlist -> IntFactor1
      IntFactor1 %>% paste0( collapse = "x" ) -> IntFactorName
      IntFactor1[1] %>% paste("^" , . , "$" , sep = "") %>% grep( sif %>% colnames ) %>% sif[ , . ] -> F1
      IntFactor1[2] %>% paste("^" , . , "$" , sep = "") %>% grep( sif %>% colnames ) %>% sif[ , . ] -> F2
      expand.grid( F1 %>% levels , F2 %>% levels ) -> LevelsIntFactor0
      paste( LevelsIntFactor0[,1] , LevelsIntFactor0[,2] , sep = "-") -> LevelsIntFactor1
      paste( F1 , F2 , sep = "-" ) -> IntFactor2
      IntFactor2 %>% ordered( LevelsIntFactor1 ) -> IntFactor3
      sif %>% cbind( IntFactor3 ) -> sif
      colnames( sif )[ ncol(sif)  ] <- IntFactorName
      # change interaction name in CompFactors
      CompFactors[ grep("\\*" , CompFactors , invert = F ) ] <- IntFactorName
      INTERACTION <- T }
    # nested factor
    NESTED <- F
    if( !is.null(nested) )
    {
      nested %>% strsplit("in") %>% unlist -> NestedFactor0
      nested %>% sub("in","%in%",.) -> nested
      # create nested factor
      NestedFactor0[1] %>% paste("^" , . , "$" , sep = "") %>% grep( sif %>% colnames ) %>% sif[ , . ] -> F1
      NestedFactor0[2] %>% paste("^" , . , "$" , sep = "") %>% grep( sif %>% colnames ) %>% sif[ , . ] -> F2
      expand.grid( F1 %>% levels , F2 %>% levels ) -> LevelsNestedFactor0
      paste( LevelsNestedFactor0[,1],LevelsNestedFactor0[,2],sep="-") -> LevelsNestedFactor1
      paste(F1,F2,sep=":") -> NestedFactor1
      NestedFactor1 %>% as.factor -> NestedFactor2
      sif %>% cbind( NestedFactor2 ) -> sif
      colnames(sif)[ncol(sif)] <- "NESTED"
      NESTED <- T
      paired <- "NESTED"
    }
    ### additional factor
    if( !is.null( addfactor ) )
    {
      addfactor %>% gsub( " " , "" , . ) -> addfactor
      addfactor %>% strsplit( "\\+" ) %>% unlist %>% unique -> AddFactors
      AddFactors %>% c( AnovaFactors ) -> AnovaFactors
      AddFactors %>% c( CompFactors ) -> RbeFactors
      paste0( model , "+" , addfactor ) -> model
      paste0( model , "+" , addfactor ) -> RbeModel
    }
    # create model name for log
    paste( length(AnovaFactors)," way(s) anova", sep="" ) -> ModelName
    # interaction
    if( INTERACTION ){ paste( ModelName,"with interaction between",gsub("x"," and ",IntFactorName), sep=" " ) -> ModelName }
    # nested factor
    if( NESTED ){paste( "nested paired",ModelName,sep=" " ) -> ModelName ; paste0( model,"+","NESTED" ) -> model }
    # paired
    if( !is.null(paired) & !NESTED ){ ModelName %>% paste( "paired " , . , sep = "" ) -> ModelName ; paste0( model , "+" , paired ) -> model }
    # add batch
    if( !is.null(batch) ){ paste(ModelName,"and batch effect for factor",batch,sep=" ") -> ModelName ; paste0(model,"+",batch ) -> model }
    paste("model :","Y = mu + ",model," + E \n",ModelName,"\n") %>% cat
    #
    # RBE : remove batch effect
    #
    # for batch
    RBE <- F
    if( !is.null(batch) & is.null(paired) )
    {
      dat -> DatRbe
      batch %>% paste( "^" , . , "$" , sep = ""  ) %>% paste0( collapse = "|" ) -> grep
      sif %>% colnames %>% grep( grep , . ) %>% sif[ , . ] -> Batch0
      model.matrix( eval( expr = parse( text = paste0( "~" , RbeModel ) ) ), sif )  -> Design
      dat %>% dplyr::select(-1) %>% removeBatchEffect( design = Design , batch = Batch0 ) %>%
        data.frame( "rowID" = dat[,1] %>% as.character, . , stringsAsFactors = F ) -> dat
      # output
      paste( "datanorm_RBE_" , batch ,"_" , ( ncol( dat ) - 1 ) , "_" , nrow( dat ) , ".tsv" , sep = "" ) %>%
        file.path( Path , "inputdata", . ) -> FileName
      dat %>% output( FileName )
      # QC RBE
      paste( "QC_all_RBE_",ncol(dat)-1,"_",nrow(dat),sep="" ) -> DirName
      dat %>% qc( sif = sif , path = Path , dirname = DirName , inputdata = F )
      # add factor for anova
      AnovaFactors %>% c(batch) -> AnovaFactors
      RBE <- T
    }
    # paired (or nested paired)
    if( !is.null(paired) & is.null(batch) )
    {
      dat -> DatRbe
      paired %>% paste( "^" , . , "$" , sep = ""  ) %>% paste0( collapse = "|" ) -> grep
      sif %>% colnames %>% grep( grep , . ) %>% sif[ , . ] -> Paired0
      model.matrix( eval( expr = parse( text = paste0( "~" , RbeModel ) ) ), sif )  -> Design
      dat %>% dplyr::select(-1) %>%
        removeBatchEffect( . , design = Design , batch = Paired0 ) %>%
        data.frame( "rowID" = dat[,1] %>% as.character, . , stringsAsFactors = F ) -> dat
      # output
      paste( "datanorm_RBE_" , paired ,"_" , ( ncol( dat ) - 1 ) , "_" , nrow( dat ) , ".tsv" , sep = "" ) %>%
        file.path( Path , "inputdata", . ) -> FileName
      dat %>% output( FileName )
      # QC RBE
      paste( "QC_all_RBE_",ncol(dat)-1,"_",nrow(dat),sep="") -> DirName
      dat %>% qc( sif = sif , path = Path , dirname = DirName , inputdata = F )
      AnovaFactors %>% c( paired ) -> AnovaFactors
      RBE <- T
    }
    # both batch and paired (or nested paired)
    if( !is.null(paired) & !is.null(batch) )
    {
      dat -> DatRbe
      paired %>% paste( "^" , . , "$" , sep = ""  ) %>% paste0( collapse = "|" ) -> grep
      sif %>% colnames %>% grep( grep , . ) %>% sif[ , . ] -> Paired0
      batch %>% paste( "^" , . , "$" , sep = ""  ) %>% paste0( collapse = "|" ) -> grep
      sif %>% colnames %>% grep( grep , . ) %>% sif[ , . ] -> Batch0
      model.matrix( eval( expr = parse( text = paste0( "~" , RbeModel ) ) ), sif )  -> Design
      dat %>% dplyr::select(-1) %>%
        removeBatchEffect( . , design = Design , batch = Batch0 , batch2 = Paired0 ) %>%
        data.frame( "rowID" = dat[,1] %>% as.character, . , stringsAsFactors = F ) -> dat
      # output
      paste( "datanorm_RBE_" , paired ,"_" , ( ncol( dat ) - 1 ) , "_" , nrow( dat ) , ".tsv" , sep = "" ) %>%
        file.path( Path , "inputdata", . ) -> FileName
      dat %>% output( FileName )
      # QC RBE
      paste( "QC_RBE_all_" , paired , "_" , ncol( dat )-1 , "_" , nrow( dat ) , sep = "" ) -> DirName
      dat %>% qc( sif = sif , path = Path , dirname = DirName , inputdata = F )
      AnovaFactors %>% c( paired , batch ) -> AnovaFactors
      RBE <- T
    }
    ### log infos
    log[[1]] <- ModelFactors
    log[[2]] <- model
    log[[3]] <- ModelName
    #
    # anova
    #
    paste( "anova processing...\n", sep ="" ) %>% cat
    # dopar
    if( is.null( dopar ) ){ detectCores() -> NbCore }else{ parallel::detectCores()/dopar -> NbCore }
    makeCluster( NbCore ) -> cl ; registerDoParallel( cl )
    i=1
    All0 <- foreach( i=1:nrow(dat), .packages=c("uem915","magrittr","dplyr","broom"), .combine=rbind ) %dopar%
      { anova( dat = sif %>% dplyr::select( AnovaFactors ) %>% data.frame( y=dat[i,-1] %>% as.numeric, . ) , model=model ) %>% unlist }
    # anovastat data
    All0 %>% colnames %>% grep( "(Sumsq_|Meansq_|Fratio)",. ,value=F ) -> sel
    All0 %>% as.data.frame %>% dplyr::select( all_of( sel ) ) %>%
      data.frame( "rowID" = dat[,1] %>% as.character,. ,stringsAsFactors=F ) -> Allstat0
    All0 %>% as.data.frame %>% dplyr::select( -all_of( sel ) ) -> All
    # remove comp for batch factor
    if( !is.null( batch ) )
    {
      batch %>% paste( "^" , . , "$" , sep = ""  ) %>% paste0( collapse = "|" ) -> grep
      sif %>% colnames %>% grep( grep , . ) %>% sif[ , . ] -> Batch0
      Batch0 %>% levels %>% rev %>% as.character %>% combn( 2 , simplify = F ) %>%
        lapply( paste0 , collapse = "vs" ) %>% unlist %>% paste0( collapse = "|" ) -> grep
      All %>% colnames %>% grep( grep , . , value = F ) -> sel
      All %>% as.data.frame %>% dplyr::select( -all_of( sel ) ) -> All
    }
    # remove comparisons for paired factor
    if( !is.null( paired ) )
    {
      paired %>% paste( "^" , . , "$" , sep = ""  ) %>% paste0( collapse = "|" ) -> grep
      sif %>% colnames %>% grep( grep , . ) %>% sif[ , . ] -> Paired0
      Paired0 %>% as.factor %>% levels %>% rev %>% as.character %>% combn( 2 , simplify = F ) %>%
        lapply(gsub,pattern=":",replacement="-") %>% lapply( paste0 , collapse = "vs" ) %>% unlist %>% paste0( collapse = "|" ) -> grep
      All %>% colnames %>% grep( grep , . , value = F ) -> sel
      All %>% as.data.frame %>% dplyr::select( -all_of( sel ) ) -> All
    }
    # add sd and mean to annotation
    dat %>% dplyr::select(-1) %>% "^"( 2 , . ) %>% apply( 1 , mean ) -> Mean0
    dat %>% dplyr::select(-1) %>% "^"( 2 , . ) %>% apply( 1 , sd ) -> Sd0
    data.frame(
      "rowID" = Annot0[,1] , "mean" = Mean0 , "sd" = Sd0,
      Annot0[,-1] %>% as.data.frame %>% setNames( colnames( Annot0 )[-1] ) , stringsAsFactors = F ) -> Annot1
    # output anova all
    cbind( "rowID" = dat[,1] %>% as.character, All , Annot1[,-1], stringsAsFactors = F ) -> r0
    paste( "anova_all_" , dim(dat)[2]-1 , "_", dim(dat)[1] , ".tsv" , sep = "" ) -> FileName
    if( padj == "none" ){ PVALADJ <- FALSE }else{ PVALADJ <- TRUE ; Padj0 <- padj }
    if( PVALADJ )
    { 
      r0 %>% colnames %>% grep("^p_",.) -> sel
      foreach( i=1:length(sel) ) %do%
        { r0[,sel[i]] %>% unlist %>% stats::p.adjust( . , method = Padj0) -> r0[,sel[i]] }
    }
    r0 %>% output( Path %>% file.path( FileName ) )
    paste( "anova done.\n" ) %>% cat
    # anovastat all
    paste("anovastat all processing... \n",sep="") %>% cat
    paste("anovastat_all_" ,ncol(dat)-1,"_",nrow(dat),sep="") -> DirName
    Allstat0 %>% anovastat(dat=.,path=Path,dirname=DirName)
    paste("anovastat all done. \n",sep="") %>% cat
    # anovastat all RBE
    if( RBE )
    {
      paste("anovastat all RBE processing... \n",sep="") %>% cat
      AllRbe <- foreach(i=1:nrow(dat),.packages=c("uem915","magrittr","dplyr","broom"),.combine=rbind ) %dopar%
        { anova( dat = sif %>% dplyr::select(AnovaFactors) %>% data.frame( y=DatRbe[i,-1] %>% as.numeric,.) , model=model ) %>% unlist }
      AllRbe %>% colnames %>% grep( "(Sumsq_|Meansq_|Fratio)", . , value = F ) -> sel
      AllRbe %>% as.data.frame %>% dplyr::select(all_of(sel)) %>% data.frame("rowID"=DatRbe[,1] %>% as.character,.,stringsAsFactors=F) -> AllRbestat
      paste( "anovastat_all_RBE_" ,ncol(dat)-1,"_",nrow(dat),sep="") -> DirName
      AllRbestat %>% anovastat(dat = . , path = Path , dirname = DirName )
      paste("anovastat all RBE done. \n",sep="") %>% cat
    }
  }
  #
  # filter
  #
  if( is.null(threshold) )
  {
    venn <- F ; cluster <- F
    heatmap <- F ; lineplot <- F ; volcanoplot <- F
    ena <- F ; network <- F ; topgo <- F
  }
  if( !is.null( threshold ) )
  {
    ifelse( threshold == "all", 1:length(thresholdlist) -> Threshold0, threshold -> Threshold0  ) 
    paste( "Filter processing...\n", sep="" ) %>% cat
    # anova filter for model factor
    foreach( i=1:length( CompFactors ) ) %do%
      {
        paste( "anova_",CompFactors[i], sep="" ) -> DirName
        anovafilter( dat=r0, sif=sif, annot=Annot1, comp=CompFactors[i],
                     threshold=Threshold0, path=Path, dirname=DirName )
      }
    # anova filter for additionnal factor
    if( !is.null( addfactor ) )
    {
      foreach( i=1:length( AddFactors ) ) %do%
        {
          paste( "anova_",AddFactors[i], sep="" ) -> DirName
          anovafilter( dat=r0, sif=sif, annot=Annot1, comp=AddFactors[i],
                       threshold=Threshold0, path=Path, dirname=DirName )
        }
    }
    paste( "anova factor filter done. \n", sep="" ) %>% cat
    # comp filter for model factor
    foreach( i=1:length( CompFactors ) ) %do%
      {
        paste( "DiffLists_",CompFactors[i], sep="" ) -> DirName
        anovacompfilter( dat=r0, sif=sif, annot=Annot1, comp=CompFactors[i],
                         threshold=Threshold0, path=Path, dirname=DirName )
      }
    ### comp filter for additionnal factor
    if( !is.null( addfactor ) )
    {
      foreach( i=1:length( AddFactors ) ) %do%
        {
          paste( "DiffLists_",AddFactors[i], sep="" ) -> DirName
          anovacompfilter( dat=r0, sif=sif, annot=Annot1, comp=AddFactors[i],
                           threshold=Threshold0, path=Path, dirname=DirName )
        }
    }
    paste( "DiffLists filter done. \n", sep="" ) %>% cat
    # pattern filter for model factor
    if( pattern )
    {
      foreach( i=1:length( CompFactors ) ) %do%
        {
          paste( "pattern_",CompFactors[i], sep="" ) -> DirName
          anovapatternfilter( dat=r0, sif=sif, annot=Annot1, comp=CompFactors[i],
                              threshold=Threshold0, path=Path, dirname=DirName )
        }
    }
    paste( "pattern filter done. \n", sep="" ) %>% cat
    paste( "Filter done.\n", sep="" ) %>% cat
    # add sum row to directory name
    Path %>% list.dirs( recursive=F ) -> d0
    d0 %>% basename %>% grep( "^DiffLists_|^anova_|^pattern_" ,. ) %>% d0[.] -> d1
    d2 <- foreach( i=1:length( d1 ) ) %do% { d1[i] %>% list.dirs( recursive=F ) }
    d2 %>% unlist -> d3
    d3 %>% lapply( list.files, recursive=T ) %>% lapply( length ) %>% "=="(.,0) %>% which -> sel
    if(length( sel ) > 0)
    {
      d3[sel] -> Path0
      Path0 %>% basename %>% paste( "_0",sep="" ) %>% file.path( Path0 %>% dirname, . ) -> Path1
      file.rename( from=Path0, to=Path1 )
      d3[-sel] -> d3
    }
    # create statList0 for rename file, anovastat, venn and cluster
    statList0 <- foreach( i=1:length(d3) ) %do%
      {
        d3[i] %>% list.files( recursive=T, full.names=T ) -> l0
        l1 <- foreach( j=1:length( l0 ) ) %do% { l0[j] %>% input %>% "["(.,1) }
        list( d3[i],l1 %>% unlist %>% as.character %>% unique,
              l1 %>% unlist %>% as.character %>% unique %>% length )
      }
    paste( "stats list done.\n", sep="" ) %>% cat
    # rename
    d3 -> Path0
    Path0 %>% basename %>% paste( "_",statList0 %>% sapply("[[",3),sep="") %>% file.path( Path0 %>% dirname,. ) -> Path1
    file.rename( from=Path0, to=Path1 )
    paste( "rename done.\n", sep="" ) %>% cat
    # anovastat threshold
    paste( "anovastat threshold processing...\n", sep="" ) %>% cat
    foreach( i=1:length( statList0 ) ) %do% { statList0[[i]][[1]] <- Path1[i] }
    statList0 %>% sapply( "[[", 3 ) %>% ">="(.,5) %>% which %>% statList0[.] -> statList1
    if( length( statList1 ) > 0 ){ DoAnovastatLists <- T }else{ paste( "Not enough rows to make anovastat \n" ) %>% cat ; DoAnovastatLists <- F }
    if( DoAnovastatLists )
    {
      file.path( Path , "anovastat_threshold" ) %>% dir.create
      foreach( i=1:length(statList1), .packages=c("magrittr","uem915","dplyr","foreach","ggplot2") ) %dopar%
        {
          statList1[[i]][2] %>% unlist %>% as.character %>% data.frame( rowID = .,stringsAsFactors = F ) %>% inner_join( Allstat0,by="rowID" ) -> Allstat1
          statList1[[i]][[1]] %>% strsplit("\\/") %>% unlist %>% length -> Nb0
          statList1[[i]][[1]] %>% strsplit("\\/") %>% unlist %>% "["(Nb0-1) %>% sub( "(.*)_.*" , "\\1" , . ) %>%
            paste("_",statList1[[i]][[1]] %>% strsplit("\\/") %>% unlist %>% "["(Nb0) , sep = "" ) -> DirName
          Allstat1 %>% anovastat(path = Path %>% file.path("anovastat_threshold") , dirname = DirName)
        }
      paste("anovastat threshold done.\n") %>% cat
    }
    # hc threshold
    if( DoAnovastatLists )
    {
      paste("hierarchical clustering threshold preprocessing...\n", sep="") %>% cat
      file.path( Path,"hc_threshold" ) %>% dir.create
      foreach( i=1:length(statList1), .packages=c("magrittr","uem915","dplyr","foreach") ) %dopar%
        {
          statList1[[i]][2] %>% unlist %>% as.character %>% data.frame( rowID=., stringsAsFactors=F ) %>%
            inner_join( dat, by="rowID" ) -> Allstat1
          statList1[[i]][[1]] %>% strsplit("\\/") %>% unlist %>% length -> Nb0
          statList1[[i]][[1]] %>% strsplit("\\/") %>% unlist %>% "["(Nb0-1) %>% sub( "(.*)_.*", "\\1",. ) %>%
            paste( "_",statList1[[i]][[1]] %>% strsplit("\\/") %>% unlist %>% "["(Nb0), sep="" ) -> DirName
          Allstat1 %>% qc( sif=sif, inputdata = F, histo = F , boxplot = F , hc = T , acp = F,
                           path=Path %>% file.path("hc_threshold"), dirname=DirName )
        }
      paste( "hc threshold done.\n" , sep="") %>% cat
    }
    # acp threshold
    if( DoAnovastatLists )
    {
      paste("pca threshold processing...\n",sep="") %>% cat
      file.path( Path , "pca_threshold" ) %>% dir.create
      foreach( i=1:length(statList1), .packages=c("magrittr","uem915","dplyr","foreach") ) %dopar%
        {
          statList1[[i]][2] %>% unlist %>% as.character %>% data.frame( rowID=., stringsAsFactors=F ) %>%
            inner_join( dat, by="rowID" ) -> Allstat1
          statList1[[i]][[1]] %>% strsplit("\\/") %>% unlist %>% length -> Nb0
          statList1[[i]][[1]] %>% strsplit("\\/") %>% unlist %>% "["(Nb0-1) %>% sub( "(.*)_.*" , "\\1" , . ) %>%
            paste("_",statList1[[i]][[1]] %>% strsplit("\\/") %>% unlist %>% "["(Nb0) , sep = "" ) -> DirName
          Allstat1 %>% qc( sif=sif, inputdata=F, histo=F, boxplot=F, hc=F, acp=T,
                           path=Path %>% file.path("pca_threshold"), dirname=DirName )
        }
    }
    paste( "pca threshold done.\n" ) %>% cat
  }
  #
  # Venn
  #
  # venn preprocessing
  if( venn )
  {
    paste( "Venn preprocessing...\n", sep="" ) %>% cat
    Path %>% list.dirs( recursive=F ) -> d0
    d0 %>% basename %>% grep( "^DiffLists" , . ) %>% d0[.] -> d1
    d1 %>% basename %>% sub( "^.*_(.*)", "\\1" , . ) -> FactorVenn0
    FactorVenn0 %>% paste("^",.,"$",sep="") %>% paste0( collapse="|" ) %>% grep( colnames(sif) ) %>% dplyr::select(sif,.) -> FactorVenn1
    # check if more than 2 levels in factors
    FactorVenn1 %>% lapply(levels) %>% lapply(length) %>% ">"(.,2) %>% which -> sel
    if( length(sel) > 0 ){ FactorVenn1[sel] -> FactorVenn2 ; venn <- T }else{ venn <- F }
    paste( "Venn preprocessing done.\n", sep="" ) %>% cat
  }
  # venn processing
  if( venn )
  {
    paste( "Venn processing...\n", sep="" ) %>% cat
    foreach( i=1:ncol(FactorVenn2) ) %do%
      {
        VENN3 <- F ; VENN4 <- F ; VENN5 <- F ; VENNSUP5 <- F ; VENN4INT <- F
        if( FactorVenn2[i] %>% lapply(levels) %>% lapply(length) %>% "=="(.,3) ){ VENN3 <- T }
        if( FactorVenn2[i] %>% lapply(levels) %>% lapply(length) %>% "=="(4) & !INTERACTION ){ VENN4 <- T }
        if( FactorVenn2[i] %>% lapply(levels) %>% lapply(length) %>% "=="(5) & !INTERACTION ){ VENN5 <- T }
        if( FactorVenn2[i] %>% lapply(levels) %>% lapply(length) %>% ">"(5) & !INTERACTION ){ VENNSUP5 <- T }
        if( FactorVenn2[i] %>% lapply(levels) %>% lapply(length) %>% ">="(4) & INTERACTION ){ VENN4INT <- T }
        # VENN3
        if( VENN3 )
        {
          paste( "venn_",colnames(FactorVenn2)[i], sep="" ) -> DirName0
          Path %>% file.path( DirName0 ) %>% dir.create
          d1 %>% basename %>% sub(".*_(.*)","\\1",.) %>% "=="(.,colnames(FactorVenn2)[i]) %>%
            which %>% d1[.] %>% list.files( full.names=T ) -> d2
          foreach( j=1:length(d2) ) %do%
            {
              d2[j] %>% list.files( full.names= T) %>% grep( ".*.tsv$", . ,value=T ) -> ListFileVenn0
              # up
              ListFileVenn0 %>% grep( "_u_",. , value=T ) -> ListFileVenn1
              ListFileVenn1 %>% lapply(input) %>% lapply("[",1) %>% lapply(unlist) %>% lapply(as.character) -> ListFileVenn2
              ListFileVenn1 %>% basename %>% strsplit("_") %>% sapply("[",6) %>% paste0(collapse="x") -> DirNameVenn0
              venn3( list = ListFileVenn2,
                     listnames = ListFileVenn1 %>% basename %>% strsplit("_") %>% sapply("[",6),
                     title = ListFileVenn1 %>% basename %>% strsplit("_") %>% sapply("[",6) %>% paste0(collapse="x"),
                     export=T, plot=F,
                     dirname = d2[j] %>% basename %>% sub("(.*_.*_.*)_.*","\\1",.) %>% paste( "_u_",DirNameVenn0, sep = ""),
                     path=file.path( Path, DirName0 ) )
              # down
              ListFileVenn0 %>% grep( "_d_",. , value=T ) -> ListFileVenn1
              ListFileVenn1 %>% lapply(input) %>% lapply("[",1) %>% lapply(unlist) %>% lapply(as.character) -> ListFileVenn2
              ListFileVenn1 %>% basename %>% strsplit("_") %>% sapply("[",6) %>% paste0(collapse="x") -> DirNameVenn0
              venn3( list = ListFileVenn2,
                     listnames = ListFileVenn1 %>% basename %>% strsplit("_") %>% sapply("[",6),
                     title = ListFileVenn1 %>% basename %>% strsplit("_") %>% sapply("[",6) %>% paste0(collapse="x"),
                     export=T, plot=F,
                     dirname = d2[j] %>% basename %>% sub("(.*_.*_.*)_.*","\\1",.) %>% paste( "_d_",DirNameVenn0, sep = ""),
                     path=file.path( Path, DirName0 ) )
              # up + down
              ListFileVenn0 %>% grep( "_ud_",. , value=T ) -> ListFileVenn1
              ListFileVenn1 %>% lapply(input) %>% lapply("[",1) %>% lapply(unlist) %>% lapply(as.character) -> ListFileVenn2
              ListFileVenn1 %>% basename %>% strsplit("_") %>% sapply("[",6) %>% paste0(collapse="x") -> DirNameVenn0
              venn3( list = ListFileVenn2,
                     listnames = ListFileVenn1 %>% basename %>% strsplit("_") %>% sapply("[",6),
                     title = ListFileVenn1 %>% basename %>% strsplit("_") %>% sapply("[",6) %>% paste0(collapse="x"),
                     export=T, plot=F,
                     dirname = d2[j] %>% basename %>% sub("(.*_.*_.*)_.*","\\1",.) %>% paste( "_ud_",DirNameVenn0, sep = ""),
                     path=file.path( Path, DirName0 ) )
            }
        }
        # VENN4 sans interaction
        if( VENN4 )
        {
          paste( "venn_",colnames(FactorVenn2)[i], sep="" ) -> DirName0
          Path %>% file.path( DirName0 ) %>% dir.create
          d1 %>% basename %>% sub(".*_(.*)","\\1",.) %>% "=="(.,colnames(FactorVenn2)[i]) %>%
            which %>% d1[.] %>% list.files( full.names=T ) -> d2
          foreach( j=1:length(d2) ) %do%
            {
              d2[j] %>% list.files( full.names= T) %>% grep( ".*.tsv$", . ,value=T ) -> ListFileVenn0
              # venn3
              FactorVenn2[,i] %>% levels %>% rev %>% combn( 2, simplify=F ) %>% lapply(paste0,collapse="vs") -> VennComp0
              VennComp0 %>% grep( FactorVenn2[,i] %>% levels %>% "["(1),.,value=T) -> VennComp1
              VennComp1 %>% combn( 3, simplify=F ) -> VennComp2
              foreach( k=1:length(VennComp2) ) %do%
                {
                  VennComp2[k] %>% unlist %>% paste0(collapse="|") -> grep
                  ListFileVenn0 %>% grep( grep, .,value=T ) -> ListFileVenn1
                  # up
                  ListFileVenn1 %>% grep( "_u_",. , value=T ) -> ListFileVenn2
                  ListFileVenn2 %>% lapply(input) %>% lapply("[",1) %>% lapply(unlist) %>% lapply(as.character) -> ListFileVenn3
                  ListFileVenn2 %>% basename %>% strsplit("_") %>% sapply("[",6) %>% paste0(collapse="x") -> DirNameVenn0
                  venn3( list = ListFileVenn3,
                         listnames = ListFileVenn2 %>% basename %>% strsplit("_") %>% sapply("[",6),
                         title = ListFileVenn2 %>% basename %>% strsplit("_") %>% sapply("[",6) %>% paste0(collapse="x"),
                         export=T, plot=F,
                         dirname = d2[j] %>% basename %>% sub("(.*_.*_.*)_.*","\\1",.) %>% paste( "_u_",DirNameVenn0, sep = ""),
                         path=file.path( Path, DirName0 ) )
                  # down
                  ListFileVenn1 %>% grep( "_d_",. , value=T ) -> ListFileVenn2
                  ListFileVenn2 %>% lapply(input) %>% lapply("[",1) %>% lapply(unlist) %>% lapply(as.character) -> ListFileVenn3
                  ListFileVenn2 %>% basename %>% strsplit("_") %>% sapply("[",6) %>% paste0(collapse="x") -> DirNameVenn0
                  venn3( list = ListFileVenn3,
                         listnames = ListFileVenn2 %>% basename %>% strsplit("_") %>% sapply("[",6),
                         title = ListFileVenn2 %>% basename %>% strsplit("_") %>% sapply("[",6) %>% paste0(collapse="x"),
                         export=T, plot=F,
                         dirname = d2[j] %>% basename %>% sub("(.*_.*_.*)_.*","\\1",.) %>% paste( "_d_",DirNameVenn0, sep = ""),
                         path=file.path( Path, DirName0 ) )
                  # up + down
                  ListFileVenn1 %>% grep( "_ud_",. , value=T ) -> ListFileVenn2
                  ListFileVenn2 %>% lapply(input) %>% lapply("[",1) %>% lapply(unlist) %>% lapply(as.character) -> ListFileVenn3
                  ListFileVenn2 %>% basename %>% strsplit("_") %>% sapply("[",6) %>% paste0(collapse="x") -> DirNameVenn0
                  venn3( list = ListFileVenn3,
                         listnames = ListFileVenn2 %>% basename %>% strsplit("_") %>% sapply("[",6),
                         title = ListFileVenn2 %>% basename %>% strsplit("_") %>% sapply("[",6) %>% paste0(collapse="x"),
                         export=T, plot=F,
                         dirname = d2[j] %>% basename %>% sub("(.*_.*_.*)_.*","\\1",.) %>% paste( "_ud_",DirNameVenn0, sep = ""),
                         path=file.path( Path, DirName0 ) )
                }
            }
        }# fi VENN4 == 4 without interaction
        # VENN5 sans interaction
        if( VENN5 )
        {
          paste( "venn_",colnames(FactorVenn2)[i], sep="" ) -> DirName0
          Path %>% file.path( DirName0 ) %>% dir.create
          d1 %>% basename %>% sub(".*_(.*)","\\1",.) %>% "=="(.,colnames(FactorVenn2)[i]) %>%
            which %>% d1[.] %>% list.files( full.names=T ) -> d2
          foreach( j=1:length(d2) ) %do%
            {
              d2[j] %>% list.files( full.names= T) %>% grep( ".*.tsv$", . ,value=T ) -> ListFileVenn0
              # venn4
              FactorVenn2[,i] %>% levels %>% rev %>% combn( 2, simplify=F ) %>% lapply(paste0,collapse="vs") -> VennComp0
              VennComp0 %>% grep( FactorVenn2[,i] %>% levels %>% "["(1),.,value=T) -> VennComp1
              VennComp1 %>% combn( 4, simplify=F ) -> VennComp2
              foreach( k=1:length(VennComp2) ) %do%
                {
                  VennComp2[k] %>% unlist %>% paste0(collapse="|") -> grep
                  ListFileVenn0 %>% grep( grep, .,value=T ) -> ListFileVenn1
                  # up
                  ListFileVenn1 %>% grep( "_u_",. , value=T ) -> ListFileVenn2
                  if( ListFileVenn2 %>% basename %>% sub(".*_(.*).tsv","\\1",.) %>% as.numeric %>% "=="(.,0) %>% "!"(.) %>% all )
                  {
                    ListFileVenn2 %>% lapply(input) %>% lapply("[",1) %>% lapply(unlist) %>% lapply(as.character) -> ListFileVenn3
                    ListFileVenn2 %>% basename %>% strsplit("_") %>% sapply("[",6) %>% paste0(collapse="x") -> DirNameVenn0
                    venn4( list = ListFileVenn3,
                           listnames = ListFileVenn2 %>% basename %>% strsplit("_") %>% sapply("[",6),
                           title = ListFileVenn2 %>% basename %>% strsplit("_") %>% sapply("[",6) %>% paste0(collapse="x"),
                           export=T, plot=F,
                           dirname = d2[j] %>% basename %>% sub("(.*_.*_.*)_.*","\\1",.) %>% paste( "_u_",DirNameVenn0, sep = ""),
                           path=file.path( Path, DirName0 ) )
                  }
                  # down
                  ListFileVenn1 %>% grep( "_d_",. , value=T ) -> ListFileVenn2
                  if( ListFileVenn2 %>% basename %>% sub(".*_(.*).tsv","\\1",.) %>% as.numeric %>% "=="(.,0) %>% "!"(.) %>% all )
                  {
                    ListFileVenn2 %>% lapply(input) %>% lapply("[",1) %>% lapply(unlist) %>% lapply(as.character) -> ListFileVenn3
                    ListFileVenn2 %>% basename %>% strsplit("_") %>% sapply("[",6) %>% paste0(collapse="x") -> DirNameVenn0
                    venn4( list = ListFileVenn3,
                           listnames = ListFileVenn2 %>% basename %>% strsplit("_") %>% sapply("[",6),
                           title = ListFileVenn2 %>% basename %>% strsplit("_") %>% sapply("[",6) %>% paste0(collapse="x"),
                           export=T, plot=F,
                           dirname = d2[j] %>% basename %>% sub("(.*_.*_.*)_.*","\\1",.) %>% paste( "_d_",DirNameVenn0, sep = ""),
                           path=file.path( Path, DirName0 ) )
                  }
                  # up + down
                  ListFileVenn1 %>% grep( "_ud_",. , value=T ) -> ListFileVenn2
                  if( ListFileVenn2 %>% basename %>% sub(".*_(.*).tsv","\\1",.) %>% as.numeric %>% "=="(.,0) %>% "!"(.) %>% all )
                  {
                    ListFileVenn2 %>% lapply(input) %>% lapply("[",1) %>% lapply(unlist) %>% lapply(as.character) -> ListFileVenn3
                    ListFileVenn2 %>% basename %>% strsplit("_") %>% sapply("[",6) %>% paste0(collapse="x") -> DirNameVenn0
                    venn4( list = ListFileVenn3,
                           listnames = ListFileVenn2 %>% basename %>% strsplit("_") %>% sapply("[",6),
                           title = ListFileVenn2 %>% basename %>% strsplit("_") %>% sapply("[",6) %>% paste0(collapse="x"),
                           export=T, plot=F,
                           dirname = d2[j] %>% basename %>% sub("(.*_.*_.*)_.*","\\1",.) %>% paste( "_ud_",DirNameVenn0, sep = ""),
                           path=file.path( Path, DirName0 ) )

                  }
                }
            }
        }# fi VENN5 (levels == 5 without interaction)
        # VENN4 avec interaction
        if( VENN4INT )
        {
          paste( "venn_",IntFactorName, sep="" ) -> DirName0
          Path %>% file.path( DirName0 ) %>% dir.create
          d1 %>% basename %>% sub(".*_(.*)","\\1",.) %>% "=="(.,colnames(FactorVenn2)[i]) %>%
            which %>% d1[.] %>% list.files( full.names=T ) -> d2
          colnames(FactorVenn2)[i] %>% strsplit("x") %>% unlist %>%
            paste("^",.,"$",sep="") %>% paste0( collapse="|" ) %>% grep( colnames(FactorVenn1) ) %>% FactorVenn1[,.] -> FactorVennInt0
          FactorVennInt0 %>% lapply(levels) %>% lapply(length) %>% unlist %>% which.max %>% FactorVennInt0[.] -> FactorVennIntMax0
          FactorVennIntMax0 %>% lapply(levels) %>% unlist %>% combn( 2, simplify=F ) -> FactorVennIntMax1
          #
          foreach( j=1:length(FactorVennIntMax1) ) %do%
            {
              FactorVennIntMax1[j] %>% unlist -> FactorVennIntMax2
              foreach( k=1:length(d2) ) %do%
                {
                  d2[k] %>% list.files( full.names= T) %>% grep( ".*.tsv$", . ,value=T ) -> ListFileVenn0
                  FactorVennIntMax2 %>% paste0(collapse="|") -> grep
                  IntFactor3 %>% levels %>% grep( grep, . , value=T ) -> FactorVennIntMax3
                  FactorVennIntMax3 %>% rev %>% combn( 2 , simplify = F ) %>% lapply(gsub,pattern=":",replacement="-") %>%
                    lapply( paste0 , collapse = "vs" ) %>% unlist -> FactorVennIntMax4
                  FactorVennIntMax4 %>% sub("(.*)-(.*)vs(.*)-(.*)","\\1 \\2 \\3 \\4",.) %>% strsplit(" ") %>%
                    lapply(table) %>% lapply(length) %>% "=="(4) %>% "!"(.) %>% which %>% FactorVennIntMax4[.] -> FactorVennIntMax5
                  FactorVennIntMax5 %>% paste0(collapse="|") -> grep
                  ListFileVenn0 %>% basename %>% grep( grep , . , value=F ) %>% ListFileVenn0[.] -> ListFileVenn1
                  # up
                  ListFileVenn1 %>% grep( "_u_",. , value=T ) -> ListFileVenn2
                  if( ListFileVenn2 %>% basename %>% sub(".*_(.*).tsv","\\1",.) %>% as.numeric %>% "=="(.,0) %>% "!"(.) %>% all )
                  {
                    ListFileVenn2 %>% lapply(input) %>% lapply("[",1) %>% lapply(unlist) %>% lapply(as.character) -> ListFileVenn3
                    ListFileVenn2 %>% basename %>% strsplit("_") %>% sapply("[",6) %>% paste0(collapse="x") -> DirNameVenn0
                    venn4( list = ListFileVenn3,
                           listnames = ListFileVenn2 %>% basename %>% strsplit("_") %>% sapply("[",6),
                           title = ListFileVenn2 %>% basename %>% strsplit("_") %>% sapply("[",6) %>% paste0(collapse="x"),
                           export=T, plot=F,
                           dirname = d2[k] %>% basename %>% sub("(.*_.*_.*)_.*","\\1",.) %>% paste( "_u_",DirNameVenn0, sep = ""),
                           path=file.path( Path, DirName0 ) )
                  }
                  # down
                  ListFileVenn1 %>% grep( "_d_",. , value=T ) -> ListFileVenn2
                  if( ListFileVenn2 %>% basename %>% sub(".*_(.*).tsv","\\1",.) %>% as.numeric %>% "=="(.,0) %>% "!"(.) %>% all )
                  {
                    ListFileVenn2 %>% lapply(input) %>% lapply("[",1) %>% lapply(unlist) %>% lapply(as.character) -> ListFileVenn3
                    ListFileVenn2 %>% basename %>% strsplit("_") %>% sapply("[",6) %>% paste0(collapse="x") -> DirNameVenn0
                    venn4( list = ListFileVenn3,
                           listnames = ListFileVenn2 %>% basename %>% strsplit("_") %>% sapply("[",6),
                           title = ListFileVenn2 %>% basename %>% strsplit("_") %>% sapply("[",6) %>% paste0(collapse="x"),
                           export=T, plot=F,
                           dirname = d2[k] %>% basename %>% sub("(.*_.*_.*)_.*","\\1",.) %>% paste( "_d_",DirNameVenn0, sep = ""),
                           path=file.path( Path, DirName0 ) )
                  }

                  # up + down
                  ListFileVenn1 %>% grep( "_ud_",. , value=T ) -> ListFileVenn2
                  if( ListFileVenn2 %>% basename %>% sub(".*_(.*).tsv","\\1",.) %>% as.numeric %>% "=="(.,0) %>% "!"(.) %>% all )
                  {
                    ListFileVenn2 %>% lapply(input) %>% lapply("[",1) %>% lapply(unlist) %>% lapply(as.character) -> ListFileVenn3
                    ListFileVenn2 %>% basename %>% strsplit("_") %>% sapply("[",6) %>% paste0(collapse="x") -> DirNameVenn0
                    venn4( list = ListFileVenn3,
                           listnames = ListFileVenn2 %>% basename %>% strsplit("_") %>% sapply("[",6),
                           title = ListFileVenn2 %>% basename %>% strsplit("_") %>% sapply("[",6) %>% paste0(collapse="x"),
                           export=T, plot=F,
                           dirname = d2[k] %>% basename %>% sub("(.*_.*_.*)_.*","\\1",.) %>% paste( "_ud_",DirNameVenn0, sep = ""),
                           path=file.path( Path, DirName0 ) )
                  }

                }
            }
        }
        # VENN2 levels > 4 without INTERACTION
        if( VENNSUP5 )
        {
          paste( "venn_",colnames(FactorVenn2)[i], sep="" ) -> DirName0
          Path %>% file.path( DirName0 ) %>% dir.create
          d1 %>% basename %>% sub(".*_(.*)","\\1",.) %>% "=="(.,colnames(FactorVenn2)[i]) %>%
            which %>% d1[.] %>% list.files( full.names=T ) -> d2
          foreach( j=1:length(d2) ) %do%
            {
              d2[j] %>% list.files( full.names= T) %>% grep( ".*.tsv$", . ,value=T ) -> ListFileVenn0
              # venn2
              FactorVenn2[,i] %>% levels %>% rev %>% combn( 2, simplify=F ) %>% lapply(paste0,collapse="vs") -> VennComp0
              VennComp0 %>% grep( FactorVenn2[,i] %>% levels %>% "["(1),.,value=T) -> VennComp1
              VennComp1 %>% combn( 2, simplify=F ) -> VennComp2
              foreach( k=1:length(VennComp2) ) %do%
                {
                  VennComp2[k] %>% unlist %>% paste0(collapse="|") -> grep
                  ListFileVenn0 %>% grep( grep, .,value=T ) -> ListFileVenn1
                  # up
                  ListFileVenn1 %>% grep( "_u_",. , value=T ) -> ListFileVenn2
                  ListFileVenn2 %>% lapply(input) %>% lapply("[",1) %>% lapply(unlist) %>% lapply(as.character) -> ListFileVenn3
                  ListFileVenn2 %>% basename %>% strsplit("_") %>% sapply("[",6) %>% paste0(collapse="x") -> DirNameVenn0
                  venn2( list = ListFileVenn3,
                         listnames = ListFileVenn2 %>% basename %>% strsplit("_") %>% sapply("[",6),
                         title = ListFileVenn2 %>% basename %>% strsplit("_") %>% sapply("[",6) %>% paste0(collapse="x"),
                         export=T, plot=F,
                         dirname = d2[j] %>% basename %>% sub("(.*_.*_.*)_.*","\\1",.) %>% paste( "_u_",DirNameVenn0, sep = ""),
                         path=file.path( Path, DirName0 ) )
                  # down
                  ListFileVenn1 %>% grep( "_d_",. , value=T ) -> ListFileVenn2
                  ListFileVenn2 %>% lapply(input) %>% lapply("[",1) %>% lapply(unlist) %>% lapply(as.character) -> ListFileVenn3
                  ListFileVenn2 %>% basename %>% strsplit("_") %>% sapply("[",6) %>% paste0(collapse="x") -> DirNameVenn0
                  venn2( list = ListFileVenn3,
                         listnames = ListFileVenn2 %>% basename %>% strsplit("_") %>% sapply("[",6),
                         title = ListFileVenn2 %>% basename %>% strsplit("_") %>% sapply("[",6) %>% paste0(collapse="x"),
                         export=T, plot=F,
                         dirname = d2[j] %>% basename %>% sub("(.*_.*_.*)_.*","\\1",.) %>% paste( "_d_",DirNameVenn0, sep = ""),
                         path=file.path( Path, DirName0 ) )
                  # up + down
                  ListFileVenn1 %>% grep( "_ud_",. , value=T ) -> ListFileVenn2
                  ListFileVenn2 %>% lapply(input) %>% lapply("[",1) %>% lapply(unlist) %>% lapply(as.character) -> ListFileVenn3
                  ListFileVenn2 %>% basename %>% strsplit("_") %>% sapply("[",6) %>% paste0(collapse="x") -> DirNameVenn0
                  venn2( list = ListFileVenn3,
                         listnames = ListFileVenn2 %>% basename %>% strsplit("_") %>% sapply("[",6),
                         title = ListFileVenn2 %>% basename %>% strsplit("_") %>% sapply("[",6) %>% paste0(collapse="x"),
                         export=T, plot=F,
                         dirname = d2[j] %>% basename %>% sub("(.*_.*_.*)_.*","\\1",.) %>% paste( "_ud_",DirNameVenn0, sep = ""),
                         path=file.path( Path, DirName0 ) )
                }
            }
        }# fi VENN2 > 4 without interaction
    }# rof factors
    # add annotation to list
    Path %>% list.files(full.names = T) %>% grep("venn",.,value=T) -> d0
    d0 %>% list.files( full.names=T, recursive=T ) %>% grep( "List_.*.tsv$" , . , value = T ) -> l0
    if( length(l0) > 1 )
    {
      foreach( i=1:length(l0), .packages=c("magrittr","dplyr","uem915") ) %dopar%
        {
          l0[i] %>% dirname %>% basename %>% strsplit("_") %>% unlist %>% "["(5) %>% gsub("x","\\|",.) -> grep
          r0 %>% dplyr::select( 1, r0 %>% colnames %>% grep(grep,.) ) %>% cbind( Annot0[,-1] ) -> r1
          l0[i] %>% input %>% dplyr::select(1) %>% unlist %>% as.character %>%
            data.frame( "rowID" = . , stringsAsFactors = F ) %>% inner_join( r1  , by = "rowID" ) -> r2
          r2 %>% output( l0[i] )
        }
    }
    # rename lists (for heatmap)
    l0 %>% basename %>% strsplit("_") %>% sapply("[",2) -> temp1
    l0 %>% basename %>% sub("List_.*_(.*).tsv","\\1",.) -> temp2
    l0 %>% dirname %>% basename %>% strsplit("_") %>% lapply("[",1:3) %>% lapply(paste0,collapse="_") %>% unlist -> temp3
    file.path( l0 %>% dirname, paste("List_",temp3,"_",temp1,"_",temp2,".tsv",sep="") ) -> to
    file.rename( from = l0 , to = to )
    paste( "Venn done.\n") %>% cat
  }
  #
  # Cluster
  #
  # cluster preprocessing
  if( cluster )
  {
    statList1 %>% sapply("[[",3) %>% ">"(.,15) %>% which %>% statList1[.] -> statList1
    if(length(statList1) < 1)
    { cluster <- FALSE ; paste("Not enough rows to make cluster analysis \n") %>% cat }
  }
  # cluster processing
  if( cluster )
  {
    paste("Cluster analysis processing...\n",sep="") %>% cat
    Path %>% file.path("cluster") %>% dir.create
    foreach( i=1:length(statList1),.packages=c("magrittr","uem915","dplyr","foreach") ) %dopar%
      {
        statList1[[i]][[1]] %>% strsplit("\\/") %>% unlist %>% length -> Nb0
        statList1[[i]][[1]] %>% strsplit("\\/") %>% unlist %>% "["(Nb0-1) %>% sub("(.*)_.*" , "\\1" , . ) %>%
          paste("_",statList1[[i]][[1]] %>% strsplit("\\/") %>% unlist %>% "["(Nb0) , sep = "" ) -> DirName
        # create list for Symbol and output
        statList1[[i]][[2]] %>%
          data.frame( "rowID" = . , stringsAsFactors = F ) %>%
          inner_join( r0 , by = "rowID" ) -> r1
        # factor for heatmap
        DirName %>% strsplit("_") %>% unlist %>% length -> Nb0
        DirName %>% strsplit("_") %>% unlist %>% "["(Nb0-1) -> ClusterFactor0
        ClusterFactor0 %>% paste("^" , . , "$" , sep = "" ) -> grep
        sif %>% colnames %>% grep( grep , . ) %>% sif[ , . ] -> ClusterFactor1
        # matrix data
        statList1[[i]][[2]] %>%
          data.frame( "rowID" = . , stringsAsFactors = F ) %>%
          inner_join( dat , by = "rowID" ) -> Dat2
        # clustering
        Dat2 %>% hcluster( Symbol=r1$Symbol, factor=ClusterFactor1 , nc=nc,
                           path=Path %>% file.path("cluster"), dirname=DirName )
        # output
        DirName %>% strsplit("_") %>% unlist %>% "["(c(2:4)) %>% paste0(collapse = "_") %>%
          paste( "List_",.,"_cluster_",nrow(r1),".tsv" , sep = "" ) -> FileName
        Path %>% file.path( "cluster" , DirName , FileName ) %>% output( r1 , . )
      }
    # annotations
    Path %>% file.path("cluster") %>% list.files(recursive = T ,full.names = T) -> l0
    l0 %>% basename %>% grep("^List_.*.tsv" , . ) %>% l0[.] -> l1
    foreach(i=1:length(l1),.packages=c("magrittr","uem915","dplyr","foreach")) %dopar%
      {
        l1[i] %>% input -> clm0
        clm0[,1] %>% as.character %>% data.frame(rowID=.) %>% inner_join(r0) %>% output(l1[i])
      }
    paste( "Clusters done.\n" ) %>% cat
  }
  #
  #  Heatmap
  #
  # Heatmap preprocessing
  if( heatmap )
  {
    paste("Heatmap preprocessing...\n",sep = "") %>% cat
    Path %>% list.dirs(recursive = F) -> d0
    d0 %>% basename %>% grep("^anova_|^DiffLists|^pattern|^venn|^cluster",.) %>% d0[.] -> d1
    d1 %>% list.files(recursive=T,full.names=T) -> l0
    l0 %>% basename %>% grep("^List_.*.tsv",.) %>% l0[.] -> l0
    l0 %>% basename %>% sub("^List_.*_(.*).tsv$","\\1",.) %>% as.numeric %>% ">"(.,5) %>% which -> sel
    if( length(sel) > 0 ){ l0[sel] -> l1 }else{ heatmap <- FALSE }
    if( heatmap )
    {
      if( is.null(maxheatmap) ){ maxheatmap <- nrow(dat )}
      l1 %>% basename %>% sub( "^List_.*_(.*).tsv$", "\\1",. ) %>% as.numeric %>% "<"(.,maxheatmap) %>% which -> sel
      if( length(sel)>0 ){ l1[sel] -> l2 }else
        { heatmap <- FALSE ; paste("Not enough rows to make heatmap \n",sep="") %>% cat }
    }
    paste("Heatmap preprocessing done.\n",sep = "") %>% cat
  }
  # heatmap processing
  if( heatmap )
  {
    paste("Heatmap both clustering processing...\n",sep = "") %>% cat
    # both clustering
    foreach( i=1:length(l2), .packages=c("magrittr","uem915","dplyr") ) %dopar%
      {
        l2[i] %>% input -> ll0
        # create factor
        l2[i] %>% basename %>% strsplit("_") %>% unlist %>% "["(.,4) -> FactorHeatmap0
        ll0 %>% colnames %>% grep("^fc_" , . , value = T ) %>% sub("^fc_(.*)" , "\\1" , . ) %>%
          strsplit("vs") %>% unlist %>% unique  -> LevelsFactorHeatmap0
        FactorHeatmap0 %>% paste("^" , . , "$" , sep = "") %>% paste0( collapse = "|") -> grep
        sif %>% colnames %>% grep( grep , . ) %>% sif[ , . ] -> FactorHeatmap1
        LevelsFactorHeatmap0 %>% paste( "^" , . , "$" , sep = "" ) %>% paste0( collapse = "|" ) -> grep
        FactorHeatmap1 %>% levels %>% grep( grep , . , value = T ) -> LevelsFactorHeatmap0
        FactorHeatmap1 %>% grep( grep , . ) -> sel
        FactorHeatmap1[sel]  -> FactorHeatmap2
        FactorHeatmap2 %>% ordered( LevelsFactorHeatmap0 ) -> FactorHeatmap3
        # create matrix data
        ll0[,1] %>% as.character %>% data.frame("rowID"=.,stringsAsFactors=F) %>% inner_join(dat,by="rowID") -> DatHeatmap0
        DatHeatmap0 %>% dplyr::select(all_of(sel+1)) -> DatHeatmap1
        heatmapsd <- TRUE
        DatHeatmap1 %>% apply(1,sd) %>% "=="(.,0) %>% which -> sel
        if( length(sel) > 0){ heatmapsd <- FALSE }
        DatHeatmap1 %>% apply(2,sd) %>% "=="(.,0) %>% which -> sel
        if( length(sel) > 0){ heatmapsd <- FALSE }
        if( heatmapsd )
        {
          # Symbol list
          ll0[,1] %>% as.character %>% data.frame("rowID"=.,stringsAsFactors=F) %>%
            inner_join(Annot0,by="rowID") %>% dplyr::select(Symbol) %>% unlist %>% as.character -> Symbol
          l2[i] %>% basename %>% sub("(.*).tsv","\\1",.) %>% paste("Heatmap_",.,".pdf",sep="") -> FileName
          l2[i] %>% dirname %>% file.path(FileName) -> FileName
          pdf(FileName)
          DatHeatmap1 %>% heatmap(factor=FactorHeatmap3,labRow=Symbol,cexCol = 0.45, labCol=colnames(DatHeatmap1) )
          graphics.off()
        }
      }
    paste("Heatmap both clustering done.\n",sep = "") %>% cat
    paste("Heatmap row clustering processing...\n",sep = "") %>% cat
    # row clustering
    foreach( i=1:length(l2), .packages=c("magrittr","uem915","dplyr") ) %dopar%
      {
        # create factor
        l2[i] %>% basename %>% strsplit("_") %>% unlist %>% "["(.,4) -> FactorHeatmap0
        l2[i] %>% input -> ll0
        ll0 %>% colnames %>% grep("^fc_" , . , value = T ) %>% sub("^fc_(.*)" , "\\1" , . ) %>%
          strsplit("vs") %>% unlist %>% unique  -> LevelsFactorHeatmap0
        FactorHeatmap0 %>% paste("^" , . , "$" , sep = "") %>% paste0( collapse = "|") -> grep
        sif %>% colnames %>% grep( grep , . ) %>% sif[ , . ] -> FactorHeatmap1
        LevelsFactorHeatmap0 %>% paste( "^" , . , "$" , sep = "" ) %>% paste0( collapse = "|" ) -> grep
        FactorHeatmap1 %>% levels %>% grep( grep , . , value = T ) -> LevelsFactorHeatmap0
        FactorHeatmap1 %>% grep( grep , . ) -> sel
        FactorHeatmap1[sel]  -> FactorHeatmap2
        FactorHeatmap2 %>% ordered( LevelsFactorHeatmap0 ) -> FactorHeatmap3
        # create matrix data
        ll0[,1] %>% as.character %>% data.frame("rowID"=., stringsAsFactors=F) %>% inner_join(dat, by="rowID") -> DatHeatmap0
        DatHeatmap0 %>% dplyr::select(., all_of(sel+1)) -> DatHeatmap1
        heatmapsd <- TRUE
        DatHeatmap1 %>% apply(1,sd) %>% "=="(.,0) %>% which -> sel
        if( length(sel) > 0){ heatmapsd <- FALSE}
        DatHeatmap1 %>% apply(2,sd) %>% "=="(.,0) %>% which -> sel
        if( heatmapsd )
        {
          # Symbol list
          ll0[,1] %>% as.character %>% data.frame( "rowID" = . , stringsAsFactors = F ) %>%
            inner_join( Annot0 , by = "rowID" ) %>% dplyr::select(Symbol) %>% unlist %>% as.character -> Symbol
          l2[i] %>% basename %>% sub( "(.*).tsv", "\\1" , .  ) %>% paste("Heatmap_row_" , . , ".pdf" , sep = "") -> FileName
          l2[i] %>% dirname %>% file.path( FileName ) -> FileName
          pdf( FileName )
          DatHeatmap1 %>% heatmap( factor=FactorHeatmap3, labRow=Symbol,
                                   labCol=colnames(DatHeatmap1), cexCol = 0.45, dendrogram="row" )
          graphics.off()
        }
      }
    paste("Heatmap done.\n",sep="") %>% cat
  }
  #
  # volcanoplot
  #
  if( volcanoplot )
  {
    paste("Volcanoplot processing...\n") %>% cat
    if( is.character(threshold) ){1:length(thresholdlist) -> threshold }
    thresholdlist[threshold] -> Threshold0
    Path %>% file.path("volcanoplot") %>% dir.create
    r0 %>% colnames %>% grep("^p_.*vs.*$",.,value=T) %>% sub("p_(.*vs.*)","\\1",.) -> CompVolcano
    foreach(i=1:length(CompVolcano), .packages=c("magrittr","dplyr","uem915","ggplot2")) %dopar%
      {
        r0 %>% dplyr::select(dplyr::contains(CompVolcano[i])) -> DatVolcanoplot0
        data.frame("rowID"=r0[,1] %>% as.character,DatVolcanoplot0,"Symbol"=Annot0$Symbol,stringsAsFactors=F) -> DatVolcanoplot1
        foreach(j=1:length(Threshold0), .packages=c("magrittr","dplyr","uem915")) %do%
          {
            DatVolcanoplot1 %>% volcanoplot(pval=Threshold0[[j]][[1]], fc=Threshold0[[j]][[2]], GeneName=T) -> Vplot
            which(
              DatVolcanoplot1[,2] < Threshold0[[j]][[1]] &
              abs(DatVolcanoplot1[,3]) > Threshold0[[j]][[2]] ) %>% length -> NbGeneDiff
            paste( "Volcanoplot_p" , sub( "0\\.(.*)" , "\\1" , Threshold0[[j]][[1]] ),"_","fc",
                   gsub( "\\." , "" , Threshold0[[j]][[2]] ),"_", CompVolcano[i] , "_" , NbGeneDiff, ".pdf" , sep = "" ) -> FileName
            Path %>% file.path("volcanoplot", FileName) -> FileName
            ggsave(plot=Vplot, filename=FileName)
          }
      }
    paste("Volcanoplot done.\n") %>% cat
  }
  #
  # lineplot
  #
  if( lineplot )
  {
    paste( "Line plot processing...\n" , sep = "" ) %>% cat
    Path %>% list.dirs(recursive = F) -> d0
    d0 %>% basename %>% grep("^DiffLists_" , . ) %>% d0[.] -> d1
    d1 %>% list.files( recursive = T , full.names = T ) -> l0
    l0 %>% basename %>% grep("^List_.*.tsv" , . , value = F ) %>% l0[.] -> l1
    l1 %>% basename %>% sub( "List_.*_(.*).tsv" , "\\1" , .  ) %>%
      as.numeric %>% ">"(.,0) %>% which %>% l1[.] -> l2
    if( length(l2)==0 ){ lineplot <- FALSE ; paste("No significant feature for lineplot done.\n",sep="") %>% cat }
    if( lineplot )
    {
      ll1 <- foreach( i=1:length(l2), .packages=c("magrittr","dplyr") ) %dopar%
        { l2[i] %>% input %>% dplyr::select(1) %>% unlist %>% unique }
      ll1 %>% unlist %>% as.character %>% unique -> ll2
      # only for symbol in line plot
      ll2 %>% data.frame("rowID"=., stringsAsFactors=FALSE) %>% inner_join(r0, by="rowID") -> r1
      ll2 %>% data.frame("rowID"=., stringsAsFactors=FALSE) %>% inner_join(dat, by="rowID") -> Dat2
      paste("lineplot_",nrow(Dat2),sep="") -> DirName
      Path %>% file.path(DirName) %>% dir.create
      foreach( i=1:length(CompFactors) ) %do%
        {
          Path %>% file.path( DirName , CompFactors[i] ) %>% dir.create
          CompFactors[i] -> CompFactors0
          if( CompFactors[i] %>% grepl("x",.) ){ IntFactor1 -> CompFactors0 }
          foreach( j=1:nrow(Dat2), .packages = c("uem915","magrittr","dplyr","broom","ggplot2") ) %dopar%
            {
              r1$Symbol[j] %>% gsub("/"," ",.) %>% gsub("\'","",.) %>% gsub("-","",.) -> Symbol
              paste("lineplot_",CompFactors[i],"_",Symbol,".pdf",sep="") -> FileName
              Path %>% file.path(DirName, CompFactors[i], FileName) -> FileName
              CompFactors0 %>% paste("^",.,"$",sep="") %>% paste0(collapse="|") -> grep
              sif %>% colnames %>% grep(grep, ., value=FALSE) %>% sif[,.] -> Sif0
              Sif0 %>% data.frame( y=Dat2[j,-1] %>% as.numeric, . ) -> Sif1
              dplot( dat=Sif1, title=r1$Symbol[j]  ) -> p0
              ggsave( plot=p0, filename=FileName )
            }
        }
      paste("Lineplot done.\n",sep="") %>% cat
    }
  }
  #
  # boxplot
  #
  if( boxplotrow )
  {
    paste( "boxplot processing...\n" , sep = "" ) %>% cat
    Path %>% list.dirs(recursive = F) -> d0
    d0 %>% basename %>% grep("^DiffLists_" , . ) %>% d0[.] -> d1
    d1 %>% list.files( recursive = T , full.names = T ) -> l0
    l0 %>% basename %>% grep("^List_.*.tsv" , . , value = F ) %>% l0[.] -> l1
    l1 %>% basename %>% sub( "List_.*_(.*).tsv" , "\\1" , .  ) %>%
      as.numeric %>% ">"(.,0) %>% which %>% l1[.] -> l2
    if( length(l2)==0 ){ boxplotrow <- FALSE ; paste("No significant feature for boxplot done.\n",sep="") %>% cat }
    if( boxplotrow )
    {
      ll1 <- foreach( i=1:length(l2), .packages=c("magrittr","dplyr") ) %dopar%
        { l2[i] %>% input %>% dplyr::select(1) %>% unlist %>% unique }
      ll1 %>% unlist %>% as.character %>% unique -> ll2
      # only for symbol in line plot
      ll2 %>% data.frame("rowID"=., stringsAsFactors=FALSE) %>% inner_join(r0, by="rowID") -> r1
      ll2 %>% data.frame("rowID"=., stringsAsFactors=FALSE) %>% inner_join(dat, by="rowID") -> Dat2
      paste("boxplot_",nrow(Dat2),sep="") -> DirName
      Path %>% file.path(DirName) %>% dir.create
      foreach( i=1:length(CompFactors) ) %do%
        {
          Path %>% file.path( DirName , CompFactors[i] ) %>% dir.create
          CompFactors[i] -> CompFactors0
          # if( CompFactors[i] %>% grepl("x",.) ){ IntFactor1 -> CompFactors0 }
          foreach( j=1:nrow(Dat2), .packages = c("uem915","magrittr","dplyr","broom","ggplot2","ggpubr") ) %dopar%
            {
              r1$Symbol[j] %>% gsub("/"," ",.) %>% gsub("\'","",.) %>% gsub("-","",.) -> Symbol
              paste("boxplot_",CompFactors[i],"_",Symbol,".pdf",sep="") -> FileName
              Path %>% file.path(DirName, CompFactors[i], FileName) -> FileName
              CompFactors0 %>% paste("^",.,"$",sep="") %>% paste0(collapse="|") -> grep
              sif %>% colnames %>% grep(grep, ., value=FALSE) %>% sif[,.] -> Sif0
              Sif0 %>% data.frame( y=Dat2[j,-1] %>% as.numeric, . ) -> Sif1
              boxplot1( dat=Sif1 , xlab =CompFactors[i], ylab = r1$Symbol[j] ) -> p0
              ggsave( plot=p0, filename=FileName )
            }
        }
      paste("boxplot done.\n",sep="") %>% cat
    }
  }
  #
  # Functional analysis :
  #
  # enrichment analysis
  #
  if( ena  )
  {
    paste("#\n") %>% cat
    paste("Functional analysis :\n") %>% cat
    paste("ena preprocessing...\n") %>% cat
    #
    Path %>% list.dirs(recursive = F) -> d0
    # Lists
    d0 %>% basename %>% grep("^DiffLists_|^pattern|^venn|^cluster", .) %>% d0[.] -> d1
    d1 %>% lapply(list.files,recursive=T, full.names=T ) %>% unlist -> l0
    l0 %>% basename %>% grep("^List_p.*_fc", .) %>% l0[.] -> l1
    l1 %>% basename %>% sub("^List_.*_(.*).tsv$", "\\1", .) %>% as.numeric %>% ">="(.,enamin) %>% which -> sel
    if( length(sel) > 0 ){ l1[sel]-> l1 }
    l1 %>% basename %>% sub("^List_.*_(.*).tsv$", "\\1", .) %>% as.numeric %>% "<"(.,3000) %>% which -> sel
    if( length(sel) > 0 ){ l1[sel]-> l1 }
    if( length(l1) > 0 ){ ena <- TRUE }else{ ena <- FALSE }
    if( ena )
    {
      ll0 <- foreach(i=1:length(l1), .packages=c("uem915","magrittr","dplyr")) %dopar%
        { l1[i] %>% input %>% dplyr::select(.data$Symbol) %>% unlist %>% as.character }
      ll0 %>% unlist %>% unique -> SymbolListEna
      if( species != "hs" )
      {
        SymbolListEna %>% annot(species=species, ortholog=T) -> AnnotEnaOrtho
        nrow(AnnotEnaSpecies)+nrow(AnnotEnaOrtho) -> NbMatch
      }else
        {
          SymbolListEna %>% annot(species="hs", ortholog=F) -> AnnotEnaHs
          nrow(AnnotEnaHs) -> NbMatch
        }
      if( NbMatch < 3 )
      {
        paste("ena preprocessing done \n") %>% cat
        paste("No match symbol \n") %>% cat
        ena <- FALSE
      }else{ paste("ena preprocessing done \n") %>% cat }
    }
  }
  #
  # ena processing
  #
  if( ena )
  {
    paste("Enrichment analysis processing for ",length(l1)," lists.\n",sep = "") %>% cat
    Path %>% file.path("ena") %>% dir.create
    l1 %>% lapply(strsplit,Path) %>% lapply(unlist) %>% sapply("[",2) %>%
      dirname %>% file.path(Path , "ena" , . ) -> PathDirNameEna
    PathDirNameEna %>% unique %>% lapply(dir.create,recursive=T)
    if( species!="hs" )
    {
      ll1 <- foreach(i=1:length(ll0), .packages=c("uem915","magrittr","dplyr")) %dopar%
        {
          ll0[[i]] %>% as.character %>% data.frame(Symbol=., stringsAsFactors=FALSE) %>% inner_join(AnnotEnaSpecies) %>%
            inner_join(AnnotEnaOrtho, c("GeneID"="OtherGeneID")) %>% dplyr::select(.data$GeneID.y,.data$Symbol.y) %>%
            setNames(c("GeneID","Symbol"))
        }
    }else
      {
        ll1 <- foreach(i=1:length(ll0), .packages=c("magrittr","dplyr")) %dopar%
          { ll0[[i]] %>% data.frame(Symbol=., stringsAsFactors=FALSE) %>% inner_join(AnnotEnaHs) %>%
              dplyr::select(.data$GeneID,.data$Symbol) }
      }
    foreach(i=1:length(ll1),.packages=c("magrittr","dplyr","uem915")) %dopar%
      {
        l1[i] %>% basename %>% sub("(.*_.*).tsv", "\\1", .) -> DirName
        ena( geneannot=ll1[[i]], filtergeneset=filtergeneset, bg=bg,
             overlapmin=enamin, enaScoremin=enaScoremin, path=PathDirNameEna[i], dirname=DirName )
      }
    paste( "Enrichment analysis done.\n" ) %>% cat
  }
  #
  # network
  #
  #
  # ena summary DiffList
  #
  #
  # ena summary clustering
  #
  #
  # dopar
  #
  parallel::stopCluster(cl) ; doParallel::stopImplicitCluster()
  #
  # log
  #
  LogOutput <- list()
  paste("##############################\n",sep="") %>%
  paste("####### omic analysis  #######\n",sep="") %>%
  paste("##############################\n",sep="") %>% c(LogOutput,.) -> LogOutput
  #
  paste("Factor(s) : ",AnovaFactors,sep="") %>% c(LogOutput,.) -> LogOutput
  paste("anova model : Y = mu + ",model," + E",sep="") %>% c(LogOutput,.) -> LogOutput
  paste(ModelName,sep="") %>% c(LogOutput,.) -> LogOutput
  paste("Remove batch effect : ",ifelse(is.null(batch), "NO", "YES"),sep="") %>% c(LogOutput,.) -> LogOutput
  paste("adjusted p-value : ",padj,sep="") %>% c(LogOutput,.) -> LogOutput
  if( !is.null(threshold) ){
    paste("Thresholds used (pval|fold-change) :",
          Threshold0 %>% lapply("[",c(1,2)) %>% lapply(unlist) %>%
            lapply(paste0,collapse="|") %>% unlist %>% paste0(collapse = ",") ) %>% c(LogOutput,.) -> LogOutput }
  list(
    list( pattern, "pattern done"),
    list( T , "anovastat done"),
    list( venn , "Venn done"),
    list( cluster , paste("Cluster analysis done with cut value = ",nc,sep="")),
    list( heatmap , "Heatmap done"),
    list( volcanoplot , "Volcanoplot done"),
    list( lineplot , "Lineplot done")) -> par0
  par0 %>% sapply("[[", 1) %>% which %>% par0[.] %>% lapply("[[",2) %>% c(LogOutput,.) -> LogOutput
  Path %>% list.files(recursive=TRUE) %>% basename %>% grep("^List_.*.tsv$", .) %>% length -> NbList
  paste("Number of list diff created by all analysis : ",NbList,sep="") %>% c(LogOutput,.) -> LogOutput
  if(ena)
  {
    paste("# Functionnal analysis : ",sep="") %>% c(LogOutput,.) -> LogOutput
    paste("Enrichment analysis done.",sep="") %>% c(LogOutput,.) -> LogOutput
    paste("species : ",species,sep="") %>% c(LogOutput,.) -> LogOutput
    paste("Background used : ",bg,sep="") %>% c(LogOutput,.) -> LogOutput
    paste("Filter geneset : ",filtergeneset,sep="") %>% c(LogOutput,.) -> LogOutput
  }
  paste("Date : ",format(Sys.time(), "%a %b %d %X %Y"),sep = "") %>% c(LogOutput,.) -> LogOutput
  capture.output(Sys.time() - start) %>% unlist %>% strsplit(" ") %>% unlist %>% "["(.,-c(1:3)) %>% paste0(collapse=" ") -> Time
  paste("analysis time : " , Time, sep = "") %>% c(LogOutput,.) -> LogOutput
  paste("version : ",packageVersion("uem915") %>% as.character, sep = "") %>% c(LogOutput,.) -> LogOutput
  #
  Path %>% file.path("log.txt") %>% file("w") -> f
  foreach(i=1:length(LogOutput)) %do% { LogOutput[[i]] %>% paste(.,"\n",sep="") %>% writeLines(f) }
  close(f)
  # zip
  if( zip )
  {
    utils::zip( zipfile=Path %>% basename, files=Path %>% basename )
    paste( "zip done.\n" ) %>% cat
  }
  if( zip & remove ){ base::unlink( Path, recursive=T ) ; paste( "remove done.\n" ) %>% cat }
  paste( "###########################\nOmic analysis done \n###########################\n",sep="") %>% cat
  # Rplots.pdf
  if( file.exists( file.path(path,"Rplots.pdf") ) ){ file.remove( file.path(path,"Rplots.pdf") ) }
}