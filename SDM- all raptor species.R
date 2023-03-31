## biomod2 video: Multi species modelling ----

## load the required packages
library(biomod2)
library(raster)
library(rasterVis)
library(gridExtra)
library(reshape2)

## Obtaining species data from GBIF

setwd("E:/species/")
mysp <- read.csv("data/Finalraptdata.csv")
head (myspecies)
myspecies <- mysp%>%
  select(species, Longitude, Latitude)%>%
  group_by(myspecies, species)
  
unique(myspecies$species)

#Loading predictor variables
lst <- list.files(path="data_kenya/",pattern='asc$',all.files = TRUE, full.names = T) 
preds<-stack(lst)

plot(preds[[1]])

#####Select the variables to use#######################################################################
spg<- myspecies%>%select(Longitude,Latitude)
head(spg)
spg$species<-1
spg<- spg%>%drop_na()
nrow(spg)
coordinates(spg)<-c('Longitude', 'Latitude')
class(spg)
head(spg)
plot(spg,add=T)

ex<-raster ::extract (preds,spg)
head(ex)

v<-usdm::vifstep(ex)
preds<-usdm::exclude(preds, v)

#######################################################################################################
## curent climatic variables
stk_current <- 
  raster::stack(
    c(preds),
    RAT = FALSE
  )
plot(stk_current)

###################################FUTURE CLIMATIC VARIABLES##################################
## 2050 climatic variables
##stk_2050_BC_45 <- 
  #raster::stack(
   # c(preds),
    #RAT = FALSE
  #)

## 2070 climatic variables
##stk_2070_BC_45 <- 
  #raster::stack(
    #c(preds),
    #RAT = FALSE
  #)
###################################################################################################

## build species modelling wrapper ----
biomod2_wrapper <- function(sp){
  cat("\n> species : ", sp)
  
  ## get occurrences points
  sp_dat <- data[data$species == sp, ]
  
  ## formating the data
  sp_format <- 
    BIOMOD_FormatingData(
      resp.var = rep(1, nrow(sp_dat)), 
      expl.var = stk_current,
      resp.xy = sp_dat[, c("long", "lat")],
      resp.name = sp,
      PA.strategy = "random", 
      PA.nb.rep = 1, 
      PA.nb.absences = 10000
    )
  
  ## print formatting summary
  sp_format
  
  ## save image of input data summary
  if(!exists(sp)) dir.create(sp)
  pdf(paste(sp, "/", sp ,"_data_formated.pdf", sep="" ))
  try(plot(sp_format))
  dev.off()
  
  ## define models options
  sp_opt <- BIOMOD_ModelingOptions(GLM = list(type = 'quadratic', interaction.level = 1),
                                   GBM = list(n.trees = 5000),
                                   GAM = list(algo = 'GAM_mgcv'),
                                   MARS = list(interaction.level =2),
                                   RF = list(ntree =750),
                                   MAXENT.Phillips = list(maximumiterations = 1000))
  BIOMOD_ModelingOptions() 
  ## model species
  sp_model <- BIOMOD_Modeling( 
    sp_format, 
    models = c("GLM", "GBM", "RF", "GAM", "ANN", "MARS", "MAXENT.Phillips.2"), 
    models.options = sp_opt, 
    NbRunEval = 5, 
    DataSplit = 70, 
    Yweights = NULL, 
    VarImport = 3, 
    models.eval.meth = c('TSS', 'ROC'),
    SaveObj = TRUE,
    rescal.all.models = FALSE,
    do.full.models = TRUE,
    modeling.id = "demo1"
  )
  
  ## save some graphical outputs
  #### models scores
  pdf(paste0(sp, "/", sp , "_models_scores.pdf"))
  try(gg1 <- models_scores_graph(sp_model, metrics = c("TSS", "ROC"), by = 'models', plot = TRUE))
  try(gg2 <- models_scores_graph(sp_model, metrics = c("TSS", "ROC"), by = 'data_set', plot = TRUE))
  try(gg3 <- models_scores_graph(sp_model, metrics = c("TSS", "ROC"), by = 'cv_run', plot = TRUE))
  try(grid.arrange(gg1, gg2, gg3))
  dev.off()
  
  ## build ensemble models
  sp_ens_model <- 
    BIOMOD_EnsembleModeling(
      modeling.output = sp_model,
      em.by = 'all',
      eval.metric = 'TSS',
      eval.metric.quality.threshold = 0.7,
      models.eval.meth = c('TSS','ROC'),
      prob.mean = FALSE,
      prob.mean.weight = TRUE,
      VarImport = 0
    )
  
  ## do projections
  proj_scen <- ("current")
  
  for(scen in proj_scen){
    cat("\n> projections of ", scen)
    
    ## single model projections
    sp_proj <- 
      BIOMOD_Projection(
        modeling.output = sp_model,
        new.env = get(paste0("stk_", scen)),
        proj.name = scen,
        selected.models = 'all',
        binary.meth = "TSS",
        filtered.meth = NULL,
        compress = TRUE,
        build.clamping.mask = FALSE,
        do.stack = FALSE,
        output.format = ".img"
      )
    
    ## ensemble model projections
    sp_ens_proj <- 
      BIOMOD_EnsembleForecasting(
        EM.output = sp_ens_model,
        projection.output = sp_proj,
        binary.meth = "TSS",
        compress = TRUE,
        do.stack = FALSE,
        output.format = ".img"
      )
  }
  
  return(paste0(sp," modelling completed !"))
}


## launch the spiecies modelling wrapper over species list ----
if(require(snowfall)){ ## parallel computation
  ## start the cluster
  sfInit(parallel = TRUE, cpus = 5) ## here we only require 4 cpus
  sfExportAll()
  sfLibrary(biomod2)
  ## launch our wrapper in parallel
  sf_out <- sfLapply(spp_to_model, biomod2_wrapper)
  ## stop the cluster
  sfStop()
} else { ## sequencial computation
  for (sp in spp_to_model){
    biomod2_wrapper(sp)
  }
  ## or with a lapply function in sequential model
  ## all_species_bm <- lapply(spp_to_model, biomod2_wrapper)
}

