# #########################################################################
# DESCRIPTION:
#   Convert cifti (1) functional time series and (2) parcellation file into 
#   correlation matrix
#
# Inputs: 
#           func_file,          Functional time series in CIFTI (.cii)
#
#           parcel_file,        Pracellation file in CIFTI (.nii)
#
#           tmask,              Tmask indicating whether volumes should be retained.
#                               If a tmask is not included, all volumes will be retained.
#
#           surface_only,       A boolean on whether only cortical surface data is needed.
#                               If the parcellation file does not contain subcortical data, 
#                               select TRUE (default=FALSE)
#
# Ouput:    
#           m,                  A parcel x parcel correlation matrix
#
# #########################################################################
# myc, UTD 9/25/2019 - Initial
# #########################################################################
### **** Need to test with HCP, MSC and in house data! - 9/25/2019

cii2mat <- function(func_file, parcel_file, tmask=NULL, surface_only=FALSE){
  library(cifti)
  
  # Load CIFTI files
  cii <- read_cifti(func_file, drop_data = FALSE, trans_data = T) 
  parcel <- as.matrix(read_cifti(parcel_file)$data)
  
  if(is.null(tmask)){
    tmask <- matrix(data = 1, nrow = dim(cii$data)[,2])
  }else{
    if(is.character(tmask)){
      tmask <- read.table(tmask)$V1
    }else if(is.logical(tmask)){
      tmask <- tmask
    }else{
      stop("tmask variable not character (file name) or logical (vector of 1 or 0)")
    }
  }
  
  u_parcel <- unique(parcel)
  u_parcel <- u_parcel[u_parcel!=0] # Remove parcel 0 and order parcel by their number
  
  # == Make brainstructure index ==
  # The brainstrucure index lets us filter out anatomical structures based 
  # on an index (this mirrors the cifti packages in MATLAB). As of June 2019, 
  # the cifti package on CRAN would ignore some subcortical labels, so make 
  # sure to use the development version on github
  cii$brainstructureindex <- as.matrix(NA, dim(cii$data)[1])
  for(i in 1:length(cii$BrainModel)){
    startindx <- attributes(cii$BrainModel[[i]])$IndexOffset + 1
    endindx <- attributes(cii$BrainModel[[i]])$IndexOffset + 
      attributes(cii$BrainModel[[i]])$IndexCount
    
    cii$brainstructureindex[startindx:endindx] <- i
  }
  
  # Check Dimension of data
  d_cii <- dim(cii$data) # ~ 64k vertices, includes subcortical volumes
  d_parcel <- dim(parcel)   # surface only, excluded medial wall
  
  # Exclude subcortical data if surface_only==TRUE
  if(surface_only==TRUE){
    cdata <- as.matrix(cii$data[cii$brainstructureindex==1 | cii$brainstructureindex==2,,])
  }else{
    cdata <- cii$data
  }
    
  
  ctmask <- cdata[,as.logical(tmask)]
  
  # Create parcel x mean time series matrix
  tp <- matrix(0, length(u_parcel), sum(tmask)) # initialize empty matrix
  
  for(i in 1:length(u_parcel)){               
    tp[i,]<- colMeans(ctmask[which(parcel==u_parcel[i]),])
  }
  
  tp <- tp[order(u_parcel),] # Order matrix by parcel number
  
  r <- cor(t(tp))         # Correlation matrix between all nodes
  z <- psych::fisherz(r)  # Fisher's z-transform: 0.5 * log((1+r)/(1-r))
  diag(z) <- 0       
  
  return(r, z)
}