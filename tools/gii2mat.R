# #########################################################################
# DESCRIPTION:
#   Convert GIFTI (1) functional time series and (2) parcellation file into 
#   correlation matrix
#
# Inputs: 
#           gl_file | gr_file,        Functional time series in CIFTI (.cii)
#
#           roi_l_file | roi_r_file, Pracellation file in CIFTI (.nii)
#
#           tmask,              Tmask indicating whether volumes should be retained.
#                               If a tmask is not included, all volumes will be retained.
#

#
# Ouput:    
#           r,                  A parcel x parcel correlation matrix in r
#           z,                  A parcel x parcel correlation matrix in fisher's transformed z
# #########################################################################
gii2mat <- function(gl_file, gr_file, roi_l_file, roi_r_file, tmask=NULL){

library(gifti)
  
##################
## Read in files #
##################
gL <- read_gifti(gL_file)
gR <- read_gifti(gR_file)

# nodes (L/R): Chan_RSFC_Nodes/gifti_multiple_columns/ROI_L/R_dis8_fwhm0_limit3_overlapEXCLUDE.func.gii
# metadata: Chan_RSFC_Nodes/Chan_RSFC_Nodes_PNAS2014_metadata.txt
node_L <- read_gifti(nodeL_file)    
node_R <- read_gifti(nodeR_file)    
node_order <- read.table(node_meta_table_file, sep="\t", header=T)  


# Load motion/tmask files
fd <- read.table(fd_file, col.names = "FD")
fd$vol <- 1:nrow(fd)

motion <- read.table(motion_file, col.names = c("X","Y","Z","P","Y","R"))
dv <- read.table(dv_file, col.names = "DVARS")
dv$vol <- 1:nrow(dv)

tmask <- read.table(tmask_file, col.names = "tmask")

# Tmask the motion files
fd <- data.frame(fd[as.logical(tmask$tmask),])
motion <- motion[as.logical(tmask$tmask),]
dv <- dv[as.logical(tmask$tmask),]

######################################################
## Extract Nodes' mean time series from surface data #
######################################################
gL <- as.matrix(data.frame(gL$data))
gR <- as.matrix(data.frame(gR$data))
node_L <- as.matrix(data.frame(node_L$data))
node_R <- as.matrix(data.frame(node_R$data))

if(ncol(gL)!=ncol(gR)){
  stop("Column size (# volumes) of left & right hemisphere should be equal. Check input data.")
}

tp_L <- matrix(0, ncol(node_L), ncol(gL))
tp_R <- matrix(0, ncol(node_R), ncol(gR))

for(i in 1:ncol(node_L)){                 # Left Hemipshere
  tp_L[i,]<- colMeans(gL[node_L[,i]==1,])
}

for(i in 1:ncol(node_R)){                 # Right Hemipshere
  tp_R[i,]<- colMeans(gR[node_R[,i]==1,])
}

tp <- rbind(tp_L, tp_R)                   # combine L and R
rm(tp_L, tp_R) # cleanup

r <- cor(t(tp))         # Correlation matrix between all nodes
z <- psych::fisherz(r)  # Fisher's z-transform: 0.5 * log((1+r)/(1-r))


}