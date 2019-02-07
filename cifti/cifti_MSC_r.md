CIFTI in R, MSC single-subject
================
Micalea Chan
2/5/2019

## Read in cifti files

  - MSC-01’s data are used
      - Individual specific parcellation and community (network) are
        loaded.

<!-- end list -->

``` r
cii <- read_cifti(cii_file, drop_data = FALSE, trans_data = T)         # Could read in just the data using cifti_data(cii_file)

# Make brainstructure index
cii$brainstructureindex <- as.matrix(NA, dim(cii$data)[1])
for(i in 1:length(cii$BrainModel)){
  startindx <- as.numeric(attributes(cii$BrainModel[[i]])$IndexOffset + 1)
  endindx <- as.numeric(attributes(cii$BrainModel[[i]])$IndexOffset+attributes(cii$BrainModel[[i]])$IndexCount)
  
  cii$brainstructureindex[startindx:endindx] <- i
  
}


parcel <- read_cifti(parcel_file)
parcel <- as.matrix(parcel$data)

comm <- read_cifti(comm_file)
comm <- as.matrix(comm$data)
  
# Check dimension of cifti data (volume/frame x vertices)
dim(cii$data) # ~ 64k vertices, includes subcortical volumes
```

    ## [1] 65890   818     1

``` r
dim(parcel)   # surface only, excluded medial wall
```

    ## [1] 59412     1

``` r
dim(comm)     # surface only, excluded medial wall
```

    ## [1] 59412     1

``` r
# What are the labeled brain structures in the cii file? 
cifti_brain_structs(cii)
```

    ##  [1] "CIFTI_STRUCTURE_CORTEX_LEFT"      
    ##  [2] "CIFTI_STRUCTURE_CORTEX_RIGHT"     
    ##  [3] "CIFTI_STRUCTURE_ACCUMBENS_LEFT"   
    ##  [4] "CIFTI_STRUCTURE_ACCUMBENS_RIGHT"  
    ##  [5] "CIFTI_STRUCTURE_AMYGDALA_LEFT"    
    ##  [6] "CIFTI_STRUCTURE_AMYGDALA_RIGHT"   
    ##  [7] "CIFTI_STRUCTURE_CAUDATE_LEFT"     
    ##  [8] "CIFTI_STRUCTURE_CAUDATE_RIGHT"    
    ##  [9] "CIFTI_STRUCTURE_CEREBELLUM_LEFT"  
    ## [10] "CIFTI_STRUCTURE_CEREBELLUM_RIGHT" 
    ## [11] "CIFTI_STRUCTURE_HIPPOCAMPUS_LEFT" 
    ## [12] "CIFTI_STRUCTURE_HIPPOCAMPUS_RIGHT"
    ## [13] "CIFTI_STRUCTURE_PALLIDUM_LEFT"    
    ## [14] "CIFTI_STRUCTURE_PALLIDUM_RIGHT"   
    ## [15] "CIFTI_STRUCTURE_PUTAMEN_LEFT"     
    ## [16] "CIFTI_STRUCTURE_PUTAMEN_RIGHT"    
    ## [17] "CIFTI_STRUCTURE_THALAMUS_LEFT"    
    ## [18] "CIFTI_STRUCTURE_THALAMUS_RIGHT"

``` r
# ==== Take the Left & Right cortex only
cdata <- as.matrix(cii$data[cii$brainstructureindex==1 | cii$brainstructureindex==2,,])
# cdata <- cii$data[,,]
dim(cdata)
```

    ## [1] 59412   818

``` r
u_parcel <- unique(parcel)
u_parcel <- u_parcel[u_parcel!=0] # Remove parcel 0 and order parcel by #
```

## Extract Nodes’ mean time series from surface data

``` r
# ==== Mask out bad volumes from data
tmask <- read.table(tmask_file)$V1
ctmask <- cdata[,as.logical(tmask)]

# ==== Extract mean time series from each parcel -> Parcel x Vol matrix
tp <- matrix(0, length(u_parcel), sum(tmask))   

for(i in 1:length(u_parcel)){               
  tp[i,]<- colMeans(ctmask[which(parcel==u_parcel[i]),])
}
```

## Plot processed mean time series of each node

  - The heatmaps here are generated using a customized version of the
    [superheat (github)](https://github.com/mychan24/superheat) package.

<!-- end list -->

``` r
superheat::superheat(tp,
                     heat.lim = c(-20, 20), 
                     heat.pal = c("black","white"),
                     grid.hline = FALSE,
                     grid.vline = FALSE,
                     title="Mean Time series of each parcel")
```

![](cifti_MSC_r_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

## Correlation Matrix (z-transformed)

``` r
r <- cor(t(tp))         # Correlation matrix between all nodes
z <- psych::fisherz(r)  # Fisher's z-transform: 0.5 * log((1+r)/(1-r))

diag(z) <- 0            # Set diagonal to '0'; not informative

superheat::superheat(z, 
                     y.axis.reverse = TRUE, # Used to make origin (0,0) on top left corner
                     heat.lim = c(-.1, .5), 
                     heat.pal = parula(100),# rev(brewer.rdylbu(100)), 
                     heat.pal.values = c(0, 0.5,1),
                     grid.hline = FALSE,
                     grid.vline = FALSE,
                     title="Node x Node Correlation Matrix (z-transformed)")
```

![](cifti_MSC_r_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

## Correlation Matrix, nodes ordered by systems

### Setup System Color for Plot

``` r
parlabel <- data.frame(parcel_num=u_parcel, community=NA)
plotlabel <- read.csv("../data/systemlabel.txt", header=F,
                          col.names = c("community","comm_label","color","comm_shortlabel"))

for(i in 1:length(u_parcel)){
  parlabel$community[i] <- unique(comm[which(parcel==u_parcel[i])])
}
```

``` r
superheat::superheat(X = z, 
                     y.axis.reverse = TRUE,
                     membership.rows = parlabel$community,
                     membership.cols = parlabel$community,
                     left.label.col=plotlabel$color,
                     bottom.label.col=plotlabel$color,
                     extreme.values.na = FALSE,
                     heat.lim = c(-.2, .2), 
                     heat.pal = rev(brewer.rdylbu(100)),
                     heat.pal.values = c(0, 0.15, 0.25, 0.75,1),
                     title="Parcel x Parcel Correlation Matrix (z-transformed)")
```

![](cifti_MSC_r_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

## Splitting Negative and Positive

``` r
# ==== Setup positive matrix plot
z_pos <- z
z_pos[z<0] <- 0
ss_pos <- superheat::superheat(X = z_pos, 
                     y.axis.reverse = TRUE,
                     membership.rows = parlabel$community,
                     membership.cols = parlabel$community,
                     left.label.col=plotlabel$color,
                     bottom.label.col=plotlabel$color,
                     extreme.values.na = FALSE,
                     heat.lim = c(0, .3), 
                     heat.pal = parula(20),
                     heat.pal.values = c(0, 0.5, 1),
                     title="Node x Node Positive Correlation Matrix (z-transformed")
```

``` r
# ==== Setup negative matrix plot
z_neg <- z
z_neg[z>0] <- 0
ss_neg <- superheat::superheat(X = z_neg, 
                     y.axis.reverse = TRUE,
                     membership.rows = parlabel$community,
                     membership.cols = parlabel$community,
                     left.label.col=plotlabel$color,
                     bottom.label.col=plotlabel$color,
                     extreme.values.na = FALSE,
                     heat.lim = c(-.3, 0), 
                     heat.pal = rev(parula(20)),
                     heat.pal.values = c(0, 0.5, 1),
                     title="Node x Node Negative Correlation Matrix (z-transformed")
```

``` r
gridExtra::grid.arrange(ggplotify::as.grob(ss_pos$plot), ggplotify::as.grob(ss_neg$plot), 
                        nrow=1)
```

![](cifti_MSC_r_files/figure-gfm/pn_matrices-1.png)<!-- -->
