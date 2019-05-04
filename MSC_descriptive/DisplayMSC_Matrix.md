MSC All Subjects RSFC Matrices
================

``` r
# rename communities
#--- read the file with community information: label number & community name & color & abbreviation for community name
colnames(CommName) <- c("Comm", "CommLabel","CommColor","CommLabel.short")
CommName[,'Community'] <- c("UnAssign", "Default","lateral Visual","Frontoparietal",
                            "medial Visual","dorsal Attention","Premotor",
                            "ventral Attention","Salience","Cingular opercular",
                            "Sensorimotor - hand","Sensorimotor - face","Auditory",
                            "anterior Medial temporal","posterior Medial temporal",
                            "Parietal memory","Context","Sensorimotor - foot", 
                            "Unknown","Unknown","Unknown","UnKnown","UnKnown")
# pander::pander(CommName[,c("Comm", "Community","CommLabel.short")])
```

``` r
# read parcel labels for each subject
parcel.list <- lapply(1:length(parcelfile2read), function(x){
                      parcel <- read.table(paste0(parcel.comm.path, parcelfile2read[x]),sep = ",")
                      getVoxDes(parcel,CommName)
                      })
names(parcel.list) <- subj.name

#-- Create colors for heatmap 
labelcol <- list()
textcol <- list()
for(i in 1:length(subj.name)){
  labelcol[[i]] <- parcel.list[[i]]$Comm.col$gc[order(rownames(parcel.list[[i]]$Comm.col$gc))]
  names(labelcol)[i] <- subj.name[i]
  
  textcol[[i]] <- rep("black", length(labelcol[[i]]))
  textcol[[i]][as(colorspace::hex2RGB(labelcol[[i]]), "polarLUV")@coords[,1] < 35] <- "white"  # Convert hex2RGB to lum
}
```

``` r
### Make these not objects, but part of a list? ### MYC NEED TO CHANGE 20190404
#--- plot heatmap
hmaplist <- list()
hmaplist_i <- 1
for(subj.count in 1:length(parcel.list)){
  load(sprintf("%s/sub-MSC%02d_zcube_rcube.RData",zmat.path,subj.count))
  hmaplist[[subj.count]] <- list()
  for(session.count in 1:dim(cubes$rcube)[3]){
    # heat map of correlation matrices
    hmap.name <- sprintf("sub%02d_hmap%s", subj.count, session.count)
    hmaplist[[subj.count]][[session.count]] <- superheat(cubes$rcube[,,session.count],
                                membership.cols = parcel.list[[subj.count]]$vox.des$Comm.rcd,
                                membership.rows = parcel.list[[subj.count]]$vox.des$Comm.rcd,
                                clustering.method = NULL,
                                heat.col.scheme = "viridis",
                                left.label.size = 0.08,
                                bottom.label.size = 0.05,
                                y.axis.reverse = TRUE,
                                left.label.col = labelcol[[subj.count]], # order by community name
                                bottom.label.col = labelcol[[subj.count]],
                                left.label.text.size = 3,
                                bottom.label.text.size = 2,
                                left.label.text.col = textcol[[subj.count]],
                                bottom.label.text.col = textcol[[subj.count]],
                                left.label.text.alignment = "left",
                                title = sprintf("Correlation matrix of #%s session",session.count)
                                )$plot
  }
}
```

``` r
for(i in 1:length(hmaplist)){
  print(sprintf("%s", subj.name[i]))
  do.call(grid.arrange, c(hmaplist[[i]], ncol=5))
}
```

    ## [1] "sub01"

![](DisplayMSC_Matrix_files/figure-gfm/plot_all_sub_session_hmap-1.png)<!-- -->

    ## [1] "sub02"

![](DisplayMSC_Matrix_files/figure-gfm/plot_all_sub_session_hmap-2.png)<!-- -->

    ## [1] "sub03"

![](DisplayMSC_Matrix_files/figure-gfm/plot_all_sub_session_hmap-3.png)<!-- -->

    ## [1] "sub04"

![](DisplayMSC_Matrix_files/figure-gfm/plot_all_sub_session_hmap-4.png)<!-- -->

    ## [1] "sub05"

![](DisplayMSC_Matrix_files/figure-gfm/plot_all_sub_session_hmap-5.png)<!-- -->

    ## [1] "sub06"

![](DisplayMSC_Matrix_files/figure-gfm/plot_all_sub_session_hmap-6.png)<!-- -->

    ## [1] "sub07"

![](DisplayMSC_Matrix_files/figure-gfm/plot_all_sub_session_hmap-7.png)<!-- -->

    ## [1] "sub08"

![](DisplayMSC_Matrix_files/figure-gfm/plot_all_sub_session_hmap-8.png)<!-- -->

    ## [1] "sub09"

![](DisplayMSC_Matrix_files/figure-gfm/plot_all_sub_session_hmap-9.png)<!-- -->

    ## [1] "sub010"

![](DisplayMSC_Matrix_files/figure-gfm/plot_all_sub_session_hmap-10.png)<!-- -->
