MSC All Subjects RSFC Matrices
================

| Comm |         Community         | CommLabel.short |
| :--: | :-----------------------: | :-------------: |
|  0   |         UnAssign          |      00Bd       |
|  1   |          Default          |      01DMN      |
|  2   |      lateral Visual       |     02lVis      |
|  3   |      Frontoparietal       |      03FPN      |
|  4   |       medial Visual       |     04mVis      |
|  5   |     dorsal Attention      |      05DAN      |
|  6   |         Premotor          |      06PMo      |
|  7   |     ventral Attention     |      07VAN      |
|  8   |         Salience          |      08SLC      |
|  9   |    Cingular opercular     |      09CON      |
|  10  |    Sensorimotor - hand    |     10hSMN      |
|  11  |    Sensorimotor - face    |     11fSMN      |
|  12  |         Auditory          |      12Aud      |
|  13  | anterior Medial temporal  |     13aMTL      |
|  14  | posterior Medial temporal |     14pMTL      |
|  15  |      Parietal memory      |      15PMN      |
|  16  |          Context          |      16CAN      |
|  17  |    Sensorimotor - foot    |     17fSMN      |
|  21  |          Unknown          |      21UN       |
|  25  |          Unknown          |      25UN       |
|  29  |          Unknown          |      29UN       |
|  34  |          UnKnown          |      34UN       |
|  37  |          UnKnown          |      37UN       |

``` r
# read parcel labels for each subject
parcel.comm.path <- "../data/parcel_community"
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
