# R for brain network

* Documentation of how to 
  * Manipulate gifti (surface), nifti (volume), cifti (surf+vol) in R to generate network matrices
  * Calculate network metrics/measures
  * Analyze data 
  * Visualize results.

### Prerequisites

* The [gifti](https://github.com/muschellij2/gifti) package is required to read in gifti files. 
* A forked version of the r package [cifti](https://github.com/mychan24/cifti) is used for loading in cifti files from MSC data. The cifti files from MSC (at least the .dtseries storing functional time series) are in non-standard arrangement. Standard cifti would have rows x dense (vertex), whereas MSC's ciftis have dense (vertex) x column. Using the original `cifti` package would arrange the cifti from MSC incorrectly (as of 2/7/2019).
* `gridExtra` & `ggplotify` for arranging the plots. `pals` for color pallete. 

* A forked version of the r package [superheat](https://github.com/mychan24/superheat) is used for most of the visualization. Example:  
<img src="https://raw.githubusercontent.com/mychan24/r-for-brain-network/master/gifti/gifti_in_r_files/figure-gfm/unnamed-chunk-6-1.png" width="331" height="370">

## Acknowledgments

* The `superheat` package was super convenient for my work. Check out the original author's [vignette ](https://rlbarter.github.io/superheat/) for examples. 
* The `gifti` and `cifti` packages from [John Muschelli](https://github.com/muschellij2) are obviously critical for interacting with these data format. 
