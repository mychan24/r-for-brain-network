# R for brain network

### Why R for brain network?  
Human brain network derived from functional magnetic resonance imaging (fMRI) data are often analyzed and visualize using a combination of MATLAB or network-visualization software (e.g., cytoscape, gephi). 

However, much of graph analysis can be accomplished in R. Here, I document a graph-analysis workflow of fMRI resting-state network in R. This allows easy integration with other statistical analysis conducted on the targeted graph measures (e.g., degree, participation coefficient, modularity, system segregation).

**This repository includes documentation of the following:**

  * Manipulate gifti (surface), nifti (volume), cifti (surface + volume) in R to generate network matrices
    * Checkout the [cifti markdown file](https://github.com/mychan24/r-for-brain-network/blob/master/cifti_single_sub_sample/cifti_MSC_r.md) for most updated workflow (3/5/2019). 
	* Most future addition would focus on cifti, as the combination of surface and volumetric data (subcortical) are preferrable. 
  * Calculate network metrics/measures
  * Visualize network matrices and graphs
  * Analyze data 
  * Visualize results

### Prerequisites

* The [gifti](https://github.com/muschellij2/gifti) and [cifti](https://github.com/muschellij2/cifti) packages are required to read in gifti or cifti files. 
* `gridExtra` & `ggplotify` for arranging the plots. `pals` for color pallete. 
* A forked version of the r package [superheat](https://github.com/mychan24/superheat) is used for most of the visualization. Example:  
<img src="https://raw.githubusercontent.com/mychan24/r-for-brain-network/master/docs/matrix.png" width="331" height="370">

* The `igraph` is used to make network visualizaiton (i.e., hairball/spring graphs)
<img src="https://raw.githubusercontent.com/mychan24/r-for-brain-network/master/docs/spring.png" width="500" height="370">

## Acknowledgments

* The `superheat` package was super convenient for my work. Check out the original author's [vignette ](https://rlbarter.github.io/superheat/) for examples. 
* The `gifti` and `cifti` packages from [John Muschelli](https://github.com/muschellij2) are obviously critical for interacting with these data format. 
