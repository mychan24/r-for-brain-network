## A trick to loop superheat() and show it in R markdown

### Generate heatmap object using superheat in an rchunk that hides the figure
* Use `fig.show='hide'` 
* Create listA
  + `listA <- list()`
* Creat a list inside listA
  + `listA[[1]] <- list()`
* Assign a superheat()$plot into the second list
  + `listA[[1]][[1]] <- superheat(mat)$plot`
* Do that in a loop that oragnizes it (probably can use `apply` to do this...)

### Show heatmaps using grid.arrange()
* In a for loop, use `grid.arrange()` in `do.call()`
```r
`for(i in 1:length(listA)){
  do.call(grid.arrange, c(list(A)[[i]], ncol=3))
  }`
```
