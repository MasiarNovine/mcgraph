# mcgraph

Analysis of synthetic networks based on Monte Carlo simulations

## Installing

Download the current tarball `mcgraph_0.6.0.tar.gz` containing the compiled mcgraph package. Move to the directory, where you stored the file and do 

```
R CMD INSTALL mcgraph_0.6.0.tar.gz 
```

The package should now be available after opening R and running 

```
library(mcgraph)
```

## Build from source

Download the `mcgraph` folder of this repository.

Go to the directory where the folder has been saved to.

Open R and load the `roxygen2` package (install the package if necessary) and run

```
library(roxygen2)
roxygenize("mcgraph")	# needed for package building
``` 

Afterwards do 

```
R CMD build --no-build-vignettes mcgraph
```

You can check the resulting `mcgraph_0.6.0.tar.gz` file by 

```
R CMD check --cran mcgraph_0.6.0.tar.gz
```

Install the package by

```
R CMD INSTALL mcgraph_0.6.0.tar.gz
```

## License

MIT
