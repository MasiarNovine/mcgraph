# mcgraph

Analysis of synthetic networks based on Monte Carlo simulations

## Common installation

Download the current tarball `mcgraph_0.6.0.tar.gz` containing the compiled mcgraph package. Move to the directory, where you stored the file and do

```
R CMD INSTALL mcgraph_0.6.0.tar.gz
```

The package should now be available after opening R and running

```
library(mcgraph)
```

## Build from source (advanced)

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

### Build vignette (not mandatory)

By the above step, the package is workable, but the vignette might be not up-to-data.
To update the vignette, you need the R package `knitr` and the `pandoc` tool.

Before building the core source, change into the `vignette` subfolder of the directory you downloaded.

Open R and load the `knitr` package

```
library(knitr)
knit('mcgraph.Rmd')
```
Close R and run on the terminal

```
pandoc -i mcgraph.md --citeproc --bibliography bibliography.bib -o mcgraph.pdf
```

You can remove the subfolder `figure` and copy the vignette `mcgraph.pdf` into the folder `inst/doc/`.

## License

MIT
