# mcgraph

An R package for the analysis of synthetic large networks based on Monte Carlo resampling. *mcgraph* allows for creation and visualisation of different types of random networks based on a given graph structure and test the ability of algorithms to rediscover the initial graph structure.

Several algorithms (e.g. lasso, step-wise linear regression, GLM, a.o.) and evaluation metrics (e.g. accuracy, sensitivity, MCC, ROC/PROC curves, a.o.) are currently implemented.

## Common installation

Download the current tarball `mcgraph_0.6.1.tar.gz` containing the compiled mcgraph package, either by right click and saving it to your system or by

```bash
curl -L https://github.com/MasiarNovine/mcgraph/releases/download/v0.6.1/mcgraph_0.6.1.tar.gz mcgraph_0.6.1.tar.gz
```

Move to the directory, where you stored the file and do

```bash
R CMD INSTALL mcgraph_0.6.1.tar.gz
```

The package should now be available after opening R and running

```bash
library(mcgraph)
```

## Build from source (advanced)

Download the `mcgraph` folder of this repository.

Go to the directory where the folder has been saved to.

Open R and load the `roxygen2` package (install the package if necessary) and run

```r
library(roxygen2)
roxygenize("mcgraph")   # needed for package building
```

Afterwards do

```bash
R CMD build --no-build-vignettes mcgraph
```

You can check the resulting `mcgraph_0.6.1.tar.gz` file by

```bash
R CMD check --cran mcgraph_0.6.1.tar.gz
```

Install the package by

```bash
R CMD INSTALL mcgraph_0.6.1.tar.gz
```

### Build vignette

By the above step, the package is workable, but the vignette might be not up-to-data.
To update the vignette, you need the R package `knitr` and the `pandoc` tool.

Before building the core source, change into the `vignette` subfolder of the directory you downloaded.

Open R and load the `knitr` package

```r
library(knitr)
knit('mcgraph.Rmd')
```

Close R and run on the terminal

```bash
pandoc -i mcgraph.md --citeproc --bibliography bibliography.bib -o mcgraph.pdf
```

You can remove the subfolder `figure` and copy the vignette `mcgraph.pdf` into the folder `inst/doc/`.

## License

MIT
