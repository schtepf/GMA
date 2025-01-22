
# gmatools: Geometric multivariate analysis in R

<!-- badges: start -->
<!-- badges: end -->

**TODO** short description

See [NEWS](NEWS.md) for the current version of `gmatools` and recent changes.

## Installation

You can install the development version of `gmatools` from [GitHub](https://github.com/) with:

```r
# install.packages("devtools")
devtools::install_github("schtepf/GMA/pkg/gmatools")
```

## Example

**TODO** show an illustrative example

```r
library(gmatools)
## basic example code
```

## Future work

- Package documentation
  - expand package documentation, adding code examples where possible
  - write a tutorial vignette based on data sets from `corpora` package

- Software tests
  - write software tests for various GMA functions (to the extent possible)

- Integrate further existing functions into `gmatools`
  - `discriminant.plot()` as `gma.discriminant.plot()`
  - `gma.plot.features()`
  - export `.scaleMargins()` as `gma.scale.margins()` or provide a dedicated method for feature contributions
  - `correlation.plot()` ?
  - `parcoord.plot()` ?
  - `lda.cv()` as `gma.lda.cv()`, but needs to be rethought with additional parameters
  - `gma.clust()`, possibly also utilities `majorityLabels()` and `adjustedRandIndex()`
  - interactive 3D visualisations (depending on `rgl` package)
    - `gma.3d()` and associated functions (need to add `Suggest: rgl` and conditional checks in function)
    - `make.movie()` is probably too specialised, interactive WebGL displays are better for end users

- Implement `GMA$drop()` method to remove basis dimensions

- Possibly implement rotations from factor analysis? (`varimax()`, packages `GPArotation` and `psych`)
  - use same interface as `factanal()` rotations, so an arbitrary function can be supplied by user
    (which means we must pass something that can be interpreted as a loadings matrix)
  - only orthogonal rotations will be allowed in GMA, of course
