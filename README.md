# GMA – Geometric Multivariate Analysis

This repository collects materials and software for the method of Geometric Multivariate Analysis (GMA), which has been developed by Diwersy et al. (2014), Evert & Neumann (2017), and Neumann & Evert (2021).

## Mathematical background

A work in progress report on the mathematical background of GMA can be found in [`doc/`](doc/). So far it contains a detailed exposition of the key technique of linear discriminant analysis (LDA) as well as some novel extensions. Further aspects of GMA will be described in future.

## Software

An experimental R package `gmatools` can be found in [`pkg/gmatools/`](pkg/gmatools/), currently at version 0.1. While not available on CRAN yet, it can easily be installed from R with the help of `devtools`.

```r
devtools::install_github("schtepf/GMA/pkg/gmatools")
```

## References

Diwersy, S., Evert, S., and Neumann, S. (2014).
A weakly supervised multivariate approach to the study of language variation.
In Szmrecsanyi, B. and Wälchli, B., editors, _Aggregating Dialectology, Typology, and Register Analysis. Linguistic Variation in Text and Speech_, Linguae et Litterae: Publications of the School of Language and Literature, Freiburg Institute for Advanced Studies, pages 174--204. De Gruyter, Berlin, Boston.

Evert, S. and Neumann, S. (2017).
The impact of translation direction on characteristics of translated texts. A multivariate analysis for English and German.
In De Sutter, G., Lefer, M.-A., and Delaere, I., editors, _Empirical Translation Studies. New Theoretical and Methodological Traditions_, number 300 in Trends in Linguistics. Studies and Monographs (TiLSM), pages 47--80. Mouton de Gruyter, Berlin.
Online supplement: [https://www.stephanie-evert.de/PUB/EvertNeumann2017/](https://www.stephanie-evert.de/PUB/EvertNeumann2017/).

Neumann, S. and Evert, S. (2021).
A register variation perspective on varieties of English.
In Seoane, E. and Biber, D., editors, _Corpus based approaches to register variation_, chapter 6, pages 143--178. Benjamins, Amsterdam.
Online supplement: [https://www.stephanie-evert.de/PUB/NeumannEvert2021/](https://www.stephanie-evert.de/PUB/NeumannEvert2021/).

