## <img src = "docs/figures/bendr_hex.png" align="right" width="200" height = "200"> `bendr`: Built Environment Nested Dirichlet Processes in R 
<!-- badges: start -->
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R build status](https://github.com/apeterson91/rndpp/workflows/R-CMD-check/badge.svg)](https://github.com/apeterson91/rndpp/actions)
<!-- badges: end -->

## About

This is an R package that fits the Nested Dirichlet Process to grouped distance data according to 
an Inhomogenous Poisson Process model. The primary target audience is researchers interested in the effect of built environment features (BEFs) on human health, 
though other applications are possible. See the package's [website](https://apeterson91.github.io/bendr/) for an [introduction](https://apeterson91.github.io/bendr/articles/Introduction.html).
 Currently both normal and beta base measures are implemented. See the documentation for more information.


## Installation

### Development Version

 Currently this package is only available via Github. In order to install the software use the following 
 lines of R code

 ```r
 if(!require(devtools)){
	install.packages("devtools")
	library(devtools)
 }

install_github("apeterson91/bendr",dependencies = TRUE)
 ```

## Contributing

 Examples and code contributions are welcome. Feel free to start/address a feature in the issue tracker and I'll be notified shortly. 

#### Code of Conduct

Please note that `bendr` is released with a [Contributor Code of Conduct](https://www.contributor-covenant.org/). By contributing to this project, you agree to abide by its terms.


## How to cite this package

The software can be cited using the first citation. If looking for an application of the software see the second citation.
```
@manual{Peterson2013bendr,
		author = {Peterson, Adam},
		title = { {bendr}: Built Environment Nested Dirichlet Processes},
		year = {2020},
		howpublished = {\url{https://github.com/apeterson91/bendr}},
		note = {{R} package version 1.0.4}
}
```


```
@misc{peterson2020close,
	  title={How Close and How Much? Linking Health Outcomes to Built Environment Spatial Distributions}, 
	  author={Adam Peterson and Veronica Berrocal and Emma Sanchez-Vaznaugh and Brisa Sanchez},
	  year={2020},
	  eprint={2010.07348},
	  archivePrefix={arXiv},
	  primaryClass={stat.AP}
}
```

## Acknowledgments

This work was developed with support from NIH grant R01-HL131610 (PI: Sanchez).

Special thanks to Emily Hector and Andrew Whiteman for help with the package name.

