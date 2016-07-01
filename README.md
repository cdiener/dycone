[![Build Status](https://travis-ci.org/cdiener/dycone.svg?branch=master)](https://travis-ci.org/cdiener/dycone)
[![wercker status](https://app.wercker.com/status/b8b682f3278fd21b27afab12fb050db7/s "wercker status")](https://app.wercker.com/project/bykey/b8b682f3278fd21b27afab12fb050db7)
[![codecov](https://codecov.io/gh/cdiener/dycone/branch/master/graph/badge.svg)](https://codecov.io/gh/cdiener/dycone)
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.49987.svg)](http://dx.doi.org/10.5281/zenodo.49987)

![dycone](stuff/logo.png)

...analyze enzyme regulation with metabolome data

What does it do?
--------------

Dycone ("dynamic cone") allows you infer enzymatic regulation from 
metabolome mesurements. It employs formalisms based on flux and k-cone 
analysis to connect metabolome data to distinct regulations of enzyme activity. 
Most of the analysis methods can be applied to genome-scale data. 

Attribution
-----------

If you use dycone in your academic work you may use the following paper for a citation

> Diener C, Muñoz-Gonzalez F, Encarnación S, Resendis-Antonio O.     
> The space of enzyme regulation in HeLa cells can be inferred from its intracellular metabolome.      
> Sci Rep. 2016 Jun 23;6:28415. doi: 10.1038/srep28415. PubMed PMID: 27335086.

Installation
-----------

## External dependencies

Even though dycone can be installed directly via R as described in the following 
section, some of its dependencies require the installation of packages via
your package manager. Under Debian and Ubuntu you can install those dependencies
via

```bash
sudo apt-get install libcurl4-openssl-dev libssl-dev libgmp-dev
```

In case you do not have R installed you will also need to install it via

```bash
sudo apt-get install r-base r-base-dev
```

## R package

The package can be installed directly from github by using `devtools`:

```R
install.packages("devtools")
devtools::install_github("cdiener/dycone")
```

This will also install all required dependencies and you are good to go :)

You can start using the package by importing it:
```R
library(dycone)
```

Docs
----

Documentation can be found using the included help in R (`?function_name` for 
instance). An example using the Human Red Blood Cell model can be found in the 
package vignette "eryth": 

```R
vignette("eryth")
```

For an exhaustive example you can look at the Rmarkdown generated protocol at 
https://github.com/cdiener/kcone-paper.

Roadmap
-------

To see planned changes to dycone please consult the [project milestones](https://github.com/cdiener/dycone/milestones).
