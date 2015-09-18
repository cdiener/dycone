[![Build Status](https://travis-ci.org/cdiener/dycone.svg?branch=master)](https://travis-ci.org/cdiener/dycone)
[![codecov.io](http://codecov.io/github/cdiener/dycone/coverage.svg?branch=master)](http://codecov.io/github/cdiener/dycone?branch=master)

![dycone](stuff/logo.png)

...analyze enzyme regulation with metabolome data

What does it do?
--------------

Dycone ("dynamic cone") allows you infer enzymatic regulation from 
metabolome mesurements. It employs formalisms based on flux and k-cone 
analysis to connect metabolome data to distinct regulations of enzyme activity. 
Most of the analysis methods can be applied to genome-scale data. 

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
instance). An example can be found in the package vignette. For an exhaustive
example you can look at the Rmarkdown generated protocol at
https://github.com/cdiener/kcone-paper.
