[![Build Status](https://travis-ci.org/cdiener/dycone.svg?branch=master)](https://travis-ci.org/cdiener/dycone)
![dycone](stuff/logo.svg)

...analyze enzyme regulation with metabolome data

What does it do?
--------------

Dycone allows you to analyze the dynamic landscape of a metabolic system 
without knowing its kinetic parameters by using steady state and dynamic 
metabolome data (for instance from proteomics).

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
