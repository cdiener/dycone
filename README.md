[![Travis-CI Build Status](https://travis-ci.org/cdiener/dycone.svg?branch=master)](https://travis-ci.org/cdiener/dycone)

dycone
======

...create the dynamic landscape from metabolome data

What does it do?
--------------

Dycone allows you to analyze the dynamic landscape of a metabolic system 
without knowing its kinetic parameters by using steady state and dynamic 
metabolome data (for instance from proteomics).

Installation
-----------

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
instance). Complete docuemntation can be found in the package vignette. The
github repo also contains a scripts folder with some sample analysis. 
