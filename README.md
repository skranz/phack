This R package implements the tests from Elliot et al. (2022) for detecting p-Hacking. 

For documentation and examples, take a look at the vignette here:

[https://github.com/skranz/phack/blob/main/documentation/vignette_phack.pdf](https://github.com/skranz/phack/blob/main/documentation/vignette_phack.pdf)


The package is essentially a simple wrapper to the code provided
in the code and data supplement of the article, with some cosmetical changes.
The original code can be found in the code and data supplement of the article.
  
You can install the package from my r-universe repository by running the following code:

```r
options(repos = c(
  skranz = 'https://skranz.r-universe.dev',
  CRAN = 'https://cloud.r-project.org'))
install.packages('phack')
```

  
  
Reference:

Elliott, G., Kudrin, N., & Wüthrich, K. (2022). Detecting p‐Hacking. Econometrica, 90(2), 887-906.
