SIRFindeR: Stable Estimation of Intron Retention levels from RNA-seq data
============================================================================

SIRFindeR is an R package for detecting and quantifying Intron Retention levels using Second and Third Generation RNA-Sequencing data.

## Package installation 

So as to install SIRFindeR in a R session from GitHub, you will need the R package *devtools*. 
You can install it using:

```

install.packages("devtools")
   
```
Then, install SIRFindeR with:

```

devtools::install_github("lbroseus/SIRFindeR", build_vignette = TRUE)
   
```

## Using SIRFindeR

Guidelines to run SIRFindeR on short or long RNA-seq data can be found in the accompanying [Tutorial](https://github.com/lbroseus/SIRFindeR/blob/master/vignettes/SIRFindeR.pdf); or within a R session:

```

require(SIRFindeR)

vignette("SIRFindeR")
 
```

## References

1. _Broseus and Ritchie. Challenges in detecting and quantifying intron retention from RNA-seq data. CSBJ. (2020)_
2. _Broseus and Ritchie. S-IRFindeR: stable and accurate measurement of intron retention. bioRXiv. (2020)_
