# AP-SKAT
[AP-SKAT: highly-efficient genome-wide rare variant association test](https://pubmed.ncbi.nlm.nih.gov/27654840/)

## Availability and requirements  

Project name:AP-SKAT  

Project home page: http://nagasakilab.csml.org/data/ aSKAT.zip  

Operating system(s): Platform independent Programming language:R  

Any restrictions to use by non-academics:Pleasecontact authors for commercial use. 



## Disadvantages
An inefficient "for" loop is used.

## Optimization
```R
library(parallel) # Support for parallel computation, including random-number generation.
```

## Notes

The original article had no contribution from me, I just optimized it.
