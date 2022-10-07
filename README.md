# PRS-CS

**PRS-CS** is a Python based command line tool that infers posterior SNP effect sizes under continuous shrinkage (CS) priors
using GWAS summary statistics and an external LD reference panel. Details of the method are described in the article:

T Ge, CY Chen, Y Ni, YCA Feng, JW Smoller. Polygenic Prediction via Bayesian Regression and Continuous Shrinkage Priors. *Nature Communications*, 10:1776, 2019.

**PRScs.jl** is a Julia implementation of the above.

Current status:

faster than PRScs just by virtue of the language (the command line arg parsing is much slower, but you can also just run the function)

Future additions:
more flexibility with the CS priors
