# PRS-CS

**PRS-CS** is a Python based command line tool that infers posterior SNP effect sizes under continuous shrinkage (CS) priors
using GWAS summary statistics and an external LD reference panel. Details of the method are described in the article:

T Ge, CY Chen, Y Ni, YCA Feng, JW Smoller. Polygenic Prediction via Bayesian Regression and Continuous Shrinkage Priors. *Nature Communications*, 10:1776, 2019.

**PRScs.jl** is a Julia implementation of the above.

Current status:

* still needs to be tested/verified on genetic data
* DL prior isn't super stable 

Future additions:
* more flexibility with the CS priors ?
* better paralellization

# What is PRS-cs?

SNP associations with phenotypes are hard, because p>>(>)n. We often use GWAS to derive
SNP measures of association, which are beta coefficients from regressions of a
genotype at a given locus on a phenotype of interest. Obviously there's a problem here:
you can't include all the genotypes in a regression--so each coefficient is picking up loads of other noise.
To improve the accuracy, we need to regularize this. One common and intuitive way is to
impose a bayesian "spike and slab" prior, a mix of a mass at 0 a distribution for nonzero effects.
The problem with this is that you need to search over $2^p$ models, which is hard to do even
before you incorporate linkage disequilibrium information. A different class of prior,
called "continuous shrinkage" priors, can closely approximate this without introducing
a wacky discrete space grid into the problem, which really simplifies the problem.

# What continuous shrinkage priors do we use?

Right now, PRScs.jl implements:
* Strawderman-Berger (a=1, b = 1/2)
* Horseshoe (a=1/2, b = 1/2)
* Dirichlet-Laplace (note, good candidates for a here are 1/n, 1/p, 0.5)


# How do you derive the posterior mean for the beta coefficients?

At some point, I'll do a more detailed writeup, but the key point is that the LD
matrix is $D = \frac{X'X}{N}$, and that the GWAS $\hat{\beta}$ we calculate are $\hat{\beta} = \frac{X'Y}{N}$.
The true $\beta$ should be $\beta = (X'X)^{-1}X'Y$. This gives $\beta = D^{-1}\hat{\beta}$.
