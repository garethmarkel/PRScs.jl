# PRS-CS

**PRS-CS** is a Python based command line tool that infers posterior SNP effect sizes under continuous shrinkage (CS) priors
using GWAS summary statistics and an external LD reference panel. Details of the method are described in the article:

T Ge, CY Chen, Y Ni, YCA Feng, JW Smoller. Polygenic Prediction via Bayesian Regression and Continuous Shrinkage Priors. *Nature Communications*, 10:1776, 2019.

**PRScs.jl** is a Julia implementation of the above, with an additional prior added.


# What is PRS-cs?

SNP associations with phenotypes are hard, because p>>(>)n. We often use GWAS to derive
SNP measures of association, which are beta coefficients from regressions of a
genotype at a given locus on a phenotype of interest. Obviously there's a problem here:
you can't include the genotype at all locii in a regression--so each coefficient is picking up loads of other noise.
To improve the accuracy, we need to correct and regularize this. One common and intuitive way is to
impose a bayesian "spike and slab" prior, a mix of a mass at 0 a distribution for nonzero effects.
The problem with this is that you need to search over $2^p$ models, which is hard to do even
before you incorporate linkage disequilibrium information. A different class of prior,
called "continuous shrinkage" priors, can closely approximate this without introducing
a wacky discrete space grid into the problem, which really simplifies the problem.

# What continuous shrinkage priors do we use?

Right now, PRScs.jl implements:
* Strawderman-Berger (*a*=1, *b* = 1/2)
* Horseshoe (*a*=1/2, *b* = 1/2)
* Dirichlet-Laplace (note, good candidates for *a* here are 1/n, 1/p, 0.5). *b* doesn't matter

# Dirichlet-Laplace

The basic model is:

$y = Z\beta + \epsilon$

$\epsilon \sim N(0, \sigma^2 I)$

$p(\sigma^2) \propto \sigma^{-2}$

$\beta_j \sim N(0, \frac{\sigma^2}{N}\tau_j\phi_j^2\lambda^2)$

$\tau_j \sim Exp(1/2)$

$\phi_j \sim Dir(a,...a)$

$\lambda \sim G(pa, 1/2)$

where both $y$ and $z$ have been standardized. Let $\hat{\beta} = Z'y/N$ represent the GWAS estimates for $M$ SNPs. Let $\Psi = diag(\frac{\sigma^2}{N}\tau_j\phi_j^2\lambda^2)_{\{j\}}$, and $D = Z'Z/N$ (the LD matrix).

The Gibbs sampler is as follows:

1. Update $\beta \sim MVN(\frac{N}{\sigma^2}{\Sigma\hat{\beta}}, \Sigma = \frac{\sigma^2}{N}(D + \Psi^{-1})^{-1})$
2. Update $\sigma^2 \sim iG(\frac{N+M}{2}, \frac{N}{2}(1 - 2\beta^T \hat{\beta} + \beta^T(D + \Psi^{-1})\beta))$
3. Update $\tau_j$ by drawing $\tau_j^{-1} \sim iG(\lambda*\phi_j*\sigma/|\beta_j|,1)$ and inverting
4. Update $\lambda \sim gIG(pa-p, 1.0, (2/\sigma)*(\Sigma(|\beta_j|/\phi_j)))$
5. Update $\phi_j \sim gIG(a-1, 1.0, 2|\beta_j|/\sigma)$ then normalize $\phi_j = \phi_j/sum(\phi_j)$

# Strawderman Berger

The basic model is:

$y = Z\beta + \epsilon$

$\epsilon \sim N(0, \sigma^2 I)$

$p(\sigma^2) \propto \sigma^{-2}$

$\beta_j \sim N(0, \frac{\sigma^2}{N}\psi_j)$

$\psi_j \sim G(a, \delta_j)$

$\delta_j \sim G(b, \phi)$

where both $y$ and $z$ have been standardized. Let $\hat{\beta} = Z'y/N$ represent the GWAS estimates for $M$ SNPs. Let $\Psi = diag(\psi_{\{j\}})$, and $D = Z'Z/N$ (the LD matrix).

The Gibbs sampler is as follows:

1. Update $\beta \sim MVN(\frac{N}{\sigma^2}{\Sigma\hat{\beta}}, \Sigma = \frac{\sigma^2}{N}(D + \Psi^{-1})^{-1})$
2. Update $\sigma^2 \sim iG(\frac{N+M}{2}, \frac{N}{2}(1 - 2\beta^T \hat{\beta} + \beta^T(D + \Psi^{-1})\beta))$
3. update $\psi_j \sim giG(a - 0.5, 2\delta_j, \frac{N}{\sigma^2}\beta_j^2)$
4. Update $\delta_j \sim G(a+b, \psi_j + \phi)$