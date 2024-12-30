using DataFrames, CSV, Distributions, HDF5, LinearAlgebra
using Random, Distributions, GenInvGaussian

include("./parse_genet.jl")
include("./mcmc_gtb.jl")


# PRS-CS: a polygenic prediction method that infers posterior SNP effect sizes under continuous shrinkage (CS) priors
# using GWAS summary statistics and an external LD reference panel.
#
# Reference: T Ge, CY Chen, Y Ni, YCA Feng, JW Smoller. Polygenic Prediction via Bayesian Regression and Continuous Shrinkage Priors.
#            Nature Communications, 10:1776, 2019.
