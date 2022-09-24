
# PRS-CS: a polygenic prediction method that infers posterior SNP effect sizes under continuous shrinkage (CS) priors
# using GWAS summary statistics and an external LD reference panel.
#
# Reference: T Ge, CY Chen, Y Ni, YCA Feng, JW Smoller. Polygenic Prediction via Bayesian Regression and Continuous Shrinkage Priors.
#            Nature Communications, 10:1776, 2019.



using ArgParse
#include("parse_genet.jl")
#include("mcmc_gtb.jl")

function parse_param()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--ref_dir"
            help = "Full path (including folder name) to the directory that contains information on the LD reference panel (the snpinfo file and hdf5 files)."
            arg_type = String
            required = true
        "--bim_prefix"
            help = "Full path and the prefix of the bim file for the target (validation/testing) dataset."
            arg_type = String
        "--sst_file"
            help = "Full path and the file name of the GWAS summary statistics."
            arg_type = String
            required = true
        "--n_gwas"
            help = "GWAS sampel size"
            arg_type = Int
        "--out_dir"
            help = "Output directory"
            arg_type = String
            required = true
        "--a"
            help = "Param a in the gamma-gamma prior. Defualts to 1"
            arg_type=Float64
            default=1.0
        "--b"
            help = "Param b in the gamma-gamma prior. Defualts to 0.5"
            arg_type=Float64
            default=0.5
        "--phi"
            help = """Global shrinkage parameter phi. If phi is not specified, it will be learnt from the data using a fully Bayesian approach.
                                    This usually works well for polygenic traits with large GWAS sample sizes (hundreds of thousands of subjects).
                                    For GWAS with limited sample sizes (including most of the current disease GWAS), fixing phi to 1e-4 or 1e-2,
                                    or doing a small-scale grid search (e.g., phi=1e-6, 1e-4, 1e-2, 1) to find the optimal phi value often improves perdictive performance."""
            arg_type = Float64
            default = -1000.0
        "--n_iter"
            help = "Total number of MCMC iterations. Default is 1,000."
            arg_type = Int
            default = 1000
        "--n_burnin"
            help = "Total number of burnin iterations. Default is 5000."
            arg_type = Int
            default = 500
        "--thin"
            help = "Thinning of the Markov chain. Default is 5."
            arg_type = Int
            default = 5
        "--chrom"
            help = "The chromosome on which the model is fitted, separated by comma, e.g., --chrom=\"1,3,5.\""
            arg_type=String
            default = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22"
        "--seed"
            help = "Non-negative integer which seeds the random number generator."
            arg_type = Int
            default = 3
    end

    return parse_args(s)

end

function main()
    parsed_args = parse_param()

    println(parsed_args["chrom"])

    chrom_split = parse.(Int, split(parsed_args["chrom"],","))

    for i in 1:lenth(chrom_split)
        chrom = chrom_split[i]
        print("process chromosome $chrom")

        ref_di
    end

    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end
end

main()


#
# ldblk_dir = "/Users/garethmarkel/Documents/ldblk_1kg_eur"
# chrom = 22
# n_subj = 200000
#
#
#
# ref_dict = parse_ref("/Users/garethmarkel/Documents/ldblk_1kg_eur/snpinfo_1kg_hm3", chrom)
#
# vld_dict = parse_bim("/Users/garethmarkel/Documents/PRScs/test_data/test.bim",chrom)
#
# sst_dict = parse_sumstats(ref_dict, vld_dict, "/Users/garethmarkel/Documents/PRScs/test_data/sumstats.txt", n_subj)
#
# tst = parse_ldblk(ldblk_dir, sst_dict, chrom)
