include("parse_genet.jl")


ldblk_dir = "/Users/garethmarkel/Documents/ldblk_1kg_eur"
chrom = 22
n_subj = 200000



ref_dict = parse_ref("/Users/garethmarkel/Documents/ldblk_1kg_eur/snpinfo_1kg_hm3", chrom)

vld_dict = parse_bim("/Users/garethmarkel/Documents/PRScs/test_data/test.bim",chrom)

sst_dict = parse_sumstats(ref_dict, vld_dict, "/Users/garethmarkel/Documents/PRScs/test_data/sumstats.txt", n_subj)

tst = parse_ldblk(ldblk_dir, sst_dict, chrom)
