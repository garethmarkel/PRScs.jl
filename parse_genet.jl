using DataFrames, CSV, Distributions, HDF5, LinearAlgebra

##read in the reference file and subset to a specific chromosome
function parse_ref(ref_file, chrom)
    print("... parse reference file: $ref_file ...")

    df = DataFrame(CSV.File(ref_file))

    ref_dict = filter(:CHR => n -> n == chrom, df)

    n_snp = nrow(ref_dict)

    print("... $n_snp SNPs on chromosome $chrom read from $ref_file ...")

    return ref_dict

end


function parse_bim(bim_file, chrom)
    print("... parse bim file: $bim_file.bim ...")

    df = DataFrame(CSV.File("$bim_file.bim", header = false))

    #Assign appropriate names to the columns
    #like, yes, this doesn't check to see if your column order is valid but also:
    #   1. this is a script, not a production grade system
    #   2. If i were a user i would simply have a correctly formatted bim file
    rename!(df, [:CHR, :SNP, :POS ,:BP, :A1 , :A2])

    #select SNPs on our chromosome of interest
    vld_dict = select(filter(:CHR => n -> n == chrom, df), :SNP, :A1,:A2)

    n_snp = nrow(vld_dict)

    print("... $n_snp SNPs on chromosome $chrom read from $bim_file.bim ...")
    return vld_dict
end

#function to plug into filter to make sure both allelles are valid
function valid_allelles(a1, a2)::Bool
    valid_a1 = a1 in ["A","C","G","T"]
    valid_a2 = a2 in ["A","C","G","T"]
    valid_a1 && valid_a2
end

#function to generate valid pairs for allelles
#TODO: find a clever way to refactor this
function pair_allelles(a)

    a2 = Vector{String}()

    for i in 1:length(a)
        if a[i] == "T"
            push!(a2, "A")
        end
        if a[i] == "A"
            push!(a2, "T")
        end
        if a[i] == "G"
            push!(a2, "C")
        end
        if a[i] == "C"
            push!(a2, "G")
        end
    end
    return a2
end


function parse_sumstats(ref_dict, vld_dict, sst_file, n_subj)
    print("... parse sumstats file: $sst_file ...")

    #filter for valid allelles--no real need to split this across two lines
    df = filter([:A1, :A2] => valid_allelles, DataFrame(CSV.File(sst_file)))

    #SNPs and allelles in our GWAS summary statistics
    sst_dict = select(df, :SNP, :A1, :A2)

    ##our "truth" is the validation dict, which has what we think is the right SNP/a1/a2
    ##here, for each snp we create dictionaries for A1/A2, A2/A1, A1/(correct pairing for A1),A2/correct pairing for A2
    ##We could just merge on "SNP" and take the .bim ordering--by not doing so,
    ##  we imply that some of the snp/a1/a2 orderings in the bim file don't
    ##  appear in the 4 sets of combinations we see above. I don't know how realistic this is
    ref_snp = [select(ref_dict, :SNP,:A1,:A2);
        select(rename(ref_dict, Dict(:A1 => :A2, :A2 => :A1)), :SNP,:A1,:A2);
        DataFrame(SNP = ref_dict[:,:SNP], A1 = ref_dict[:,:A1], A2 = pair_allelles(ref_dict[:,:A1]));
        DataFrame(SNP = ref_dict[:,:SNP],A1 = pair_allelles(ref_dict[:,:A2]), A2 = ref_dict[:,:A2])]

    sst_snp = [select(sst_dict, :SNP,:A1,:A2);
        select(rename(sst_dict, Dict(:A1 => :A2, :A2 => :A1)), :SNP,:A1,:A2);
        DataFrame(SNP = sst_dict[:,:SNP], A1 = sst_dict[:,:A1], A2 = pair_allelles(sst_dict[:,:A1]));
        DataFrame(SNP = sst_dict[:,:SNP],A1 = pair_allelles(sst_dict[:,:A2]), A2 = sst_dict[:,:A2])]

    comm_snp = innerjoin(innerjoin(vld_dict, ref_snp, on = [:SNP, :A1, :A2]),sst_snp, on = [:SNP, :A1, :A2])


    n_snp = nrow(comm_snp)

    #convert the common snp dict to a vector of tuples
    comm_snp = Tuple.(eachrow(comm_snp))

    print("... $n_snp common SNPs in the reference, sumstats, and validation set ...")

    #N subjects in GWAS
    n_sqrt = sqrt(n_subj)

    df[!,:A1_map] = pair_allelles(df[:,:A1])
    df[!,:A2_map] = pair_allelles(df[:,:A2])

    df[!, :BETA_STD] .= 0.0

    Ndist = Normal(0,1)

    if "OR" in names(df)
        df[!, :BETA] = log.(df[:,:OR])
    end

    #for each SNP/(A1,A2) pairing in the sumstats file, see if it appears in our list of
    # valid SNP choices
    # would just be much faster to see if the sst snps are in the valiation dict
    # this is to get flipped SNPs (probably not the most efficient way)
    in_comm1 = Tuple.(eachrow(DataFrame(SNP = df[:,:SNP], A1 = df[:,:A1], A2 = df[:,:A2]))) .∈ Ref(comm_snp)
    in_comm2 = Tuple.(eachrow(DataFrame(SNP = df[:,:SNP], A1 = df[:,:A1_map], A2 = df[:,:A2_map]))) .∈ Ref(comm_snp)

    ##standardize beta coefficients
    df[in_comm1 .| in_comm2, :BETA_STD] = sign.(df[in_comm1 .| in_comm2, :BETA]) .* abs.(quantile.(Ndist, max.(df[in_comm1 .| in_comm2, :P]./2, 1e-323))) ./ n_sqrt

    #standardize beta coefficients for flipped paits
    in_comm3 = Tuple.(eachrow(DataFrame(SNP = df[:,:SNP], A1 = df[:,:A2], A2 = df[:,:A1]))) .∈ Ref(comm_snp)
    in_comm4 = Tuple.(eachrow(DataFrame(SNP = df[:,:SNP], A1 = df[:,:A2_map], A2 = df[:,:A1_map]))) .∈ Ref(comm_snp)
    df[in_comm3 .| in_comm4, :BETA_STD] = -1 .* sign.(df[in_comm3 .| in_comm4, :BETA]) .* abs.(quantile.(Ndist, max.(df[in_comm3 .| in_comm4, :P]./2, 1e-323))) ./ n_sqrt

    #subset to valid SNPs
    df = df[in_comm1 .| in_comm4 .| in_comm2 .| in_comm3,:]

    #subset the reference LD dict to valid SNPs
    in_sst = ref_dict[:,:SNP] .∈ Ref(df[:,:SNP])

    #get the information from the reference dict and join it to the effect sizes
    sst_dict = innerjoin(select(ref_dict, :CHR, :SNP, :BP,:A1,:A2, :MAF), select(df, :SNP, :BETA_STD), on = :SNP)
    sst_dict[!,:A1_map] = pair_allelles(sst_dict[:,:A1])
    sst_dict[!,:A2_map] = pair_allelles(sst_dict[:,:A2])
    sst_dict[!,:tmp] = copy(sst_dict[:,:A2])
    sst_dict[!, :FLP] .= 0

    #detect flipped allelle sets
    in_comm1 = Tuple.(eachrow(DataFrame(SNP = sst_dict[:,:SNP], A1 = sst_dict[:,:A1], A2 = sst_dict[:,:A2]))) .∈ Ref(comm_snp)
    in_comm2 = Tuple.(eachrow(DataFrame(SNP = sst_dict[:,:SNP], A1 = sst_dict[:,:A1_map], A2 = sst_dict[:,:A2_map]))) .∈ Ref(comm_snp)
    in_comm3 = Tuple.(eachrow(DataFrame(SNP = sst_dict[:,:SNP], A1 = sst_dict[:,:A2], A2 = sst_dict[:,:A1]))) .∈ Ref(comm_snp)
    in_comm4 = Tuple.(eachrow(DataFrame(SNP = sst_dict[:,:SNP], A1 = sst_dict[:,:A2_map], A2 = sst_dict[:,:A1_map]))) .∈ Ref(comm_snp)

    sst_dict[in_comm2,:A1] = sst_dict[in_comm2,:A1_map]
    sst_dict[in_comm2,:A2] = sst_dict[in_comm2,:A2_map]
    sst_dict[in_comm1 .| in_comm2, :FLP] .= 1

    sst_dict[in_comm3,:tmp] = sst_dict[in_comm3,:A1]
    sst_dict[in_comm3,:A1] = sst_dict[in_comm3,:A2]
    sst_dict[in_comm3,:A2] = sst_dict[in_comm3,:tmp]

    sst_dict[in_comm4,:A1] = sst_dict[in_comm2,:A2_map]
    sst_dict[in_comm4,:A2] = sst_dict[in_comm2,:A1_map]

    sst_dict[in_comm3 .| in_comm4, :FLP] .= -1
    sst_dict[in_comm3 .| in_comm4, :MAF] = 1. .- sst_dict[in_comm3 .| in_comm4, :MAF]

    sst_dict[!,:BETA] = sst_dict[:, :BETA_STD]


    ret_df = select(sst_dict, :CHR, :SNP, :BP, :A1, :A2, :MAF, :BETA, :FLP)
    return ret_df
end

function parse_ldblk(ldblk_dir, sst_dict, chrom)

    print("... parse reference LD on chromosome $chrom ...")

    if occursin("1kg", ldblk_dir)
        chr_name = "$ldblk_dir/ldblk_1kg_chr$chrom.hdf5"
    end

    if occursin("ukbb", ldblk_dir)
        chr_name = "$ldblk_dir/ldblk_ukbb_chr$chrom.hdf5"
    end

    hdf_chr = h5open(chr_name, "r")

    n_blk = length(hdf_chr)

    ###NOTE: THIS LD BLOCK IS NOT THE SAME AS THE PYTHON VERSION--check block 1
    ld_blk = [read(hdf_chr["blk_$i"]["ldblk"]) for i in 1:n_blk]

    snp_blk = [read(hdf_chr["blk_$i"]["snplist"]) for i in 1:n_blk]

    idx_blks = [[j for j in 1:length(i) if i[j] in sst_dict[:,:SNP]] for i in snp_blk]

    blk_size = length.(idx_blks)

    mm = 0

    for i in 1:n_blk
        if idx_blks[i] != []
            flip = sst_dict[(mm+1):(mm + blk_size[i]), :FLP]

            #negative correlations for flipped alleles
            ld_blk[i] = ld_blk[i][idx_blks[i], idx_blks[i]] .* (flip * flip')

            U, S, V = svd(ld_blk[i])

            ##reconstruct the matrix, take the average of the reconstructed values and
            ##original
            h = V*Diagonal(S)*V'
            ld_blk[i] = (ld_blk[i] .+ h) ./ 2

            mm = mm + blk_size[i]
        end

        #empty blocks
        if idx_blks[i] == []
            ld_blk[i] = Array{Float64}(undef, 0,0)
        end
    end

    return ld_blk, blk_size

end
