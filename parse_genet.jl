"""
    ParseRef(ref_file, chrom)

This function reads in the snp info LD reference file for a given chromosome
...
# Arguments
- `ref_file`: Reference file path
- `chrom`: chromosome number
...
# Returns
- `ref_dict`: dataframe of SNPs appearing on chromosome chrom
...
"""
function ParseRef(ref_file, chrom)

    print("... parse reference file: $ref_file ...")

    df = DataFrame(CSV.File(ref_file))

    ref_dict = filter(:CHR => n -> n == chrom, df)

    n_snp = nrow(ref_dict)

    print("... $n_snp SNPs on chromosome $chrom read from $ref_file ...")

    return ref_dict

end


"""
    ParseBim(ref_file, chrom)

This function reads in the sample bim file for a given chromosome
...
# Arguments
- `bim_file`: bim file path, with columns chr, snp, pos, bp, a1, a2 in that order
- `chrom`: chromosome number
...
# Returns
- `valid_snp_dict`: dataframe of SNPs in the sample on the given chromosome 
...
"""
function ParseBim(bim_file, chrom)
    print("... parse bim file: $bim_file.bim ...")

    raw_bim = DataFrame(CSV.File("$bim_file.bim", header = false))

    rename!(raw_bim, [:CHR, :SNP, :POS ,:BP, :A1 , :A2])

    #select SNPs on our chromosome of interest
    valid_snp_dict = select(filter(:CHR => n -> n == chrom, raw_bim), :SNP, :A1,:A2)

    n_snp = nrow(valid_snp_dict)

    print("... $n_snp SNPs on chromosome $chrom read from $bim_file.bim ...")
    return valid_snp_dict
end

"""
    CheckValidAllelles(a1, a2)

This function checks to make sure a pair of allelles is valid
...
# Arguments
- `a1`: chr for major SNP allelle
- `a2`: chr for minor SNP allelle
...
# Returns
- boolean for whether allelle pair is valid
...
"""
function CheckValidAllelles(a1, a2)::Bool
    valid_a1 = a1 in ["A","C","G","T"]
    valid_a2 = a2 in ["A","C","G","T"]
    valid_a1 && valid_a2
end

"""
    GetAllelePair(a)

This function gets the appropriate complementary allele for A,C,G,T. Each allele
    can also pair with itself (So A can be paired as AA or AT).
...
# Arguments
- `a`: vector of chr for major SNP allelle
...
# Returns
- `a2`: vector of chr for minor SNP allelle corresponding to a
...
"""
function GetAllelePair(a)

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


"""
    ParseSumstats(ref_dict, vld_dict, sst_file, n_subj)

This function parses a GWAS summary statistics file.
...
# Arguments
- `ref_dict`: dataframe of SNPs appearing on chromosome chrom
- `vld_dict`: dataframe of SNPs in the sample on a given chromosome, from ParseBim
- `sst_file`: GWAS file path, with columns :SNP, :A1, :A2 appearing
- `n_subj`: number of GWAS subjects for that SNP (e.g. MTAGed statistics)
...
# Returns
- `ret_df`: valid SNPs from summary statistic files, with columns
     :CHR, :SNP, :BP, :A1, :A2, :MAF, :BETA, :FLP
...
"""
function ParseSumstats(ref_dict, vld_dict, sst_file, n_subj)
    print("... parse sumstats file: $sst_file ...")

    #filter for valid allelles--no real need to split this across two lines
    df = filter(
        [:A1, :A2] => CheckValidAllelles, 
        DataFrame(CSV.File(sst_file))
        )

    #SNPs and allelles in our GWAS summary statistics
    sst_dict = select(df, :SNP, :A1, :A2)


    # Our "truth" is the validation dict, which has what we think is the right SNP/a1/a2.
    # Here, for each snp we create dictionaries for A1/A2, A2/A1, 
    #     A1/(correct pairing for A1), A2/correct pairing for A2.
    # We could just merge on "SNP" and take the .bim ordering--by not doing so,
    #   we imply that some of the snp/a1/a2 orderings in the bim file don't
    #   appear in the 4 sets of combinations we see above. This is probably 
    #   overkill to worry about, but it's not a bottleneck.
    ref_snp = [
        select(
            ref_dict, 
            :SNP, :A1, :A2
            );
        select(
            rename(ref_dict, Dict(:A1 => :A2, :A2 => :A1)), 
            :SNP, :A1, :A2
            );
        DataFrame(
            SNP = ref_dict[:,:SNP], 
            A1 = ref_dict[:,:A1], 
            A2 = GetAllelePair(ref_dict[:,:A1])
            );
        DataFrame(
            SNP = ref_dict[:,:SNP],
            A1 = GetAllelePair(ref_dict[:,:A2]), 
            A2 = ref_dict[:,:A2]
            )
        ]

    sst_snp = [
        select(
            sst_dict, 
            :SNP, :A1, :A2
            );
        select(
            rename(sst_dict, Dict(:A1 => :A2, :A2 => :A1)), 
            :SNP, :A1, :A2
            );
        DataFrame(
            SNP = sst_dict[:,:SNP], 
            A1 = sst_dict[:,:A1], 
            A2 = GetAllelePair(sst_dict[:,:A1])
            );
        DataFrame(
            SNP = sst_dict[:,:SNP],
            A1 = GetAllelePair(sst_dict[:,:A2]), 
            A2 = sst_dict[:,:A2]
            )
        ]

    # common SNPs betweent he bim file, the reference panel, and the sumstats file
    comm_snp = innerjoin(
        innerjoin(
            vld_dict, 
            ref_snp, 
            on = [:SNP, :A1, :A2]),
        sst_snp, on = [:SNP, :A1, :A2]
        )

    n_snp = nrow(comm_snp)

    #convert the common snp dict to a vector of tuples
    comm_snp = Tuple.(eachrow(comm_snp))

    print("... $n_snp common SNPs in the reference, sumstats, and validation set ...")

    #N subjects in GWAS
    n_sqrt = sqrt(n_subj)

    # query the mapped pairing for each SNP (Which can map to itself or a friend)
    df[!,:A1_map] = GetAllelePair(df[:,:A1])
    df[!,:A2_map] = GetAllelePair(df[:,:A2])

    df[!, :BETA_STD] .= 0.0

    # create a normal distribution object to be used later

    Ndist = Normal(0,1)

    # convert betas from logistic GWA
    if "OR" in names(df)
        df[!, :BETA] = log.(df[:,:OR])
    end

    #for each SNP/(A1,A2) pairing in the sumstats file, see if it appears in our list of
    # valid SNP choices
    # would just be much faster to see if the sst snps are in the valiation dict
    # this is to get flipped SNPs (probably not the most efficient way)
    in_comm1 = Tuple.(
        eachrow(DataFrame(SNP = df[:,:SNP], A1 = df[:,:A1], A2 = df[:,:A2]))
        ) .∈ Ref(comm_snp)

    in_comm2 = Tuple.(
        eachrow(DataFrame(SNP = df[:,:SNP], A1 = df[:,:A1_map], A2 = df[:,:A2_map]))
        ) .∈ Ref(comm_snp)

    ##standardize beta coefficients--z score/sqrt(n_subjects)

    common_snp_not_flipped = in_comm1 .| in_comm2

    df[common_snp_not_flipped, :BETA_STD] = (
        sign.(df[common_snp_not_flipped, :BETA])
        ) .* abs.(
            quantile.(Ndist, max.(df[common_snp_not_flipped, :P]./2, 1e-323))
            ) ./ n_sqrt

    #standardize beta coefficients for flipped paits-- (-1) *z score/sqrt(n_subjects)
    in_comm3 = Tuple.(
        eachrow(DataFrame(SNP = df[:,:SNP], A1 = df[:,:A2], A2 = df[:,:A1]))
        ) .∈ Ref(comm_snp)

    in_comm4 = Tuple.(
        eachrow(DataFrame(SNP = df[:,:SNP], A1 = df[:,:A2_map], A2 = df[:,:A1_map]))
        ) .∈ Ref(comm_snp)

    common_snp_flipped = in_comm3 .| in_comm4

    df[common_snp_flipped, :BETA_STD] = (
        -1 .* sign.(df[common_snp_flipped, :BETA])
        ) .* abs.(
            quantile.(Ndist, max.(df[common_snp_flipped, :P]./2, 1e-323))
        ) ./ n_sqrt

    #subset to valid SNPs
    df = df[in_comm1 .| in_comm4 .| in_comm2 .| in_comm3, :]

    #subset the reference LD dict to valid SNPs
    in_sst = ref_dict[:,:SNP] .∈ Ref(df[:,:SNP])

    #get the information from the reference dict and join it to the effect sizes
    sst_dict = innerjoin(
        select(
            ref_dict, 
            :CHR, :SNP, :BP,:A1,:A2, :MAF
            ), 
        select(
            df, 
            :SNP, :BETA_STD
            ), 
        on = :SNP
        )

    sst_dict[!,:A1_map] = GetAllelePair(sst_dict[:,:A1])
    sst_dict[!,:A2_map] = GetAllelePair(sst_dict[:,:A2])

    sst_dict[!,:tmp] = copy(sst_dict[:,:A2])
    sst_dict[!, :FLP] .= 0

    #detect flipped allelle sets
    in_comm1 = Tuple.(
        eachrow(
            DataFrame(
                SNP = sst_dict[:,:SNP], 
                A1 = sst_dict[:,:A1], 
                A2 = sst_dict[:,:A2]
                )
            )
        ) .∈ Ref(comm_snp)
    
    in_comm2 = Tuple.(
        eachrow(
            DataFrame(
                SNP = sst_dict[:,:SNP], 
                A1 = sst_dict[:,:A1_map], 
                A2 = sst_dict[:,:A2_map]
            )
        )
    ) .∈ Ref(comm_snp)

    in_comm3 = Tuple.(
        eachrow(
            DataFrame(
            SNP = sst_dict[:,:SNP], 
                A1 = sst_dict[:,:A2], 
                A2 = sst_dict[:,:A1]
                )
            )
        ) .∈ Ref(comm_snp)

    in_comm4 = Tuple.(
        eachrow(
            DataFrame(
                SNP = sst_dict[:,:SNP], 
                A1 = sst_dict[:,:A2_map], 
                A2 = sst_dict[:,:A1_map]
            )
        )
    ) .∈ Ref(comm_snp)


    # summary statistics with non-flipped alleles
    sst_dict[in_comm2,:A1] = sst_dict[in_comm2,:A1_map]
    sst_dict[in_comm2,:A2] = sst_dict[in_comm2,:A2_map]
    
    sst_dict[in_comm1 .| in_comm2, :FLP] .= 1

    # summary statistics with flipped alleles
    sst_dict[in_comm3,:tmp] = sst_dict[in_comm3,:A1]
    sst_dict[in_comm3,:A1] = sst_dict[in_comm3,:A2]
    sst_dict[in_comm3,:A2] = sst_dict[in_comm3,:tmp]

    sst_dict[in_comm4,:A1] = sst_dict[in_comm2,:A2_map]
    sst_dict[in_comm4,:A2] = sst_dict[in_comm2,:A1_map]

    sst_dict[in_comm3 .| in_comm4, :FLP] .= -1

    # set MAF info for flupped alleles
    sst_dict[in_comm3 .| in_comm4, :MAF] = 1. .- sst_dict[in_comm3 .| in_comm4, :MAF]

    # output standardized beta
    sst_dict[!,:BETA] = sst_dict[:, :BETA_STD]

    ret_df = select(sst_dict, :CHR, :SNP, :BP, :A1, :A2, :MAF, :BETA, :FLP)
    return ret_df
end

"""
    ParseLDBlock(ldblk_dir, sst_dict, chrom)

This function parses a GWAS summary statistics file.
...
# Arguments
- `ldblk_dir`: LD Reference file 
- `sst_dict`: output of ParseSumstats
- `chrom`: which chromosome
...
# Returns
- `ld_blk`: LD matrix for SNPs appearing in the GWAS sumstat file
- `blk_size`: number of SNPs in each block
...
"""
function ParseLDBlock(ldblk_dir, sst_dict, chrom)

    print("... parse reference LD on chromosome $chrom ...")

    if occursin("1kg", ldblk_dir)
        chr_name = "$ldblk_dir/ldblk_1kg_chr$chrom.hdf5"
    end

    if occursin("ukbb", ldblk_dir)
        chr_name = "$ldblk_dir/ldblk_ukbb_chr$chrom.hdf5"
    end


    hdf_chr = h5open(chr_name, "r")

    n_blk = length(hdf_chr)

    # read LD blocks
    ld_blk = [
        read(hdf_chr["blk_$i"]["ldblk"]) for i in 1:n_blk
        ]

        # read the snplist for ldblocks
    snp_blk = [
        read(hdf_chr["blk_$i"]["snplist"]) for i in 1:n_blk
        ]

    # get indexes of snps in the snplist that appear in our sumstats 
    idx_blks = [
            [
                j for j in 1:length(i) if i[j] in sst_dict[:,:SNP]
            ] for i in snp_blk
        ]

    blk_size = length.(idx_blks)

    mm = 0

    for i in 1:n_blk
        if idx_blks[i] != []
            
            flip = sst_dict[(mm+1):(mm + blk_size[i]), :FLP]

            #subset block to SNPs that appear, and set negative correlations for flipped alleles
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
