using DataFrames, CSV

##read in the reference file and subset to a specific chromosome
function parse_ref(ref_file, chrom)
    print("... parse reference file: $ref_file ...")

    df = DataFrame(CSV.File(ref_file))

    ref_dict = filter(:CHR => n -> n == chrom, df)

    n_snp = nrow(ref_dict)

    print("... $n_snp SNPs on chromosome $chrom read from $ref_file ...")

    return ref_dict

end

ref_dict = parse_ref("/Users/garethmarkel/Documents/ldblk_1kg_eur/snpinfo_1kg_hm3", 22)

function parse_bim(bim_file, chrom)
    print("... parse bim file: $bim_file.bim ...")

    df = DataFrame(CSV.File(bim_file, header = false))
    rename!(df, [:CHR, :SNP, :POS ,:BP, :A1 , :A2])

    vld_dict = select(filter(:CHR => n -> n == chrom, df), :SNP, :A1,:A2)

    n_snp = nrow(vld_dict)

    print("... $n_snp SNPs on chromosome $chrom read from $bim_file.bim ...")
    return vld_dict
end


vld_dict = parse_bim("/Users/garethmarkel/Documents/PRScs/test_data/test.bim",22)


df = DataFrame(CSV.File("/Users/garethmarkel/Documents/PRScs/test_data/sumstats.txt"))

sst_dict = subset(filter([:A1, :A2] => valid_allelles, df), :SNP, :A1, :A2)

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

print("... $n_snp common SNPs in the reference, sumstats, and validation set ...")



function valid_allelles(a1, a2)::Bool
    valid_a1 = a1 in ["A","C","G","T"]
    valid_a2 = a2 in ["A","C","G","T"]
    valid_a1 && valid_a2
end

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


end
