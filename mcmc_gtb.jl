using Random, Distributions, GenInvGaussian

#steps:
# 1. strawderman
#   a. beta
#   b. sigma
#   c. delta
#   d. psi
#   e. phi (sometimes)
# 2. DL
#   a. beta
#   b. sigma

"""
    SetParamDict(model, p,a,b,phi)

This function creates a dictionary to store parameters for the model
...
# Arguments
- `model`: Name of model, one of strawderman or dirichlet-laplace
- `p`: dimensionality
- `a`: a parameter for strawderman berger
- `b`: b parameter for strawderman berger
- `phi`: phi parameter for strawderman berger
...
# Returns
- `prior_params`: dictionary of parameters for the appropriate prior
...
"""
function SetParamDict(model, p,a,b,phi)
    if model == "strawderman"
        prior_params = Dict(
            "model" => "strawderman",
            "beta" => zeros(Float64,p,1),
            "beta_est" => zeros(Float64,p,1),
            "psi" => ones(Float64,p,1),
            "psi_est" => zeros(Float64,p,1),
            "delta" => ones(Float64,p,1),
            "delta_est" => zeros(Float64,p,1),
            "phi" => phi,
            "w" => 1.0,
            "sigma" => 1.0,
            "a" => a, #TODO: fix this placeholder
            "b" => b #TODO: fix this placeholder
        )
    elseif model =="dirichlet-laplace"
        prior_params = Dict(
            "model" => "dirichlet-laplace",
            "beta" => zeros(Float64,p,1),
            "beta_est" => zeros(Float64,p,1),
            "psi" => ones(Float64,p,1),
            "psi_est" => zeros(Float64,p,1),
            "lambda" => 1.0,
            "lambda_est" => 0.0,
            "tau" => ones(Float64,p,1),
            "tau_est" => zeros(Float64,p,1),
            "sigma" => 1.0,
            "sigma_est" => 0.0,
            "a" => 1.0/p #TODO: fix this placeholder
        )
    else
         println("placeholder")
    end

    return prior_params
end


"""
    UpdateHyperprior!(param_dict, hyperprior_update)

This function updates the hyperpriors for the strawderman berger model
...
# Arguments
- `param_dict`: dictionary of parameters for the appropriate prior
- `update_hyperprior`: boolean to run update
...
"""
function UpdateHyperprior!(param_dict, hyperprior_update)
    p = length(param_dict["beta"])

    if param_dict["model"] == "strawderman"
        if hyperprior_update == true
            param_dict["w"] = rand(
                Gamma(
                    1.0, 
                    1.0/(param_dict["phi"]+1.0)
                    )
                )
            
            param_dict["phi"] = rand(
                Gamma(
                    p*param_dict["b"]+0.5, 
                    1.0/(sum(param_dict["delta"])+param_dict["w"])
                    )
                )
        end
    end
end

"""
    UpdateBeta(param_dict, dinvt, beta_mrg, idx_blk, n)

This function updates the coefficients for the global-local cs models
...
# Arguments
- `param_dict`: dictionary of parameters for the appropriate prior
- `dinvt`: shrunk correlation matrix of snp locations
- `beta_mrg`: fitted betas from GWAS (normalized)
- `idx_blk`: indices of betas to be used for the model
- `n`: gwas size 
...
# Returns
- `beta_tmp`: newly smapled betas
...
"""
function UpdateBeta(param_dict, dinvt, beta_mrg, idx_blk, n)

    dinvt_chol = cholesky(Hermitian(dinvt))

    block_size = length(beta_mrg[idx_blk])

    # sample β = MVN(μ = beta_mrg*(N/σ^2)(D+ψ⁻¹)⁻¹, Σ = ((σ^2)/N)*(D+ψ⁻¹)⁻¹
    beta_tmp = dinvt_chol.L\beta_mrg[idx_blk] .+  (sqrt(param_dict["sigma"]/n) .* rand(Normal(0,1),block_size))
    beta_tmp = dinvt_chol.U\beta_tmp

    #in the DL formulation, there's a step where we have to divide by beta
    beta_tmp[beta_tmp .== 0.0] .= 1e-16

    return beta_tmp

end

"""
    UpdateGlobalScale!(param_dict, beta_mrg, n, quad)

This function updates the shared variance of beta parameters
...
# Arguments
- `param_dict`: dictionary of parameters for the appropriate prior
- `beta_mrg`: fitted betas from GWAS (normalized)
- `n`: gwas size 
- `quad`: quadratic term β'(D+ψ⁻¹)β
...
"""
function UpdateGlobalScale!(param_dict, beta_mrg, n, quad)

    # #generate the variance for sigma^2
    p = length(param_dict["beta"])

    if param_dict["model"] == "strawderman"

        #this is to deal with numerical issues
        if isnan(quad) | isinf(quad)
            err = n/2.0*sum((param_dict["beta"].^2)./param_dict["psi"])
        else
            err = max(
                n/2.0 * (
                    1.0 - 2.0 * sum(param_dict["beta"] .* beta_mrg) + quad
                    ), 
                n/2.0 * sum((param_dict["beta"].^2)./param_dict["psi"])
                )
        end


        #σ^2 ∼ iG((n+p)/2, (n/2)*(1 - 2*β'*beta_mrg + β'(D+ψ⁻¹)β))

        param_dict["sigma"] = rand(InverseGamma((n+p)/2.0, err))

    elseif param_dict["model"] =="dirichlet-laplace"

        #this is to deal with numerical issues
        if isnan(quad) | isinf(quad)
            err = n / 2.0 * sum((param_dict["beta"].^2)./param_dict["psi"])
        else
            err = max(
                n / 2.0 * (
                    1.0 - 2.0*sum(param_dict["beta"] .* beta_mrg) + quad
                    ), 
                n / 2.0 * sum((param_dict["beta"].^2) ./ param_dict["psi"])
                )
        end

        #σ^2 ∼ iG((n+p)/2, (n/2)*(1 - 2*β'*beta_mrg + β'(D+ψ⁻¹)β))
        param_dict["sigma"] = rand(InverseGamma((n+p)/2.0, err))

    else
         println("placeholder")
    end
end

"""
    UpdateLocalScale!(param_dict, n)

This function updates the local part of the parameter variance
...
# Arguments
- `param_dict`: dictionary of parameters for the appropriate prior
- `n`: gwas size 
...
# Returns
- `beta_tmp`: newly smapled betas
...
"""
function UpdateLocalScale!(param_dict, n)

    p = length(param_dict["beta"])

    if param_dict["model"] == "strawderman"
        
        delta = rand.(
            Gamma.(
                param_dict["a"]+param_dict["b"], 
                1.0 ./ (param_dict["psi"] .+ param_dict["phi"])
                ),
            1
            )

        # δ_i ∼ G(a+b, ψ_i + ϕ)

        param_dict["delta"] = [i[1] for i in delta]

        # ψ_i = giG(a - 0.5, 2δ_i, (n/σ^2)*β_i)

        psi_gig = GeneralizedInverseGaussian.(
            param_dict["a"]-0.5, 
            2 .* param_dict["delta"], 
            n .* (param_dict["beta"].^2) ./ param_dict["sigma"]
            )

        for jj in 1:length(param_dict["psi"])
            param_dict["psi"][jj] = rand(psi_gig[jj])
        end

        param_dict["psi"][param_dict["psi"] .> 1] .= 1.0

    elseif param_dict["model"] == "dirichlet-laplace"


        # τ_i ∼ iG(λ_i * ψ_i * σ / β_i, 1.0)

        tau_inv = rand.(
            InverseGaussian.(
                param_dict["lambda"] .* param_dict["psi"] .* (
                    sqrt(param_dict["sigma"]) ./ abs.(param_dict["beta"])
                    ),
                1.0),
            1
            )

        #tau_inv = ones(Float64,p,1)
        param_dict["tau"] = 1.0 ./ [i[1] for i in tau_inv]


        # λ = giG(pa - p, 1.0, 2*Sum(abs(β_i/ψ_i))/σ

        param_dict["lambda"] = rand(
            GeneralizedInverseGaussian(
                p*param_dict["a"]-p,
                1.0,
                2*sum(abs.(param_dict["beta"])./ param_dict["psi"])/sqrt(param_dict["sigma"])
                )
            )

        # ψ_i = giG(a - 1, 1.0, 2*|β_i|/σ)

        psi_gig = GeneralizedInverseGaussian.(
            param_dict["a"] - 1, 
            1.0, 
            2 .* abs.(param_dict["beta"]) ./ sqrt(param_dict["sigma"])  
            )

        for jj in 1:length(param_dict["psi"])
            param_dict["psi"][jj] = rand(psi_gig[jj])
        end

        param_dict["psi"] = param_dict["psi"] ./ sum(param_dict["psi"])

        #deal with numerical issues
        param_dict["psi"][param_dict["psi"] .== 0.0] .= 1e-16

        

    else
         println("placeholder")
    end
end



"""
    GetShrinkageMat(idk_blk, param_dict)

This function gets the shrinkagemat for each LD block
...
# Arguments
- `idx_blk`: indices of betas to be used for the model
- `param_dict`: dictionary of parameters for the appropriate prior
...
# Returns
- shrinkage matrix for LD block
...
"""
function GetShrinkageMat(idx_blk, param_dict)
    
    if param_dict["model"] == "strawderman"

        # diagonal values are given 1/ψ_i
        
        return Diagonal(1 ./ param_dict["psi"][idx_blk])
    
    elseif param_dict["model"] == "dirichlet-laplace"
        
        # diagonal values are given 1/((ψ_i^2)*(λ_i^2)*(τ_i))

        return Diagonal(
            1 ./ (
                (param_dict["psi"][idx_blk].^2) .* (param_dict["lambda"]^2) .* param_dict["tau"][idx_blk]
                )
            )
    
    else
         println("placeholder")
    end
end



"""
    ContinuousShrinkageGibbsMCMC(
    model, a, b, phi, sst_dict, n, ld_blk, 
    blk_size, n_iter, n_burnin, thin, 
    chrom, out_dir, beta_std, seed)

This runs the MCMC and writes the output to a file
...
# Arguments
- `model`: string wiht model name
- `a`: a parameter for SB
- `b`: b parameter for SB
- `phi` phi parameter for SB/DL
- `sst_dict`: summary statistics from GWAS, after being cleaned by ParseSumstats
- `n`: gwas size
- `ld_blk`: vector of LD blocks parsed to only have snps that appear in our reference
    bim file and gwas
- `blk_size`: size of ld blocks
- `n_iter`: how many iterations?
- `n_burnin`: how many initial iterations to discard
- `thin`: how often to throw out iterations
- `chrom`: which chromosome
- `out_dir`: where should output be written
- `beta_std` boolean to standardize betas
- `seed`: random seed
...
"""
function ContinuousShrinkageGibbsMCMC(
    model, a, b, phi, sst_dict, n, ld_blk, 
    blk_size, n_iter, n_burnin, thin, 
    chrom, out_dir, beta_std, seed)

    print("... MCMC ...")

    Random.seed!(seed)

    beta_mrg = sst_dict[:,:BETA]
    maf = sst_dict[:,:MAF]

    n_pst = (n_iter-n_burnin)/thin

    p = nrow(sst_dict)
    n_blk = length(ld_blk)

    param_dict = SetParamDict(model, p, a,b,phi)

    phi_updt = true

    #MCMC
    for itr in 1:n_iter

        mm = 0
        quad = 0.0

        ##here we update beta
        for kk in 1:n_blk

            if blk_size[kk] == 0
                continue
            end
            if blk_size[kk] != 0
                idx_blk = (mm+1):(mm + blk_size[kk])

                shrinkage_mat = GetShrinkageMat(idx_blk, param_dict)

                dinvt = ld_blk[kk] .+ shrinkage_mat

                param_dict["beta"][idx_blk] = UpdateBeta(param_dict, dinvt, beta_mrg, idx_blk, n)

                #assume the LD blocks are independent, and iterate on the sum
                quad = quad + param_dict["beta"][idx_blk]'*dinvt*param_dict["beta"][idx_blk]

                mm = mm + blk_size[kk]
            end
        end

        UpdateGlobalScale!(param_dict, beta_mrg, n, quad)

        UpdateLocalScale!(param_dict, n)

        UpdateHyperprior!(param_dict, phi_updt)

        # posterior
        if itr>n_burnin & itr % thin == 0
            param_dict["beta_est"] = param_dict["beta_est"] .+ param_dict["beta"]./n_pst
        end

    end

    if beta_std == false
        param_dict["beta_est"] = param_dict["beta_est"] ./ sqrt.(2.0 .* maf .* (1.0 .- maf))
    end


    aw = round(a, digits = 1)
    bw = round(b, digits = 1)

    if phi_updt == true
        eff_file = "$out_dir-$model-pst_eff_a$aw-b$bw-phiauto-chr$chrom.txt"
    else
        phiw - round(param_dict["phi"], digits = 1)
        eff_file = "$out_dir-$model-pst_eff_a$aw-b$bw-phi$phiw-chr$chrom.txt" % (a, b, phi, chrom)
    end

    touch(eff_file)
    open(eff_file, "w") do file
        for i in 1:length(param_dict["beta_est"])
            snpw = sst_dict[:,:SNP][i]
            bpw = sst_dict[:,:BP][i]
            a1w = sst_dict[:,:A1][i]
            a2w = sst_dict[:,:A2][i]
            betaw = round(param_dict["beta_est"][i], digits = 6)
            write(file, "$chrom\t$snpw\t$bpw\t$a1w\t$a2w\t$betaw\n")
        end
    end

    print("... Done ...")
end
