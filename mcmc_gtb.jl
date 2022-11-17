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

function set_dict(model, p,a,b,phi)
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

function update_hyperprior!(param_dict, hyperprior_update)
    p = length(param_dict["beta"])
    if param_dict["model"] == "strawderman"
        if hyperprior_update == true
            param_dict["w"] = rand(Gamma(1.0, 1.0/(param_dict["phi"]+1.0)))
            param_dict["phi"] = rand(Gamma(p*param_dict["b"]+0.5, 1.0/(sum(param_dict["delta"])+param_dict["w"])))
        end
    end
end

#could also implement this: https://arxiv.org/pdf/1506.04778.pdf
function update_beta(param_dict, dinvt, beta_mrg, idx_blk, n)

    dinvt_chol = cholesky(Hermitian(dinvt))

    block_size = length(beta_mrg[idx_blk])

    #Here, we use the cholesky decomposition to (Relatively) efficiently invert
    #the matrix (D+ψ⁻¹) to get Σ = (D+ψ⁻¹)⁻¹ and sample in one fell swoop
    #this gets more useful with big LD blocks (or if you want to fit this on like a whole chromosome)
    beta_tmp = dinvt_chol.L\beta_mrg[idx_blk] .+  (sqrt(param_dict["sigma"]/n) .* rand(Normal(0,1),block_size))
    beta_tmp = dinvt_chol.U\beta_tmp

    #in the DL formulation, there's a step where we have to divide by beta
    #TODO: reparametrize the tau^-2 sampling so this isn't an issue
    beta_tmp[beta_tmp .== 0.0] .= 1e-16
    return beta_tmp

end

#note--max not NaN tolerant
function update_global_scale!(param_dict, beta_mrg, n, quad)

    # #generate the variance for sigma^2
    p = length(param_dict["beta"])

    if param_dict["model"] == "strawderman"

        #this is to deal with numerical issues
        if isnan(quad) | isinf(quad)
            err = n/2.0*sum((param_dict["beta"].^2)./param_dict["psi"])
        else
            err = max(n/2.0 * (1.0 - 2.0*sum(param_dict["beta"].*beta_mrg) + quad), n/2.0*sum((param_dict["beta"].^2)./param_dict["psi"]))
        end

        #update sigma^2
        param_dict["sigma"] = rand(InverseGamma((n+p)/2.0, err))

    elseif param_dict["model"] =="dirichlet-laplace"

        #this is to deal with numerical issues
        if isnan(quad) | isinf(quad)
            err = n/2.0*sum((param_dict["beta"].^2)./param_dict["psi"])
        else
            err = max(n/2.0 * (1.0 - 2.0*sum(param_dict["beta"].*beta_mrg) + quad), n/2.0*sum((param_dict["beta"].^2)./param_dict["psi"]))
        end

        #update sigma^2
        param_dict["sigma"] = rand(InverseGamma((n+p)/2.0, err))

    else
         println("placeholder")
    end
end

#problem: this does all the modification *inside* the function--should be able to run update_local_scale!(model, param_dict, n) but need to test
function update_local_scale!(param_dict, n)

    p = length(param_dict["beta"])

    if param_dict["model"] == "strawderman"
        delta = rand.(Gamma.(param_dict["a"]+param_dict["b"], 1.0 ./ (param_dict["psi"] .+ param_dict["phi"])),1)
        param_dict["delta"] = [i[1] for i in delta]

        psi_gig = GeneralizedInverseGaussian.(param_dict["a"]-0.5, 2 .* param_dict["delta"], n .* (param_dict["beta"].^2) ./ param_dict["sigma"])

        ##TODO: fix this--should be able to broadcast it but the rand function in the implementation
        ##of the generlaized inverse gaussian function does one at a time
        for jj in 1:length(param_dict["psi"])
            param_dict["psi"][jj] = rand(psi_gig[jj])
        end
        param_dict["psi"][param_dict["psi"] .> 1] .= 1.0

    elseif param_dict["model"] == "dirichlet-laplace"

        psi_gig = GeneralizedInverseGaussian.(param_dict["a"] - 1, 1.0, 2 .* abs.(param_dict["beta"]) ./ sqrt(param_dict["sigma"])  )

        for jj in 1:length(param_dict["psi"])
            param_dict["psi"][jj] = rand(psi_gig[jj])
        end

        param_dict["psi"] = param_dict["psi"] ./ sum(param_dict["psi"])

        #deal with numerical issues
        param_dict["psi"][param_dict["psi"] .== 0.0] .= 0.0001

        param_dict["lambda"] = rand(GeneralizedInverseGaussian(p*param_dict["a"]-p,1.0,2*sum(abs.(param_dict["beta"])./ param_dict["psi"])/sqrt(param_dict["sigma"])))

        #TODO: figure out how to reparametrize this
        tau_inv = rand.(InverseGaussian.(param_dict["lambda"] .* param_dict["psi"] .* sqrt(param_dict["sigma"]) ./ abs.(param_dict["beta"]),1.0),1)

        #tau_inv = ones(Float64,p,1)
        param_dict["tau"] = 1.0 ./ [i[1] for i in tau_inv]

    else
         println("placeholder")
    end
end

function get_shrinkage_mat(idx_blk, param_dict)
    if param_dict["model"] == "strawderman"
        return Diagonal(1 ./ param_dict["psi"][idx_blk])
    elseif param_dict["model"] == "dirichlet-laplace"
        return Diagonal(1 ./ ((param_dict["psi"][idx_blk].^2) .* (param_dict["lambda"]^2) .* param_dict["tau"][idx_blk]))
    else
         println("placeholder")
    end
end



function mcmc(model, a, b, phi, sst_dict, n, ld_blk, blk_size, n_iter, n_burnin, thin, chrom, out_dir, beta_std, seed)

    print("... MCMC ...")

    Random.seed!(seed)

    beta_mrg = sst_dict[:,:BETA]
    maf = sst_dict[:,:MAF]

    n_pst = (n_iter-n_burnin)/thin

    p = nrow(sst_dict)
    n_blk = length(ld_blk)

    param_dict = set_dict(model, p, a,b,phi)

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

                shrinkage_mat = get_shrinkage_mat(idx_blk, param_dict)

                dinvt = ld_blk[kk] .+ shrinkage_mat

                param_dict["beta"][idx_blk] = update_beta(param_dict, dinvt, beta_mrg, idx_blk, n)

                #assume the LD blocks are independent, and iterate on the sum
                quad = quad + param_dict["beta"][idx_blk]'*dinvt*param_dict["beta"][idx_blk]

                mm = mm + blk_size[kk]
            end
        end

        update_global_scale!(param_dict, beta_mrg, n, quad)

        update_local_scale!(param_dict, n)

        update_hyperprior!(param_dict, phi_updt)

        # posterior
        if itr>n_burnin & itr % thin == 0
            # beta_est = beta_est + beta/n_pst
            param_dict["beta_est"] = param_dict["beta_est"] .+ param_dict["beta"]./n_pst
            # psi_est = psi_est + psi/n_pst
            # sigma_est = sigma_est + sigma/n_pst
            # phi_est = phi_est + phi/n_pst
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
