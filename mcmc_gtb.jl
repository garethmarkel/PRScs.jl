using Random, Distributions, GenInvGaussian

function mcmc(a, b, phi, sst_dict, n, ld_blk, blk_size, n_iter, n_burnin, thin, chrom, out_dir, beta_std, seed)
    print("... MCMC ...")

    ##un hard code this! TODO
    Random.seed!(seed)

    beta_mrg = sst_dict[:,:BETA]
    maf = sst_dict[:,:MAF]

    n_pst = (n_iter-n_burnin)/thin

    p = nrow(sst_dict)
    n_blk = length(ld_blk)

    # initialization
    beta = zeros(Float64,p,1)
    psi = ones(Float64,p,1)
    sigma = 1.0

    # if phi == None:
    #     phi = 1.0; phi_updt = True
    # else:
    #     phi_updt = False

    phi_updt = true

    beta_est = zeros(Float64,p,1)
    psi_est = ones(Float64,p,1)
    sigma_est = 0.0
    phi_est = 0.0

    #create normal dist
    ndist = Normal(0,1)

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

                dinvt = ld_blk[kk] .+ Diagonal(1 ./ psi[idx_blk])
                dinvt_chol = cholesky(Hermitian(dinvt))

                #Here, we use the cholesky decomposition to (Relatively) efficiently invert
                #the matrix (D+ψ⁻¹) to get Σ = (D+ψ⁻¹)⁻¹ and sample in one fell swoop
                #this gets more useful with big LD blocks (or if you want to fit this on like a whole chromosome)
                beta_tmp = dinvt_chol.L\beta_mrg[idx_blk] .+  (sqrt(sigma/n) .* rand(ndist,blk_size[kk]))
                beta[idx_blk] = dinvt_chol.U\beta_tmp

                #assume the LD blocks are independent, and iterate on the sum
                quad = quad + beta[idx_blk]'*dinvt*beta[idx_blk]
                mm = mm + blk_size[kk]
            end
        end

        #generate the variance for sigma^2
        err = max(n/2.0 * (1.0 - 2.0*sum(beta.*beta_mrg) + quad), n/2.0*sum((beta.^2)./psi))

        #update sigma^2
        #prior versions sampled from a gamma(a,b^-1) but let's not do that anymore
        sigma = rand(InverseGamma((n+p)/2.0, err))

        delta = rand.(Gamma.(a+b, 1.0 ./ (psi .+ phi)),1)
        delta = [i[1] for i in delta]

        #to run this, you may have to pull the fork https://github.com/garethmarkel/GenInvGaussian.jl
        #again, this is not my original work--the original repo was for an older version of julia
        #so I just updated the "paperwork" to get back to working on reimplementing PRScs
        psi_gig = GeneralizedInverseGaussian.(a-0.5, 2 .* delta, n .* (beta.^2) ./ sigma)

        ##TODO: fix this--should be able to broadcast it but the rand function in the implementation
        ##of the generlaized inverse gaussian function does one at a time
        for jj in 1:length(psi)
            psi[jj] = rand(psi_gig[jj])
        end
        psi[psi .> 1] .= 1.0

        if phi_updt == true
            w = rand(Gamma(1.0, 1.0/(phi+1.0)))
            phi = rand(Gamma(p*b+0.5, 1.0/(sum(delta)+w)))
        end

        # posterior
        if itr>n_burnin & itr % thin == 0
            beta_est = beta_est + beta/n_pst
            psi_est = psi_est + psi/n_pst
            sigma_est = sigma_est + sigma/n_pst
            phi_est = phi_est + phi/n_pst
        end

    end

    if beta_std == false
        beta_est = beta_est ./ sqrt.(2.0 .* maf .* (1.0 .- maf))
    end


    aw = round(a, digits = 1)
    bw = round(b, digits = 1)

    if phi_updt == true

        eff_file = "$out_dir-pst_eff_a$aw-b$bw-phiauto-chr$chrom.txt"
    else
        phiw - round(phi, digits = 1)
        eff_file = "$out_dir-pst_eff_a$aw-b$bw-phi$phiw-chr$chrom.txt" % (a, b, phi, chrom)
    end

    open(eff_file, "w") do file
        for i in 1:length(beta_est)
            snpw = sst_dict[:,:SNP][i]
            bpw = sst_dict[:,:BP][i]
            a1w = sst_dict[:,:A1][i]
            a2w = sst_dict[:,:A2][i]
            betaw = round(beta_est[i], digits = 6)
            write(file, "$chrom\t$snpw\t$bpw\t$a1w\t$a2w\t$betaw\n")
        end
    end

    if phi_updt == true
        print("... Estimated global shrinkage parameter: $phi_est ..." )
    end

    print("... Done ...")
end
