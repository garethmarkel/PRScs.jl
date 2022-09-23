using Random, Distributions, GenInvGaussian

#include("gigrnd.jl")


a = 1
b = 0.5
n = 200000
n_iter = 1000
n_burnin = 500
thin = 5
phi = 0.01

ld_blk = tst[1]
blk_size = tst[2]

Random.seed!(3)

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

            ##very not sure this is right--here's the original python:
            # dinvt = ld_blk[kk]+sp.diag(1.0/psi[(mm+1):(mm + blk_size[i])].T[0])
            dinvt = ld_blk[kk] .+ Diagonal(1 ./ psi[idx_blk])
            dinvt_chol = cholesky(Hermitian(dinvt))

            #we add in the randomness to create some sampling variation
            #TODO: determine why this is in the first command instead of the second
            beta_tmp = dinvt_chol.L\beta_mrg[idx_blk] .+  (sqrt(sigma/n) .* rand(ndist,blk_size[kk]))
            beta[idx_blk] = dinvt_chol.U\beta_tmp

            quad = quad + beta[idx_blk]'*dinvt*beta[idx_blk]
            mm = mm + blk_size[kk]
        end
    end

    #generate the variance for sigma^2
    err = max(n/2.0 * (1.0 - 2.0*sum(beta.*beta_mrg) + quad), n/2.0*sum((beta.^2)./psi))

    #update sigma^2
    sigma = 1.0 ./ rand(Gamma((n+p)/2.0, 1/err),1)
    sigma = sigma[1]

    delta = rand.(Gamma.(a+b, 1.0 ./ (psi .+ phi)),1)
    delta = [i[1] for i in delta]

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



function mcmc(a, b, phi, sst_dict, n, ld_blk, blk_size, n_iter, n_burnin, thin, chrom, out_dir, beta_std, seed)
    print("... MCMC ...")

    ##un hard code this! TODO
    Random.seed!(3)

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

                ##very not sure this is right--here's the original python:
                # dinvt = ld_blk[kk]+sp.diag(1.0/psi[(mm+1):(mm + blk_size[i])].T[0])
                dinvt = ld_blk[kk] .+ Diagonal(1 ./ psi[idx_blk])
                dinvt_chol = cholesky(Hermitian(dinvt))

                #we add in the randomness to create some sampling variation
                #TODO: determine why this is in the first command instead of the second
                beta_tmp = dinvt_chol.L\beta_mrg[idx_blk] .+  (sqrt(sigma/n) .* rand(ndist,blk_size[kk]))
                beta[idx_blk] = dinvt_chol.U\beta_tmp

                quad = quad + beta[idx_blk]'*dinvt*beta[idx_blk]
                mm = mm + blk_size[kk]
            end
        end

        #generate the variance for sigma^2
        err = max(n/2.0 * (1.0 - 2.0*sum(beta.*beta_mrg) + quad), n/2.0*sum((beta.^2)./psi))

        #update sigma^2
        sigma = 1.0 ./ rand(Gamma((n+p)/2.0, 1/err),1)
        sigma = sigma[1]

        delta = rand.(Gamma.(a+b, 1.0 ./ (psi .+ phi)),1)
        delta = [i[1] for i in delta]

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
