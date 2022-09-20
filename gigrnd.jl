# """
# Random variate generator for the generalized inverse Gaussian distribution.
# Reference: L Devroye. Random variate generation for the generalized inverse Gaussian distribution.
#            Statistics and Computing, 24(2):239-246, 2014.
# TODO: convert these to Distributions.jl
# """

using Distributions

function fnpsi(x, alpha, lam)
    f = -alpha*(cosh(x)-1.0)-lam*(exp(x)-x-1.0)
    return f
end


function dpsi(x, alpha, lam)
    f = -alpha*sinh(x)-lam*(exp(x)-1.0)
    return f
end

function g(x, sd, td, f1, f2)
    if x >= -sd & x <= td
        f = 1.0
    elseif x > td
        f = f1
    else x < -sd
        f = f2
    end

    return f
end


function gigrnd(p, a, b)
    # setup -- sample from the two-parameter version gig(lam,omega)

    lam = p

    omega = sqrt(a*b)

    if lam < 0
        lam = -lam
        swap = true
    else
        swap = false
    end

    alpha = sqrt(omega^2 + lam^2)-lam

    # find t
    x = -psi(1.0, alpha, lam)
    if x >= 0.5 & x <= 2.0
        t = 1.0
        s = 1.0
    elseif x > 2.0
        if alpha == 0 & lam == 0
            t = 1.0
            s = 1.0
        else
            t = sqrt(2.0/(alpha+lam))
            s = sqrt(4.0/(alpha*cosh(1)+lam))
        end
    elseif x < 0.5
        if alpha == 0 & lam == 0
            t = 1.0
            s = 1.0
        else
            t = log(4.0/(alpha+2.0*lam))
            if alpha == 0
                s = 1.0/lam
            elseif lam == 0
                s = log(1.0+1.0/alpha+sqrt(1.0/(alpha^2)+2.0/alpha))
            else
                s = min(1.0/lam, log(1.0+1.0/alpha+sqrt(1.0/(alpha^2)+2.0/alpha)))
            end
        end
    end

    # find auxiliary parameters
    eta = -fnpsi(t, alpha, lam)
    zeta = -dpsi(t, alpha, lam)
    theta = -fnpsi(-s, alpha, lam)
    xi = dpsi(-s, alpha, lam)

    p = 1.0/xi
    r = 1.0/zeta

    td = t-r*eta
    sd = s-p*theta
    q = td+sd

    # random variate generation
    while True
        U = rand(Uniform(),1)
        V = rand(Uniform(),1)
        W = rand(Uniform(),1)
        if U < q/(p+q+r)
            rnd = -sd+q*V
        elseif U < (q+r)/(p+q+r)
            rnd = td-r*log(V)
        else
            rnd = -sd+p*log(V)
        end

        f1 = exp(-eta-zeta*(rnd-t))
        f2 = exp(-theta+xi*(rnd+s))
        if W*g(rnd, sd, td, f1, f2) <= exp(fnpsi(rnd, alpha, lam))
            break
        end
    end

    # transform back to the three-parameter version gig(p,a,b)
    rnd = exp(rnd)*(lam/omega+sqrt(1.0+(lam^2)/(omega^2)))
    if swap
        rnd = 1.0/rnd
    end

    rnd = rnd/sqrt(a/b)
    return rnd
end
