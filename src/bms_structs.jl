# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
@concrete struct BMSBase
    n_chains        # number of chains
    n_samples       # number of samples per chain
    n_thinning      # thinning period
    n_burnin        # lenght of burnout
    n_change        # number of changes in M
    n_scale         # n_scale of P_r|m
    ϵ               # ϵ for P_r|m
    uniform_initial # whether initial sample is unofrm
end
function BMSBase(; n_chains = 40, n_samples = 100000,
        n_thinning = 50, n_burnin = n_thinning, n_change = 1, n_scale = 1., ϵ = 1., uniform_initial = false)
    BMSBase(n_chains, n_samples, n_thinning, n_burnin,  n_change, n_scale, ϵ, uniform_initial)
end
export BMSBase

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
@concrete struct BMSBorSampNaive
    n_samples       # number of samples
end
function BMSBorSampNaive(; n_samples = 1000000)
    BMSBorSampNaive(n_samples)
end
export BMSBorSampNaive

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
@concrete struct BMSObj
    l_matrix        # size = n_subject * n_model
    n_subject       # number of subjects
    n_model         # number of models
    α               # prior parameters
    bms_basechain   # base chain information
    bms_borsamp     # bor sampling information
end
function BMSObj(l_matrix; 
        n_subject=size(l_matrix)[1], n_model=size(l_matrix)[2], α=1/n_model, 
        bms_basechain=BMSBase(), bms_borsamp=BMSBorSampNaive())
    BMSObj(l_matrix,n_subject,n_model,α,bms_basechain,bms_borsamp)
end
export BMSObj

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
@concrete struct BMSStats
    exp_M           # P(m_n | data)
    exp_r           # E[r | data]
    d_exp_r         # std[r | data]
    xp              # excedence prob
    bor             # BOR
    pxp             # protected excedence prob
end
export BMSStats

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
@concrete struct BMSResults
    object
    stats
    posterior_results
    bor_results
end
export BMSResults
