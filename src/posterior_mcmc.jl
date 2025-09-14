# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
function sample_posterior(BMS::BMSObj; base::BMSBase=BMS.bms_basechain)
    l_matrix = BMS.l_matrix
    n_model = BMS.n_model
    n_subject = BMS.n_subject

    n_chains = base.n_chains
    n_samples = base.n_samples
    
    ϵ=base.ϵ
    n_scale=base.n_scale
    n_change=base.n_change

    r_matrix = zeros(n_samples, n_model, n_chains);
    M_matrix = zeros(n_samples, n_subject, n_chains);
    π_matrix = zeros(n_samples, n_chains);

    P_r = Dirichlet(n_model, BMS.α)

    M0_star = zeros(n_subject)
    for i = 1:n_subject
            M0_star[i] = findmax(l_matrix[i,:])[2]
    end
    r0_star = [count(==(i),M0_star) for i = 1:n_model] ./ n_subject

    for i_chain = 1:n_chains
            if base.uniform_initial
                    r0,M0 = sample_prior(1.;n_model=n_model,n_subject=n_subject)
            else
                    r0,M0 = sample_basechain(r0_star,M0_star;
                                    n_scale=n_scale,ϵ=ϵ,n_change=n_change)
            end
            r_matrix[1,:,i_chain] = r0
            M_matrix[1,:,i_chain] = M0
            π_matrix[1,i_chain]   = logpdf_joint(r0,M0,l_matrix;P_r=P_r)
            for t = 2:n_samples
                    r1 = r_matrix[t-1,:,i_chain]
                    M1 = M_matrix[t-1,:,i_chain]
                    r2,M2 = sample_basechain(r1,M1;n_scale=n_scale,ϵ=ϵ,
                                                    n_change=n_change)
                    a12 = acceptprob(r1,M1,r2,M2,l_matrix;
                                    n_scale=n_scale,ϵ=ϵ,P_r=P_r)
                    if rand() < a12
                            r_matrix[t,:,i_chain] = r2
                            M_matrix[t,:,i_chain] = M2
                            π_matrix[t,i_chain] = logpdf_joint(r2,M2,l_matrix;P_r=P_r)
                    else
                            r_matrix[t,:,i_chain] = r1
                            M_matrix[t,:,i_chain] = M1
                            π_matrix[t,i_chain] = π_matrix[t-1,i_chain]
                    end
            end
    end

    sample_ids = (base.n_burnin):base.n_thinning:n_samples
    r_matrix = r_matrix[sample_ids,:,:];
    M_matrix = M_matrix[sample_ids,:,:];
    π_matrix = π_matrix[sample_ids,:];

    return (;r_matrix=r_matrix, M_matrix=M_matrix, π_matrix=π_matrix)
end
export sample_posterior

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
function sample_prior(α::Real;n_model=3,n_subject=10,H0=false)
        if H0
                r = ones(n_model) ./ n_model
        else
                P_r = Dirichlet(n_model, α)
                r = rand(P_r)
        end
        M = rand(Categorical(r),n_subject)
        return r,M
end
function sample_prior(r, α::Real; n_subject=10)
        M = rand(Categorical(r),n_subject)
        return r,M
end
export sample_prior

function sample_basechain(r, M; n_scale=10, ϵ=1., n_change = 1)
        M_new = deepcopy(M)
        P_unif = DiscreteUniform(1, length(M))
        for i = shuffle(1:length(M))[1:n_change]
                j = rand(P_unif)
                M_new[i] = M[j]
        end
        r_new = rand(Dirichlet(ϵ .+
                     [count(==(i),M) for i = 1:length(r)] ./ n_scale))
        return r_new, M_new
end
export sample_basechain

function logpdf_basechain(r1,M1,r2,M2;n_scale=10,ϵ=1.)
        if sum(abs.(M2 .- M1) .> 0) == 0
                log_Prob_m = 0 # This probability doesn't matter if M2 and M1 are equal
        else
                log_Prob_m = [0.]
                for i = (1:length(M1))[abs.(M2 .- M1) .> 0]
                        log_Prob_m[1] = log_Prob_m[1] +
                                        log(count(==(M2[i]),M1)) -
                                        log(length(M1))
                end
                log_Prob_m = log_Prob_m[1]
        end
        Prob_r = pdf(Dirichlet(ϵ .+
                     [count(==(i),M1) for i = 1:length(r2)] ./ n_scale),
                     r2)
        return log_Prob_m + log(Prob_r)
end
export logpdf_basechain

function logpdf_joint(r,M,l_matrix;P_r=Dirichlet(3, 1.))
        n_subject = length(M)
        log_prior = log(pdf(P_r,r))
        log_M_r = sum([count(==(i),M) for i = 1:length(r)] .* log.(r))
        log_Y_M = sum([l_matrix[i,Int64(M[i])] for i=1:n_subject])
        return log_prior + log_M_r + log_Y_M
end
export logpdf_joint

function acceptprob(r1,M1,r2,M2,l_matrix;
                            n_scale=10,ϵ=1., P_r=Dirichlet(3, 1))
        log_π1 = logpdf_joint(r1,M1,l_matrix;P_r=P_r)
        log_π2 = logpdf_joint(r2,M2,l_matrix;P_r=P_r)
        log_ϕ12 = logpdf_basechain(r1,M1,r2,M2;n_scale=n_scale,ϵ=ϵ)
        log_ϕ21 = logpdf_basechain(r2,M2,r1,M1;n_scale=n_scale,ϵ=ϵ)
        log_a12 = log_π2 + log_ϕ21 - log_π1 - log_ϕ12
        return min(1,exp(log_a12))
end
export acceptprob
