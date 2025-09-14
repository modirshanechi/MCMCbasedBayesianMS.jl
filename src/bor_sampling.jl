# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
function sample_bor(BMS::BMSObj; borobj::BMSBorSampNaive=BMS.bms_borsamp)
    l_matrix = BMS.l_matrix
    n_model = BMS.n_model
    n_subject = BMS.n_subject

    n_samples = borobj.n_samples

    logP0 = logp_datagivenr(l_matrix, ones(n_model) ./ n_model)

    R1_matrix = zeros(n_samples, n_model);
    π1_matrix = zeros(n_samples);
    for i_sample = 1:n_samples
        r1,M1 = sample_prior(BMS.α;n_model=n_model,n_subject=n_subject)
        R1_matrix[i_sample,:] = r1
        π1_matrix[i_sample] = logp_datagivenr(l_matrix,r1)
    end
    π1_max = findmax(π1_matrix)[1]
    logP1 = π1_max + log(sum(exp.(π1_matrix .- π1_max))) - log(n_samples)
    
    BOR = 1 / (1 + exp(logP1 - logP0))

    return (BOR = BOR, logP0=logP0, logP1=logP1, 
            R1_matrix=R1_matrix, π1_matrix=π1_matrix)
end
export sample_posterior

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
function logp_datagivenr(l_matrix,r)
        n_subject = size(l_matrix)[1]
        l_max_vec = [findmax(l_matrix[n,:])[1] for n = 1:n_subject]
        l_xr_vec = [(l_max_vec[n] +
                    log(sum(r .* exp.(l_matrix[n,:] .- l_max_vec[n]))))
                        for n = 1:n_subject]
        return sum(l_xr_vec)
end
export logp_datagivenr