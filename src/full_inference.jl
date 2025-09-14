function BMSStats(posterior_results,bor_results)
    BOR = bor_results.BOR
    r_samples = posterior_results.r_matrix
    M_samples = posterior_results.M_matrix
    
    # evaluating exp_r, d_exp_r, and xp
    r_samples_concat = zeros(   size(r_samples)[1]*size(r_samples)[3],
                                size(r_samples)[2]);
    for i = 1:size(r_samples)[2]
            r_samples_concat[:,i] = r_samples[:,i,:][:]
    end
    exp_r   = mean(r_samples_concat, dims = 1)[:]
    d_exp_r = std(r_samples_concat, dims = 1)[:]

    best_model_samples = zeros(size(r_samples_concat));
    for i = 1:size(r_samples_concat)[1]
            best_model_samples[i,findmax(r_samples_concat[i,:])[2]] = 1
    end
    xp = mean(best_model_samples, dims = 1)[:]

    # evaluating pxp
    pxp = xp .* (1 - BOR) .+ ones(length(xp)) ./ length(xp) .* BOR

    # evaluating exp_M
    M_samples_all = zeros(  size(M_samples)[1]*size(M_samples)[3],
                            size(M_samples)[2]);
    for i = 1:size(M_samples)[2]
            M_samples_all[:,i] = M_samples[:,i,:][:]
    end
    exp_M = zeros(size(M_samples_all)[2],size(r_samples)[2])
    for i = 1:size(r_samples)[2]
            exp_M[:,i] = mean(M_samples_all .== i,dims=1)[:]
    end
    return BMSStats(exp_M, exp_r, d_exp_r, xp, BOR, pxp)
end

function doBMSinference(BMS::BMSObj)
    posterior_results = sample_posterior(BMS)
    bor_results = sample_bor(BMS)
    stats = BMSStats(posterior_results,bor_results)
    BMSResults(BMS,stats,posterior_results,bor_results)
end
export doBMSinference