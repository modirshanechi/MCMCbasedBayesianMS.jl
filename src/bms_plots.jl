function plotstats(BMS::BMSResults; figsize=(1200,400))
    # Create a 1x2 layout
    plt = plot(layout = (1,3), size=figsize)
    
    # Left subplot: estimates with error bars
    bar!(
        plt[1],
        1:BMS.object.n_model,
        BMS.stats.exp_r,
        yerr = BMS.stats.d_exp_r,
        label = "exp_r",
        ylim = (0,1),
        legend = :topleft,
    )

    # Right subplot: ground truth
    bar!(
        plt[2],
        1:BMS.object.n_model,
        BMS.stats.pxp,
        label = "pxp",
        ylim = (0,1),
        legend = :topleft,
    )

    # Right subplot: ground truth
    heatmap!(
        plt[3],
        BMS.stats.exp_M,
        label = "exp_M"
    )

    display(plt) 
    return plt
end
export plotstats

function plotrsamples(BMS::BMSResults; figsize=(300,400))
    R = BMS.posterior_results.r_matrix
    N = size(R)[2]; M = size(R)[3]

    plt = plot(layout=(1, N), size=(figsize[1]*N, figsize[2]))

    for i in 1:N
        # pick subplot i
        for j in 1:M
            density!(
                plt[i],
                R[:, i, j],
                label = "chain $j",
                lw = 2,
                alpha = 0.6
            )
        end
        xlabel!(plt[i], "x")
        ylabel!(plt[i], "pdf")
        title!(plt[i], "r[$i]")
    end

    display(plt)
    return plt
end
export plotrsamples