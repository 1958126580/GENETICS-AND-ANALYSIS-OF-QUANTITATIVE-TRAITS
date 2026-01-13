# =============================================================================
# Visualization.jl - GWAS可视化
# =============================================================================

"""
    manhattan_plot(result::GWASResult) -> NamedTuple

生成Manhattan图数据。
"""
function manhattan_plot(result::GWASResult)
    neg_log_p = -log10.(max.(result.pvalues, 1e-300))
    
    return (positions=result.positions, 
            chromosomes=result.chromosomes,
            neg_log_p=neg_log_p,
            threshold_genome=(-log10(GWAS_GENOME_WIDE_ALPHA)),
            threshold_suggestive=(-log10(GWAS_SUGGESTIVE_ALPHA)))
end

"""
    qq_plot(pvalues::Vector{Float64}) -> NamedTuple

生成QQ图数据。
"""
function qq_plot(pvalues::Vector{Float64})
    valid_p = sort(filter(p -> 0 < p <= 1, pvalues))
    n = length(valid_p)
    
    expected = [(i - 0.5) / n for i in 1:n]
    observed = valid_p
    
    return (expected=-log10.(expected), 
            observed=-log10.(observed),
            identity_line=[0.0, maximum(-log10.(expected))])
end

manhattan_qq_viz(stats) = (manhattan=manhattan_plot(stats), qq=qq_plot(stats.pvalues))
