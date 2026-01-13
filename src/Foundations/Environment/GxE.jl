# =============================================================================
# GxE.jl - 基因型×环境互作
# =============================================================================
# 对应章节: Walsh 2nd Ed, Chapter 6
# =============================================================================

"""
    reaction_norm_slope(phenotypes::Vector{Float64}, environments::Vector{Float64}) -> NamedTuple

计算反应规范斜率（环境敏感性）。

# 中文说明
反应规范描述基因型在不同环境中的表型表达。
斜率反映基因型对环境变化的敏感性。
"""
function reaction_norm_slope(phenotypes::Vector{Float64}, environments::Vector{Float64})
    result = simple_regression(environments, phenotypes)
    
    return (slope=result.b, intercept=result.a, r_squared=result.r2,
            sensitivity=abs(result.b))
end

"""
    environmental_variance_partition(within_env::Vector{Float64}, 
                                    between_env::Vector{Float64}) -> NamedTuple

分解环境方差为组内和组间成分。

# 参数
- within_env: 同一环境内的方差估计
- between_env: 不同环境间的方差估计

# 中文说明
V_E = V_Ew + V_Eb (特殊环境方差 + 一般环境方差)
"""
function environmental_variance_partition(within_env::Vector{Float64},
                                          between_env::Vector{Float64})
    V_Ew = mean(within_env)   # 特殊环境方差
    V_Eb = var(between_env)    # 一般环境方差
    V_E = V_Ew + V_Eb
    
    return (V_E=V_E, V_Ew=V_Ew, V_Eb=V_Eb, 
            prop_within=V_Ew/V_E, prop_between=V_Eb/V_E)
end

"""
    genotype_environment_covariance(G::AbstractVector, E::AbstractVector) -> Float64

计算遗传与环境的协方差。

# 中文说明
当基因型与环境非随机关联时（如生态位选择），Cov(G,E)≠0。
"""
function genotype_environment_covariance(G::AbstractVector, E::AbstractVector)
    @assert length(G) == length(E)
    return cov(Float64.(G), Float64.(E))
end

"""
    gxe_anova(phenotypes::Matrix{Float64}, genotypes::Vector{Int}, 
              environments::Vector{Int}) -> NamedTuple

GxE互作方差分析。

phenotypes: 表型矩阵 (基因型 × 环境)
"""
function gxe_anova(phenotypes::Matrix{Float64})
    n_g, n_e = size(phenotypes)
    grand_mean = mean(phenotypes)
    
    # 基因型均值和环境均值
    g_means = [mean(phenotypes[i, :]) for i in 1:n_g]
    e_means = [mean(phenotypes[:, j]) for j in 1:n_e]
    
    # 平方和
    SS_G = n_e * sum((g_means .- grand_mean).^2)
    SS_E = n_g * sum((e_means .- grand_mean).^2)
    
    SS_total = sum((phenotypes .- grand_mean).^2)
    
    # GxE交互作用
    SS_GxE = 0.0
    for i in 1:n_g, j in 1:n_e
        expected = g_means[i] + e_means[j] - grand_mean
        SS_GxE += (phenotypes[i,j] - expected)^2
    end
    
    # 自由度
    df_G = n_g - 1
    df_E = n_e - 1
    df_GxE = df_G * df_E
    
    # 均方
    MS_G = SS_G / df_G
    MS_E = SS_E / df_E
    MS_GxE = SS_GxE / df_GxE
    
    return (SS_G=SS_G, SS_E=SS_E, SS_GxE=SS_GxE,
            MS_G=MS_G, MS_E=MS_E, MS_GxE=MS_GxE,
            df_G=df_G, df_E=df_E, df_GxE=df_GxE)
end
