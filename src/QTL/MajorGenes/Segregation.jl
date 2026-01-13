# =============================================================================
# Segregation.jl - 主效基因检测
# =============================================================================
# 对应章节: Walsh 2nd Ed, Chapters 15-16
# =============================================================================

"""
    segregation_ratio_testing(observed::Vector{Int}, 
                             expected_ratio::Vector{Float64}) -> NamedTuple

分离比卡方检验。
"""
function segregation_ratio_testing(observed::Vector{Int},
                                   expected_ratio::Vector{Float64})
    n_total = sum(observed)
    expected = expected_ratio ./ sum(expected_ratio) .* n_total
    
    chi_sq = sum((observed .- expected).^2 ./ expected)
    df = length(observed) - 1
    pvalue = ccdf(Chisq(df), chi_sq)
    
    return (chi_square=chi_sq, df=df, pvalue=pvalue,
            observed=observed, expected=expected)
end

"""
    bartlett_test_homogeneity(variances::Vector{Float64}, 
                              sample_sizes::Vector{Int}) -> NamedTuple

Bartlett方差齐性检验。
"""
function bartlett_test_homogeneity(variances::Vector{Float64},
                                   sample_sizes::Vector{Int})
    k = length(variances)
    @assert length(sample_sizes) == k
    
    n_total = sum(sample_sizes)
    df = sample_sizes .- 1
    
    # 合并方差
    s2_pooled = sum(df .* variances) / sum(df)
    
    # Bartlett统计量
    numerator = sum(df) * log(s2_pooled) - sum(df .* log.(variances))
    C = 1 + (sum(1 ./ df) - 1/sum(df)) / (3 * (k - 1))
    
    B = numerator / C
    pvalue = ccdf(Chisq(k - 1), B)
    
    return (B=B, df=k-1, pvalue=pvalue)
end

"""
    mixture_normal_analysis(phenotypes::AbstractVector, k::Int=2) -> NamedTuple

混合正态分析检测主效基因。

返回各成分的均值、方差和混合比例。
"""
function mixture_normal_analysis(phenotypes::AbstractVector, k::Int=2;
                                 max_iter::Int=100, tol::Float64=1e-6)
    valid = filter(v -> !ismissing(v) && !(v isa Number && isnan(v)), phenotypes)
    y = Float64.(valid)
    n = length(y)
    
    # 初始化
    pi_k = fill(1.0/k, k)
    mu_k = quantile(y, range(0.2, 0.8, length=k))
    sigma2_k = fill(var(y)/k, k)
    
    # EM算法
    for _ in 1:max_iter
        # E步: 计算责任
        gamma = zeros(n, k)
        for i in 1:n
            for j in 1:k
                gamma[i, j] = pi_k[j] * pdf(Normal(mu_k[j], sqrt(sigma2_k[j])), y[i])
            end
            gamma[i, :] ./= sum(gamma[i, :])
        end
        
        # M步: 更新参数
        N_k = sum(gamma, dims=1)[:]
        
        pi_k_new = N_k ./ n
        mu_k_new = [sum(gamma[:, j] .* y) / N_k[j] for j in 1:k]
        sigma2_k_new = [sum(gamma[:, j] .* (y .- mu_k_new[j]).^2) / N_k[j] for j in 1:k]
        
        # 检查收敛
        if maximum(abs.(mu_k_new .- mu_k)) < tol
            break
        end
        
        pi_k, mu_k, sigma2_k = pi_k_new, mu_k_new, sigma2_k_new
    end
    
    return (proportions=pi_k, means=mu_k, variances=sigma2_k, k=k)
end
