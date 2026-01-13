# =============================================================================
# GenerationMeans.jl - 世代均值分析与品系杂交
# =============================================================================
# 对应章节: Walsh 2nd Ed, Chapter 11
# =============================================================================

"""
    generation_mean_analysis(P1::Float64, P2::Float64, F1::Float64, 
                            F2::Float64, BC1::Float64, BC2::Float64) -> LineCrossResult

Mather-Jinks世代均值分析。

# 模型
μ = m + [a]x_a + [d]x_d

各世代期望:
- P1: m + [a]
- P2: m - [a]
- F1: m + [d]
- F2: m + [d]/2
- BC1: m + [a]/2 + [d]/2
- BC2: m - [a]/2 + [d]/2

# 中文说明
用于估计品系杂交的加性[a]和显性[d]效应。
"""
function generation_mean_analysis(P1::Float64, P2::Float64, F1::Float64,
                                  F2::Float64, BC1::Float64, BC2::Float64)
    # 设计矩阵 (6个观测，3个参数)
    X = [1.0  1.0  0.0;   # P1
         1.0 -1.0  0.0;   # P2
         1.0  0.0  1.0;   # F1
         1.0  0.0  0.5;   # F2
         1.0  0.5  0.5;   # BC1
         1.0 -0.5  0.5]   # BC2
    
    y = [P1, P2, F1, F2, BC1, BC2]
    
    # 最小二乘估计
    params = X \ y
    m, a, d = params
    
    # 拟合值和残差
    y_hat = X * params
    resid = y - y_hat
    chi_sq = sum(resid.^2 ./ abs.(y_hat))  # 简化的拟合检验
    
    df = 3  # 6 - 3 = 3 剩余自由度
    pvalue = ccdf(Chisq(df), chi_sq)
    
    # Castle-Wright基因数估计
    n_genes = castle_wright_n_factors(P1, P2, var(y), var(resid))
    
    return LineCrossResult(m, a, d, chi_sq, pvalue, n_genes, df)
end

"""
    joint_scaling_test_m_a_d(means::Vector{Float64}, 
                            design::Matrix{Float64}) -> NamedTuple

Mather-Jinks联合尺度检验。

# 参数
- means: 各世代观测均值
- design: 设计矩阵 [1, a系数, d系数]
"""
function joint_scaling_test_m_a_d(means::Vector{Float64}, design::Matrix{Float64})
    n = length(means)
    @assert size(design, 1) == n
    
    params = design \ means
    y_hat = design * params
    resid = means - y_hat
    
    SS_resid = sum(resid.^2)
    df = n - 3
    chi_sq = SS_resid / var(means)
    pvalue = df > 0 ? ccdf(Chisq(df), chi_sq) : NaN
    
    return (m=params[1], a=params[2], d=params[3],
            chi_square=chi_sq, pvalue=pvalue, fitted=y_hat)
end

"""
    castle_wright_n_factors(P1::Float64, P2::Float64, 
                           V_seg::Float64, V_env::Float64) -> Float64

Castle-Wright有效基因数估计。

# 公式
n_e = (P1 - P2)² / (8 × V_seg)

其中 V_seg 是分离群体的遗传方差。

# 中文说明
估计控制两品系差异的最小基因数。
假设所有基因效应相等且加性。
"""
function castle_wright_n_factors(P1::Float64, P2::Float64,
                                 V_seg::Float64, V_env::Float64)
    V_G = V_seg - V_env
    V_G = max(V_G, 0.0)
    
    if V_G ≈ 0
        return Inf
    end
    
    n_e = (P1 - P2)^2 / (8 * V_G)
    return n_e
end
