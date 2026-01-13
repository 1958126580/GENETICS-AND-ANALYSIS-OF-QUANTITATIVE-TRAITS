# =============================================================================
# Transformations.jl - 尺度与阈值模型
# =============================================================================
# 对应章节: Walsh 2nd Ed, Chapter 14
# =============================================================================

"""
    mather_jinks_scaling_test(P1::Float64, P2::Float64, F1::Float64,
                              F2::Float64, BC1::Float64, BC2::Float64) -> NamedTuple

Mather-Jinks A, B, C尺度检验。

# 检验统计量
- A = 2BC1 - P1 - F1 (应接近0)
- B = 2BC2 - P2 - F1 (应接近0)  
- C = 4F2 - 2F1 - P1 - P2 (应接近0)

# 中文说明
用于检验是否需要尺度变换或是否存在上位性。
"""
function mather_jinks_scaling_test(P1::Float64, P2::Float64, F1::Float64,
                                   F2::Float64, BC1::Float64, BC2::Float64)
    A = 2*BC1 - P1 - F1
    B = 2*BC2 - P2 - F1
    C = 4*F2 - 2*F1 - P1 - P2
    
    return (A=A, B=B, C=C,
            A_significant=abs(A) > 0.1 * abs(P1 - P2),
            B_significant=abs(B) > 0.1 * abs(P1 - P2),
            C_significant=abs(C) > 0.1 * abs(P1 - P2))
end

"""
    threshold_model_liability(proportions::Vector{Float64}) -> NamedTuple

阈值模型责任估计。

# 模型
假设存在连续的责任(liability)分布，超过阈值T则表现疾病。
P(affected) = 1 - Φ(T)

# 参数
- proportions: 各类别的受影响比例
"""
function threshold_model_liability(proportions::Vector{Float64})
    thresholds = [quantile(Normal(), 1 - p) for p in proportions]
    
    return (thresholds=thresholds, proportions=proportions)
end

"""
    polychoric_correlation(contingency::Matrix{Int}) -> Float64

估计潜在连续变量间的多分相关。

# 中文说明
对于分类数据，假设潜在连续责任变量，估计其相关。
"""
function polychoric_correlation(contingency::Matrix{Int})
    # 简化实现：使用四分相关近似
    n = sum(contingency)
    
    if size(contingency) == (2, 2)
        a, b, c, d = contingency[1,1], contingency[1,2], contingency[2,1], contingency[2,2]
        # Tetrachoric correlation approximation (Digby)
        ad = a * d
        bc = b * c
        
        if ad == bc
            return 0.0
        end
        
        r = (sqrt(ad) - sqrt(bc)) / (sqrt(ad) + sqrt(bc))
        return r * π / 2  # 近似转换
    else
        # 更复杂的多分情况，返回简单Pearson近似
        r, c = size(contingency)
        x = repeat(1:r, 1, c)
        y = repeat((1:c)', r, 1)
        
        weights = contingency ./ n
        mean_x = sum(x .* weights)
        mean_y = sum(y .* weights)
        
        cov_xy = sum((x .- mean_x) .* (y .- mean_y) .* weights)
        var_x = sum((x .- mean_x).^2 .* weights)
        var_y = sum((y .- mean_y).^2 .* weights)
        
        return cov_xy / sqrt(var_x * var_y)
    end
end

"""
    liability_heritability(prevalence::Float64, relative_risk_siblings::Float64) -> Float64

从患病率和同胞相对风险估计责任遗传力。

# 公式 (Falconer)
h² = 2r / i² 

其中 r 是同胞间责任相关，i 是选择强度。
"""
function liability_heritability(prevalence::Float64, relative_risk_siblings::Float64)
    T = quantile(Normal(), 1 - prevalence)
    i = pdf(Normal(), T) / prevalence  # 选择强度
    
    # 同胞患病率
    K_s = prevalence * relative_risk_siblings
    T_s = quantile(Normal(), 1 - K_s)
    
    # 责任相关近似
    r = (T - T_s) / i
    
    # 遗传力（同胞相关×2因为共享1/2基因）
    h2 = min(1.0, 2 * r)
    
    return h2
end
