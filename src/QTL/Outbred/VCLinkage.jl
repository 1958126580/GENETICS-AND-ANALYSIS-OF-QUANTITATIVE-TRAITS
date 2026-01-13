# =============================================================================
# VCLinkage.jl - 远交系方差组分连锁分析
# =============================================================================
# 对应章节: Walsh 2nd Ed, Chapter 19
# =============================================================================

"""
    haseman_elston_regression(Y_diff::AbstractVector, pi_hat::AbstractVector) -> NamedTuple

Haseman-Elston回归用于检测连锁。

# 模型
(Y₁ - Y₂)² = α + β×π̂ + ε

负的β值表示连锁（共享IBD等位基因导致表型相似）。

# 参数
- Y_diff: 同胞对的表型差平方
- pi_hat: 估计的IBD共享比例
"""
function haseman_elston_regression(Y_diff::AbstractVector, pi_hat::AbstractVector)
    @assert length(Y_diff) == length(pi_hat)
    
    result = simple_regression(pi_hat, Y_diff)
    
    # 负斜率表示连锁
    linkage_evidence = result.b < 0
    
    return (intercept=result.a, slope=result.b, 
            se_slope=result.se_b, r_squared=result.r2,
            linkage_evidence=linkage_evidence)
end

"""
    outbred_variance_comp_qtl(y::AbstractVector, K::AbstractMatrix, 
                              pi_matrix::AbstractMatrix) -> VarianceComponents

方差组分连锁分析。

# 模型
y = μ + a + q + e

其中:
- a ~ N(0, Kσ²_a) (多基因效应)
- q ~ N(0, Πσ²_q) (QTL效应)

# 参数
- K: 亲缘关系矩阵
- pi_matrix: 位点特异IBD共享矩阵
"""
function outbred_variance_comp_qtl(y::AbstractVector, K::AbstractMatrix,
                                   pi_matrix::AbstractMatrix)
    n = length(y)
    X = ones(n, 1)
    
    # 同时估计多基因和QTL方差
    var_y = var(y)
    sigma2_a = var_y * 0.3
    sigma2_q = var_y * 0.2
    sigma2_e = var_y * 0.5
    
    # 简化的方差组分估计
    for _ in 1:20
        V = K * sigma2_a + pi_matrix * sigma2_q + I(n) * sigma2_e
        
        V_inv = try
            inv(V)
        catch
            inv(V + 1e-6 * I(n))
        end
        
        P = V_inv - V_inv * X * inv(X' * V_inv * X) * X' * V_inv
        Py = P * y
        
        # 更新
        sigma2_a = max(1e-6, (Py' * K * Py) / tr(P * K))
        sigma2_q = max(1e-6, (Py' * pi_matrix * Py) / tr(P * pi_matrix))
        sigma2_e = max(1e-6, (Py' * Py) / tr(P))
    end
    
    h2 = sigma2_a / (sigma2_a + sigma2_q + sigma2_e)
    
    return VarianceComponents(sigma2_a + sigma2_q, sigma2_e; method=:VC_Linkage)
end
