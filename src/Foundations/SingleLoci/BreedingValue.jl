# =============================================================================
# BreedingValue.jl - 育种值与遗传方差分解
# =============================================================================
# 对应章节: Walsh 2nd Ed, Chapter 4 (Fisher分解)
# =============================================================================

"""
    fisher_decomposition(p::Float64, q::Float64, a::Float64, d::Float64) -> NamedTuple

Fisher加性-显性分解。

# 基因型值
- G_AA = +a
- G_Aa = d  
- G_aa = -a

# 返回
- alpha: 平均效应 α = a + d(q-p)
- V_A: 加性方差 = 2pqα²
- V_D: 显性方差 = (2pqd)²

# 中文说明
将遗传方差分解为加性和显性成分。
"""
function fisher_decomposition(p::Float64, q::Float64, a::Float64, d::Float64)
    @assert 0 <= p <= 1 && isapprox(p + q, 1.0)
    
    # 平均效应(替代效应)
    alpha = a + d * (q - p)
    
    # 方差分解
    V_A = 2 * p * q * alpha^2
    V_D = (2 * p * q * d)^2
    V_G = V_A + V_D
    
    return (alpha=alpha, V_A=V_A, V_D=V_D, V_G=V_G, a=a, d=d, p=p, q=q)
end

"""
    average_excesses(p::Float64, q::Float64, G_AA::Float64, G_Aa::Float64, G_aa::Float64) -> NamedTuple

计算等位基因平均超值。

α*_A = (pG_AA + qG_Aa) - μ
α*_a = (pG_Aa + qG_aa) - μ

# 中文说明
平均超值是携带某等位基因的个体与群体均值的偏差。
"""
function average_excesses(p::Float64, q::Float64, G_AA::Float64, G_Aa::Float64, G_aa::Float64)
    # 群体均值
    mu = p^2 * G_AA + 2*p*q * G_Aa + q^2 * G_aa
    
    # 平均超值
    alpha_star_A = (p * G_AA + q * G_Aa) - mu
    alpha_star_a = (p * G_Aa + q * G_aa) - mu
    
    return (alpha_star_A=alpha_star_A, alpha_star_a=alpha_star_a, mu=mu)
end

"""
    average_effects_alpha(p::Float64, q::Float64, a::Float64, d::Float64) -> Float64

计算平均效应(基因替代效应)。

α = a + d(q - p)

# 中文说明
用a替换A等位基因时，预期表型变化量。
"""
average_effects_alpha(p::Float64, q::Float64, a::Float64, d::Float64) = a + d * (q - p)

"""
    breeding_value_partition(p::Float64, q::Float64, alpha::Float64) -> NamedTuple

计算各基因型的育种值。

A_AA = 2qα, A_Aa = (q-p)α, A_aa = -2pα

# 中文说明
育种值是个体作为亲本的遗传价值，等于其等位基因平均效应之和。
"""
function breeding_value_partition(p::Float64, q::Float64, alpha::Float64)
    A_AA = 2 * q * alpha
    A_Aa = (q - p) * alpha
    A_aa = -2 * p * alpha
    
    return (A_AA=A_AA, A_Aa=A_Aa, A_aa=A_aa)
end

"""
    dominance_deviation_calc(p::Float64, q::Float64, d::Float64) -> NamedTuple

计算各基因型的显性偏差。

D_AA = -2q²d, D_Aa = 2pqd, D_aa = -2p²d

# 中文说明
显性偏差是基因型值与育种值的差异。
"""
function dominance_deviation_calc(p::Float64, q::Float64, d::Float64)
    D_AA = -2 * q^2 * d
    D_Aa = 2 * p * q * d
    D_aa = -2 * p^2 * d
    
    return (D_AA=D_AA, D_Aa=D_Aa, D_aa=D_aa)
end

"""
    variance_partition_single_locus(p::Float64, q::Float64, a::Float64, d::Float64) -> BreedingValueResult

单位点遗传方差的完整分解。

V_A = 2pqα², V_D = (2pqd)²
"""
function variance_partition_single_locus(p::Float64, q::Float64, a::Float64, d::Float64)
    decomp = fisher_decomposition(p, q, a, d)
    
    return BreedingValueResult(
        p, q, a, d, decomp.alpha,
        decomp.V_A, decomp.V_D, decomp.V_G
    )
end

"""
    dominance_ratio(a::Float64, d::Float64) -> Float64

显性度 = d/|a|

# 解释
- 0: 无显性(加性)
- 1: 完全显性
- >1: 超显性
"""
dominance_ratio(a::Float64, d::Float64) = abs(a) ≈ 0 ? Inf : d / abs(a)
