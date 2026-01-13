# =============================================================================
# TwoLocus.jl - 双位点上位性分析
# =============================================================================
# 对应章节: Walsh 2nd Ed, Chapter 5
# =============================================================================

"""
    two_locus_epistasis(p1::Float64, p2::Float64, effects::Dict{Symbol,Float64}) -> NamedTuple

双位点上位性模型分析。

# effects字典键
- :a1, :a2: 加性效应
- :d1, :d2: 显性效应  
- :aa, :ad, :da, :dd: 上位性效应

# 中文说明
Cockerham正交分解双位点遗传方差。
"""
function two_locus_epistasis(p1::Float64, p2::Float64, effects::Dict{Symbol,Float64})
    q1, q2 = 1-p1, 1-p2
    
    a1 = get(effects, :a1, 0.0)
    a2 = get(effects, :a2, 0.0)
    d1 = get(effects, :d1, 0.0)
    d2 = get(effects, :d2, 0.0)
    i_aa = get(effects, :aa, 0.0)
    i_ad = get(effects, :ad, 0.0)
    i_da = get(effects, :da, 0.0)
    i_dd = get(effects, :dd, 0.0)
    
    # 各方差成分
    V_A1 = 2*p1*q1 * (a1 + d1*(q1-p1))^2
    V_A2 = 2*p2*q2 * (a2 + d2*(q2-p2))^2
    V_D1 = (2*p1*q1*d1)^2
    V_D2 = (2*p2*q2*d2)^2
    
    V_AA = (2*p1*q1)*(2*p2*q2) * i_aa^2
    V_AD = (2*p1*q1)*(2*p2*q2*d2)^2 * i_ad^2 / (p2*q2)
    V_DA = (2*p1*q1*d1)^2*(2*p2*q2) * i_da^2 / (p1*q1)
    V_DD = (2*p1*q1*d1)^2 * (2*p2*q2*d2)^2 * i_dd^2 / (4*p1*q1*p2*q2)
    
    V_A = V_A1 + V_A2
    V_D = V_D1 + V_D2
    V_I = V_AA + V_AD + V_DA + V_DD
    V_G = V_A + V_D + V_I
    
    return (V_A=V_A, V_D=V_D, V_I=V_I, V_G=V_G,
            V_A1=V_A1, V_A2=V_A2, V_D1=V_D1, V_D2=V_D2,
            V_AA=V_AA, V_AD=V_AD, V_DA=V_DA, V_DD=V_DD)
end

"""
    cockerham_orthogonal_partition(genotype_means::Matrix{Float64}, 
                                   p1::Float64, p2::Float64) -> NamedTuple

Cockerham正交方差分解。

genotype_means: 3×3矩阵，[AA,Aa,aa] × [BB,Bb,bb]
"""
function cockerham_orthogonal_partition(genotype_means::Matrix{Float64},
                                        p1::Float64, p2::Float64)
    @assert size(genotype_means) == (3, 3)
    
    q1, q2 = 1-p1, 1-p2
    
    # 加权均值
    w1 = [p1^2, 2*p1*q1, q1^2]
    w2 = [p2^2, 2*p2*q2, q2^2]
    
    mu = sum(w1[i] * w2[j] * genotype_means[i,j] for i in 1:3, j in 1:3)
    
    # 边际基因型均值
    G1 = [sum(w2[j] * genotype_means[i,j] for j in 1:3) for i in 1:3]
    G2 = [sum(w1[i] * genotype_means[i,j] for i in 1:3) for j in 1:3]
    
    # 使用边际计算加性和显性效应
    a1 = (G1[1] - G1[3]) / 2
    d1 = G1[2] - (G1[1] + G1[3]) / 2
    a2 = (G2[1] - G2[3]) / 2
    d2 = G2[2] - (G2[1] + G2[3]) / 2
    
    return (mu=mu, a1=a1, a2=a2, d1=d1, d2=d2,
            marginal_1=G1, marginal_2=G2)
end
