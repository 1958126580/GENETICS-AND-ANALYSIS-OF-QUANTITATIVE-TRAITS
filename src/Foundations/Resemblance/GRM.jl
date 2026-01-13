# =============================================================================
# GRM.jl - 基因组关系矩阵
# =============================================================================
# 对应章节: Walsh 2nd Ed, Chapter 8
# =============================================================================

"""
    grm_vanraden_yang(Z::AbstractMatrix, p::AbstractVector; 
                      method::Symbol=:vanraden) -> RelationshipMatrix

计算基因组关系矩阵G。

# 方法
- :vanraden: G = ZZ' / 2Σp(1-p) (VanRaden 2008)
- :yang: G_ij = Σ(x_i - 2p)(x_j - 2p) / Σ2p(1-p) (Yang et al 2010)

# 参数
- Z: 中心化基因型矩阵或原始基因型(0/1/2)
- p: 等位基因频率向量
"""
function grm_vanraden_yang(Z::AbstractMatrix{<:Real}, p::AbstractVector{Float64};
                           method::Symbol=:vanraden)
    n, m = size(Z)
    @assert length(p) == m "等位基因频率数量与标记数不匹配"
    
    # 中心化
    Z_centered = Z .- 2 .* p'
    
    # 标准化因子
    scale = 2 * sum(p .* (1 .- p))
    
    if scale <= 0
        @warn "标准化因子为零或负，使用单位矩阵"
        return RelationshipMatrix(Matrix{Float64}(I, n, n), 
                                  ["ind_$i" for i in 1:n], :genomic)
    end
    
    # G矩阵
    G = (Z_centered * Z_centered') ./ scale
    
    return RelationshipMatrix(G, ["ind_$i" for i in 1:n], :genomic)
end

"""
    grm_vanraden_yang(pop::Population) -> RelationshipMatrix

从Population对象计算GRM。
"""
function grm_vanraden_yang(pop::Population)
    geno = Float64.(pop.genotypes)
    n_markers = n_markers(pop)
    
    # 估计等位基因频率
    p = [mean(geno[:, j]) / 2 for j in 1:n_markers]
    
    G = grm_vanraden_yang(geno, p)
    
    return RelationshipMatrix(G.matrix, pop.individuals, :genomic)
end

"""
    realized_relatedness(genotypes::AbstractMatrix) -> RelationshipMatrix

从基因型计算实现的亲缘关系（简化版本）。
"""
function realized_relatedness(genotypes::AbstractMatrix{<:Real})
    n, m = size(genotypes)
    
    # 标准化基因型
    geno = Float64.(genotypes)
    geno_mean = mean(geno, dims=1)
    geno_std = std(geno, dims=1)
    geno_std[geno_std .== 0] .= 1.0
    
    Z = (geno .- geno_mean) ./ geno_std
    
    # IBS相似矩阵
    G = (Z * Z') ./ m
    
    return RelationshipMatrix(G, ["ind_$i" for i in 1:n], :genomic)
end

"""
    assortative_mating_models(h2::Float64, r::Float64) -> NamedTuple

选型交配对方差的影响。

# 参数
- h2: 遗传力
- r: 配偶表型相关
"""
function assortative_mating_models(h2::Float64, r::Float64)
    # 平衡状态下加性方差的增加
    V_A_increase = (1 + r * h2) / (1 - r * h2)
    
    # 表型方差变化
    mu = r * h2  # 配偶育种值相关
    
    return (V_A_multiplier=V_A_increase, mate_bv_correlation=mu)
end

"""
    blend_grm_arm(G::AbstractMatrix, A::AbstractMatrix; 
                  alpha::Float64=0.95) -> Matrix{Float64}

混合G矩阵和A矩阵（处理奇异性）。

G_blend = α×G + (1-α)×A
"""
function blend_grm_arm(G::AbstractMatrix, A::AbstractMatrix;
                       alpha::Float64=0.95)
    @assert size(G) == size(A)
    return alpha .* G .+ (1 - alpha) .* A
end
