# =============================================================================
# LD.jl - 连锁不平衡分析
# =============================================================================
# 对应章节: Walsh 2nd Ed, Chapters 5-6
# =============================================================================

"""
    linkage_disequilibrium_metrics(P_AB::Float64, p_A::Float64, p_B::Float64) -> NamedTuple

计算LD指标: D, D', r²

# 参数
- P_AB: AB单倍型频率
- p_A, p_B: 边际等位基因频率

# 返回
- D: 配子不平衡系数
- D_prime: 标准化D (Lewontin)
- r_squared: 相关系数平方
"""
function linkage_disequilibrium_metrics(P_AB::Float64, p_A::Float64, p_B::Float64)
    # D = P_AB - p_A × p_B
    D = P_AB - p_A * p_B
    
    # D的边界
    p_a = 1 - p_A
    p_b = 1 - p_B
    
    if D >= 0
        D_max = min(p_A * p_b, p_a * p_B)
    else
        D_max = max(-p_A * p_B, -p_a * p_b)
    end
    
    # D' (Lewontin's standardization)
    D_prime = abs(D_max) > 0 ? D / D_max : 0.0
    
    # r² (correlation squared)
    denom = p_A * p_a * p_B * p_b
    r_squared = denom > 0 ? D^2 / denom : 0.0
    
    return (D=D, D_prime=D_prime, r_squared=r_squared, D_max=D_max)
end

"""
    ld_decay_model(D0::Float64, c::Float64, t::Int) -> Float64

LD衰减模型: D_t = D_0 × (1-c)^t

# 参数
- D0: 初始D值
- c: 重组率
- t: 代数
"""
ld_decay_model(D0::Float64, c::Float64, t::Int) = D0 * (1 - c)^t

"""
    ld_from_genotypes(geno1::AbstractVector, geno2::AbstractVector) -> NamedTuple

从基因型数据估计LD。
"""
function ld_from_genotypes(geno1::AbstractVector{<:Real}, geno2::AbstractVector{<:Real})
    @assert length(geno1) == length(geno2)
    
    # 过滤缺失
    valid_idx = [i for i in 1:length(geno1) 
                 if !isnan(geno1[i]) && !isnan(geno2[i]) && 
                    geno1[i] >= 0 && geno2[i] >= 0]
    
    n = length(valid_idx)
    n < 10 && return (D=NaN, r_squared=NaN, n=n)
    
    g1 = Float64[geno1[i] for i in valid_idx]
    g2 = Float64[geno2[i] for i in valid_idx]
    
    # r² from correlation
    r = cor(g1, g2)
    r_squared = r^2
    
    # 估计等位基因频率
    p1 = mean(g1) / 2
    p2 = mean(g2) / 2
    
    # 估计D
    D = r * sqrt(p1 * (1-p1) * p2 * (1-p2))
    
    return (D=D, r_squared=r_squared, r=r, p1=p1, p2=p2, n=n)
end

"""
    ld_matrix(genotypes::AbstractMatrix; min_maf::Float64=0.01) -> Matrix{Float64}

计算标记间的LD矩阵 (r²)。
"""
function ld_matrix(genotypes::AbstractMatrix{<:Real}; min_maf::Float64=0.01)
    n_ind, n_markers = size(genotypes)
    ld_mat = zeros(n_markers, n_markers)
    
    for i in 1:n_markers
        ld_mat[i, i] = 1.0
        for j in (i+1):n_markers
            result = ld_from_genotypes(genotypes[:, i], genotypes[:, j])
            ld_mat[i, j] = result.r_squared
            ld_mat[j, i] = result.r_squared
        end
    end
    
    return ld_mat
end
