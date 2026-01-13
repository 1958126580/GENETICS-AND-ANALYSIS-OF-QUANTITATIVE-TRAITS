# =============================================================================
# IBD.jl - 同源同源分析
# =============================================================================
# 对应章节: Walsh 2nd Ed, Chapter 7
# =============================================================================

"""
    identity_by_descent_vs_state(g1::Genotype, g2::Genotype) -> NamedTuple

区分IBD和IBS（同源同状态）。

# 中文说明
- IBD: 等位基因来自共同祖先
- IBS: 等位基因相同但不一定IBD
"""
function identity_by_descent_vs_state(g1::Genotype, g2::Genotype)
    # IBS计数 (等位基因匹配数)
    ibs = 0
    if g1.allele1 == g2.allele1 || g1.allele1 == g2.allele2
        ibs += 1
    end
    if g1.allele2 == g2.allele1 || g1.allele2 == g2.allele2
        ibs += 1
    end
    
    return (ibs_count=ibs, ibs_proportion=ibs/2)
end

"""
    recursive_kinship(ped::Pedigree) -> RelationshipMatrix

Henderson递归算法计算亲缘系数矩阵。

# 公式
Φ_ij = 0.5 × (Φ_i,sire_j + Φ_i,dam_j)
Φ_ii = 0.5 × (1 + Φ_sire,dam) = 0.5 × (1 + F_i的亲本亲缘)

# 中文说明
从谱系递归计算Malécot亲缘系数矩阵。
A矩阵 = 2 × Φ矩阵
"""
function recursive_kinship(ped::Pedigree)
    n = length(ped.ids)
    id_to_idx = Dict(id => i for (i, id) in enumerate(ped.ids))
    
    # 亲缘系数矩阵 (Φ)
    Phi = zeros(n, n)
    
    for i in 1:n
        # 自身亲缘 (对角线)
        sire_i = ped.sires[i]
        dam_i = ped.dams[i]
        
        if sire_i === nothing || dam_i === nothing ||
           !haskey(id_to_idx, sire_i) || !haskey(id_to_idx, dam_i)
            Phi[i, i] = 0.5  # 奠基者假设无近交
        else
            s_idx = id_to_idx[sire_i]
            d_idx = id_to_idx[dam_i]
            Phi[i, i] = 0.5 * (1 + Phi[s_idx, d_idx])
        end
        
        # 与其他个体的亲缘
        for j in (i+1):n
            sire_j = ped.sires[j]
            dam_j = ped.dams[j]
            
            phi_ij = 0.0
            count = 0
            
            if sire_j !== nothing && haskey(id_to_idx, sire_j)
                phi_ij += Phi[i, id_to_idx[sire_j]]
                count += 1
            end
            if dam_j !== nothing && haskey(id_to_idx, dam_j)
                phi_ij += Phi[i, id_to_idx[dam_j]]
                count += 1
            end
            
            Phi[i, j] = count > 0 ? 0.5 * phi_ij : 0.0
            Phi[j, i] = Phi[i, j]
        end
    end
    
    return RelationshipMatrix(Phi, ped.ids, :kinship)
end

"""
    coefficient_of_fraternity(ped::Pedigree, id1::String, id2::String) -> Float64

计算同胞显性系数Δ₇（两个体共享两个IBD等位基因的概率）。

# 中文说明
Δ系数用于估计显性遗传方差的贡献。
"""
function coefficient_of_fraternity(ped::Pedigree, id1::String, id2::String)
    sire1, dam1 = get_parents(ped, id1)
    sire2, dam2 = get_parents(ped, id2)
    
    # 全同胞: Δ₇ = 0.25
    if sire1 !== nothing && dam1 !== nothing &&
       sire1 == sire2 && dam1 == dam2
        return 0.25
    end
    
    # 半同胞: Δ₇ = 0
    if (sire1 == sire2 && sire1 !== nothing) ||
       (dam1 == dam2 && dam1 !== nothing)
        return 0.0
    end
    
    return 0.0
end

"""
    inbreeding_coef_f(ped::Pedigree, id::String) -> Float64

从谱系计算个体近交系数F。

F = Φ_sire,dam (父母的亲缘系数)
"""
function inbreeding_coef_f(ped::Pedigree, id::String)
    K = recursive_kinship(ped)
    idx = findfirst(==(id), ped.ids)
    idx === nothing && error("个体不在谱系中")
    
    # F = 2Φ_ii - 1 = A_ii - 1
    return 2 * K.matrix[idx, idx] - 1.0
end

"""
    additive_relationship_matrix(ped::Pedigree) -> RelationshipMatrix

计算加性遗传关系矩阵A = 2Φ。
"""
function additive_relationship_matrix(ped::Pedigree)
    K = recursive_kinship(ped)
    A = 2.0 .* K.matrix
    return RelationshipMatrix(A, K.ids, :additive)
end
