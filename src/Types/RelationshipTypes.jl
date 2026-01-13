# =============================================================================
# RelationshipTypes.jl - 谱系与亲缘关系类型定义
# =============================================================================
#
# 定义谱系结构和亲缘关系矩阵
#
# 对应章节: Walsh 2nd Ed, Chapters 7-8, 10
# =============================================================================

"""
    Pedigree

谱系结构，记录个体间的亲子关系。

# 字段
- `ids::Vector{String}`: 个体ID列表
- `sires::Vector{Union{String, Nothing}}`: 父本ID列表
- `dams::Vector{Union{String, Nothing}}`: 母本ID列表
- `founder_ids::Set{String}`: 奠基者（无已知亲本）ID集合
- `n_generations::Int`: 估计的世代数

# 中文说明
谱系（Pedigree）记录个体间的血缘关系，是计算亲缘系数矩阵的基础。

奠基者（Founder）是指父母均未知的个体，通常假定：
- 奠基者间不相关
- 奠基者不近交（F = 0）

谱系必须按时间顺序排列：后代总是出现在亲本之后。

# 示例
```julia
ped = Pedigree(
    ["A", "B", "C", "D", "E"],           # 个体ID
    [nothing, nothing, "A", "A", "C"],    # 父本
    [nothing, nothing, "B", "B", "D"]     # 母本
)

# A和B是奠基者
# C和D是A×B的后代
# E是C×D的后代（近交个体）
```

参见: [`recursive_kinship`](@ref), [`inbreeding_coef_f`](@ref)
"""
struct Pedigree
    ids::Vector{String}
    sires::Vector{Union{String, Nothing}}
    dams::Vector{Union{String, Nothing}}
    founder_ids::Set{String}
    n_generations::Int
    
    function Pedigree(ids::Vector{String}, 
                     sires::Vector{<:Union{String, Nothing}},
                     dams::Vector{<:Union{String, Nothing}})
        n = length(ids)
        if length(sires) != n || length(dams) != n
            throw(ArgumentError("ids、sires、dams长度必须相同"))
        end
        
        # 转换类型
        sires_conv = Vector{Union{String, Nothing}}(sires)
        dams_conv = Vector{Union{String, Nothing}}(dams)
        
        # 识别奠基者
        id_set = Set(ids)
        founders = Set{String}()
        
        for i in 1:n
            is_founder = true
            
            # 如果父本已知且在谱系中，则不是奠基者
            if sires_conv[i] !== nothing && sires_conv[i] in id_set
                is_founder = false
            end
            # 如果母本已知且在谱系中，则不是奠基者
            if dams_conv[i] !== nothing && dams_conv[i] in id_set
                is_founder = false
            end
            
            if is_founder
                push!(founders, ids[i])
            end
        end
        
        # 估计世代数（简化方法）
        n_gen = _estimate_generations(ids, sires_conv, dams_conv, founders)
        
        new(ids, sires_conv, dams_conv, founders, n_gen)
    end
end

function _estimate_generations(ids, sires, dams, founders)
    id_to_idx = Dict(id => i for (i, id) in enumerate(ids))
    generation = Dict{String, Int}()
    
    # 奠基者为第0代
    for f in founders
        if f in keys(id_to_idx)
            generation[f] = 0
        end
    end
    
    # 迭代计算后代的世代
    changed = true
    max_iter = length(ids)
    iter = 0
    
    while changed && iter < max_iter
        changed = false
        iter += 1
        
        for (i, id) in enumerate(ids)
            if haskey(generation, id)
                continue
            end
            
            sire_gen = nothing
            dam_gen = nothing
            
            if sires[i] !== nothing && haskey(generation, sires[i])
                sire_gen = generation[sires[i]]
            end
            if dams[i] !== nothing && haskey(generation, dams[i])
                dam_gen = generation[dams[i]]
            end
            
            if sire_gen !== nothing || dam_gen !== nothing
                parent_gen = max(
                    sire_gen === nothing ? 0 : sire_gen,
                    dam_gen === nothing ? 0 : dam_gen
                )
                generation[id] = parent_gen + 1
                changed = true
            end
        end
    end
    
    return isempty(generation) ? 1 : maximum(values(generation)) + 1
end

# 显示
function Base.show(io::IO, ped::Pedigree)
    n = length(ped.ids)
    n_founders = length(ped.founder_ids)
    print(io, "Pedigree($n individuals, $n_founders founders, $(ped.n_generations) generations)")
end

"""
    n_individuals(ped::Pedigree) -> Int

返回谱系中的个体数量。
"""
n_individuals(ped::Pedigree) = length(ped.ids)

"""
    n_founders(ped::Pedigree) -> Int

返回谱系中的奠基者数量。
"""
n_founders(ped::Pedigree) = length(ped.founder_ids)

"""
    is_founder(ped::Pedigree, id::String) -> Bool

判断个体是否为奠基者。
"""
is_founder(ped::Pedigree, id::String) = id in ped.founder_ids

"""
    get_parents(ped::Pedigree, id::String) -> Tuple

获取个体的父母ID。

# 返回值
`(sire, dam)` 元组，如果亲本未知则为 `nothing`
"""
function get_parents(ped::Pedigree, id::String)
    idx = findfirst(==(id), ped.ids)
    if idx === nothing
        throw(ArgumentError("个体 '$id' 不在谱系中"))
    end
    return (ped.sires[idx], ped.dams[idx])
end

"""
    get_offspring(ped::Pedigree, id::String) -> Vector{String}

获取个体的所有后代ID。
"""
function get_offspring(ped::Pedigree, id::String)
    offspring = String[]
    for (i, ind_id) in enumerate(ped.ids)
        if ped.sires[i] == id || ped.dams[i] == id
            push!(offspring, ind_id)
        end
    end
    return offspring
end

"""
    RelationshipMatrix

亲缘关系矩阵的封装。

# 字段
- `matrix::Matrix{Float64}`: 关系矩阵
- `ids::Vector{String}`: 对应的个体ID
- `type::Symbol`: 矩阵类型

# 支持的类型
- `:kinship`: Malécot亲缘系数矩阵 Φ（IBD概率）
- `:additive`: 加性遗传关系矩阵 A = 2Φ（Henderson）
- `:genomic`: 基因组关系矩阵 G（VanRaden/Yang）
- `:dominance`: 显性关系矩阵 D

# 中文说明
亲缘关系矩阵（Relationship Matrix）量化个体间的遗传相似程度。

亲缘系数（Kinship Coefficient）Φ_ij：
两个体随机各取一个等位基因，这两个等位基因IBD的概率。

加性遗传关系矩阵 A：
A_ij = 2 × Φ_ij
对角线元素 A_ii = 1 + F_i（其中F是近交系数）

# 示例
```julia
# 从谱系计算A矩阵
A = recursive_kinship(pedigree)

# 从基因型计算G矩阵
G = grm_vanraden_yang(genotypes, allele_freqs)
```

参见: [`recursive_kinship`](@ref), [`grm_vanraden_yang`](@ref)
"""
struct RelationshipMatrix
    matrix::Matrix{Float64}
    ids::Vector{String}
    type::Symbol
    
    function RelationshipMatrix(matrix::AbstractMatrix{<:Real}, 
                               ids::Vector{String},
                               type::Symbol)
        mat = Matrix{Float64}(matrix)
        n = length(ids)
        
        if size(mat, 1) != n || size(mat, 2) != n
            throw(ArgumentError("矩阵维度 $(size(mat)) 与个体数 $n 不匹配"))
        end
        
        if !(type in [:kinship, :additive, :genomic, :dominance, :custom])
            @warn "未知的关系矩阵类型: $type"
        end
        
        new(mat, ids, type)
    end
end

# 显示
function Base.show(io::IO, rm::RelationshipMatrix)
    n = length(rm.ids)
    print(io, "RelationshipMatrix(:$(rm.type), $n × $n)")
end

# 索引操作
Base.size(rm::RelationshipMatrix) = size(rm.matrix)
Base.getindex(rm::RelationshipMatrix, i::Int, j::Int) = rm.matrix[i, j]

"""
    get_relationship(rm::RelationshipMatrix, id1::String, id2::String) -> Float64

获取两个个体间的关系系数。
"""
function get_relationship(rm::RelationshipMatrix, id1::String, id2::String)
    i = findfirst(==(id1), rm.ids)
    j = findfirst(==(id2), rm.ids)
    
    if i === nothing
        throw(ArgumentError("个体 '$id1' 不在矩阵中"))
    end
    if j === nothing
        throw(ArgumentError("个体 '$id2' 不在矩阵中"))
    end
    
    return rm.matrix[i, j]
end

"""
    inbreeding_from_A(rm::RelationshipMatrix) -> Vector{Float64}

从加性关系矩阵提取近交系数。

F_i = A_ii - 1

# 中文说明
近交系数等于个体与自身的亲缘系数减1。
"""
function inbreeding_from_A(rm::RelationshipMatrix)
    if rm.type != :additive
        @warn "该函数设计用于加性关系矩阵，当前类型为 $(rm.type)"
    end
    return [rm.matrix[i, i] - 1.0 for i in 1:length(rm.ids)]
end

"""
    is_positive_definite(rm::RelationshipMatrix) -> Bool

检查关系矩阵是否正定。

# 中文说明
正定性是混合模型方程求解的必要条件。
如果矩阵不正定，可能存在线性相关的个体或数据错误。
"""
function is_positive_definite(rm::RelationshipMatrix)
    try
        cholesky(Symmetric(rm.matrix))
        return true
    catch
        return false
    end
end

"""
    IBDProbabilities

同源同源（IBD）概率结构，用于远交系连锁分析。

# 字段
- `pi_0::Float64`: 共享0个IBD等位基因的概率
- `pi_1::Float64`: 共享1个IBD等位基因的概率
- `pi_2::Float64`: 共享2个IBD等位基因的概率

# 中文说明
对于有亲缘关系的个体对，在特定位点上的IBD共享状态：
- π₀: 0个等位基因IBD
- π₁: 1个等位基因IBD
- π₂: 2个等位基因IBD（对于二倍体）

这些概率满足：π₀ + π₁ + π₂ = 1

# 亲缘关系对应的期望IBD
- 同卵双生：π₂ = 1
- 同胞兄妹：π₀ = 0.25, π₁ = 0.50, π₂ = 0.25
- 半同胞：π₀ = 0.50, π₁ = 0.50, π₂ = 0
- 亲子：π₀ = 0, π₁ = 1, π₂ = 0

参见: [`haseman_elston_regression`](@ref)
"""
struct IBDProbabilities
    pi_0::Float64
    pi_1::Float64
    pi_2::Float64
    
    function IBDProbabilities(pi_0::Real, pi_1::Real, pi_2::Real)
        total = pi_0 + pi_1 + pi_2
        if !isapprox(total, 1.0, atol=1e-6)
            throw(ArgumentError("IBD概率之和必须等于1，当前为 $total"))
        end
        if pi_0 < 0 || pi_1 < 0 || pi_2 < 0
            throw(ArgumentError("IBD概率不能为负"))
        end
        new(Float64(pi_0), Float64(pi_1), Float64(pi_2))
    end
end

"""
    expected_ibd(ibd::IBDProbabilities) -> Float64

计算期望IBD共享数（π̂）。

π̂ = π₁/2 + π₂

# 中文说明
期望IBD共享数，取值范围[0, 1]。
"""
expected_ibd(ibd::IBDProbabilities) = ibd.pi_1 / 2 + ibd.pi_2

"""
    Δ_coefficients

Jacquard的9个身份系数，用于精确描述任意亲缘关系。

# 中文说明
Δ系数提供了比简单IBD概率更精确的亲缘关系描述：
- Δ₁-Δ₉: 9种可能的同源状态概率

对于大多数应用，简化的Φ（亲缘系数）和Δ₇（同胞共享显性效应）已足够。
"""
struct Δ_coefficients
    delta::NTuple{9, Float64}
    
    function Δ_coefficients(delta::NTuple{9, Float64})
        total = sum(delta)
        if !isapprox(total, 1.0, atol=1e-6)
            throw(ArgumentError("Δ系数之和必须等于1"))
        end
        new(delta)
    end
end

function Base.show(io::IO, d::Δ_coefficients)
    print(io, "Δ_coefficients(Δ₇=$(round(d.delta[7], digits=4)), ...)")
end

"""
    kinship_from_delta(d::Δ_coefficients) -> Float64

从Δ系数计算亲缘系数Φ。

Φ = Δ₁ + Δ₃/2 + Δ₅/2 + Δ₇/4 + Δ₈/2
"""
function kinship_from_delta(d::Δ_coefficients)
    return d.delta[1] + d.delta[3]/2 + d.delta[5]/2 + d.delta[7]/4 + d.delta[8]/2
end
