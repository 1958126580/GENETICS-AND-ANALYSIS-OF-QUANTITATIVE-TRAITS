# =============================================================================
# PopulationTypes.jl - 群体与个体类型定义
# =============================================================================
#
# 定义群体遗传学中的核心数据结构：个体、群体、性状数据
#
# 对应章节: Walsh 2nd Ed, Chapters 1, 4
# =============================================================================

"""
    AbstractPopulation

群体的抽象基类型。

# 中文说明
所有群体类型的基类，支持不同的实验设计和群体结构。
"""
abstract type AbstractPopulation end

"""
    TraitData

表型性状数据的容器。

# 字段
- `name::String`: 性状名称
- `values::Vector{Union{Float64, Missing}}`: 各个体的表型值

# 中文说明
存储单个数量性状的测量值。缺失值用 `missing` 表示。

# 示例
```julia
trait = TraitData("plant_height", [125.3, 118.7, missing, 132.1])
```
"""
struct TraitData
    name::String
    values::Vector{Union{Float64, Missing}}
    
    function TraitData(name::String, values::AbstractVector)
        # 转换为允许missing的向量
        converted = Vector{Union{Float64, Missing}}(undef, length(values))
        for (i, v) in enumerate(values)
            if ismissing(v) || (v isa Number && isnan(v))
                converted[i] = missing
            else
                converted[i] = Float64(v)
            end
        end
        new(name, converted)
    end
end

Base.length(t::TraitData) = length(t.values)
Base.getindex(t::TraitData, i::Int) = t.values[i]

function Base.show(io::IO, t::TraitData)
    n = length(t.values)
    n_missing = count(ismissing, t.values)
    print(io, "TraitData(\"$(t.name)\", n=$n, missing=$n_missing)")
end

"""
    n_observed(trait::TraitData) -> Int

返回有有效观测值的个体数。

# 中文说明
统计非缺失表型值的数量。
"""
n_observed(trait::TraitData) = count(!ismissing, trait.values)

"""
    trait_mean(trait::TraitData) -> Float64

计算性状的均值（忽略缺失值）。

# 中文说明
计算性状均值，自动排除缺失值。
"""
function trait_mean(trait::TraitData)
    valid = skipmissing(trait.values)
    return mean(valid)
end

"""
    trait_var(trait::TraitData) -> Float64

计算性状的方差（忽略缺失值）。

# 中文说明
计算性状方差，自动排除缺失值。
"""
function trait_var(trait::TraitData)
    valid = skipmissing(trait.values)
    return var(valid)
end

"""
    Individual

单个个体的完整遗传和表型信息。

# 字段
- `id::String`: 个体唯一标识符
- `sire::Union{String, Nothing}`: 父本ID（如未知则为nothing）
- `dam::Union{String, Nothing}`: 母本ID（如未知则为nothing）
- `sex::Union{Char, Nothing}`: 性别（'M'/'F'/nothing）
- `genotypes::Dict{String, Int}`: 各标记的基因型（标记名 => 0/1/2编码）
- `phenotypes::Dict{String, Union{Float64, Missing}}`: 各性状的表型值

# 中文说明
个体（Individual）是群体遗传分析的基本单位。
每个个体包含其系谱信息、基因型数据和表型测量。

基因型编码约定：
- 0: 纯合子 AA
- 1: 杂合子 Aa
- 2: 纯合子 aa
- -9: 缺失

# 示例
```julia
ind = Individual(
    "Ind001",
    "Sire01", "Dam01",
    'M',
    Dict("M1" => 0, "M2" => 1, "M3" => 2),
    Dict("height" => 125.3, "weight" => 45.2)
)
```
"""
struct Individual
    id::String
    sire::Union{String, Nothing}
    dam::Union{String, Nothing}
    sex::Union{Char, Nothing}
    genotypes::Dict{String, Int}
    phenotypes::Dict{String, Union{Float64, Missing}}
    
    function Individual(id::String, 
                       sire::Union{String, Nothing}=nothing,
                       dam::Union{String, Nothing}=nothing,
                       sex::Union{Char, Nothing}=nothing,
                       genotypes::Dict{String, Int}=Dict{String, Int}(),
                       phenotypes::Dict=Dict{String, Union{Float64, Missing}}())
        new(id, sire, dam, sex, genotypes, phenotypes)
    end
end

function Base.show(io::IO, ind::Individual)
    n_geno = length(ind.genotypes)
    n_pheno = length(ind.phenotypes)
    print(io, "Individual(\"$(ind.id)\", geno=$n_geno, pheno=$n_pheno)")
end

"""
    Population{T} <: AbstractPopulation

群体的主数据结构，优化用于列式访问（矩阵形式存储基因型）。

# 字段
- `name::String`: 群体名称
- `individuals::Vector{String}`: 个体ID列表
- `markers::Vector{String}`: 标记名称列表
- `genotypes::Matrix{T}`: 基因型矩阵（行=个体，列=标记）
- `traits::Dict{String, Vector{Union{Float64, Missing}}}`: 表型数据字典
- `design::Symbol`: 实验设计类型（:F2, :BC, :RIL, :DH, :Outbred等）
- `metadata::Dict{Symbol, Any}`: 额外元数据

# 类型参数
- `T`: 基因型编码类型，通常为 `Int` 或 `Float64`

# 中文说明
群体（Population）是一组具有共同遗传背景的个体集合。
基因型数据以矩阵形式存储，便于高效的向量化计算。

支持的实验设计：
- `:F2`: F2代杂交群体
- `:BC`: 回交群体
- `:RIL`: 重组近交系
- `:DH`: 双单倍体
- `:Outbred`: 远交群体

# 示例
```julia
# 创建F2群体
pop = Population(
    "F2_Cross",
    ["Ind1", "Ind2", "Ind3"],
    ["M1", "M2", "M3"],
    [0 1 2; 1 1 0; 2 0 1],  # 3个个体 × 3个标记
    Dict("height" => [125.3, 118.7, 132.1]),
    :F2
)

n_individuals(pop)  # 3
n_markers(pop)      # 3
```

参见: [`Individual`](@ref), [`TraitData`](@ref)
"""
struct Population{T} <: AbstractPopulation
    name::String
    individuals::Vector{String}
    markers::Vector{String}
    genotypes::Matrix{T}
    traits::Dict{String, Vector{Union{Float64, Missing}}}
    design::Symbol
    metadata::Dict{Symbol, Any}
    
    function Population(name::String,
                       individuals::Vector{String},
                       markers::Vector{String},
                       genotypes::Matrix{T},
                       traits::Dict{String, <:AbstractVector}=Dict{String, Vector{Union{Float64, Missing}}}(),
                       design::Symbol=:Unknown;
                       metadata::Dict{Symbol, Any}=Dict{Symbol, Any}()) where T
        # 验证维度
        n_ind = length(individuals)
        n_mark = length(markers)
        
        if size(genotypes, 1) != n_ind
            throw(ArgumentError("基因型矩阵行数 ($(size(genotypes, 1))) 与个体数 ($n_ind) 不匹配"))
        end
        if size(genotypes, 2) != n_mark
            throw(ArgumentError("基因型矩阵列数 ($(size(genotypes, 2))) 与标记数 ($n_mark) 不匹配"))
        end
        
        # 验证表型数据维度
        for (trait_name, values) in traits
            if length(values) != n_ind
                throw(ArgumentError("性状 '$trait_name' 的观测数 ($(length(values))) 与个体数 ($n_ind) 不匹配"))
            end
        end
        
        # 转换表型数据类型
        converted_traits = Dict{String, Vector{Union{Float64, Missing}}}()
        for (trait_name, values) in traits
            converted = Vector{Union{Float64, Missing}}(undef, length(values))
            for (i, v) in enumerate(values)
                if ismissing(v) || (v isa Number && isnan(v))
                    converted[i] = missing
                else
                    converted[i] = Float64(v)
                end
            end
            converted_traits[trait_name] = converted
        end
        
        new{T}(name, individuals, markers, genotypes, converted_traits, design, metadata)
    end
end

# 显示
function Base.show(io::IO, pop::Population)
    n_ind = n_individuals(pop)
    n_mark = n_markers(pop)
    n_trait = n_traits(pop)
    print(io, "Population(\"$(pop.name)\", $n_ind individuals, $n_mark markers, $n_trait traits, design=$(pop.design))")
end

"""
    n_individuals(pop::Population) -> Int

返回群体中的个体数量。

# 中文说明
获取群体样本量。
"""
n_individuals(pop::Population) = length(pop.individuals)

"""
    n_markers(pop::Population) -> Int

返回群体中的标记数量。

# 中文说明
获取基因型标记数目。
"""
n_markers(pop::Population) = length(pop.markers)

"""
    n_traits(pop::Population) -> Int

返回群体中的性状数量。

# 中文说明
获取表型性状数目。
"""
n_traits(pop::Population) = length(pop.traits)

"""
    get_genotypes(pop::Population, marker_idx::Int) -> Vector

获取特定标记的所有个体基因型。

# 中文说明
提取单个标记位点的全部基因型数据。
"""
get_genotypes(pop::Population, marker_idx::Int) = pop.genotypes[:, marker_idx]

"""
    get_genotypes(pop::Population, marker_name::String) -> Vector

通过标记名称获取基因型。
"""
function get_genotypes(pop::Population, marker_name::String)
    idx = findfirst(==(marker_name), pop.markers)
    if idx === nothing
        throw(ArgumentError("标记 '$marker_name' 不存在于群体中"))
    end
    return pop.genotypes[:, idx]
end

"""
    get_phenotypes(pop::Population, trait_name::String) -> Vector

获取特定性状的所有个体表型值。

# 中文说明
提取单个性状的全部表型数据。
"""
function get_phenotypes(pop::Population, trait_name::String)
    if !haskey(pop.traits, trait_name)
        throw(ArgumentError("性状 '$trait_name' 不存在于群体中"))
    end
    return pop.traits[trait_name]
end

"""
    get_phenotypes(pop::Population, trait_name::Symbol) -> Vector

通过Symbol获取表型数据。
"""
get_phenotypes(pop::Population, trait_name::Symbol) = get_phenotypes(pop, String(trait_name))

"""
    get_individual_genotypes(pop::Population, ind_idx::Int) -> Vector

获取特定个体的所有标记基因型。

# 中文说明
提取单个个体的全部基因型数据。
"""
get_individual_genotypes(pop::Population, ind_idx::Int) = pop.genotypes[ind_idx, :]

"""
    trait_names(pop::Population) -> Vector{String}

返回所有性状名称。

# 中文说明
列出群体中所有表型性状的名称。
"""
trait_names(pop::Population) = collect(keys(pop.traits))

"""
    marker_names(pop::Population) -> Vector{String}

返回所有标记名称。

# 中文说明
列出群体中所有遗传标记的名称。
"""
marker_names(pop::Population) = pop.markers

"""
    subset_individuals(pop::Population, indices::Vector{Int}) -> Population

创建包含指定个体的子群体。

# 中文说明
根据索引提取部分个体形成新群体。
"""
function subset_individuals(pop::Population{T}, indices::Vector{Int}) where T
    new_inds = pop.individuals[indices]
    new_genos = pop.genotypes[indices, :]
    new_traits = Dict{String, Vector{Union{Float64, Missing}}}()
    for (name, values) in pop.traits
        new_traits[name] = values[indices]
    end
    return Population(pop.name, new_inds, pop.markers, new_genos, new_traits, pop.design; 
                     metadata=pop.metadata)
end

"""
    subset_markers(pop::Population, indices::Vector{Int}) -> Population

创建包含指定标记的子群体。

# 中文说明
根据索引提取部分标记形成新群体。
"""
function subset_markers(pop::Population{T}, indices::Vector{Int}) where T
    new_markers = pop.markers[indices]
    new_genos = pop.genotypes[:, indices]
    return Population(pop.name, pop.individuals, new_markers, new_genos, pop.traits, pop.design;
                     metadata=pop.metadata)
end

"""
    CrossDesign

杂交设计的详细描述。

# 字段
- `design_type::Symbol`: 设计类型
- `n_parents::Int`: 亲本数量
- `n_generations::Int`: 杂交代数
- `description::String`: 设计描述

# 中文说明
描述产生群体的杂交方案，影响遗传模型的选择。
"""
struct CrossDesign
    design_type::Symbol
    n_parents::Int
    n_generations::Int
    description::String
end

# 预定义的杂交设计
const BACKCROSS_DESIGN = CrossDesign(:BC, 2, 1, "Backcross to one parent")
const F2_DESIGN = CrossDesign(:F2, 2, 2, "F2 intercross")
const RIL_DESIGN = CrossDesign(:RIL, 2, 8, "Recombinant Inbred Lines (SSD)")
const DH_DESIGN = CrossDesign(:DH, 2, 1, "Doubled Haploid")
