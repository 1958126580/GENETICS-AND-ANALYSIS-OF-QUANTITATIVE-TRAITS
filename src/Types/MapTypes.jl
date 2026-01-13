# =============================================================================
# MapTypes.jl - 遗传图谱类型定义
# =============================================================================
#
# 定义遗传图谱相关的数据结构：标记、遗传图谱
#
# 对应章节: Walsh 2nd Ed, Chapter 17 (Mapping Functions and Linkage Maps)
# =============================================================================

"""
    AbstractMarker

遗传标记的抽象基类型。

# 中文说明
所有遗传标记类型的基类，用于多态性分析和连锁作图。
"""
abstract type AbstractMarker end

"""
    AbstractMap

遗传图谱的抽象基类型。

# 中文说明
所有遗传图谱类型的基类，用于存储标记位置信息。
"""
abstract type AbstractMap end

"""
    Marker{T} <: AbstractMarker

遗传标记的完整表示。

# 字段
- `name::String`: 标记的唯一名称/ID
- `chromosome::String`: 所在染色体
- `position::Float64`: 在染色体上的位置（通常以cM为单位）
- `alleles::Vector{Allele{T}}`: 该标记的所有已知等位基因
- `physical_position::Union{Int, Nothing}`: 物理位置（bp，可选）

# 中文说明
遗传标记是基因组上可识别的多态性位点，用于连锁分析和QTL作图。
常见类型包括：
- SSR（简单重复序列）
- SNP（单核苷酸多态性）
- RFLP（限制性片段长度多态性）

位置通常以厘摩（centiMorgan, cM）为单位，反映标记间的重组率。

# 示例
```julia
# 创建一个SNP标记
marker = Marker(
    "rs12345",           # 标记名称
    "Chr1",              # 染色体
    25.5,                # 位置 (cM)
    [Allele('A'), Allele('G')],  # 等位基因
    1234567              # 物理位置 (bp)
)
```

参见: [`GeneticMap`](@ref), [`Allele`](@ref)
"""
struct Marker{T} <: AbstractMarker
    name::String
    chromosome::String
    position::Float64        # cM
    alleles::Vector{Allele{T}}
    physical_position::Union{Int, Nothing}
    
    function Marker(name::String, chromosome::String, position::Float64,
                   alleles::Vector{Allele{T}}, 
                   physical_pos::Union{Int, Nothing}=nothing) where T
        if position < 0
            throw(ArgumentError("标记位置不能为负数: $position"))
        end
        if isempty(alleles)
            throw(ArgumentError("标记必须至少有一个等位基因"))
        end
        new{T}(name, chromosome, position, alleles, physical_pos)
    end
    
    # 简化构造函数
    function Marker(name::String, chromosome::String, position::Float64)
        new{Char}(name, chromosome, position, Allele{Char}[], nothing)
    end
end

Base.show(io::IO, m::Marker) = print(io, "Marker($(m.name), $(m.chromosome):$(m.position)cM)")

# 标记比较（按位置排序）
function Base.isless(m1::Marker, m2::Marker)
    if m1.chromosome != m2.chromosome
        return m1.chromosome < m2.chromosome
    end
    return m1.position < m2.position
end

"""
    n_alleles(marker::Marker) -> Int

返回标记的等位基因数量。

# 中文说明
获取该标记位点的等位基因数目。
"""
n_alleles(marker::Marker) = length(marker.alleles)

"""
    is_biallelic(marker::Marker) -> Bool

判断标记是否为双等位基因标记（如SNP）。

# 中文说明
判断是否为双等位基因标记，SNP通常是双等位的。
"""
is_biallelic(marker::Marker) = n_alleles(marker) == 2

"""
    GeneticMap <: AbstractMap

遗传图谱，包含一组有序排列的遗传标记。

# 字段
- `name::String`: 图谱名称
- `markers::Vector{Marker}`: 按位置排序的标记列表
- `map_function::Symbol`: 使用的映射函数类型（:haldane 或 :kosambi）
- `chromosomes::Vector{String}`: 图谱中包含的染色体列表

# 中文说明
遗传图谱是通过连锁分析确定的标记线性排列顺序和间距。
图谱距离以厘摩（cM）为单位，1 cM对应约1%的重组率。

常用的映射函数：
- Haldane: 假设无干扰（no interference）
- Kosambi: 考虑正干扰（positive interference）

# 示例
```julia
# 创建遗传图谱
markers = [
    Marker("M1", "Chr1", 0.0),
    Marker("M2", "Chr1", 10.5),
    Marker("M3", "Chr1", 25.0)
]
gmap = GeneticMap("MyMap", markers, :haldane)

# 获取图谱长度
total_length(gmap)  # 返回 25.0
```

参见: [`Marker`](@ref), [`haldane_mapping`](@ref), [`kosambi_mapping`](@ref)
"""
struct GeneticMap <: AbstractMap
    name::String
    markers::Vector{Marker}
    map_function::Symbol
    chromosomes::Vector{String}
    
    function GeneticMap(name::String, markers::Vector{<:Marker}; 
                       map_function::Symbol=:haldane)
        # 按染色体和位置排序
        sorted_markers = sort(markers)
        # 获取唯一染色体列表
        chroms = unique([m.chromosome for m in sorted_markers])
        
        if !(map_function in [:haldane, :kosambi, :carter_falconer])
            throw(ArgumentError("不支持的映射函数: $map_function"))
        end
        
        new(name, sorted_markers, map_function, chroms)
    end
    
    # 无名称的简化构造
    function GeneticMap(markers::Vector{<:Marker}; map_function::Symbol=:haldane)
        GeneticMap("Unnamed", markers; map_function=map_function)
    end
end

# 迭代协议
Base.length(gm::GeneticMap) = length(gm.markers)
Base.getindex(gm::GeneticMap, i::Int) = gm.markers[i]
Base.iterate(gm::GeneticMap) = iterate(gm.markers)
Base.iterate(gm::GeneticMap, state) = iterate(gm.markers, state)
Base.firstindex(gm::GeneticMap) = 1
Base.lastindex(gm::GeneticMap) = length(gm.markers)

function Base.show(io::IO, gm::GeneticMap)
    n_markers = length(gm.markers)
    n_chroms = length(gm.chromosomes)
    print(io, "GeneticMap(\"$(gm.name)\", $n_markers markers, $n_chroms chromosomes)")
end

"""
    n_markers(gmap::GeneticMap) -> Int

返回遗传图谱中的标记数量。

# 中文说明
获取图谱中标记总数。
"""
n_markers(gmap::GeneticMap) = length(gmap.markers)

"""
    n_chromosomes(gmap::GeneticMap) -> Int

返回遗传图谱中的染色体数量。

# 中文说明
获取图谱覆盖的染色体数目。
"""
n_chromosomes(gmap::GeneticMap) = length(gmap.chromosomes)

"""
    total_length(gmap::GeneticMap) -> Float64

计算遗传图谱的总长度（所有染色体之和，cM）。

# 中文说明
计算图谱总遗传长度，单位为厘摩。
"""
function total_length(gmap::GeneticMap)
    if isempty(gmap.markers)
        return 0.0
    end
    
    total = 0.0
    for chrom in gmap.chromosomes
        chrom_markers = filter(m -> m.chromosome == chrom, gmap.markers)
        if length(chrom_markers) >= 2
            total += chrom_markers[end].position - chrom_markers[1].position
        end
    end
    return total
end

"""
    chromosome_map(gmap::GeneticMap, chrom::String) -> GeneticMap

提取指定染色体的子图谱。

# 参数
- `gmap::GeneticMap`: 完整遗传图谱
- `chrom::String`: 目标染色体名称

# 返回值
- `GeneticMap`: 仅包含指定染色体标记的子图谱

# 中文说明
从完整图谱中提取单条染色体的标记。
"""
function chromosome_map(gmap::GeneticMap, chrom::String)
    chrom_markers = filter(m -> m.chromosome == chrom, gmap.markers)
    if isempty(chrom_markers)
        throw(ArgumentError("染色体 $chrom 不存在于图谱中"))
    end
    GeneticMap("$(gmap.name)_$chrom", chrom_markers; map_function=gmap.map_function)
end

"""
    get_marker_interval(gmap::GeneticMap, position::Float64, chrom::String) -> Tuple{Int, Int}

找到包含给定位置的标记区间。

# 返回值
返回一个元组 (left_idx, right_idx)，表示包含该位置的相邻标记索引。

# 中文说明
在区间作图中，需要找到测试位置两侧的侧翼标记。
"""
function get_marker_interval(gmap::GeneticMap, position::Float64, chrom::String)
    chrom_markers = [(i, m) for (i, m) in enumerate(gmap.markers) if m.chromosome == chrom]
    
    if isempty(chrom_markers)
        throw(ArgumentError("染色体 $chrom 不存在于图谱中"))
    end
    
    # 找到position所在的区间
    for i in 1:(length(chrom_markers)-1)
        left_idx, left_marker = chrom_markers[i]
        right_idx, right_marker = chrom_markers[i+1]
        
        if left_marker.position <= position <= right_marker.position
            return (left_idx, right_idx)
        end
    end
    
    # 如果位置超出范围
    if position < chrom_markers[1][2].position
        return (chrom_markers[1][1], chrom_markers[1][1])
    else
        return (chrom_markers[end][1], chrom_markers[end][1])
    end
end

"""
    marker_distance(m1::Marker, m2::Marker) -> Float64

计算两个标记之间的遗传距离（cM）。

# 中文说明
返回两个标记间的图距。标记必须在同一条染色体上。
"""
function marker_distance(m1::Marker, m2::Marker)
    if m1.chromosome != m2.chromosome
        throw(ArgumentError("标记位于不同染色体，无法计算距离"))
    end
    return abs(m1.position - m2.position)
end

"""
    IntervalPosition

表示遗传图谱上的一个连续位置（用于QTL扫描）。

# 字段
- `chromosome::String`: 染色体
- `position::Float64`: 位置（cM）
- `left_marker_idx::Int`: 左侧标记索引
- `right_marker_idx::Int`: 右侧标记索引
- `left_distance::Float64`: 到左侧标记的距离
- `right_distance::Float64`: 到右侧标记的距离

# 中文说明
在区间作图中，假设的QTL位置通常在两个标记之间。
该结构记录QTL位置及其与侧翼标记的关系。
"""
struct IntervalPosition
    chromosome::String
    position::Float64
    left_marker_idx::Int
    right_marker_idx::Int
    left_distance::Float64
    right_distance::Float64
end

function Base.show(io::IO, ip::IntervalPosition)
    print(io, "IntervalPos($(ip.chromosome):$(ip.position)cM)")
end
