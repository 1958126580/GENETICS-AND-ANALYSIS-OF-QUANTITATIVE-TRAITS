# =============================================================================
# GeneticsTypes.jl - 基础遗传学类型定义
# =============================================================================
#
# 定义遗传学中最基本的数据结构：等位基因、基因型、单倍型
# 
# 对应章节: Walsh 2nd Ed, Chapter 4 (Properties of Single Loci)
# =============================================================================

"""
    Allele{T}

表示一个遗传等位基因的结构体。

等位基因可以用整数、字符或字符串编码。类型参数T决定了编码方式。

# 字段
- `id::T`: 等位基因的唯一标识符

# 中文说明
等位基因（Allele）是位于染色体特定位置（位点）的一个基因的不同形式。
例如，控制豌豆种子颜色的基因可能有黄色(Y)和绿色(y)两个等位基因。

# 示例
```julia
# 整数编码
a1 = Allele(1)
a2 = Allele(2)

# 字符编码
aA = Allele('A')
aa = Allele('a')

# 字符串编码（用于复杂等位基因名称）
a_del = Allele("DEL")
a_ins = Allele("INS_5bp")
```
"""
struct Allele{T}
    id::T
    
    function Allele(id::T) where T
        new{T}(id)
    end
end

# 显示格式化
Base.show(io::IO, a::Allele) = print(io, a.id)

# 相等性比较
Base.isequal(a1::Allele, a2::Allele) = isequal(a1.id, a2.id)
Base.:(==)(a1::Allele, a2::Allele) = a1.id == a2.id
Base.hash(a::Allele, h::UInt) = hash(a.id, h)

"""
    Genotype{T}

表示特定位点的基因型，由两个等位基因组成（二倍体生物）。

# 字段
- `allele1::Allele{T}`: 第一个等位基因（通常来自父本）
- `allele2::Allele{T}`: 第二个等位基因（通常来自母本）

# 中文说明
基因型（Genotype）是个体在特定位点上拥有的等位基因组合。
对于二倍体生物，每个位点有两个等位基因，一个来自父本，一个来自母本。

常见基因型类型：
- 纯合子（Homozygote）: 两个等位基因相同，如AA或aa
- 杂合子（Heterozygote）: 两个等位基因不同，如Aa

# 示例
```julia
# 创建纯合子
geno_AA = Genotype('A', 'A')
geno_aa = Genotype('a', 'a')

# 创建杂合子
geno_Aa = Genotype('A', 'a')

# 检查杂合性
is_heterozygous(geno_Aa)  # true
is_homozygous(geno_AA)    # true
```

参见: [`Allele`](@ref), [`is_homozygous`](@ref), [`is_heterozygous`](@ref)
"""
struct Genotype{T}
    allele1::Allele{T}
    allele2::Allele{T}
    
    # 从两个等位基因ID直接构造
    function Genotype(a1::T, a2::T) where T
        new{T}(Allele(a1), Allele(a2))
    end
    
    # 从两个Allele对象构造
    function Genotype(a1::Allele{T}, a2::Allele{T}) where T
        new{T}(a1, a2)
    end
end

# 显示格式化
Base.show(io::IO, g::Genotype) = print(io, "$(g.allele1)/$(g.allele2)")

# 相等性（考虑等位基因顺序无关）
function Base.isequal(g1::Genotype, g2::Genotype)
    (isequal(g1.allele1, g2.allele1) && isequal(g1.allele2, g2.allele2)) ||
    (isequal(g1.allele1, g2.allele2) && isequal(g1.allele2, g2.allele1))
end

function Base.:(==)(g1::Genotype, g2::Genotype)
    (g1.allele1 == g2.allele1 && g1.allele2 == g2.allele2) ||
    (g1.allele1 == g2.allele2 && g1.allele2 == g2.allele1)
end

"""
    is_homozygous(g::Genotype) -> Bool

判断基因型是否为纯合子。

# 参数
- `g::Genotype`: 待检测的基因型

# 返回值
- `Bool`: 如果两个等位基因相同返回true

# 中文说明
纯合子指两个等位基因完全相同的基因型，如AA或aa。

# 示例
```julia
g = Genotype('A', 'A')
is_homozygous(g)  # true
```
"""
is_homozygous(g::Genotype) = isequal(g.allele1, g.allele2)

"""
    is_heterozygous(g::Genotype) -> Bool

判断基因型是否为杂合子。

# 中文说明
杂合子指两个等位基因不同的基因型，如Aa。

# 示例
```julia
g = Genotype('A', 'a')
is_heterozygous(g)  # true
```
"""
is_heterozygous(g::Genotype) = !is_homozygous(g)

"""
    Haplotype{T}

表示一条染色体上连续多个位点的等位基因组合。

# 字段
- `alleles::Vector{Allele{T}}`: 沿染色体排列的等位基因序列
- `marker_names::Vector{String}`: 对应的标记名称（可选）

# 中文说明
单倍型（Haplotype）是同一条染色体上多个连锁位点的等位基因组合。
在连锁分析和QTL作图中，追踪单倍型的遗传传递至关重要。

# 示例
```julia
# 创建一个包含3个位点的单倍型
hap = Haplotype([Allele('A'), Allele('B'), Allele('C')], 
                ["marker1", "marker2", "marker3"])
```
"""
struct Haplotype{T}
    alleles::Vector{Allele{T}}
    marker_names::Vector{String}
    
    function Haplotype(alleles::Vector{Allele{T}}, names::Vector{String}=String[]) where T
        if !isempty(names) && length(alleles) != length(names)
            throw(ArgumentError("等位基因数量与标记名称数量不匹配"))
        end
        new{T}(alleles, names)
    end
    
    # 便捷构造：从原始值创建
    function Haplotype(allele_ids::Vector{T}, names::Vector{String}=String[]) where T
        alleles = [Allele(id) for id in allele_ids]
        Haplotype(alleles, names)
    end
end

Base.length(h::Haplotype) = length(h.alleles)
Base.getindex(h::Haplotype, i::Int) = h.alleles[i]
Base.iterate(h::Haplotype) = iterate(h.alleles)
Base.iterate(h::Haplotype, state) = iterate(h.alleles, state)

function Base.show(io::IO, h::Haplotype)
    allele_str = join([string(a.id) for a in h.alleles], "-")
    print(io, "Haplotype: $allele_str")
end

"""
    DiploidGenome{T}

二倍体生物的完整基因组表示，包含两套单倍型。

# 字段
- `maternal::Haplotype{T}`: 母本来源的单倍型
- `paternal::Haplotype{T}`: 父本来源的单倍型

# 中文说明
二倍体基因组由两条同源染色体组成，分别来自父母双方。
该结构用于追踪完整的遗传信息传递。
"""
struct DiploidGenome{T}
    maternal::Haplotype{T}
    paternal::Haplotype{T}
    
    function DiploidGenome(maternal::Haplotype{T}, paternal::Haplotype{T}) where T
        if length(maternal) != length(paternal)
            throw(ArgumentError("母本和父本单倍型长度必须相同"))
        end
        new{T}(maternal, paternal)
    end
end

"""
    get_genotype(genome::DiploidGenome, locus::Int) -> Genotype

获取二倍体基因组在特定位点的基因型。

# 中文说明
从完整基因组中提取单个位点的基因型信息。
"""
function get_genotype(genome::DiploidGenome{T}, locus::Int) where T
    Genotype(genome.maternal[locus], genome.paternal[locus])
end

"""
    genotype_to_numeric(g::Genotype, ref_allele) -> Int

将基因型转换为数值编码。

编码规则（假设ref_allele为参考等位基因）：
- 0: 纯合参考 (ref/ref)
- 1: 杂合 (ref/alt)
- 2: 纯合突变 (alt/alt)

# 中文说明
将基因型转换为0/1/2编码，用于数值计算。
0表示纯合参考，1表示杂合，2表示纯合替代。
"""
function genotype_to_numeric(g::Genotype{T}, ref_allele::T) where T
    count = 0
    if g.allele1.id != ref_allele
        count += 1
    end
    if g.allele2.id != ref_allele
        count += 1
    end
    return count
end

"""
    genotype_to_numeric(g::Genotype{Int}) -> Int

简化版本：对于整数编码的基因型，直接返回0/1/2编码。
假设等位基因编码为0和1。
"""
function genotype_to_numeric(g::Genotype{Int})
    return g.allele1.id + g.allele2.id
end

"""
    numeric_to_additive(code::Int) -> Float64

将0/1/2编码转换为加性编码系数。

编码转换：
- 0 → -1.0 (纯合子，效应为 -a)
- 1 → 0.0  (杂合子，无加性贡献)
- 2 → +1.0 (纯合子，效应为 +a)

# 中文说明
加性编码用于回归分析中估计加性效应。
基于Fisher的加性-显性分解模型。
"""
function numeric_to_additive(code::Int)
    return Float64(code - 1)
end

"""
    numeric_to_dominance(code::Int) -> Float64

将0/1/2编码转换为显性编码系数。

编码转换：
- 0 → 0.0 (纯合子AA)
- 1 → 1.0 (杂合子Aa，显性效应)
- 2 → 0.0 (纯合子aa)

# 中文说明
显性编码用于估计显性效应(d)。
杂合子编码为1，纯合子编码为0。
"""
function numeric_to_dominance(code::Int)
    return code == 1 ? 1.0 : 0.0
end
