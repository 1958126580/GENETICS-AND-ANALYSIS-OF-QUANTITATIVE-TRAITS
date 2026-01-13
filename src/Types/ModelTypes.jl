# =============================================================================
# ModelTypes.jl - 遗传模型类型定义
# =============================================================================
#
# 定义QTL分析中使用的各类遗传效应模型
#
# 对应章节: Walsh 2nd Ed, Chapters 4, 18-19
# =============================================================================

"""
    AbstractModel

所有遗传模型的抽象基类型。

# 中文说明
遗传模型的基类，用于类型分派。
"""
abstract type AbstractModel end

"""
    QTLModel <: AbstractModel

QTL效应模型的抽象基类型。

# 中文说明
QTL模型基类，所有QTL效应模型都继承自此类型。
"""
abstract type QTLModel <: AbstractModel end

"""
    AdditiveModel <: QTLModel

纯加性效应模型。

# 模型方程
```math
y = \\mu + a \\cdot x_a + \\epsilon
```

其中：
- μ: 群体均值
- a: 加性效应
- x_a: 加性设计变量（-1/0/1 编码）

# 字段
- `mu::Float64`: 群体均值
- `a::Float64`: 加性效应

# 中文说明
加性模型假设基因型效应呈线性关系：
- AA: μ + a
- Aa: μ
- aa: μ - a

这是最简单的遗传效应模型，适用于无显性效应的情况。

# 公式来源
Walsh 2nd Ed, Eq 4.2 - 加性效应定义

# 示例
```julia
model = AdditiveModel(100.0, 5.0)  # 均值100，加性效应5
```
"""
struct AdditiveModel <: QTLModel
    mu::Float64
    a::Float64
    
    function AdditiveModel(mu::Real, a::Real)
        new(Float64(mu), Float64(a))
    end
end

function Base.show(io::IO, m::AdditiveModel)
    print(io, "AdditiveModel(μ=$(round(m.mu, digits=3)), a=$(round(m.a, digits=3)))")
end

"""
    DominanceModel <: QTLModel

加性-显性效应模型。

# 模型方程
```math
y = \\mu + a \\cdot x_a + d \\cdot x_d + \\epsilon
```

其中：
- μ: 群体均值（杂合子均值）
- a: 加性效应
- d: 显性效应
- x_a: 加性设计变量
- x_d: 显性设计变量（0/1 编码）

# 字段
- `mu::Float64`: 群体均值
- `a::Float64`: 加性效应
- `d::Float64`: 显性效应

# 中文说明
加性-显性模型考虑了杂合子与纯合子中值的偏离：
- AA: μ + a
- Aa: μ + d
- aa: μ - a

显性度 = d/a:
- d = 0: 无显性（加性遗传）
- d = a: 完全显性
- d > a: 超显性
- 0 < d < a: 部分显性

# 公式来源
Walsh 2nd Ed, Eq 4.8 - Fisher分解

# 示例
```julia
model = DominanceModel(100.0, 5.0, 2.0)  # 加性5，显性2
```
"""
struct DominanceModel <: QTLModel
    mu::Float64
    a::Float64
    d::Float64
    
    function DominanceModel(mu::Real, a::Real, d::Real)
        new(Float64(mu), Float64(a), Float64(d))
    end
end

function Base.show(io::IO, m::DominanceModel)
    print(io, "DominanceModel(μ=$(round(m.mu, digits=3)), a=$(round(m.a, digits=3)), d=$(round(m.d, digits=3)))")
end

"""
    dominance_ratio(model::DominanceModel) -> Float64

计算显性度 d/|a|。

# 中文说明
显性度衡量显性效应相对于加性效应的大小。
"""
function dominance_ratio(model::DominanceModel)
    if model.a ≈ 0
        return Inf
    end
    return model.d / abs(model.a)
end

"""
    EpistaticModel <: QTLModel

双位点上位性效应模型。

# 模型方程
```math
y = \\mu + a_1 x_{a1} + a_2 x_{a2} + d_1 x_{d1} + d_2 x_{d2} + i_{aa} x_{a1}x_{a2} + i_{ad} x_{a1}x_{d2} + i_{da} x_{d1}x_{a2} + i_{dd} x_{d1}x_{d2} + \\epsilon
```

# 字段
- `mu::Float64`: 群体均值
- `a1::Float64`: 位点1的加性效应
- `a2::Float64`: 位点2的加性效应
- `d1::Float64`: 位点1的显性效应
- `d2::Float64`: 位点2的显性效应
- `i_aa::Float64`: 加性×加性上位性
- `i_ad::Float64`: 加性×显性上位性
- `i_da::Float64`: 显性×加性上位性
- `i_dd::Float64`: 显性×显性上位性

# 中文说明
上位性模型描述两个位点间的非加性互作效应。

上位性类型：
- i_aa (加性×加性): 两个位点的加性效应相互影响
- i_ad (加性×显性): 位点1的加性效应受位点2杂合状态影响
- i_da (显性×加性): 位点1的杂合状态受位点2加性效应影响
- i_dd (显性×显性): 两个位点的杂合状态相互影响

# 公式来源
Walsh 2nd Ed, Chapter 5 - 上位性方差分解
Cockerham (1954) 正交模型

# 示例
```julia
model = EpistaticModel(
    100.0,      # mu
    3.0, 2.0,   # a1, a2
    1.0, 0.5,   # d1, d2
    0.5, 0.2, 0.2, 0.1  # i_aa, i_ad, i_da, i_dd
)
```
"""
struct EpistaticModel <: QTLModel
    mu::Float64
    a1::Float64
    a2::Float64
    d1::Float64
    d2::Float64
    i_aa::Float64
    i_ad::Float64
    i_da::Float64
    i_dd::Float64
    
    function EpistaticModel(mu::Real, a1::Real, a2::Real, d1::Real, d2::Real,
                           i_aa::Real, i_ad::Real, i_da::Real, i_dd::Real)
        new(Float64(mu), Float64(a1), Float64(a2), Float64(d1), Float64(d2),
            Float64(i_aa), Float64(i_ad), Float64(i_da), Float64(i_dd))
    end
end

function Base.show(io::IO, m::EpistaticModel)
    print(io, "EpistaticModel(μ=$(round(m.mu, digits=2)), a1=$(round(m.a1, digits=2)), a2=$(round(m.a2, digits=2)), ...)")
end

"""
    CompositeModel <: QTLModel

复合区间作图(CIM)模型，包含标记协因子。

# 模型方程
```math
y = \\mu + QTL效应 + \\sum_{k} \\beta_k X_k + \\epsilon
```

其中 Xk 是作为背景控制的标记协因子。

# 字段
- `qtl_model::QTLModel`: 基础QTL效应模型
- `cofactors::Vector{Int}`: 协因子标记索引
- `cofactor_effects::Vector{Float64}`: 协因子效应估计值

# 中文说明
复合区间作图通过引入标记协因子来控制背景遗传效应，
从而提高QTL检测的精度和分辨率。

# 公式来源
Walsh 2nd Ed, Chapter 18 - CIM (Zeng, 1994)

# 示例
```julia
base_model = DominanceModel(100.0, 5.0, 2.0)
cim_model = CompositeModel(base_model, [3, 7, 15], [1.2, -0.8, 2.1])
```
"""
struct CompositeModel <: QTLModel
    qtl_model::QTLModel
    cofactors::Vector{Int}
    cofactor_effects::Vector{Float64}
    
    function CompositeModel(qtl_model::QTLModel, 
                           cofactors::Vector{Int},
                           cofactor_effects::Vector{Float64})
        if length(cofactors) != length(cofactor_effects)
            throw(ArgumentError("协因子数量与效应数量不匹配"))
        end
        new(qtl_model, cofactors, cofactor_effects)
    end
end

function Base.show(io::IO, m::CompositeModel)
    n_cof = length(m.cofactors)
    print(io, "CompositeModel($(m.qtl_model), $n_cof cofactors)")
end

"""
    n_cofactors(model::CompositeModel) -> Int

返回模型中的协因子数量。
"""
n_cofactors(model::CompositeModel) = length(model.cofactors)

"""
    MixedModel <: AbstractModel

混合线性模型结构。

# 模型方程
```math
y = X\\beta + Zu + \\epsilon
```

其中：
- y: 观测向量
- X: 固定效应设计矩阵
- β: 固定效应向量
- Z: 随机效应设计矩阵
- u: 随机效应向量，u ~ N(0, Gσ²_u)
- ε: 残差，ε ~ N(0, Rσ²_e)

# 字段
- `X::Matrix{Float64}`: 固定效应设计矩阵
- `Z::Matrix{Float64}`: 随机效应设计矩阵
- `G::Matrix{Float64}`: 随机效应协方差结构
- `R::Matrix{Float64}`: 残差协方差结构

# 中文说明
混合线性模型是数量遗传学的核心工具，用于：
- 育种值预测（BLUP）
- 方差组分估计（REML）
- 基因组选择（GBLUP）

# 公式来源
Walsh 2nd Ed, Chapter 9-10 - 混合模型方程

# 示例
```julia
model = MixedModel(X, Z, A, I)  # A为亲缘矩阵
```

参见: [`mixed_model_equations`](@ref), [`reml_ai_algorithm`](@ref)
"""
struct MixedModel <: AbstractModel
    X::Matrix{Float64}
    Z::Matrix{Float64}
    G::Matrix{Float64}
    R::Matrix{Float64}
    
    function MixedModel(X::AbstractMatrix, Z::AbstractMatrix, 
                       G::AbstractMatrix, R::AbstractMatrix)
        X_f = Matrix{Float64}(X)
        Z_f = Matrix{Float64}(Z)
        G_f = Matrix{Float64}(G)
        R_f = Matrix{Float64}(R)
        
        if size(X_f, 1) != size(Z_f, 1)
            throw(ArgumentError("X和Z的行数必须相同"))
        end
        if size(Z_f, 2) != size(G_f, 1)
            throw(ArgumentError("Z的列数必须与G的维度匹配"))
        end
        if size(X_f, 1) != size(R_f, 1)
            throw(ArgumentError("X/Z的行数必须与R的维度匹配"))
        end
        
        new(X_f, Z_f, G_f, R_f)
    end
end

"""
    VarianceComponentsModel <: AbstractModel

方差组分模型。

# 字段
- `components::Vector{Symbol}`: 方差组分名称
- `variances::Vector{Float64}`: 方差估计值
- `covariances::Matrix{Float64}`: 协方差矩阵（可选）

# 中文说明
用于存储方差组分估计结果。
"""
struct VarianceComponentsModel <: AbstractModel
    components::Vector{Symbol}
    variances::Vector{Float64}
    covariances::Matrix{Float64}
    
    function VarianceComponentsModel(components::Vector{Symbol}, 
                                    variances::Vector{Float64})
        n = length(components)
        if length(variances) != n
            throw(ArgumentError("组分名称数量与方差数量不匹配"))
        end
        new(components, variances, diagm(variances))
    end
    
    function VarianceComponentsModel(components::Vector{Symbol}, 
                                    variances::Vector{Float64},
                                    covariances::Matrix{Float64})
        new(components, variances, covariances)
    end
end

function Base.show(io::IO, m::VarianceComponentsModel)
    parts = ["$(c)=$(round(v, digits=3))" for (c, v) in zip(m.components, m.variances)]
    print(io, "VarianceComponents($(join(parts, ", ")))")
end

"""
    heritability(vc::VarianceComponentsModel) -> Float64

从方差组分计算遗传力。

# 中文说明
遗传力 h² = σ²_A / σ²_P（狭义遗传力）
"""
function heritability(vc::VarianceComponentsModel)
    additive_idx = findfirst(==(:additive), vc.components)
    if additive_idx === nothing
        additive_idx = findfirst(==(:genetic), vc.components)
    end
    
    if additive_idx === nothing
        throw(ArgumentError("未找到加性或遗传方差组分"))
    end
    
    total = sum(vc.variances)
    return vc.variances[additive_idx] / total
end
