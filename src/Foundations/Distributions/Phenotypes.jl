# =============================================================================
# Phenotypes.jl - 表型分解与性状类型
# =============================================================================
#
# 实现表型的遗传与环境分解
#
# 对应章节: Walsh 2nd Ed, Chapter 1 (Overview of Quantitative Genetics)
# =============================================================================

"""
    phenotype_decomposition(G::Real, E::Real; 
                           interaction::Real=0.0, 
                           covariance::Real=0.0) -> NamedTuple

计算表型值的完整分解。

# 模型方程
```math
P = G + E + G \\times E + 2\\text{Cov}(G, E)
```

# 参数
- `G::Real`: 遗传效应值
- `E::Real`: 环境效应值
- `interaction::Real`: 基因型×环境互作效应（默认0）
- `covariance::Real`: 遗传-环境协方差（默认0）

# 返回值
NamedTuple包含：
- `P`: 表型值
- `G`: 遗传效应
- `E`: 环境效应
- `GxE`: 互作效应
- `CovGE`: 协方差项

# 中文说明
表型（Phenotype）是遗传和环境共同作用的结果。

组成成分：
- **G（遗传效应）**: 由基因型决定的部分
- **E（环境效应）**: 由环境因素决定的部分
- **G×E（互作效应）**: 特定基因型对特定环境的反应
- **Cov(G,E)**: 当基因型与环境非随机关联时出现

在随机交配群体和随机化实验中，Cov(G,E)通常为0。
但在自然群体中，由于生态位选择等原因，这一协方差可能非零。

# 公式来源
Walsh 2nd Ed, Eq 1.1

# 示例
```julia
# 基本表型分解
result = phenotype_decomposition(50.0, 10.0)
result.P  # 60.0

# 包含互作的表型
result = phenotype_decomposition(50.0, 10.0, interaction=5.0)
result.P  # 65.0
```

参见: [`variance_partition_phenotypes`](@ref)
"""
function phenotype_decomposition(G::Real, E::Real;
                                interaction::Real=0.0,
                                covariance::Real=0.0)
    P = G + E + interaction + 2 * covariance
    return (P=P, G=G, E=E, GxE=interaction, CovGE=covariance)
end

"""
    variance_partition_phenotypes(V_G::Real, V_E::Real;
                                  V_GxE::Real=0.0,
                                  Cov_GE::Real=0.0) -> NamedTuple

表型方差的分解。

# 模型方程
```math
V_P = V_G + V_E + V_{G \\times E} + 2\\text{Cov}(G, E)
```

# 返回值
NamedTuple包含：
- `V_P`: 表型方差
- `H2`: 广义遗传力 (V_G / V_P)

# 中文说明
表型方差分解是遗传力估计的基础。

广义遗传力 H² = V_G / V_P：
衡量表型变异中遗传变异所占的比例。

# 公式来源
Walsh 2nd Ed, Eq 1.2-1.3

# 示例
```julia
result = variance_partition_phenotypes(100.0, 50.0)
result.H2  # ≈ 0.667
```
"""
function variance_partition_phenotypes(V_G::Real, V_E::Real;
                                       V_GxE::Real=0.0,
                                       Cov_GE::Real=0.0)
    if V_G < 0 || V_E < 0 || V_GxE < 0
        throw(ArgumentError("方差不能为负"))
    end
    
    V_P = V_G + V_E + V_GxE + 2 * Cov_GE
    
    if V_P <= 0
        throw(ArgumentError("总表型方差必须为正"))
    end
    
    H2 = V_G / V_P  # 广义遗传力
    
    return (V_P=V_P, V_G=V_G, V_E=V_E, V_GxE=V_GxE, Cov_GE=Cov_GE, H2=H2)
end

"""
    TraitType

性状类型枚举。

# 类型
- `ContinuousTrait`: 连续性状（如身高、产量）
- `ThresholdTrait`: 阈值性状（如疾病易感性）
- `MeristicTrait`: 计数性状（如脊椎数、花瓣数）
- `BinaryTrait`: 二分性状（如存活/死亡）

# 中文说明
数量性状可分为多种类型，影响分析方法的选择。
"""
@enum TraitType begin
    ContinuousTrait   # 连续性状
    ThresholdTrait    # 阈值性状
    MeristicTrait     # 计数性状
    BinaryTrait       # 二分性状
end

"""
    continuous_vs_discrete_traits(values::AbstractVector;
                                 threshold::Float64=0.0) -> NamedTuple

区分连续性状和离散性状，并提供统计摘要。

# 参数
- `values`: 性状观测值
- `threshold`: 阈值性状的分界点

# 返回值
NamedTuple包含性状类型判断和基本统计量。

# 中文说明
根据数据特征自动判断性状类型：
- 连续性状：值域范围大，取值多样
- 计数性状：离散整数值
- 二分性状：仅有两种取值

# 公式来源
Walsh 2nd Ed, Chapter 1

# 示例
```julia
# 连续性状
result = continuous_vs_discrete_traits([1.2, 3.4, 2.1, 5.6, 4.3])
result.trait_type  # ContinuousTrait

# 计数性状
result = continuous_vs_discrete_traits([1, 2, 3, 2, 1, 3, 4])
result.trait_type  # MeristicTrait

# 二分性状
result = continuous_vs_discrete_traits([0, 1, 1, 0, 1])
result.trait_type  # BinaryTrait
```
"""
function continuous_vs_discrete_traits(values::AbstractVector;
                                       threshold::Float64=0.0)
    # 移除缺失值
    valid_values = filter(v -> !ismissing(v) && !(v isa Number && isnan(v)), values)
    valid_float = Float64.(valid_values)
    
    n = length(valid_float)
    if n == 0
        throw(ArgumentError("没有有效的观测值"))
    end
    
    # 基本统计量
    mu = mean(valid_float)
    sigma = n > 1 ? std(valid_float) : 0.0
    min_val = minimum(valid_float)
    max_val = maximum(valid_float)
    unique_vals = unique(valid_float)
    n_unique = length(unique_vals)
    
    # 判断性状类型
    trait_type = ContinuousTrait
    
    if n_unique == 2
        trait_type = BinaryTrait
    elseif all(v -> isinteger(v), valid_float) && n_unique <= min(20, n/5)
        trait_type = MeristicTrait
    elseif n_unique <= 5 && all(v -> isinteger(v), valid_float)
        # 可能是阈值性状的离散化
        trait_type = ThresholdTrait
    end
    
    # 阈值性状分析
    n_above_threshold = count(v -> v > threshold, valid_float)
    proportion_above = n_above_threshold / n
    
    return (
        trait_type = trait_type,
        n = n,
        mean = mu,
        std = sigma,
        min = min_val,
        max = max_val,
        n_unique = n_unique,
        unique_values = sort(unique_vals),
        n_above_threshold = n_above_threshold,
        proportion_above = proportion_above
    )
end

"""
    phenotypic_correlation(trait1::AbstractVector, trait2::AbstractVector) -> NamedTuple

计算两个性状间的表型相关。

# 返回值
- `r_P`: 表型相关系数
- `se`: 标准误估计
- `n`: 有效样本量

# 中文说明
表型相关是两个性状表型值之间的Pearson相关系数。
它包含遗传相关和环境相关的混合效应。

r_P = Cov(P₁, P₂) / √(V_P₁ × V_P₂)

# 示例
```julia
result = phenotypic_correlation(height, weight)
result.r_P  # 表型相关系数
```
"""
function phenotypic_correlation(trait1::AbstractVector, trait2::AbstractVector)
    if length(trait1) != length(trait2)
        throw(ArgumentError("两个性状的观测数必须相同"))
    end
    
    # 找到两者都有效的观测
    valid_idx = Int[]
    for i in 1:length(trait1)
        v1 = trait1[i]
        v2 = trait2[i]
        if !ismissing(v1) && !ismissing(v2) &&
           !(v1 isa Number && isnan(v1)) && !(v2 isa Number && isnan(v2))
            push!(valid_idx, i)
        end
    end
    
    n = length(valid_idx)
    if n < 3
        throw(ArgumentError("有效样本量不足（需要至少3个）"))
    end
    
    x = Float64[trait1[i] for i in valid_idx]
    y = Float64[trait2[i] for i in valid_idx]
    
    r_P = cor(x, y)
    
    # Fisher z变换计算标准误
    z = atanh(r_P)
    se_z = 1 / sqrt(n - 3)
    
    # 反变换得到r的标准误（近似）
    se_r = (1 - r_P^2) * se_z
    
    return (r_P=r_P, se=se_r, n=n)
end
