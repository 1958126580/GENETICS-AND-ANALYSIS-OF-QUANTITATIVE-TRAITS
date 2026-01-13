# =============================================================================
# ResultTypes.jl - 分析结果类型定义
# =============================================================================
#
# 定义各类分析方法返回的结果结构体
#
# =============================================================================

"""
    AbstractResult

所有分析结果的抽象基类型。
"""
abstract type AbstractResult end

"""
    HWEResult <: AbstractResult

Hardy-Weinberg平衡检验结果。

# 字段
- `observed::Vector{Int}`: 观测基因型计数 [AA, Aa, aa]
- `expected::Vector{Float64}`: HWE期望计数
- `chi_square::Float64`: 卡方统计量
- `pvalue::Float64`: p值
- `test_type::Symbol`: 检验类型（:chi_square 或 :exact）
- `allele_freqs::Tuple{Float64, Float64}`: 等位基因频率 (p, q)

# 中文说明
Hardy-Weinberg平衡检验用于检测群体是否处于随机交配平衡状态。

拒绝HWE的可能原因：
- 非随机交配（选型交配、近交）
- 自然选择
- 突变
- 遗传漂变
- 基因流
- 基因分型错误

# 示例
```julia
result = hwe_exact_test([25, 50, 25])
result.pvalue  # > 0.05 表示不拒绝HWE
```
"""
struct HWEResult <: AbstractResult
    observed::Vector{Int}
    expected::Vector{Float64}
    chi_square::Float64
    pvalue::Float64
    test_type::Symbol
    allele_freqs::Tuple{Float64, Float64}
end

function Base.show(io::IO, r::HWEResult)
    print(io, "HWEResult(χ²=$(round(r.chi_square, digits=3)), p=$(round(r.pvalue, digits=4)), test=:$(r.test_type))")
end

"""
    is_in_hwe(result::HWEResult; alpha::Float64=0.05) -> Bool

判断群体是否处于HWE（p > alpha）。
"""
is_in_hwe(result::HWEResult; alpha::Float64=0.05) = result.pvalue > alpha

"""
    BreedingValueResult <: AbstractResult

育种值分解结果。

# 字段
- `p::Float64`: 等位基因A的频率
- `q::Float64`: 等位基因a的频率
- `a::Float64`: 加性效应（纯合子差异的一半）
- `d::Float64`: 显性偏差
- `alpha::Float64`: 平均效应（替代效应）
- `V_A::Float64`: 加性遗传方差
- `V_D::Float64`: 显性遗传方差
- `V_G::Float64`: 总遗传方差

# 中文说明
根据Fisher的分解，单位点的遗传方差可以分解为：
- 加性方差 V_A = 2pqα²
- 显性方差 V_D = (2pqd)²

其中 α = a + d(q-p) 是平均效应。

# 公式来源
Walsh 2nd Ed, Eq 4.10-4.12
"""
struct BreedingValueResult <: AbstractResult
    p::Float64
    q::Float64
    a::Float64
    d::Float64
    alpha::Float64
    V_A::Float64
    V_D::Float64
    V_G::Float64
end

function Base.show(io::IO, r::BreedingValueResult)
    print(io, "BreedingValue(α=$(round(r.alpha, digits=3)), V_A=$(round(r.V_A, digits=3)), V_D=$(round(r.V_D, digits=3)))")
end

"""
    VarianceComponents <: AbstractResult

方差组分估计结果。

# 字段
- `sigma2_a::Float64`: 加性遗传方差
- `sigma2_e::Float64`: 残差方差
- `sigma2_d::Float64`: 显性遗传方差（可选）
- `heritability::Float64`: 狭义遗传力
- `se::Dict{Symbol, Float64}`: 标准误估计
- `method::Symbol`: 估计方法（:REML, :ML, :ANOVA）
- `converged::Bool`: 是否收敛
- `n_iterations::Int`: 迭代次数

# 中文说明
方差组分分析将表型方差分解为遗传和环境成分。

狭义遗传力 h² = σ²_A / (σ²_A + σ²_E)

# 示例
```julia
vc = reml_ai_algorithm(y, X, Z, K)
vc.heritability  # 遗传力估计
vc.sigma2_a      # 加性方差
```
"""
struct VarianceComponents <: AbstractResult
    sigma2_a::Float64
    sigma2_e::Float64
    sigma2_d::Float64
    heritability::Float64
    se::Dict{Symbol, Float64}
    method::Symbol
    converged::Bool
    n_iterations::Int
    
    function VarianceComponents(sigma2_a::Real, sigma2_e::Real;
                               sigma2_d::Real=0.0,
                               se::Dict{Symbol, Float64}=Dict{Symbol, Float64}(),
                               method::Symbol=:REML,
                               converged::Bool=true,
                               n_iterations::Int=0)
        total = sigma2_a + sigma2_d + sigma2_e
        h2 = total > 0 ? sigma2_a / total : 0.0
        new(Float64(sigma2_a), Float64(sigma2_e), Float64(sigma2_d),
            h2, se, method, converged, n_iterations)
    end
end

function Base.show(io::IO, vc::VarianceComponents)
    print(io, "VarianceComponents(σ²_A=$(round(vc.sigma2_a, digits=3)), σ²_E=$(round(vc.sigma2_e, digits=3)), h²=$(round(vc.heritability, digits=3)))")
end

"""
    QTLScanResult <: AbstractResult

QTL扫描结果。

# 字段
- `positions::Vector{Float64}`: 扫描位置（cM）
- `chromosomes::Vector{String}`: 染色体
- `lod::Vector{Float64}`: LOD分数
- `lrt::Vector{Float64}`: 似然比检验统计量
- `effects::Matrix{Float64}`: 效应估计（列：mu, a, d等）
- `threshold::Float64`: 显著性阈值
- `method::Symbol`: 扫描方法

# 中文说明
QTL扫描结果包含整个基因组的LOD曲线和效应估计。

LOD分数 = log₁₀(似然比)

常用阈值：
- LOD > 2.0: 暗示性
- LOD > 3.0: 显著（单次检验）
- LOD > 4.3: 高度显著（全基因组）

# 示例
```julia
result = haley_knott_scan(pop, :trait)
peaks = find_peaks(result, threshold=3.0)
```
"""
struct QTLScanResult <: AbstractResult
    positions::Vector{Float64}
    chromosomes::Vector{String}
    lod::Vector{Float64}
    lrt::Vector{Float64}
    effects::Matrix{Float64}
    effect_names::Vector{String}
    threshold::Float64
    method::Symbol
    trait_name::String
    
    function QTLScanResult(positions, chromosomes, lod, lrt, effects, effect_names;
                          threshold=LOD_THRESHOLD_SIGNIFICANT,
                          method=:HaleyKnott,
                          trait_name="trait")
        n = length(positions)
        if length(chromosomes) != n || length(lod) != n
            throw(ArgumentError("所有向量长度必须相同"))
        end
        new(positions, chromosomes, lod, lrt, effects, effect_names, 
            threshold, method, trait_name)
    end
end

function Base.show(io::IO, r::QTLScanResult)
    max_lod = maximum(r.lod)
    max_idx = argmax(r.lod)
    print(io, "QTLScanResult($(r.method), max_LOD=$(round(max_lod, digits=2)) at $(r.chromosomes[max_idx]):$(r.positions[max_idx])cM)")
end

"""
    max_lod(result::QTLScanResult) -> NamedTuple

获取最大LOD及其位置信息。
"""
function max_lod(result::QTLScanResult)
    idx = argmax(result.lod)
    return (lod=result.lod[idx], 
            chromosome=result.chromosomes[idx],
            position=result.positions[idx],
            effects=result.effects[idx, :])
end

"""
    significant_peaks(result::QTLScanResult; threshold=nothing) -> Vector{NamedTuple}

提取超过阈值的显著峰。
"""
function significant_peaks(result::QTLScanResult; threshold=nothing)
    thresh = threshold === nothing ? result.threshold : threshold
    peaks = NamedTuple[]
    
    for i in 1:length(result.lod)
        if result.lod[i] >= thresh
            push!(peaks, (
                lod=result.lod[i],
                chromosome=result.chromosomes[i],
                position=result.positions[i],
                effects=result.effects[i, :]
            ))
        end
    end
    
    return peaks
end

"""
    GWASResult <: AbstractResult

全基因组关联分析结果。

# 字段
- `markers::Vector{String}`: 标记ID
- `chromosomes::Vector{String}`: 染色体
- `positions::Vector{Float64}`: 位置
- `pvalues::Vector{Float64}`: p值
- `betas::Vector{Float64}`: 效应估计
- `se::Vector{Float64}`: 标准误
- `maf::Vector{Float64}`: 次等位基因频率
- `lambda_gc::Float64`: 基因组控制膨胀因子
- `method::Symbol`: 分析方法

# 中文说明
GWAS结果包含每个SNP的关联检验统计量。

基因组控制因子 λ_GC：
- λ ≈ 1.0: 无膨胀
- λ > 1.1: 可能存在群体结构

# 显著性阈值
- 全基因组显著：p < 5×10⁻⁸
- 暗示性：p < 1×10⁻⁵
"""
struct GWASResult <: AbstractResult
    markers::Vector{String}
    chromosomes::Vector{String}
    positions::Vector{Float64}
    pvalues::Vector{Float64}
    betas::Vector{Float64}
    se::Vector{Float64}
    maf::Vector{Float64}
    lambda_gc::Float64
    method::Symbol
    
    function GWASResult(markers, chromosomes, positions, pvalues, betas, se, maf;
                       lambda_gc=1.0, method=:GLM)
        new(markers, chromosomes, positions, pvalues, betas, se, maf, lambda_gc, method)
    end
end

function Base.show(io::IO, r::GWASResult)
    n_sig = count(p -> p < GWAS_GENOME_WIDE_ALPHA, r.pvalues)
    min_p = minimum(r.pvalues)
    print(io, "GWASResult($(length(r.markers)) SNPs, $n_sig significant, min_p=$(round(min_p, sigdigits=3)), λ_GC=$(round(r.lambda_gc, digits=3)))")
end

"""
    PredictionResult <: AbstractResult

基因组预测结果。

# 字段
- `gebv::Vector{Float64}`: 基因组估计育种值
- `accuracy::Float64`: 预测准确性
- `reliability::Float64`: 可靠性
- `method::Symbol`: 预测方法
- `marker_effects::Union{Vector{Float64}, Nothing}`: 标记效应（如果适用）

# 中文说明
基因组预测结果包含每个个体的GEBV估计值。

准确性 = cor(GEBV, TBV)
可靠性 = 准确性²
"""
struct PredictionResult <: AbstractResult
    gebv::Vector{Float64}
    individual_ids::Vector{String}
    accuracy::Float64
    reliability::Float64
    method::Symbol
    marker_effects::Union{Vector{Float64}, Nothing}
    
    function PredictionResult(gebv, individual_ids; 
                             accuracy=NaN, reliability=NaN,
                             method=:GBLUP, marker_effects=nothing)
        rel = isnan(reliability) && !isnan(accuracy) ? accuracy^2 : reliability
        new(gebv, individual_ids, accuracy, rel, method, marker_effects)
    end
end

function Base.show(io::IO, r::PredictionResult)
    print(io, "PredictionResult($(r.method), n=$(length(r.gebv)), accuracy=$(round(r.accuracy, digits=3)))")
end

"""
    LineCrossResult <: AbstractResult

品系杂交分析结果。

# 字段
- `m::Float64`: 中亲值
- `a::Float64`: 加性效应
- `d::Float64`: 显性效应
- `chi_square::Float64`: 拟合优度检验
- `pvalue::Float64`: p值
- `n_genes::Float64`: Castle-Wright估计基因数

# 中文说明
Mather-Jinks联合尺度检验结果。
"""
struct LineCrossResult <: AbstractResult
    m::Float64
    a::Float64
    d::Float64
    chi_square::Float64
    pvalue::Float64
    n_genes::Float64
    df::Int
end

function Base.show(io::IO, r::LineCrossResult)
    print(io, "LineCrossResult(m=$(round(r.m, digits=2)), [a]=$(round(r.a, digits=2)), [d]=$(round(r.d, digits=2)), p=$(round(r.pvalue, digits=4)))")
end
