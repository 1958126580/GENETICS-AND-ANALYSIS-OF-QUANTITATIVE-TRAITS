# =============================================================================
# Frequencies.jl - 等位基因和基因型频率计算
# =============================================================================
#
# 对应章节: Walsh 2nd Ed, Chapter 4 (Properties of Single Loci)
# =============================================================================

"""
    allele_genotype_freqs(counts::Vector{Int}) -> NamedTuple

从基因型计数估计等位基因和基因型频率。

# 参数
- `counts`: 基因型计数 [n_AA, n_Aa, n_aa]

# 公式
p = (2×n_AA + n_Aa) / (2N)

# 中文说明
从观测的基因型计数估计等位基因频率。
"""
function allele_genotype_freqs(counts::Vector{Int})
    @assert length(counts) == 3 "需要3个基因型计数"
    @assert all(c -> c >= 0, counts) "计数不能为负"
    
    n_AA, n_Aa, n_aa = counts
    N = sum(counts)
    @assert N > 0 "总样本量不能为0"
    
    p = (2n_AA + n_Aa) / (2N)
    q = 1 - p
    geno_freqs = counts ./ N
    
    return (p=p, q=q, geno_freqs=geno_freqs, n=N)
end

"""
    allele_freqs_from_matrix(genotypes; missing_code=-9) -> NamedTuple

从0/1/2编码的基因型向量估计等位基因频率。
"""
function allele_freqs_from_matrix(genotypes::AbstractVector{<:Real};
                                  missing_code::Real=-9)
    valid = [g for g in genotypes if !isnan(g) && g != missing_code]
    n = length(valid)
    n == 0 && return (p=NaN, q=NaN, maf=NaN, n_valid=0)
    
    n_0 = count(g -> round(Int, g) == 0, valid)
    n_1 = count(g -> round(Int, g) == 1, valid)
    n_2 = count(g -> round(Int, g) == 2, valid)
    
    p = (2n_0 + n_1) / (2n)
    q = 1 - p
    maf = min(p, q)
    
    return (p=p, q=q, maf=maf, n_valid=n, n_AA=n_0, n_Aa=n_1, n_aa=n_2)
end

"""
    hwe_expected_freqs(p::Float64, n::Int) -> NamedTuple

HWE期望频率: [p², 2pq, q²]
"""
function hwe_expected_freqs(p::Float64, n::Int)
    @assert 0 <= p <= 1 "频率必须在[0,1]"
    q = 1 - p
    exp_freqs = [p^2, 2p*q, q^2]
    exp_counts = exp_freqs .* n
    return (exp_freqs=exp_freqs, exp_counts=exp_counts, p=p, q=q)
end

minor_allele_frequency(p::Float64) = min(p, 1-p)

"""
    heterozygosity(p::Float64) -> NamedTuple

期望杂合度 H_exp = 2pq
"""
function heterozygosity(p::Float64)
    @assert 0 <= p <= 1
    q = 1 - p
    H_exp = 2p*q
    return (H_exp=H_exp, H_max=0.5, gene_diversity=H_exp)
end

"""
    heterozygosity_multiallelic(freqs) -> Float64

多等位基因杂合度: H = 1 - Σp²ᵢ
"""
function heterozygosity_multiallelic(freqs::AbstractVector{Float64})
    @assert isapprox(sum(freqs), 1.0, atol=1e-6)
    return 1.0 - sum(freqs.^2)
end

effective_allele_number(freqs::AbstractVector{Float64}) = 1.0 / sum(freqs.^2)

"""
    observed_heterozygosity(genotypes; missing_code=-9) -> Float64

观测杂合度 = 杂合子比例
"""
function observed_heterozygosity(genotypes::AbstractVector{<:Real}; missing_code::Real=-9)
    valid = [g for g in genotypes if !isnan(g) && g != missing_code]
    n = length(valid)
    n == 0 && return NaN
    n_het = count(g -> round(Int, g) == 1, valid)
    return n_het / n
end

inbreeding_coefficient_from_het(H_obs::Float64, H_exp::Float64) = 
    H_exp ≈ 0 ? NaN : 1.0 - H_obs/H_exp
