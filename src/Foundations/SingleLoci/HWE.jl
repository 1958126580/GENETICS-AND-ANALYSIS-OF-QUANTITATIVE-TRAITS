# =============================================================================
# HWE.jl - Hardy-Weinberg平衡检验
# =============================================================================
# 对应章节: Walsh 2nd Ed, Chapter 4
# =============================================================================

"""
    hwe_chi_square_test(observed::Vector{Int}) -> HWEResult

HWE卡方检验。

# 中文说明
使用Pearson卡方检验评估基因型频率是否符合HWE期望。
χ² = Σ(O-E)²/E, df=1
"""
function hwe_chi_square_test(observed::Vector{Int})
    @assert length(observed) == 3 "需要3个基因型计数"
    
    result = allele_genotype_freqs(observed)
    p, n = result.p, result.n
    expected = hwe_expected_freqs(p, n).exp_counts
    
    # 卡方统计量
    chi_sq = sum((observed .- expected).^2 ./ expected)
    
    # p值 (df=1 for biallelic)
    pvalue = ccdf(Chisq(1), chi_sq)
    
    return HWEResult(observed, expected, chi_sq, pvalue, :chi_square, (p, 1-p))
end

"""
    hwe_exact_test(observed::Vector{Int}) -> HWEResult

Haldane精确检验（小样本推荐）。

# 中文说明
基于完全枚举的精确检验，适用于小样本。
计算在固定等位基因数下，观测到当前或更极端杂合子数的概率。
"""
function hwe_exact_test(observed::Vector{Int})
    @assert length(observed) == 3
    
    n_AA, n_Aa, n_aa = observed
    n = sum(observed)
    
    # 等位基因计数
    n_A = 2n_AA + n_Aa
    n_a = 2n_aa + n_Aa
    
    result = allele_genotype_freqs(observed)
    p = result.p
    expected = hwe_expected_freqs(p, n).exp_counts
    
    # 计算观测配置的概率
    obs_prob = _hwe_genotype_probability(n_AA, n_Aa, n_aa, n_A, n_a)
    
    # 枚举所有可能配置，计算p值
    pvalue = 0.0
    for het in 0:min(n_A, n_a)
        if (n_A - het) % 2 != 0 || (n_a - het) % 2 != 0
            continue
        end
        aa = (n_A - het) ÷ 2
        bb = (n_a - het) ÷ 2
        if aa >= 0 && bb >= 0 && aa + het + bb == n
            prob = _hwe_genotype_probability(aa, het, bb, n_A, n_a)
            if prob <= obs_prob * (1 + 1e-10)
                pvalue += prob
            end
        end
    end
    
    pvalue = min(pvalue, 1.0)
    chi_sq = sum((observed .- expected).^2 ./ max.(expected, 1e-10))
    
    return HWEResult(observed, expected, chi_sq, pvalue, :exact, (p, 1-p))
end

function _hwe_genotype_probability(n_AA, n_Aa, n_aa, n_A, n_a)
    n = n_AA + n_Aa + n_aa
    
    # 使用对数阶乘避免溢出
    log_prob = (
        logfactorial(n) - logfactorial(n_AA) - logfactorial(n_Aa) - logfactorial(n_aa) +
        logfactorial(n_A) + logfactorial(n_a) - logfactorial(2n) +
        n_Aa * log(2)
    )
    
    return exp(log_prob)
end

# 对数阶乘
logfactorial(n::Int) = n <= 1 ? 0.0 : sum(log(i) for i in 2:n)
