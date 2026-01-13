# =============================================================================
# MapFunctions.jl - 遗传映射函数
# =============================================================================
# 对应章节: Walsh 2nd Ed, Chapter 17
# =============================================================================

"""
    haldane_mapping(r::Float64) -> Float64

Haldane映射函数：重组率→图距。

# 公式
d = -0.5 × ln(1 - 2r)

# 假设
无交叉干扰（crossing-over interference）
"""
function haldane_mapping(r::Float64)
    @assert 0 <= r <= 0.5 "重组率必须在[0, 0.5]"
    r = min(r, 0.4999)  # 避免log(0)
    return -0.5 * log(1 - 2r) * 100  # 返回cM
end

"""
    inverse_haldane(d::Float64) -> Float64

逆Haldane函数：图距(cM)→重组率。

# 公式
r = 0.5 × (1 - exp(-2d/100))
"""
function inverse_haldane(d::Float64)
    @assert d >= 0 "图距不能为负"
    return 0.5 * (1 - exp(-2 * d / 100))
end

"""
    kosambi_mapping(r::Float64) -> Float64

Kosambi映射函数（考虑正干扰）。

# 公式
d = 0.25 × ln((1+2r)/(1-2r))
"""
function kosambi_mapping(r::Float64)
    @assert 0 <= r <= 0.5
    r = min(r, 0.4999)
    return 0.25 * log((1 + 2r) / (1 - 2r)) * 100
end

"""
    inverse_kosambi(d::Float64) -> Float64

逆Kosambi函数。

r = 0.5 × tanh(2d/100)
"""
function inverse_kosambi(d::Float64)
    @assert d >= 0
    return 0.5 * tanh(2 * d / 100)
end

"""
    information_content_pic(freqs::Vector{Float64}) -> Float64

多态性信息含量PIC。

# 公式
PIC = 1 - Σpᵢ² - Σᵢ<ⱼ 2pᵢ²pⱼ²

# 中文说明
PIC衡量标记的信息量，用于连锁分析。
"""
function information_content_pic(freqs::Vector{Float64})
    @assert isapprox(sum(freqs), 1.0, atol=1e-6)
    
    n = length(freqs)
    sum_p2 = sum(freqs.^2)
    
    # 双杂合概率
    sum_pp = 0.0
    for i in 1:(n-1)
        for j in (i+1):n
            sum_pp += 2 * freqs[i]^2 * freqs[j]^2
        end
    end
    
    return 1 - sum_p2 - sum_pp
end
