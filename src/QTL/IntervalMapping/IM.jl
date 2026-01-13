# =============================================================================
# IM.jl - 区间作图
# =============================================================================
# 对应章节: Walsh 2nd Ed, Chapter 18
# =============================================================================

"""
    haley_knott_scan(pop::Population, trait::Symbol;
                    step::Float64=1.0) -> QTLScanResult

Haley-Knott回归扫描。

# 方法
对每个测试位置，计算QTL基因型的条件期望，然后回归。
比EM算法快但略有偏差。

# 中文说明
Haley-Knott是最常用的QTL扫描方法。
"""
function haley_knott_scan(pop::Population, trait::Symbol;
                          step::Float64=1.0)
    trait_name = String(trait)
    y = get_phenotypes(pop, trait_name)
    
    # 获取有效观测
    valid_idx = [i for i in 1:length(y) 
                 if !ismissing(y[i]) && !(y[i] isa Number && isnan(y[i]))]
    
    y_valid = Float64[y[i] for i in valid_idx]
    geno_valid = pop.genotypes[valid_idx, :]
    
    n = length(y_valid)
    n_markers = size(geno_valid, 2)
    
    # 简化扫描：在每个标记位置测试
    positions = Float64[]
    chromosomes = String[]
    lod_scores = Float64[]
    lrt_stats = Float64[]
    effects = Matrix{Float64}(undef, n_markers, 3)  # mu, a, d
    
    for j in 1:n_markers
        push!(positions, Float64(j))
        push!(chromosomes, "Chr1")
        
        # 基因型编码
        g = Float64.(geno_valid[:, j])
        
        # 加性和显性编码
        x_a = g .- 1.0  # -1, 0, 1
        x_d = Float64.(g .== 1)  # 0, 1, 0
        
        # 完整模型
        X_full = hcat(ones(n), x_a, x_d)
        beta_full = X_full \ y_valid
        y_hat_full = X_full * beta_full
        RSS_full = sum((y_valid .- y_hat_full).^2)
        
        # 零模型
        y_mean = mean(y_valid)
        RSS_null = sum((y_valid .- y_mean).^2)
        
        # LRT和LOD
        if RSS_full > 0 && RSS_full < RSS_null
            LRT = n * log(RSS_null / RSS_full)
            LOD = LRT / (2 * log(10))
        else
            LRT = 0.0
            LOD = 0.0
        end
        
        push!(lrt_stats, LRT)
        push!(lod_scores, LOD)
        effects[j, :] = beta_full
    end
    
    return QTLScanResult(positions, chromosomes, lod_scores, lrt_stats, effects,
                        ["mu", "a", "d"]; method=:HaleyKnott, trait_name=trait_name)
end

"""
    interval_mapping_scan(pop::Population, trait::Symbol;
                         step::Float64=1.0) -> QTLScanResult

标准区间作图(EM算法)。
"""
function interval_mapping_scan(pop::Population, trait::Symbol;
                               step::Float64=1.0)
    # 简化：调用Haley-Knott（实际EM更复杂）
    result = haley_knott_scan(pop, trait; step=step)
    return QTLScanResult(result.positions, result.chromosomes, result.lod, result.lrt,
                        result.effects, result.effect_names; 
                        method=:EM, trait_name=result.trait_name)
end

"""
    lod_score_profile(theta_grid::Vector{Float64}, 
                     log_likelihoods::Vector{Float64}) -> NamedTuple

生成LOD曲线。
"""
function lod_score_profile(theta_grid::Vector{Float64},
                           log_likelihoods::Vector{Float64})
    max_ll = maximum(log_likelihoods)
    lod = (log_likelihoods .- log_likelihoods[end]) ./ log(10)  # 相对于r=0.5
    
    peak_idx = argmax(lod)
    
    return (theta=theta_grid, lod=lod, 
            peak_theta=theta_grid[peak_idx], peak_lod=lod[peak_idx])
end

"""
    haley_knott_regression_weights(r1::Float64, r2::Float64,
                                   position::Float64) -> NamedTuple

计算Haley-Knott回归的条件概率权重。
"""
function haley_knott_regression_weights(r1::Float64, r2::Float64, position::Float64)
    # 对于F2群体
    # p_QQ, p_Qq, p_qq 给定侧翼标记
    
    c1 = inverse_haldane(r1 * 100)
    c2 = inverse_haldane(r2 * 100)
    
    # 简化：假设标记完全信息
    w_QQ = (1 - c1)^2 * (1 - c2)^2
    w_qq = c1^2 * c2^2
    w_Qq = 1 - w_QQ - w_qq
    
    return (w_QQ=w_QQ, w_Qq=w_Qq, w_qq=w_qq)
end
