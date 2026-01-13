# =============================================================================
# CIM.jl - 复合区间作图
# =============================================================================
# 对应章节: Walsh 2nd Ed, Chapter 18
# =============================================================================

"""
    cim_scan(pop::Population, trait::Symbol;
            n_cofactors::Int=5, window_size::Float64=10.0) -> QTLScanResult

复合区间作图(CIM)扫描。

# 方法
在区间作图中加入标记协因子控制背景遗传效应。
y = μ + QTL效应 + Σβₖxₖ + ε

# 参数
- n_cofactors: 协因子数量（通过逐步回归选择）
- window_size: 窗口大小（cM），窗口内标记不作为协因子
"""
function cim_scan(pop::Population, trait::Symbol;
                  n_cofactors::Int=5, window_size::Float64=10.0)
    trait_name = String(trait)
    y = get_phenotypes(pop, trait_name)
    
    valid_idx = [i for i in 1:length(y) 
                 if !ismissing(y[i]) && !(y[i] isa Number && isnan(y[i]))]
    
    y_valid = Float64[y[i] for i in valid_idx]
    geno_valid = pop.genotypes[valid_idx, :]
    
    n = length(y_valid)
    n_markers = size(geno_valid, 2)
    
    # 选择协因子（前向选择）
    cofactor_idx = select_cofactors_forward(geno_valid, y_valid, n_cofactors)
    
    # 扫描
    positions = Float64[]
    chromosomes = String[]
    lod_scores = Float64[]
    lrt_stats = Float64[]
    effects = Matrix{Float64}(undef, n_markers, 3)
    
    for j in 1:n_markers
        push!(positions, Float64(j))
        push!(chromosomes, "Chr1")
        
        # 排除窗口内的协因子
        active_cof = filter(c -> abs(c - j) > window_size, cofactor_idx)
        
        g = Float64.(geno_valid[:, j])
        x_a = g .- 1.0
        x_d = Float64.(g .== 1)
        
        # 构建设计矩阵
        X_cof = length(active_cof) > 0 ? 
                Float64.(geno_valid[:, active_cof]) .- 1.0 : 
                zeros(n, 0)
        
        X_full = hcat(ones(n), x_a, x_d, X_cof)
        X_null = hcat(ones(n), X_cof)
        
        # 回归
        beta_full = X_full \ y_valid
        y_hat_full = X_full * beta_full
        RSS_full = sum((y_valid .- y_hat_full).^2)
        
        beta_null = X_null \ y_valid
        y_hat_null = X_null * beta_null
        RSS_null = sum((y_valid .- y_hat_null).^2)
        
        if RSS_full > 0 && RSS_full < RSS_null
            LRT = n * log(RSS_null / RSS_full)
            LOD = LRT / (2 * log(10))
        else
            LRT = 0.0
            LOD = 0.0
        end
        
        push!(lrt_stats, LRT)
        push!(lod_scores, LOD)
        effects[j, :] = beta_full[1:3]
    end
    
    return QTLScanResult(positions, chromosomes, lod_scores, lrt_stats, effects,
                        ["mu", "a", "d"]; method=:CIM, trait_name=trait_name)
end

"""
    select_cofactors_forward(geno::Matrix, y::Vector, n_max::Int) -> Vector{Int}

前向逐步选择协因子。
"""
function select_cofactors_forward(geno::AbstractMatrix, y::AbstractVector, n_max::Int)
    n, m = size(geno)
    selected = Int[]
    remaining = Set(1:m)
    
    for _ in 1:n_max
        best_marker = 0
        best_r2 = -Inf
        
        for j in remaining
            test_set = [selected; j]
            X = hcat(ones(n), Float64.(geno[:, test_set]) .- 1.0)
            beta = X \ y
            y_hat = X * beta
            SS_res = sum((y .- y_hat).^2)
            SS_tot = sum((y .- mean(y)).^2)
            r2 = 1 - SS_res / SS_tot
            
            if r2 > best_r2
                best_r2 = r2
                best_marker = j
            end
        end
        
        if best_marker > 0
            push!(selected, best_marker)
            delete!(remaining, best_marker)
        else
            break
        end
    end
    
    return selected
end

"""
    cim_background_control(cofactors::Vector{Int}) -> Vector{Int}

返回CIM协因子列表。
"""
cim_background_control(cofactors::Vector{Int}) = cofactors
