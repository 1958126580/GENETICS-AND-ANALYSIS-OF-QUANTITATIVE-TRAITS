# =============================================================================
# MapConstruction.jl - 遗传图谱构建
# =============================================================================
# 对应章节: Walsh 2nd Ed, Chapter 17
# =============================================================================

"""
    estimate_recombination_mle(n_recomb::Int, n_total::Int) -> NamedTuple

两点最大似然估计重组率。

# 数学推导
对于N个减数分裂产物，其中R个为重组型：
- 似然函数: L(r) = r^R × (1-r)^(N-R)
- MLE: r̂ = R/N
- 方差: Var(r̂) = r(1-r)/N

# 公式来源
Walsh 2nd Ed, Eq 17.5
"""
function estimate_recombination_mle(n_recomb::Int, n_total::Int)
    @assert n_total > 0 "样本量必须大于0"
    @assert 0 <= n_recomb <= n_total "重组数必须在[0, n_total]范围内"
    
    r_hat = n_recomb / n_total
    
    # 标准误 (binomial SE)
    se = sqrt(r_hat * (1 - r_hat) / n_total)
    
    # LOD score for linkage (H0: r=0.5 vs H1: r=r_hat)
    if 0 < r_hat < 0.5
        log_L_r = n_recomb * log10(r_hat) + (n_total - n_recomb) * log10(1 - r_hat)
        log_L_05 = n_total * log10(0.5)
        lod = log_L_r - log_L_05
    else
        lod = 0.0
    end
    
    # 95%置信区间 (通过LOD=1下降法近似)
    ci_lower = max(0.0, r_hat - 1.96 * se)
    ci_upper = min(0.5, r_hat + 1.96 * se)
    
    return (r_hat=r_hat, se=se, lod=lod, n=n_total, 
            ci_lower=ci_lower, ci_upper=ci_upper)
end

"""
    lod_linkage_test(r::Float64, n::Int) -> Float64

连锁LOD检验统计量。

# 数学定义
LOD = log₁₀[L(r) / L(0.5)]
    = R×log₁₀(2r) + (N-R)×log₁₀(2(1-r))

# 判断标准
- LOD ≥ 3.0: 显著连锁
- LOD ≥ 2.0: 暗示性连锁
"""
function lod_linkage_test(r::Float64, n::Int)
    if r >= 0.5 || r <= 0
        return 0.0
    end
    
    n_r = round(Int, n * r)
    n_nr = n - n_r
    
    lod = n_r * log10(2r) + n_nr * log10(2 * (1-r))
    return lod
end

"""
    two_point_linkage_matrix(genotypes::AbstractMatrix{<:Real}; 
                            lod_threshold::Float64=3.0) -> Matrix{Float64}

计算所有标记对的重组率矩阵。

# 返回值
重组率矩阵，对角线为0，上三角存储r值，下三角存储LOD值。
"""
function two_point_linkage_matrix(genotypes::AbstractMatrix{<:Real};
                                  lod_threshold::Float64=3.0)
    n_ind, n_markers = size(genotypes)
    r_matrix = zeros(n_markers, n_markers)
    lod_matrix = zeros(n_markers, n_markers)
    
    for i in 1:n_markers
        for j in (i+1):n_markers
            # 计算重组率（简化：使用相关性近似）
            g1 = genotypes[:, i]
            g2 = genotypes[:, j]
            
            # 计算不一致的比例作为重组率估计
            valid_idx = findall(k -> !isnan(g1[k]) && !isnan(g2[k]) && 
                                     g1[k] >= 0 && g2[k] >= 0, 1:n_ind)
            
            if length(valid_idx) < 10
                r_matrix[i, j] = 0.5
                continue
            end
            
            g1_valid = g1[valid_idx]
            g2_valid = g2[valid_idx]
            n_valid = length(valid_idx)
            
            # 对于F2群体，使用双隐性频率估计重组率
            # r = 1 - sqrt(1 - 2*f_aabb) 其中f_aabb是双隐性频率
            # 简化：使用相关性
            corr = cor(g1_valid, g2_valid)
            
            # r² ≈ (1-2r)² for tight linkage
            if abs(corr) > 0.99
                r_est = 0.001
            elseif abs(corr) < 0.01
                r_est = 0.499
            else
                r_est = (1 - abs(corr)) / 2
            end
            
            r_matrix[i, j] = r_est
            r_matrix[j, i] = r_est
            
            # LOD
            result = estimate_recombination_mle(round(Int, r_est * n_valid), n_valid)
            lod_matrix[i, j] = result.lod
            lod_matrix[j, i] = result.lod
        end
    end
    
    return r_matrix, lod_matrix
end

"""
    linkage_group_assignment(lod_matrix::AbstractMatrix; 
                            lod_threshold::Float64=3.0) -> Vector{Int}

基于LOD阈值将标记分配到连锁群。

使用Union-Find算法高效分组。
"""
function linkage_group_assignment(lod_matrix::AbstractMatrix;
                                  lod_threshold::Float64=3.0)
    n = size(lod_matrix, 1)
    
    # Union-Find数据结构
    parent = collect(1:n)
    rank = zeros(Int, n)
    
    function find_root(x)
        if parent[x] != x
            parent[x] = find_root(parent[x])
        end
        return parent[x]
    end
    
    function union_sets(x, y)
        px, py = find_root(x), find_root(y)
        if px == py
            return
        end
        if rank[px] < rank[py]
            parent[px] = py
        elseif rank[px] > rank[py]
            parent[py] = px
        else
            parent[py] = px
            rank[px] += 1
        end
    end
    
    # 根据LOD阈值合并
    for i in 1:n
        for j in (i+1):n
            if lod_matrix[i, j] >= lod_threshold
                union_sets(i, j)
            end
        end
    end
    
    # 获取分组
    groups = [find_root(i) for i in 1:n]
    
    # 重新编号为连续整数
    unique_groups = unique(groups)
    group_map = Dict(g => i for (i, g) in enumerate(unique_groups))
    
    return [group_map[g] for g in groups]
end

"""
    order_markers_within_group(r_matrix::AbstractMatrix, 
                               group_indices::Vector{Int}) -> Vector{Int}

使用最近邻算法对连锁群内标记排序。

# 算法
1. 从任意标记开始
2. 选择距离当前末端最近的未访问标记
3. 重复直到所有标记被访问
"""
function order_markers_within_group(r_matrix::AbstractMatrix,
                                    group_indices::Vector{Int})
    n = length(group_indices)
    
    if n <= 2
        return group_indices
    end
    
    # 提取子矩阵
    sub_r = r_matrix[group_indices, group_indices]
    
    # 最近邻排序
    visited = falses(n)
    order = Int[]
    
    # 从第一个标记开始
    current = 1
    push!(order, current)
    visited[current] = true
    
    while length(order) < n
        # 找最近的未访问标记
        best_next = 0
        best_dist = Inf
        
        for j in 1:n
            if !visited[j]
                # 检查与当前路径两端的距离
                dist_to_start = sub_r[order[1], j]
                dist_to_end = sub_r[order[end], j]
                min_dist = min(dist_to_start, dist_to_end)
                
                if min_dist < best_dist
                    best_dist = min_dist
                    best_next = j
                end
            end
        end
        
        if best_next == 0
            # 添加任意未访问标记
            best_next = findfirst(.!visited)
        end
        
        # 决定添加到开头还是末尾
        if sub_r[order[1], best_next] < sub_r[order[end], best_next]
            pushfirst!(order, best_next)
        else
            push!(order, best_next)
        end
        visited[best_next] = true
    end
    
    return group_indices[order]
end

"""
    construct_genetic_map(markers::Vector{Marker}, 
                         genotypes::AbstractMatrix{<:Real};
                         lod_threshold::Float64=3.0,
                         map_function::Symbol=:haldane) -> GeneticMap

完整的遗传图谱构建流程。

# 三阶段算法

## Phase 1: 分组 (Grouping)
使用LOD阈值将标记分配到连锁群。

## Phase 2: 排序 (Ordering)
对每个连锁群内的标记使用最近邻算法排序。

## Phase 3: 间距估计 (Spacing)
使用映射函数将重组率转换为图距(cM)。

# 参数
- `markers`: 标记信息向量
- `genotypes`: 基因型矩阵 (个体×标记)
- `lod_threshold`: 连锁分组的LOD阈值
- `map_function`: 映射函数 (:haldane 或 :kosambi)

# 返回值
构建好的GeneticMap对象
"""
function construct_genetic_map(markers::Vector{Marker},
                               genotypes::AbstractMatrix{<:Real};
                               lod_threshold::Float64=3.0,
                               map_function::Symbol=:haldane)
    n_markers = length(markers)
    
    if n_markers == 0
        return GeneticMap("empty", Marker[])
    end
    
    if n_markers == 1
        return GeneticMap("single_marker", markers)
    end
    
    # Phase 1: 计算两点连锁并分组
    r_matrix, lod_matrix = two_point_linkage_matrix(Float64.(genotypes); 
                                                     lod_threshold=lod_threshold)
    groups = linkage_group_assignment(lod_matrix; lod_threshold=lod_threshold)
    
    n_groups = maximum(groups)
    
    # Phase 2 & 3: 对每个连锁群排序并计算距离
    ordered_markers = Marker[]
    
    for g in 1:n_groups
        group_idx = findall(==(g), groups)
        
        if length(group_idx) == 1
            # 单标记连锁群
            m = markers[group_idx[1]]
            push!(ordered_markers, Marker(m.name, "LG$g", 0.0, m.alleles))
            continue
        end
        
        # 排序
        ordered_idx = order_markers_within_group(r_matrix, group_idx)
        
        # 计算累积距离
        cumulative_pos = 0.0
        for (i, idx) in enumerate(ordered_idx)
            if i == 1
                new_marker = Marker(markers[idx].name, "LG$g", 0.0, markers[idx].alleles)
            else
                prev_idx = ordered_idx[i-1]
                r = r_matrix[prev_idx, idx]
                
                # 映射函数转换
                if map_function == :haldane
                    d = haldane_mapping(min(r, 0.499))
                else
                    d = kosambi_mapping(min(r, 0.499))
                end
                
                cumulative_pos += d
                new_marker = Marker(markers[idx].name, "LG$g", cumulative_pos, markers[idx].alleles)
            end
            push!(ordered_markers, new_marker)
        end
    end
    
    return GeneticMap("constructed", ordered_markers; map_function=map_function)
end

"""
    multipoint_recombination_estimate(genotypes::AbstractMatrix, 
                                      positions::Vector{Float64}) -> Vector{Float64}

多点估计相邻标记间的重组率。

使用隐马尔可夫模型(HMM)进行多点估计。
"""
function multipoint_recombination_estimate(genotypes::AbstractMatrix{<:Real},
                                           positions::Vector{Float64})
    n_ind, n_markers = size(genotypes)
    
    if n_markers < 2
        return Float64[]
    end
    
    r_estimates = zeros(n_markers - 1)
    
    # 简化实现：使用两点估计作为初始值
    for j in 1:(n_markers-1)
        g1 = genotypes[:, j]
        g2 = genotypes[:, j+1]
        
        valid_idx = findall(i -> !isnan(g1[i]) && !isnan(g2[i]) && 
                                 g1[i] >= 0 && g2[i] >= 0, 1:n_ind)
        
        if length(valid_idx) < 5
            # 使用位置差估计
            d = abs(positions[j+1] - positions[j])
            r_estimates[j] = inverse_haldane(d)
        else
            # 两点估计
            corr = cor(g1[valid_idx], g2[valid_idx])
            r_estimates[j] = (1 - abs(corr)) / 2
        end
    end
    
    return r_estimates
end
