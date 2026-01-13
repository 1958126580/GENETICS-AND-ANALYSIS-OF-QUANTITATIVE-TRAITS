# =============================================================================
# PathAnalysis.jl - 通径分析
# =============================================================================
#
# 实现Wright通径分析方法
#
# 对应章节: Walsh 2nd Ed, Chapter 3 (Path Analysis)
# =============================================================================

"""
    path_coeff_matrix(R::AbstractMatrix, structure::Dict{Tuple{Int,Int}, Symbol}) -> Matrix{Float64}

从相关矩阵和结构关系计算通径系数。

# 参数
- `R`: 观测变量间的相关矩阵
- `structure`: 因果结构字典，键为(因, 果)索引对，值为边类型

# 中文说明
通径分析用于在因果图中分解相关为直接和间接效应。

# 公式来源
Wright (1921), Walsh 2nd Ed, Chapter 3

# 示例
```julia
# 相关矩阵
R = [1.0 0.5 0.3;
     0.5 1.0 0.4;
     0.3 0.4 1.0]

# 假设变量1→3, 2→3
path_coeffs = path_coeff_matrix(R, Dict((1,3)=>:direct, (2,3)=>:direct))
```
"""
function path_coeff_matrix(R::AbstractMatrix, 
                          structure::Dict{Tuple{Int,Int}, Symbol}=Dict{Tuple{Int,Int}, Symbol}())
    n = size(R, 1)
    
    # 如果没有指定结构，假设所有变量（除最后一个）都指向最后一个变量
    if isempty(structure)
        for i in 1:(n-1)
            structure[(i, n)] = :direct
        end
    end
    
    # 识别因变量（被其他变量指向的变量）
    dependent_vars = unique([k[2] for k in keys(structure)])
    
    result = zeros(n, n)
    
    for dep_var in dependent_vars
        # 找到所有指向该变量的自变量
        predictors = [k[1] for k in keys(structure) if k[2] == dep_var]
        
        if isempty(predictors)
            continue
        end
        
        # 提取子相关矩阵（预测变量间）
        R_xx = R[predictors, predictors]
        
        # 预测变量与因变量的相关
        r_xy = R[predictors, dep_var]
        
        # 求解通径系数：p = R_xx^(-1) * r_xy
        path_coeffs = R_xx \ r_xy
        
        # 存储结果
        for (i, pred) in enumerate(predictors)
            result[pred, dep_var] = path_coeffs[i]
        end
    end
    
    return result
end

"""
    path_tracing_rules(P::AbstractMatrix, from::Int, to::Int) -> NamedTuple

使用通径追踪规则分解相关系数。

# 通径规则（Wright）
1. 可以沿箭头正向（因→果）或反向追踪
2. 同一条路径中只能改变方向一次（在"顶点"处）
3. 每条边只能通过一次

# 返回值
- `direct`: 直接通径系数
- `indirect`: 间接效应（通过共同原因）
- `total`: 总相关

# 中文说明
通径追踪规则用于分解观测相关为因果链上的各个成分。

# 公式来源
Walsh 2nd Ed, Eq 3.15

# 示例
```julia
P = path_coeff_matrix(R, structure)
decomp = path_tracing_rules(P, 1, 3)
decomp.direct   # 直接效应
decomp.indirect # 间接效应
```
"""
function path_tracing_rules(P::AbstractMatrix, from::Int, to::Int)
    n = size(P, 1)
    
    if from == to
        return (direct=1.0, indirect=0.0, total=1.0)
    end
    
    # 直接效应
    direct = P[from, to]
    
    # 计算间接效应（简化版本：仅考虑一层中介）
    indirect = 0.0
    for k in 1:n
        if k != from && k != to
            # 检查是否存在from→k→to的路径
            if P[from, k] != 0 && P[k, to] != 0
                indirect += P[from, k] * P[k, to]
            end
            # 也检查from和k相关（共同原因），然后k→to
            if P[k, from] != 0 && P[k, to] != 0
                indirect += P[k, from] * P[k, to]
            end
        end
    end
    
    total = direct + indirect
    
    return (direct=direct, indirect=indirect, total=total)
end

"""
    wright_path_decomposition(correlations::AbstractMatrix, 
                              causal_order::Vector{Int}) -> NamedTuple

Wright通径分解完整算法。

# 参数
- `correlations`: 变量间的相关矩阵
- `causal_order`: 变量的因果顺序（从原因到结果）

# 中文说明
根据假设的因果顺序，将相关系数分解为直接和间接通径系数。
假设变量按因果顺序排列，每个变量只受其前面变量的影响。

# 示例
```julia
R = [1.0 0.3 0.5;
     0.3 1.0 0.4;
     0.5 0.4 1.0]

result = wright_path_decomposition(R, [1, 2, 3])
# 变量3受变量1和2的影响
```
"""
function wright_path_decomposition(correlations::AbstractMatrix,
                                   causal_order::Vector{Int})
    n = size(correlations, 1)
    
    if length(causal_order) != n
        throw(ArgumentError("因果顺序必须包含所有变量"))
    end
    
    # 按因果顺序重排相关矩阵
    R = correlations[causal_order, causal_order]
    
    # 初始化通径系数矩阵
    P = zeros(n, n)
    
    # 对每个变量（从第二个开始），计算从前面变量的通径系数
    for j in 2:n
        # 所有可能的原因
        causes = 1:(j-1)
        
        # 子相关矩阵
        R_xx = R[causes, causes]
        r_xy = R[causes, j]
        
        # 通径系数
        p = R_xx \ r_xy
        
        for (i, cause_idx) in enumerate(causes)
            P[cause_idx, j] = p[i]
        end
    end
    
    # 计算残差通径（未解释的方差）
    residual_paths = zeros(n)
    for j in 1:n
        explained_var = 0.0
        for i in 1:(j-1)
            explained_var += P[i, j]^2
            for k in 1:(j-1)
                if k != i
                    explained_var += P[i, j] * P[k, j] * R[i, k]
                end
            end
        end
        residual_paths[j] = sqrt(max(0.0, 1.0 - explained_var))
    end
    
    return (
        path_coefficients = P,
        residual_paths = residual_paths,
        variable_order = causal_order
    )
end

"""
    standardized_coefficients(X::AbstractMatrix, y::AbstractVector) -> Vector{Float64}

计算标准化回归系数（通径系数）。

# 数学定义
```math
\\beta^* = \\beta \\frac{s_x}{s_y}
```

# 中文说明
标准化回归系数表示自变量每变化一个标准差，因变量变化的标准差数。
这些系数与通径系数等价（当没有多重共线性时）。

# 示例
```julia
std_beta = standardized_coefficients(X, y)
```
"""
function standardized_coefficients(X::AbstractMatrix, y::AbstractVector)
    n = size(X, 1)
    p = size(X, 2)
    
    # 标准化X和y
    X_sd = [std(X[:, j]) for j in 1:p]
    y_sd = std(y)
    
    X_std = (X .- mean(X, dims=1)) ./ X_sd'
    y_std = (y .- mean(y)) ./ y_sd
    
    # 在标准化数据上回归（无截距）
    beta_std = X_std \ y_std
    
    return beta_std
end

"""
    variance_inflation_factor(X::AbstractMatrix) -> Vector{Float64}

计算方差膨胀因子（VIF）。

# 数学定义
```math
VIF_j = \\frac{1}{1 - R_j^2}
```

其中 R²_j 是用其他变量预测变量j的R²。

# 中文说明
VIF用于检测多重共线性：
- VIF < 5: 无严重共线性
- VIF > 10: 严重共线性

# 示例
```julia
vif = variance_inflation_factor(X)
if any(vif .> 10)
    @warn "存在严重多重共线性"
end
```
"""
function variance_inflation_factor(X::AbstractMatrix)
    n, p = size(X)
    
    if p < 2
        return [1.0]
    end
    
    vif = zeros(p)
    
    for j in 1:p
        # 用其他变量预测第j个变量
        y_j = X[:, j]
        X_other = X[:, setdiff(1:p, j)]
        
        # 需要足够的样本
        if n <= size(X_other, 2) + 1
            vif[j] = Inf
            continue
        end
        
        # R²
        result = bols_with_intercept(X_other, y_j)
        R2_j = result.R2
        
        if R2_j >= 1.0
            vif[j] = Inf
        else
            vif[j] = 1.0 / (1.0 - R2_j)
        end
    end
    
    return vif
end

"""
    phenotypic_path_model(phenotypes::Dict{Symbol, Vector{Float64}},
                         structure::Vector{Tuple{Symbol, Symbol}}) -> NamedTuple

表型通径模型分析。

# 参数
- `phenotypes`: 表型数据字典
- `structure`: 因果关系列表，如 [(:height, :weight), (:age, :weight)]

# 中文说明
基于观测表型数据和假设的因果结构进行通径分析。

# 示例
```julia
data = Dict(:height => [...], :weight => [...], :bmi => [...])
structure = [(:height, :weight), (:height, :bmi), (:weight, :bmi)]
result = phenotypic_path_model(data, structure)
```
"""
function phenotypic_path_model(phenotypes::Dict{Symbol, Vector{Float64}},
                               structure::Vector{Tuple{Symbol, Symbol}})
    # 获取变量列表
    all_vars = unique(vcat([first(s) for s in structure], [last(s) for s in structure]))
    n_vars = length(all_vars)
    var_to_idx = Dict(v => i for (i, v) in enumerate(all_vars))
    
    # 构建相关矩阵
    R = zeros(n_vars, n_vars)
    for i in 1:n_vars
        R[i, i] = 1.0
        for j in (i+1):n_vars
            c = cor(phenotypes[all_vars[i]], phenotypes[all_vars[j]])
            R[i, j] = c
            R[j, i] = c
        end
    end
    
    # 构建结构字典
    struct_dict = Dict{Tuple{Int,Int}, Symbol}()
    for (from, to) in structure
        struct_dict[(var_to_idx[from], var_to_idx[to])] = :direct
    end
    
    # 计算通径系数
    P = path_coeff_matrix(R, struct_dict)
    
    return (
        path_coefficients = P,
        correlation_matrix = R,
        variables = all_vars,
        structure = structure
    )
end
