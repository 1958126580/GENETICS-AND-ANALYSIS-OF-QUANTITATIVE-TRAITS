# =============================================================================
# Validation.jl - 数据验证函数
# =============================================================================
#
# 提供基因型和表型数据的验证功能
#
# =============================================================================

"""
    ValidationResult

数据验证结果。

# 字段
- `valid::Bool`: 是否通过验证
- `warnings::Vector{String}`: 警告信息
- `errors::Vector{String}`: 错误信息
- `statistics::Dict{Symbol, Any}`: 统计信息
"""
struct ValidationResult
    valid::Bool
    warnings::Vector{String}
    errors::Vector{String}
    statistics::Dict{Symbol, Any}
end

function Base.show(io::IO, r::ValidationResult)
    status = r.valid ? "✓ Valid" : "✗ Invalid"
    n_warn = length(r.warnings)
    n_err = length(r.errors)
    print(io, "ValidationResult($status, $n_warn warnings, $n_err errors)")
end

"""
    validate_genotypes(genotypes::AbstractMatrix{<:Real};
                      expected_codes::Vector{Int}=[0, 1, 2],
                      missing_code::Int=-9,
                      max_missing_rate::Float64=0.2) -> ValidationResult

验证基因型矩阵。

# 检查内容
1. 是否包含非法编码
2. 个体和标记的缺失率
3. 单态标记（MAF = 0）

# 参数
- `genotypes`: 基因型矩阵（行=个体，列=标记）
- `expected_codes`: 有效基因型编码
- `missing_code`: 缺失值编码
- `max_missing_rate`: 允许的最大缺失率

# 中文说明
全面验证基因型数据质量：
- 检测非法编码值
- 统计缺失率
- 识别低变异标记

# 示例
```julia
result = validate_genotypes(geno_matrix)
if !result.valid
    println.(result.errors)
end
```
"""
function validate_genotypes(genotypes::AbstractMatrix{<:Real};
                           expected_codes::Vector{Int}=[0, 1, 2],
                           missing_code::Int=-9,
                           max_missing_rate::Float64=0.2)
    warnings = String[]
    errors = String[]
    stats = Dict{Symbol, Any}()
    
    n_ind, n_markers = size(genotypes)
    stats[:n_individuals] = n_ind
    stats[:n_markers] = n_markers
    
    # 检查非法编码
    all_valid_codes = Set([expected_codes..., missing_code])
    illegal_count = 0
    for i in 1:n_ind
        for j in 1:n_markers
            val = genotypes[i, j]
            if isnan(val)
                continue  # NaN作为缺失值处理
            end
            int_val = round(Int, val)
            if !(int_val in all_valid_codes)
                illegal_count += 1
            end
        end
    end
    
    if illegal_count > 0
        push!(errors, "发现 $illegal_count 个非法基因型编码")
    end
    stats[:illegal_codes] = illegal_count
    
    # 计算个体缺失率
    ind_missing_rates = zeros(n_ind)
    for i in 1:n_ind
        n_missing = count(j -> ismissing_geno(genotypes[i, j], missing_code), 1:n_markers)
        ind_missing_rates[i] = n_missing / n_markers
    end
    
    high_missing_inds = count(r -> r > max_missing_rate, ind_missing_rates)
    if high_missing_inds > 0
        push!(warnings, "$high_missing_inds 个个体的缺失率超过 $(max_missing_rate*100)%")
    end
    stats[:high_missing_individuals] = high_missing_inds
    stats[:mean_ind_missing_rate] = mean(ind_missing_rates)
    
    # 计算标记缺失率
    marker_missing_rates = zeros(n_markers)
    for j in 1:n_markers
        n_missing = count(i -> ismissing_geno(genotypes[i, j], missing_code), 1:n_ind)
        marker_missing_rates[j] = n_missing / n_ind
    end
    
    high_missing_markers = count(r -> r > max_missing_rate, marker_missing_rates)
    if high_missing_markers > 0
        push!(warnings, "$high_missing_markers 个标记的缺失率超过 $(max_missing_rate*100)%")
    end
    stats[:high_missing_markers] = high_missing_markers
    stats[:mean_marker_missing_rate] = mean(marker_missing_rates)
    
    # 检查单态标记
    monomorphic = 0
    for j in 1:n_markers
        col = [genotypes[i, j] for i in 1:n_ind if !ismissing_geno(genotypes[i, j], missing_code)]
        if length(unique(round.(Int, col))) <= 1
            monomorphic += 1
        end
    end
    
    if monomorphic > 0
        push!(warnings, "$monomorphic 个单态标记（无变异）")
    end
    stats[:monomorphic_markers] = monomorphic
    
    valid = isempty(errors)
    return ValidationResult(valid, warnings, errors, stats)
end

# 辅助函数
function ismissing_geno(val, missing_code)
    return isnan(val) || val == missing_code
end

"""
    validate_phenotypes(phenotypes::AbstractVector;
                       max_missing_rate::Float64=0.2,
                       outlier_threshold::Float64=4.0) -> ValidationResult

验证表型数据。

# 检查内容
1. 缺失率
2. 异常值（超过outlier_threshold个标准差）
3. 变异性（是否所有值相同）

# 中文说明
全面验证表型数据质量，识别缺失和异常值。

# 示例
```julia
result = validate_phenotypes(trait_values)
println("发现 $(result.statistics[:n_outliers]) 个异常值")
```
"""
function validate_phenotypes(phenotypes::AbstractVector;
                            max_missing_rate::Float64=0.2,
                            outlier_threshold::Float64=4.0)
    warnings = String[]
    errors = String[]
    stats = Dict{Symbol, Any}()
    
    n = length(phenotypes)
    stats[:n_observations] = n
    
    # 计算缺失率
    n_missing = count(x -> ismissing(x) || (x isa Number && isnan(x)), phenotypes)
    missing_rate = n_missing / n
    stats[:n_missing] = n_missing
    stats[:missing_rate] = missing_rate
    
    if missing_rate > max_missing_rate
        push!(errors, "缺失率 $(round(missing_rate*100, digits=1))% 超过阈值 $(max_missing_rate*100)%")
    end
    
    # 提取有效观测值
    valid_values = Float64[]
    for v in phenotypes
        if !ismissing(v) && !(v isa Number && isnan(v))
            push!(valid_values, Float64(v))
        end
    end
    
    stats[:n_valid] = length(valid_values)
    
    if length(valid_values) < 2
        push!(errors, "有效观测值不足2个，无法进行分析")
        return ValidationResult(false, warnings, errors, stats)
    end
    
    # 基本统计量
    mu = mean(valid_values)
    sigma = std(valid_values)
    stats[:mean] = mu
    stats[:std] = sigma
    stats[:min] = minimum(valid_values)
    stats[:max] = maximum(valid_values)
    
    # 检查变异性
    if sigma < 1e-10
        push!(errors, "表型无变异（标准差≈0）")
    end
    
    # 检查异常值
    if sigma > 0
        z_scores = abs.((valid_values .- mu) ./ sigma)
        n_outliers = count(z -> z > outlier_threshold, z_scores)
        stats[:n_outliers] = n_outliers
        
        if n_outliers > 0
            push!(warnings, "发现 $n_outliers 个潜在异常值（|z| > $outlier_threshold）")
        end
    else
        stats[:n_outliers] = 0
    end
    
    valid = isempty(errors)
    return ValidationResult(valid, warnings, errors, stats)
end

"""
    validate_population(pop::Population) -> ValidationResult

验证群体数据的完整性和一致性。

# 检查内容
1. 基因型矩阵维度
2. 表型数据维度
3. 基因型和表型数据质量

# 中文说明
综合验证群体数据的各个方面。
"""
function validate_population(pop::Population)
    warnings = String[]
    errors = String[]
    stats = Dict{Symbol, Any}()
    
    stats[:n_individuals] = n_individuals(pop)
    stats[:n_markers] = n_markers(pop)
    stats[:n_traits] = n_traits(pop)
    stats[:design] = pop.design
    
    # 验证基因型
    geno_result = validate_genotypes(pop.genotypes)
    if !geno_result.valid
        append!(errors, geno_result.errors)
    end
    append!(warnings, geno_result.warnings)
    stats[:genotype_stats] = geno_result.statistics
    
    # 验证每个表型
    for (trait_name, values) in pop.traits
        pheno_result = validate_phenotypes(values)
        if !pheno_result.valid
            push!(errors, "性状 '$trait_name': " * join(pheno_result.errors, "; "))
        end
        for w in pheno_result.warnings
            push!(warnings, "性状 '$trait_name': $w")
        end
    end
    
    valid = isempty(errors)
    return ValidationResult(valid, warnings, errors, stats)
end

"""
    validate_pedigree(ped::Pedigree) -> ValidationResult

验证谱系结构。

# 检查内容
1. 循环依赖（后代出现在亲本之前）
2. 自交（个体是自己的亲本）
3. 重复ID

# 中文说明
检查谱系结构的合法性。
"""
function validate_pedigree(ped::Pedigree)
    warnings = String[]
    errors = String[]
    stats = Dict{Symbol, Any}()
    
    n = length(ped.ids)
    stats[:n_individuals] = n
    stats[:n_founders] = n_founders(ped)
    stats[:n_generations] = ped.n_generations
    
    # 检查重复ID
    id_set = Set{String}()
    for id in ped.ids
        if id in id_set
            push!(errors, "重复的个体ID: $id")
        end
        push!(id_set, id)
    end
    
    # 检查自交
    for i in 1:n
        id = ped.ids[i]
        if ped.sires[i] == id || ped.dams[i] == id
            push!(errors, "个体 $id 是自己的亲本")
        end
    end
    
    # 检查顺序（简化检查：亲本应该出现在后代之前）
    id_to_idx = Dict(id => i for (i, id) in enumerate(ped.ids))
    order_issues = 0
    for i in 1:n
        if ped.sires[i] !== nothing && haskey(id_to_idx, ped.sires[i])
            if id_to_idx[ped.sires[i]] > i
                order_issues += 1
            end
        end
        if ped.dams[i] !== nothing && haskey(id_to_idx, ped.dams[i])
            if id_to_idx[ped.dams[i]] > i
                order_issues += 1
            end
        end
    end
    
    if order_issues > 0
        push!(warnings, "发现 $order_issues 处亲本出现在后代之后（谱系顺序问题）")
    end
    stats[:order_issues] = order_issues
    
    valid = isempty(errors)
    return ValidationResult(valid, warnings, errors, stats)
end

"""
    check_marker_allele_frequency(genotypes::AbstractVector{<:Real};
                                  maf_threshold::Float64=0.01) -> NamedTuple

检查标记的等位基因频率。

# 返回值
返回包含MAF、基因型计数等信息的NamedTuple。

# 中文说明
计算单个标记的等位基因频率和基因型分布。
"""
function check_marker_allele_frequency(genotypes::AbstractVector{<:Real};
                                       maf_threshold::Float64=0.01)
    valid = [g for g in genotypes if !isnan(g) && g != MISSING_GENOTYPE_CODE]
    n = length(valid)
    
    if n == 0
        return (maf=NaN, p=NaN, q=NaN, n_valid=0, low_maf=true,
                n_AA=0, n_Aa=0, n_aa=0)
    end
    
    # 计数各基因型
    n_AA = count(g -> round(Int, g) == 0, valid)
    n_Aa = count(g -> round(Int, g) == 1, valid)
    n_aa = count(g -> round(Int, g) == 2, valid)
    
    # 计算等位基因频率
    p = (2*n_AA + n_Aa) / (2*n)  # A的频率
    q = 1 - p                    # a的频率
    
    maf = min(p, q)
    low_maf = maf < maf_threshold
    
    return (maf=maf, p=p, q=q, n_valid=n, low_maf=low_maf,
            n_AA=n_AA, n_Aa=n_Aa, n_aa=n_aa)
end
