# =============================================================================
# IO.jl - 数据输入输出函数
# =============================================================================
#
# 提供数据加载和结果保存功能
#
# =============================================================================

"""
    load_population(filepath::String; 
                   format::Symbol=:auto,
                   design::Symbol=:Unknown) -> Population

从文件加载群体数据。

# 参数
- `filepath`: 数据文件路径
- `format`: 文件格式 (:csv, :tsv, :auto)
- `design`: 杂交设计类型

# 支持的文件格式
- CSV/TSV: 第一列为个体ID，后续列为标记基因型
- 表型数据可作为单独文件或合并在基因型文件中

# 列名自动映射
自动识别常见列名变体：
- "id", "ID", "ind", "individual" → 个体ID
- "chr", "chrom", "chromosome" → 染色体
- "pos", "position", "cM", "bp" → 位置

# 中文说明
从CSV/TSV文件加载群体基因型和表型数据。
支持灵活的列名识别和多种数据格式。

# 示例
```julia
pop = load_population("data/f2_cross.csv", design=:F2)
```
"""
function load_population(filepath::String;
                        format::Symbol=:auto,
                        design::Symbol=:Unknown,
                        id_col::Union{String, Int}=1,
                        marker_cols::Union{UnitRange, Vector, Symbol}=:auto,
                        trait_cols::Union{Vector{String}, Nothing}=nothing)
    
    # 确定文件格式
    if format == :auto
        if endswith(lowercase(filepath), ".tsv") || endswith(lowercase(filepath), ".txt")
            format = :tsv
        else
            format = :csv
        end
    end
    
    # 读取数据
    delim = format == :tsv ? '\t' : ','
    df = CSV.read(filepath, DataFrame; delim=delim)
    
    # 提取个体ID
    if id_col isa Int
        individual_ids = string.(df[!, id_col])
        id_col_name = names(df)[id_col]
    else
        individual_ids = string.(df[!, id_col])
        id_col_name = id_col
    end
    
    n_ind = length(individual_ids)
    
    # 识别标记列和表型列
    all_cols = names(df)
    
    # 排除ID列和已知的非标记列
    non_marker_patterns = [r"^id$"i, r"^ind"i, r"^sample"i, r"^name"i,
                          r"^sex$"i, r"^gender$"i, r"^sire$"i, r"^dam$"i]
    
    marker_names = String[]
    trait_names = String[]
    
    for col in all_cols
        if col == id_col_name
            continue
        end
        
        # 检查是否为非标记列
        is_non_marker = any(p -> occursin(p, col), non_marker_patterns)
        
        if trait_cols !== nothing && col in trait_cols
            push!(trait_names, col)
        elseif !is_non_marker
            # 尝试判断是基因型还是表型
            col_data = df[!, col]
            if is_likely_genotype(col_data)
                push!(marker_names, col)
            else
                push!(trait_names, col)
            end
        end
    end
    
    # 构建基因型矩阵
    n_markers = length(marker_names)
    genotypes = Matrix{Float64}(undef, n_ind, n_markers)
    
    for (j, marker) in enumerate(marker_names)
        col_data = df[!, marker]
        for i in 1:n_ind
            val = col_data[i]
            if ismissing(val)
                genotypes[i, j] = NaN
            elseif val isa Number
                genotypes[i, j] = Float64(val)
            else
                # 尝试解析字符串编码
                genotypes[i, j] = parse_genotype_string(string(val))
            end
        end
    end
    
    # 构建表型字典
    traits = Dict{String, Vector{Union{Float64, Missing}}}()
    for trait in trait_names
        col_data = df[!, trait]
        trait_values = Vector{Union{Float64, Missing}}(undef, n_ind)
        for i in 1:n_ind
            val = col_data[i]
            if ismissing(val) || (val isa Number && isnan(val))
                trait_values[i] = missing
            elseif val isa Number
                trait_values[i] = Float64(val)
            else
                # 尝试解析为数字
                try
                    trait_values[i] = parse(Float64, string(val))
                catch
                    trait_values[i] = missing
                end
            end
        end
        traits[trait] = trait_values
    end
    
    return Population(
        basename(filepath),
        individual_ids,
        marker_names,
        genotypes,
        traits,
        design
    )
end

# 辅助函数：判断列是否像基因型数据
function is_likely_genotype(col_data::AbstractVector)
    valid_count = 0
    geno_like_count = 0
    
    for val in col_data
        if ismissing(val)
            continue
        end
        valid_count += 1
        
        if val isa Integer || (val isa AbstractFloat && isinteger(val))
            v = round(Int, val)
            if v in [-9, 0, 1, 2]
                geno_like_count += 1
            end
        elseif val isa AbstractString
            if occursin(r"^[012]$|^-?9$|^[AaTtGgCc]{1,2}$|^[AaBb]{1,2}$", val)
                geno_like_count += 1
            end
        end
    end
    
    return valid_count > 0 && geno_like_count / valid_count > 0.8
end

# 辅助函数：解析字符串形式的基因型
function parse_genotype_string(s::AbstractString)
    s = strip(s)
    
    # 数字编码
    if occursin(r"^-?[0-9]+$", s)
        return parse(Float64, s)
    end
    
    # AA/Aa/aa 编码
    if length(s) == 2
        if uppercase(s[1]) == uppercase(s[2])
            if isuppercase(s[1])
                return 0.0  # AA
            else
                return 2.0  # aa
            end
        else
            return 1.0  # Aa
        end
    end
    
    # 缺失值
    if s in ["NA", "na", "N/A", ".", "-", ""]
        return NaN
    end
    
    return NaN
end

"""
    load_pedigree(filepath::String; format::Symbol=:auto) -> Pedigree

从文件加载谱系数据。

# 文件格式
CSV/TSV，必须包含三列：
- 个体ID
- 父本ID（未知用0或NA）
- 母本ID（未知用0或NA）

# 中文说明
从文件加载谱系结构数据。
"""
function load_pedigree(filepath::String; format::Symbol=:auto)
    if format == :auto
        if endswith(lowercase(filepath), ".tsv") || endswith(lowercase(filepath), ".txt")
            format = :tsv
        else
            format = :csv
        end
    end
    
    delim = format == :tsv ? '\t' : ','
    df = CSV.read(filepath, DataFrame; delim=delim)
    
    # 识别列
    cols = names(df)
    n_cols = length(cols)
    
    if n_cols < 3
        throw(ArgumentError("谱系文件至少需要3列：ID, Sire, Dam"))
    end
    
    # 提取数据
    ids = string.(df[!, 1])
    sires = Vector{Union{String, Nothing}}(undef, length(ids))
    dams = Vector{Union{String, Nothing}}(undef, length(ids))
    
    for i in 1:length(ids)
        # 处理父本
        sire_val = df[i, 2]
        if ismissing(sire_val) || sire_val == "0" || sire_val == 0 || 
           string(sire_val) in ["NA", "na", ".", "-", ""]
            sires[i] = nothing
        else
            sires[i] = string(sire_val)
        end
        
        # 处理母本
        dam_val = df[i, 3]
        if ismissing(dam_val) || dam_val == "0" || dam_val == 0 ||
           string(dam_val) in ["NA", "na", ".", "-", ""]
            dams[i] = nothing
        else
            dams[i] = string(dam_val)
        end
    end
    
    return Pedigree(ids, sires, dams)
end

"""
    load_genetic_map(filepath::String; format::Symbol=:auto) -> GeneticMap

从文件加载遗传图谱。

# 文件格式
CSV/TSV，至少包含：
- 标记名称
- 染色体
- 位置（cM）

# 中文说明
从文件加载遗传图谱信息。
"""
function load_genetic_map(filepath::String; 
                         format::Symbol=:auto,
                         map_function::Symbol=:haldane)
    if format == :auto
        if endswith(lowercase(filepath), ".tsv") || endswith(lowercase(filepath), ".txt")
            format = :tsv
        else
            format = :csv
        end
    end
    
    delim = format == :tsv ? '\t' : ','
    df = CSV.read(filepath, DataFrame; delim=delim)
    
    # 识别列（支持常见列名变体）
    cols = lowercase.(names(df))
    
    # 查找标记名称列
    name_idx = findfirst(c -> c in ["marker", "name", "id", "snp", "locus"], cols)
    if name_idx === nothing
        name_idx = 1
    end
    
    # 查找染色体列
    chr_idx = findfirst(c -> c in ["chr", "chrom", "chromosome", "contig"], cols)
    if chr_idx === nothing
        throw(ArgumentError("未找到染色体列"))
    end
    
    # 查找位置列
    pos_idx = findfirst(c -> c in ["pos", "position", "cm", "dist", "distance"], cols)
    if pos_idx === nothing
        throw(ArgumentError("未找到位置列"))
    end
    
    # 构建标记列表
    markers = Marker[]
    for i in 1:nrow(df)
        name = string(df[i, name_idx])
        chr = string(df[i, chr_idx])
        pos = Float64(df[i, pos_idx])
        
        push!(markers, Marker(name, chr, pos))
    end
    
    return GeneticMap(basename(filepath), markers; map_function=map_function)
end

"""
    save_results(result::AbstractResult, filepath::String; format::Symbol=:csv)

保存分析结果到文件。

# 中文说明
将分析结果导出为CSV/TSV文件。
"""
function save_results(result::QTLScanResult, filepath::String; format::Symbol=:csv)
    delim = format == :tsv ? '\t' : ','
    
    df = DataFrame(
        chromosome = result.chromosomes,
        position = result.positions,
        lod = result.lod,
        lrt = result.lrt
    )
    
    # 添加效应列
    for (i, name) in enumerate(result.effect_names)
        df[!, Symbol(name)] = result.effects[:, i]
    end
    
    CSV.write(filepath, df; delim=delim)
    return filepath
end

function save_results(result::GWASResult, filepath::String; format::Symbol=:csv)
    delim = format == :tsv ? '\t' : ','
    
    df = DataFrame(
        marker = result.markers,
        chromosome = result.chromosomes,
        position = result.positions,
        pvalue = result.pvalues,
        beta = result.betas,
        se = result.se,
        maf = result.maf
    )
    
    CSV.write(filepath, df; delim=delim)
    return filepath
end

function save_results(result::PredictionResult, filepath::String; format::Symbol=:csv)
    delim = format == :tsv ? '\t' : ','
    
    df = DataFrame(
        individual = result.individual_ids,
        gebv = result.gebv
    )
    
    CSV.write(filepath, df; delim=delim)
    return filepath
end

"""
    export_relationship_matrix(rm::RelationshipMatrix, filepath::String)

导出关系矩阵到文件。

# 中文说明
将亲缘关系矩阵导出为CSV格式。
"""
function export_relationship_matrix(rm::RelationshipMatrix, filepath::String)
    df = DataFrame(rm.matrix, Symbol.(rm.ids))
    insertcols!(df, 1, :ID => rm.ids)
    CSV.write(filepath, df)
    return filepath
end
