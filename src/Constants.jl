# =============================================================================
# Constants.jl - 全局常量定义
# =============================================================================
#
# 本文件定义了QuantitativeGenetics.jl包中使用的所有全局常量
# 包括数值计算容差、概率边界、遗传学标准值等
#
# =============================================================================

"""
    DEFAULT_TOLERANCE

数值计算的默认收敛容差。用于迭代算法（如EM、REML）的收敛判断。

默认值: 1e-8

# 中文说明
数值计算默认容差，用于判断算法是否收敛。
"""
const DEFAULT_TOLERANCE = 1e-8

"""
    STRICT_TOLERANCE

严格精度要求的容差。用于需要高精度的计算场景。

默认值: 1e-12

# 中文说明
严格容差，用于高精度计算场景。
"""
const STRICT_TOLERANCE = 1e-12

"""
    MAX_ITERATIONS

迭代算法的最大迭代次数。防止算法陷入无限循环。

默认值: 1000

# 中文说明
最大迭代次数，防止算法不收敛时无限循环。
"""
const MAX_ITERATIONS = 1000

"""
    MIN_PROBABILITY

概率计算的最小值边界。防止log(0)导致的数值问题。

默认值: 1e-300

# 中文说明
概率最小值边界，用于防止对数计算时出现-Inf。
"""
const MIN_PROBABILITY = 1e-300

"""
    MAX_PROBABILITY

概率计算的最大值边界。用于数值稳定性。

默认值: 1 - 1e-15

# 中文说明
概率最大值边界，确保概率值不超过1。
"""
const MAX_PROBABILITY = 1.0 - 1e-15

"""
    MIN_RECOMBINATION

最小重组率。重组率的理论下限。

取值: 0.0

# 中文说明
重组率下限，当两个位点完全连锁时为0。
"""
const MIN_RECOMBINATION = 0.0

"""
    MAX_RECOMBINATION

最大重组率。自由重组时的重组率上限。

取值: 0.5

# 中文说明
重组率上限，两个位点独立分离时为0.5。
"""
const MAX_RECOMBINATION = 0.5

"""
    LOD_THRESHOLD_SUGGESTIVE

暗示性LOD阈值。用于初步QTL筛选。

默认值: 2.0

# 中文说明
暗示性显著LOD阈值，用于初步筛选候选QTL。
"""
const LOD_THRESHOLD_SUGGESTIVE = 2.0

"""
    LOD_THRESHOLD_SIGNIFICANT

显著性LOD阈值。用于宣称QTL存在的标准阈值。

默认值: 3.0

# 中文说明
显著性LOD阈值，通常认为LOD>3表示存在QTL。
"""
const LOD_THRESHOLD_SIGNIFICANT = 3.0

"""
    LOD_THRESHOLD_HIGHLY_SIGNIFICANT

高度显著LOD阈值。非常强的QTL证据。

默认值: 4.3

# 中文说明
高度显著LOD阈值，对应全基因组5%显著水平（约1000个独立检验）。
"""
const LOD_THRESHOLD_HIGHLY_SIGNIFICANT = 4.3

"""
    MISSING_GENOTYPE_CODE

缺失基因型编码值。

取值: -9

# 中文说明
缺失基因型的整数编码。
"""
const MISSING_GENOTYPE_CODE = -9

"""
    MISSING_PHENOTYPE_CODE

缺失表型编码值。

取值: NaN

# 中文说明
缺失表型的浮点编码。
"""
const MISSING_PHENOTYPE_CODE = NaN

"""
    DEFAULT_SCAN_STEP

QTL扫描的默认步长（cM）。

默认值: 1.0

# 中文说明
QTL扫描时默认每1 cM计算一次LOD值。
"""
const DEFAULT_SCAN_STEP = 1.0

"""
    DEFAULT_WINDOW_SIZE

CIM中标记窗口大小（cM）。

默认值: 10.0

# 中文说明
复合区间作图中，排除窗口内标记作为协因子的距离。
"""
const DEFAULT_WINDOW_SIZE = 10.0

"""
    DEFAULT_N_PERMUTATIONS

置换检验的默认置换次数。

默认值: 1000

# 中文说明
用于确定经验显著性阈值的置换次数。
"""
const DEFAULT_N_PERMUTATIONS = 1000

"""
    ALPHA_LEVELS

常用的显著性水平集合。

取值: [0.01, 0.05, 0.10]

# 中文说明
常用的α显著性水平。
"""
const ALPHA_LEVELS = [0.01, 0.05, 0.10]

"""
    LOG2

自然对数2，用于LOD与LRT转换。

log(2) ≈ 0.693147

# 中文说明
LOD = LRT / (2 * log(10)) 转换中使用。
"""
const LOG2 = log(2)

"""
    LOG10

自然对数10，用于LOD计算。

log(10) ≈ 2.302585

# 中文说明
LOD分数计算中使用：LOD = log₁₀(LR)
"""
const LOG10 = log(10)

"""
    LRT_TO_LOD

LRT到LOD的转换因子。

取值: 1 / (2 * log(10)) ≈ 0.2171

# 中文说明
似然比检验统计量转换为LOD分数：LOD = LRT * LRT_TO_LOD
"""
const LRT_TO_LOD = 1.0 / (2.0 * LOG10)

"""
    GWAS_GENOME_WIDE_ALPHA

GWAS全基因组显著性水平。

默认值: 5e-8

# 中文说明
GWAS中常用的全基因组显著性阈值（约100万个SNP的Bonferroni校正）。
"""
const GWAS_GENOME_WIDE_ALPHA = 5e-8

"""
    GWAS_SUGGESTIVE_ALPHA

GWAS暗示性显著水平。

默认值: 1e-5

# 中文说明
GWAS暗示性关联阈值。
"""
const GWAS_SUGGESTIVE_ALPHA = 1e-5

"""
    DEFAULT_MAF_THRESHOLD

默认最小等位基因频率阈值。

默认值: 0.01

# 中文说明
过滤低频变异的MAF阈值，低于此值的SNP通常被剔除。
"""
const DEFAULT_MAF_THRESHOLD = 0.01

"""
    DEFAULT_MISSING_RATE_THRESHOLD

默认最大缺失率阈值。

默认值: 0.20

# 中文说明
单个标记或个体缺失率超过20%时发出警告或剔除。
"""
const DEFAULT_MISSING_RATE_THRESHOLD = 0.20

"""
    OUTLIER_SD_THRESHOLD

异常值检测的标准差阈值。

默认值: 4.0

# 中文说明
表型值超过4个标准差被认为是潜在异常值。
"""
const OUTLIER_SD_THRESHOLD = 4.0

"""
    DEFAULT_RANDOM_SEED

默认随机数种子。用于确保结果可重复。

默认值: 42

# 中文说明
随机数生成器默认种子，确保结果可重复。
"""
const DEFAULT_RANDOM_SEED = 42

# =============================================================================
# 遗传编码常量
# =============================================================================

"""
    GENOTYPE_CODES

标准基因型编码字典。

- AA/BB (纯合优势): 0
- Aa/Bb (杂合): 1  
- aa/bb (纯合劣势): 2

# 中文说明
数量遗传学中常用的基因型数值编码。
"""
const GENOTYPE_CODES = Dict(
    :AA => 0, :BB => 0,
    :Aa => 1, :Bb => 1, :AB => 1, :ab => 1,
    :aa => 2, :bb => 2
)

"""
    CROSS_TYPES

支持的杂交类型。

# 中文说明
不同杂交设计的编码。
"""
const CROSS_TYPES = Dict(
    :BC => "Backcross",           # 回交
    :F2 => "F2 Intercross",       # F2代杂交
    :RIL => "Recombinant Inbred", # 重组近交系
    :DH => "Doubled Haploid",     # 双单倍体
    :MAGIC => "MAGIC Population", # 多亲本高级代际杂交
    :NAM => "Nested Association"  # 巢式关联作图
)
