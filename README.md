# QuantitativeGenetics.jl

[![Build Status](https://github.com/1958126580/QuantitativeGenetics.jl/workflows/CI/badge.svg)](https://github.com/1958126580/QuantitativeGenetics.jl/actions)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://1958126580.github.io/QuantitativeGenetics.jl/stable)
[![codecov](https://codecov.io/gh/1958126580/QuantitativeGenetics.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/1958126580/QuantitativeGenetics.jl)

## 简介 | Introduction

`QuantitativeGenetics.jl` 是一个基于 Julia 语言的高性能数量遗传学软件包，完整实现了 Walsh, Visscher & Lynch 《Genetics and Analysis of Quantitative Traits: Foundations, 2nd Edition》(2025) 的全部内容。

A high-performance Julia package for quantitative genetics, implementing the complete content of Walsh, Visscher & Lynch's "Genetics and Analysis of Quantitative Traits: Foundations, 2nd Edition" (2025).

## 特性 | Features

### 基础理论 (第1-14章)
- 🧬 **群体遗传学**: Hardy-Weinberg平衡检验、等位基因频率估计
- 📊 **育种值分析**: Fisher分解、加性与显性效应
- 🔗 **连锁不平衡**: D, D', r² 计算与LD衰减建模
- 👨‍👩‍👧‍👦 **亲缘关系**: Henderson递归亲缘矩阵、VanRaden/Yang GRM
- 📐 **混合模型**: Henderson MME、AI-REML方差组分估计
- 🌾 **品系杂交**: Mather-Jinks联合尺度检验、Castle-Wright估计

### QTL作图 (第15-19章)
- 🗺️ **连锁图构建**: Haldane/Kosambi映射函数
- 📍 **区间作图**: Lander-Botstein IM、Haley-Knott回归
- 🎯 **复合区间作图 (CIM)**: Zeng算法、标记协因子选择
- 🏠 **远交系分析**: Haseman-Elston回归、方差组分连锁分析

### 基因组学 (第20-21章)
- 🔬 **GWAS**: 单标记扫描、MLM (EMMAX/GEMMA风格)
- 🧮 **基因组预测**: GBLUP、RR-BLUP、贝叶斯字母表 (A/B/Cπ)
- 📈 **多重校正**: 基因组控制λ因子、FDR

### 数学工具箱 (附录A1-A9)
- Delta方法、通径分析、矩阵微分
- MLE优化、统计功效分析
- MCMC/Gibbs采样、实验设计

## 安装 | Installation

```julia
using Pkg
Pkg.add("QuantitativeGenetics")
```

或从GitHub安装开发版本：

```julia
using Pkg
Pkg.add(url="https://github.com/1958126580/QuantitativeGenetics.jl")
```

## 快速开始 | Quick Start

### 1. Hardy-Weinberg平衡检验

```julia
using QuantitativeGenetics

# 观测基因型计数: AA=25, Aa=50, aa=25
result = hwe_exact_test([25, 50, 25])
println("p-value: $(result.pvalue)")
```

### 2. 亲缘关系矩阵计算

```julia
using QuantitativeGenetics

# 从谱系计算A矩阵
pedigree = Pedigree(
    ids = ["1", "2", "3", "4"],
    sires = [nothing, nothing, "1", "1"],
    dams = [nothing, nothing, "2", "2"]
)
A = recursive_kinship(pedigree)
```

### 3. QTL区间作图

```julia
using QuantitativeGenetics

# 加载数据
pop = load_population("data/f2_cross.csv")

# 执行Haley-Knott区间作图
result = haley_knott_scan(pop, :trait1)

# 获取LOD峰值
peak = find_qtl_peaks(result, threshold=3.0)
```

### 4. GBLUP基因组预测

```julia
using QuantitativeGenetics

# 构建基因组关系矩阵
G = grm_vanraden_yang(genotypes, allele_freqs)

# GBLUP预测
pred = gblup_prediction_solver(phenotypes, G)
println("预测准确性: $(pred.accuracy)")
```

## 数据格式 | Data Format

### 基因型编码
- `0`: 纯合子 AA
- `1`: 杂合子 Aa  
- `2`: 纯合子 aa
- `NaN`/`missing`: 缺失数据

### 输入文件
支持 CSV/TSV 格式，自动识别列名映射。

## 文档 | Documentation

- [完整API文档](https://1958126580.github.io/QuantitativeGenetics.jl/stable)
- [教程与示例](https://1958126580.github.io/QuantitativeGenetics.jl/stable/tutorials)

## 引用 | Citation

如果您在研究中使用了本软件包，请引用：

```bibtex
@software{QuantitativeGenetics_jl,
  author = {MeiBujun},
  title = {QuantitativeGenetics.jl: A Julia Package for Quantitative Genetics},
  year = {2025},
  url = {https://github.com/1958126580/QuantitativeGenetics.jl}
}

@book{Walsh2025,
  author = {Walsh, Bruce and Visscher, Peter M. and Lynch, Michael},
  title = {Genetics and Analysis of Quantitative Traits: Foundations},
  edition = {2nd},
  publisher = {Oxford University Press},
  year = {2025},
  isbn = {978-0192898180}
}
```

## 许可证 | License

MIT License - 详见 [LICENSE](LICENSE) 文件。

## 贡献 | Contributing

欢迎提交 Issues 和 Pull Requests！

## 致谢 | Acknowledgments

本项目的实现基于 Walsh, Visscher & Lynch 的经典教科书，感谢原作者对数量遗传学领域的杰出贡献。
