# AIcelltype

AIcelltype是一个基于大语言模型的单细胞RNA-seq数据细胞类型自动注释R包，可以快速、准确地识别各种组织的细胞类型。也是我第一次写的包

## 功能特点
- 支持非常多种新模型，不怕掉队，持续进化
- 支持全世界任意组织类型的任意细胞类型识别
- 支持中文和英文两种语言，更适合中国宝宝体质
- 提供分组处理功能，支持大规模细胞集群分析
- 输出结果直接整合进Seurat对象

## 安装方法

```r
install.packages("devtools")
library(devtools)
devtools::install_github("super-818/AIcelltype")

```

# 使用示例

### 1. 设置API密钥

```r
# 在使用前需要设置API密钥
Sys.setenv(API_KEY = "你的API密钥")
```

### 2. 加载示例数据

```r
# 加载AIcelltype包
library('AIcelltype')

# 示例数据为肝脏组织的单细胞RNA-seq数据（sc_KPPC_aLiver_2.rds）
data(sc_KPPC_aLiver_2)
```

### 3. 运行细胞类型注释

```r
library('Seurat')

# 查找差异基因
markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.5,
  test.use = "wilcox",
  verbose = FALSE
)

# 调用ai_celltype函数进行细胞类型注释
cluster2celltype <- ai_celltype(
  input = markers,
  api_key = "sk-yehuhgsggcrskvssfzyttwyblvuavlaxpduzittyonoasgss", #直接公开，star后就能使用
  api_url = "https://api.siliconflow.cn/v1/chat/completions", 
  model = "deepseek-ai/DeepSeek-R1-0528-Qwen3-8B",
  tissuename = "肝脏",
  topgenenumber = 10,
  group_input = FALSE,
  language = "zh"
)
```

### 4. 结果整合与可视化

```r
celltype_labels <- cluster2celltype[as.character(seurat_obj$seurat_clusters)]

names(celltype_labels) <- NULL

seurat_obj$CellType <- celltype_labels

DimPlot(sc_KPPC_aLiver_2, group.by = "CellType", label = TRUE)


```

## 函数参数说明

`ai_celltype()`函数的主要参数：

- `input`：输入数据，可以是Seurat差异基因列表（data.frame，含gene/cluster/avg_log2FC）或基因标记列表（list）
- `tissuename`：组织名称（如"肝脏"或"liver"），用于精准注释
- `api_key`：API密钥（默认从环境变量读取，需提前用Sys.setenv(API_KEY="xxx")设置）
- `api_url`：API地址（固定为"https://api.siliconflow.cn/v1/chat/completions"）
- `model`：模型名称（默认"deepseek-ai/DeepSeek-R1-0528-Qwen3-8B"）
- `topgenenumber`：每个cluster取前N个差异基因（默认10）
- `group_input`：是否分组调用API（单组>30个cluster时自动分组，默认FALSE）
- `language`：输出语言，可选"zh"（中文）或"en"（英文），默认自动识别组织名称语言

# 注意！！！ 解析偶尔有出错，解决办法是：
没有什么是运行ai_celltype这个神奇函数解决不了的，如果有就再运行一次



