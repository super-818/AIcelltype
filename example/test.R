# 加载必要的包
library('Seurat')
library('presto')

# 加载当前包
library('AIcelltype')

# 加载示例数据
data_dir <- system.file("example", package = "AIcelltype")
seurat_obj <- readRDS(file.path(data_dir, "sc_KPPC_aLiver_2.rds"))

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
# 注意：使用前请设置API密钥：Sys.setenv(API_KEY="你的API密钥")
cluster2celltype <- ai_celltype(
  input = markers,
  api_key = "sk-yehuhgsggcrskvssfzyttwyblvuavlaxpduzittyonoasgss",
  api_url = "https://api.siliconflow.cn/v1/chat/completions",  # 固定API地址
  model = "deepseek-ai/DeepSeek-R1-0528-Qwen3-8B",  # 模型
  tissuename = "肝脏",  # 组织名称
  topgenenumber = 10,
  group_input = FALSE,
  language = "zh"
)

# sc_KPPC_aLiver_2


celltype_labels <- cluster2celltype[as.character(seurat_obj$seurat_clusters)]

names(celltype_labels) <- NULL

seurat_obj$CellType <- celltype_labels

DimPlot(sc_KPPC_aLiver_2, group.by = "CellType", label = TRUE)
