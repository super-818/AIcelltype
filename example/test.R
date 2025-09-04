library('Seurat')
library('presto')

markers <- FindAllMarkers(
  sc_KPPC_aLiver_2,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.5,
  test.use = "wilcox",
  verbose = FALSE
)


# 调用ds_celltype函数
cluster2celltype <- ai_celltype(
  input = markers,
  api_key = "sk-yehuhgsggcrskvssfzyttwyblvuavlaxpduzittyonoasgss",
  api_url = "https://api.siliconflow.cn/v1/chat/completions",  # 固定API地址
  model = "deepseek-ai/DeepSeek-R1-0528-Qwen3-8B",# 模型
  tissuename = "大脑皮层",
  topgenenumber = 10,
  group_input = FALSE,
  language = "zh"
)

# sc_KPPC_aLiver_2


celltype_labels <- cluster2celltype[as.character(sc_KPPC_aLiver_2$seurat_clusters)]
names(celltype_labels) <- NULL

sc_SB1_1$CellType <- celltype_labels
DimPlot(sc_KPPC_aLiver_2, group.by = "CellType", label = TRUE)


