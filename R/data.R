#' 加载AIcelltype包的示例数据
#' 
#' 此函数用于加载包中包含的示例单细胞RNA-seq数据
#' 
#' @param package 包名，默认为"AIcelltype"
#' @param envir 环境，默认为当前环境
#' @param verbose 是否显示加载信息，默认为FALSE
#' @param overwrite 是否覆盖已存在的同名对象，默认为FALSE
#' 
#' @examples
#' # 加载示例数据
#' data(sc_KPPC_aLiver_2)
#' 
#' @export
"sc_KPPC_aLiver_2"

#' @rdname sc_KPPC_aLiver_2
"data_sc_KPPC_aLiver_2"

data_sc_KPPC_aLiver_2 <- function() {
  # 获取包的示例数据目录
  data_dir <- system.file("example", package = "AIcelltype")
  
  # 构建数据文件的完整路径
  data_path <- file.path(data_dir, "sc_KPPC_aLiver_2.rds")
  
  # 检查文件是否存在
  if (!file.exists(data_path)) {
    stop("示例数据文件不存在，请检查包的安装是否完整。")
  }
  
  # 读取RDS文件
  seurat_obj <- readRDS(data_path)
  
  # 返回读取的数据
  return(seurat_obj)
}

#' @rdname sc_KPPC_aLiver_2
#' @export
load_example_data <- function() {
  # 调用数据加载函数
  seurat_obj <- data_sc_KPPC_aLiver_2()
  
  # 将数据赋值给全局环境中的变量
  assign("sc_KPPC_aLiver_2", seurat_obj, envir = .GlobalEnv)
  
  # 显示加载成功的消息
  message("示例数据 'sc_KPPC_aLiver_2' 已成功加载到全局环境中。")
  message("这是一个肝脏组织的单细胞RNA-seq Seurat对象。")
  message("使用 ?sc_KPPC_aLiver_2 查看更多信息。")
  
  # 不可见返回数据
  invisible(seurat_obj)
}

# 设置数据文档
#' 肝脏组织单细胞RNA-seq数据示例
#' 
#' 这是一个来自肝脏组织的单细胞RNA-seq数据集，以Seurat对象的形式存储。
#' 数据来源于KPPC（Kras G12D, p53 R172H, Pdx-1-Cre）小鼠模型的肝脏组织。
#' 
#' @format Seurat对象，包含以下内容：
#' \describe{
#'   \item{assays}{包含RNA表达数据}
#'   \item{meta.data}{包含细胞的元数据信息，如seurat_clusters等}
#'   \item{reductions}{包含降维结果，如UMAP、PCA等}
#' }
#' 
#' @source 内部生成的示例数据
#' 
#' @examples
#' # 加载示例数据
#' data(sc_KPPC_aLiver_2)
#' 
#' # 查看数据基本信息
#' dim(sc_KPPC_aLiver_2)
#' head(sc_KPPC_aLiver_2@meta.data)
#' 
#' @name sc_KPPC_aLiver_2
NULL