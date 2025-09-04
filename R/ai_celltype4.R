#' @param input 输入数据：可为Seurat差异基因列表（data.frame，含gene/cluster/avg_log2FC）或基因标记列表（list）
#' @param tissuename 组织名称（如"大脑皮层"或"cerebral cortex"），用于精准注释
#' @param api_key API密钥（默认从环境变量读取，需提前用Sys.setenv(API_KEY="xxx")设置）
#' @param api_url API地址（固定为"https://api.siliconflow.cn/v1/chat/completions"）
#' @param model 模型名称（默认"deepseek-ai/DeepSeek-R1-0528-Qwen3-8B"）
#' @param topgenenumber 每个cluster取前N个差异基因（默认10）
#' @param group_input 是否分组调用API（单组>30个cluster时自动分组，默认FALSE）
#' @param language 输出语言，可选"zh"（中文）或"en"（英文），默认自动识别组织名称语言
#' @return 解析后的`cluster2celltype`命名向量（可直接用于Seurat细胞类型赋值）
#' @export
ai_celltype <- function(
    input,
    tissuename = NULL,
    api_key = Sys.getenv("API_KEY"),
    api_url = "https://api.siliconflow.cn/v1/chat/completions",
    model = "deepseek-ai/DeepSeek-R1-0528-Qwen3-8B",
    topgenenumber = 10,
    group_input = FALSE,
    language = "zh"
) {
  ###########################################
  #  Step 1：依赖包检查
  ###########################################
  required_packages <- c("jsonlite", "stringr", "httr", "curl")
  missing_pkgs <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(sprintf(
      "缺少必要的R包，请先安装：\ninstall.packages(c('%s'))",
      paste(missing_pkgs, collapse = "', '")
    ))
  }

  ###########################################
  #  Step 2：语言自动识别与验证
  ###########################################
  if (is.null(language)) {
    # 检测中文正则表达式（匹配任何中文字符）
    has_chinese <- stringr::str_detect(tissuename, "[\u4e00-\u9fa5]")
    language <- ifelse(has_chinese, "zh", "en")
    message(sprintf("自动识别语言: %s（基于组织名称）", ifelse(language=="zh", "中文", "英文")))
  } else {
    if (!language %in% c("zh", "en")) {
      stop("语言参数仅支持 'zh'（中文）或 'en'（英文）")
    }
  }

  ###########################################
  #  Step 3：输入数据预处理
  ###########################################
  if (api_key == "") {
    stop("未找到API密钥！请先设置：\nSys.setenv(API_KEY='你的API密钥')")
  }

  if (inherits(input, "list")) {
    input <- sapply(input, paste, collapse = ", ")
  } else if (inherits(input, "data.frame")) {
    required_cols <- c("gene", "cluster", "avg_log2FC")
    if (!all(required_cols %in% colnames(input))) {
      stop(sprintf("输入数据框必须包含以下列: %s", paste(required_cols, collapse = ", ")))
    }
    input <- input[input$avg_log2FC > 0, , drop = FALSE]
    input <- tapply(
      input$gene,
      list(input$cluster),
      function(i) paste0(i[1:min(length(i), topgenenumber)], collapse = ", ")
    )
  } else {
    stop("输入必须是列表（list）或包含gene/cluster/avg_log2FC的数据框（data.frame）")
  }

  ###########################################
  #  Step 4：多语言提示词模板
  ###########################################
  prompt_templates <- list(
    zh = "根据以下每行的marker基因，识别「%s」组织的细胞类型：\n1. 生成Seurat可直接使用的命名向量（格式：cluster2celltype <- c(\"0\"=\"XXX\",\"1\"=\"XXX\",...)）\n2. 混合细胞类型（如\"B细胞/肥大细胞\"）或未明确类型（如\"其他\"）直接保留\n3. 确保cluster编号连续无遗漏（从0开始）\n\n%s",
    en = "Identify cell types of %s tissue based on the following marker genes for each row:\n1. Generate a named vector usable in Seurat (format: cluster2celltype <- c(\"0\"=\"XXX\",\"1\"=\"XXX\",...))\n2. Keep mixed cell types (e.g. \"B cell/mast cell\") or unclear types (e.g. \"other\") as is\n3. Ensure cluster numbers are consecutive and complete (starting from 0)\n\n%s"
  )

  ###########################################
  #  Step 5：API调用与解析函数
  ###########################################
  call_api_and_parse <- function(prompt) {
    # 打印调试信息（根据语言显示）
    if (language == "zh") {
      cat("=== 向模型发送提示词 ===\n", prompt, "\n\n")
    } else {
      cat("=== Sending prompt to model ===\n", prompt, "\n\n")
    }

    body <- list(
      model = model,
      messages = list(list(role = "user", content = prompt)),
      temperature = 0.7,
      max_tokens = 2048,
      response_format = list(type = "text")
    )

    response <- httr::POST(
      url = api_url,
      httr::add_headers(
        `Content-Type` = "application/json",
        `Authorization` = paste0("Bearer ", api_key)
      ),
      body = jsonlite::toJSON(body, auto_unbox = TRUE),
      encode = "raw"
    )

    if (httr::status_code(response) != 200) {
      error_msg <- if (language == "zh") {
        sprintf("API请求失败（状态码: %s）\n响应内容: %s",
                httr::status_code(response),
                httr::content(response, "text", encoding = "UTF-8"))
      } else {
        sprintf("API request failed (status code: %s)\nResponse: %s",
                httr::status_code(response),
                httr::content(response, "text", encoding = "UTF-8"))
      }
      stop(error_msg)
    }

    json_response <- httr::content(response, "text", encoding = "UTF-8")
    if (language == "zh") {
      cat("=== 接收并解析响应 ===\n", json_response, "\n\n")
    } else {
      cat("=== Received and parsing response ===\n", json_response, "\n\n")
    }

    parsed_result <- parse_full_celltype_annotation(json_response, language = language)
    return(parsed_result)
  }

  ###########################################
  #  Step 6：分组处理与结果合并
  ###########################################
  if (group_input) {
    cutnum <- ceiling(length(input) / 30)
    cid <- as.numeric(cut(1:length(input), cutnum))
  } else {
    cutnum <- 1
    cid <- rep(1, length(input))
  }

  progress_msg <- if (language == "zh") {
    sprintf("共分为 %s 组进行处理...", cutnum)
  } else {
    sprintf("Processing in %s groups...", cutnum)
  }
  message(progress_msg)

  all_parsed_results <- sapply(1:cutnum, function(i) {
    if (language == "zh") {
      cat(sprintf("处理第 %s/%s 组...\n", i, cutnum))
    } else {
      cat(sprintf("Processing group %s/%s ...\n", i, cutnum))
    }

    group_idx <- which(cid == i)
    group_clusters <- names(input)[group_idx]
    group_genes <- input[group_idx]
    prompt_lines <- paste0(group_clusters, ": ", group_genes)

    prompt <- sprintf(
      prompt_templates[[language]],
      tissuename,
      paste(prompt_lines, collapse = "\n")
    )

    call_api_and_parse(prompt)
  }, simplify = FALSE)

  ###########################################
  #  Step 7：结果整理与返回
  ###########################################
  cluster2celltype <- do.call(c, all_parsed_results)
  cluster2celltype <- cluster2celltype[order(as.numeric(names(cluster2celltype)))]


  if (language == "zh") {
    message("\n=== 细胞类型注释完成 ===")
    message(sprintf("成功注释 %s 个细胞集群", length(cluster2celltype)))
    message("注意：下游分析前请检查结果，避免AI幻觉导致的错误注释！")
  } else {
    message("\n=== Annotation completed ===")
    message(sprintf("Successfully annotated %s clusters", length(cluster2celltype)))
    message("Warning: Verify results before downstream analysis to avoid AI hallucination!")
  }

  return(cluster2celltype)
}
