#' 解析细胞类型注释的JSON响应（内部函数，无需用户手动调用）
#' @param celltype_annotation_json API返回的JSON字符串
#' @param language 输出语言，可选"zh"（中文）或"en"（英文），默认"zh"
#' @return Seurat可用的`cluster2celltype`命名向量
#' @keywords internal  # 标记为内部函数，不对外导出
#'
parse_full_celltype_annotation <- function(
    celltype_annotation_json,
    language = "zh"
) {
  ###########################################
  #  Step 1：输入合法性检查（多语言提示）
  ###########################################
  if (!is.character(celltype_annotation_json) || length(celltype_annotation_json) != 1) {
    error_msg <- if (language == "zh") {
      "输入必须是长度为1的JSON字符串"
    } else {
      "Input must be a JSON string (length = 1)"
    }
    stop(error_msg)
  }

  ###########################################
  #  Step 2：解析JSON，提取content字段（多语言提示）
  ###########################################
  json_data <- tryCatch({
    jsonlite::fromJSON(celltype_annotation_json)
  }, error = function(e) {
    error_msg <- if (language == "zh") {
      sprintf("JSON解析失败: %s", e$message)
    } else {
      sprintf("JSON Parsing Failed: %s", e$message)
    }
    stop(error_msg)
  })

  if (!"choices" %in% names(json_data) || nrow(json_data$choices) == 0) {
    error_msg <- if (language == "zh") {
      "JSON缺少'choices'字段（API响应无效）"
    } else {
      "JSON lacks 'choices' field (invalid API response)"
    }
    stop(error_msg)
  }
  content <- json_data$choices$message$content[1]
  if (is.na(content) || content == "") {
    error_msg <- if (language == "zh") {
      "API响应的'content'字段为空"
    } else {
      "API response 'content' is empty"
    }
    stop(error_msg)
  }

  preview_msg <- if (language == "zh") {
    sprintf("提取的Content（前100字符）: %s", substr(content, 1, 100))
  } else {
    sprintf("Extracted Content (First 100 Chars): %s", substr(content, 1, 100))
  }
  message(preview_msg)

  ###########################################
  #  Step 3：清理content，提取cluster2celltype代码
  ###########################################
  cleaned_content <- stringr::str_remove_all(content, "```r|```|\\*\\*.*?\\*\\*")
  cleaned_content <- stringr::str_remove_all(
    cleaned_content,
    "[^a-zA-Z0-9_\\-\\+\\=\\,\\\"\\/\\(\\)\\n\\.γδ<\\s\u4e00-\u9fa5]"
  )

  start_pattern <- "cluster2celltype\\s*<-\\s*c\\("
  start_pos <- stringr::str_locate(cleaned_content, start_pattern)[1, 1]
  if (is.na(start_pos)) {
    error_msg <- if (language == "zh") {
      "在响应中未找到 'cluster2celltype <- c(' 结构"
    } else {
      "Cannot find 'cluster2celltype <- c(' in response"
    }
    stop(error_msg)
  }

  end_matches <- stringr::str_locate_all(cleaned_content, "\\)")[[1]]
  end_matches <- end_matches[end_matches[, 1] > start_pos, , drop = FALSE]
  if (nrow(end_matches) == 0) {
    error_msg <- if (language == "zh") {
      "未找到与 'cluster2celltype <- c(' 对应的闭合括号 ')'"
    } else {
      "Cannot find closing ')' for 'cluster2celltype <- c('"
    }
    stop(error_msg)
  }
  end_pos <- end_matches[1, 1]

  core_code <- substr(cleaned_content, start_pos, end_pos + 1)
  core_msg <- if (language == "zh") {
    sprintf("提取的核心代码（前150字符）: %s", substr(core_code, 1, 150))
  } else {
    sprintf("Extracted Core Code (First 150 Chars): %s", substr(core_code, 1, 150))
  }
  message(core_msg)

  ###########################################
  #  Step 4：解析核心代码，生成命名向量（主方案+备用方案）
  ###########################################
  tryCatch({
    safe_env <- new.env()
    eval(parse(text = core_code), envir = safe_env)

    if ("cluster2celltype" %in% ls(safe_env)) {
      result <- safe_env$cluster2celltype
      success_msg <- if (language == "zh") {
        sprintf("解析成功（主方案）：共 %s 个cluster", length(result))
      } else {
        sprintf("Parsed %s clusters (Main Method)", length(result))
      }
      message(success_msg)
      return(result)
    } else {
      stop_msg <- if (language == "zh") {
        "核心代码未生成 'cluster2celltype' 对象"
      } else {
        "Core code did not generate 'cluster2celltype'"
      }
      stop(stop_msg)
    }
  }, error = function(e) {
    backup_msg <- if (language == "zh") {
      sprintf("主方案失败：%s，尝试备用方案...", e$message)
    } else {
      sprintf("Main Method Failed: %s. Trying Backup Method...", e$message)
    }
    message(backup_msg)

    pairs <- stringr::str_extract_all(core_code, "\"\\d+\"\\s*=\\s*\"[^\"]+\"")[[1]]
    if (length(pairs) == 0) {
      error_msg <- if (language == "zh") {
        "备用方案失败：未找到有效的键值对"
      } else {
        "Backup Method Failed: No valid key-value pairs found"
      }
      stop(error_msg)
    }

    keys <- stringr::str_remove_all(stringr::str_extract(pairs, "\"\\d+\""), "\"")
    values <- stringr::str_remove_all(stringr::str_extract(pairs, "(?<=\\=\\s*\")([^\"]+)"), "\"")

    result <- values
    names(result) <- keys
    result <- result[order(as.numeric(names(result)))]

    backup_success_msg <- if (language == "zh") {
      sprintf("备用方案成功：共 %s 个cluster", length(result))
    } else {
      sprintf("Parsed %s clusters (Backup Method)", length(result))
    }
    message(backup_success_msg)
    return(result)
  })
}
