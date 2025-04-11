
##################################################
library(tidyverse)
library(cowplot)

# 读取数据
all_data1 <- list.files(pattern = "*.merged_100windows.txt") %>% 
  map_df(~read_tsv(.x) %>% 
           group_by(species) %>% 
           mutate(window_id = row_number()))

# 为每个物种创建图形
plot_list <- map(unique(all_data1$species), ~{
  current_species <- .x
  sp_data <- filter(all_data1, species == current_species)
  
  # 计算比例因子
  max_gene <- max(sp_data$gene_count)
  max_te <- max(sp_data$te_count)
  ratio <- max_gene / max_te
  
  ggplot(sp_data, aes(x = window_id)) +
    # 趋势线
    geom_line(aes(y = gene_count, color = "Gene Count"), size = 0.6) +
    geom_line(aes(y = te_count * ratio, color = "TE Count"), size = 0.6) +
    
    # 坐标轴设置
    scale_x_continuous(
      name = "Window Sequence",
      breaks = seq(0, max(sp_data$window_id), by = 20),  # 每20窗口显示一个刻度
      expand = c(0.02, 0.02)  # 减少边缘空白
    ) +
    scale_y_continuous(
      name = "Gene Count",
      sec.axis = sec_axis(~./ratio, name = "TE Count"),
      breaks = scales::pretty_breaks(n = 4),  # 自动生成4个主要刻度
      expand = c(0.05, 0.05)
    ) +
    
    # 颜色设置
    scale_color_manual(
      values = c("Gene Count" = "#1b9e77", "TE Count" = "#d95f02")
    ) +
    
    # 主题设置
    theme_minimal(base_size = 9) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
      panel.border = element_rect(color = "grey70", fill = NA, size = 0.3),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(margin = margin(t = 5)),  # X轴上边距
      axis.title.y.left = element_text(color = "#1b9e77", margin = margin(r = 6)),
      axis.title.y.right = element_text(color = "#d95f02", margin = margin(l = 6)),
      axis.ticks = element_line(color = "grey60", size = 0.3),
      axis.text = element_text(color = "grey40", size = 7),
      legend.position = "none",
      plot.margin = unit(c(5,5,5,5), "mm")
    ) +
    labs(title = current_species)
})

# 组合图形
final_plot <- plot_grid(
  plotlist = plot_list,
  ncol = 2,
  align = "hv",  # 确保坐标轴对齐
  axis = "tblr",
  labels = "AUTO",
  label_size = 10,
  label_x = 0.1,
  label_y = 0.95
) 
final_plot
# 保存高清图片
ggsave("gene_te_final_plot.pdf", final_plot, width = 28, height = 24, units = "cm")









#################################
library(tidyverse)
library(cowplot)

# --------------------------
# 数据读取与预处理模块
# --------------------------

# 验证文件存在性
data_files <- list.files(pattern = "*merged_100windows.txt$")
if(length(data_files) == 0) stop("未找到输入文件，请确认：
1. 工作目录是否正确（当前目录：", getwd(), "）
2. 文件名是否包含'merged_100windows.txt'后缀")

# 安全读取数据
all_data <- map_df(data_files, ~{
  tryCatch({
    df <- read_tsv(.x, show_col_types = FALSE) %>%
      group_by(species) %>%
      mutate(window_id = row_number())  # 添加顺序窗口编号
    
    # 数据完整性验证
    required_cols <- c("chr", "start", "end", "gene_count", "te_count", "species")
    if(!all(required_cols %in% names(df))) {
      stop("文件", .x, "缺少必要列：", 
           paste(setdiff(required_cols, names(df)), collapse = ", "))
    }
    df
  }, error = function(e) {
    message("处理文件", .x, "时出错：", e$message)
    NULL
  })
})

# 验证数据内容
if(nrow(all_data) == 0) stop("所有文件读取失败，请检查文件格式")
print(str(all_data))  # 打印数据结构

# --------------------------
# 图形生成模块
# --------------------------

# 创建输出目录
output_dir <- "Species_Plots"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 定义标准化绘图函数
create_species_plot <- function(sp) {
  sp_data <- all_data %>%
    filter(species == sp) %>%
    arrange(window_id)  # 确保窗口顺序
  
  # 计算比例因子
  gene_max <- max(sp_data$gene_count, na.rm = TRUE)
  te_max <- max(sp_data$te_count, na.rm = TRUE)
  ratio <- ifelse(te_max == 0, 1, gene_max / te_max)  # 防止除零错误
  
  # 核心绘图逻辑
  ggplot(sp_data, aes(x = window_id)) +
    geom_line(aes(y = gene_count, color = "Gene"), linewidth = 0.8) +
    geom_line(aes(y = te_count * ratio, color = "TE"), linewidth = 0.8) +
    scale_color_manual(
      name = "",
      values = c(Gene = "#d95f02", TE = "#1b9e77"),
      labels = c(
        Gene = paste("Gene Count (max", gene_max, ")"), 
        TE = paste("TE Count (max", te_max, ")")
      )
    ) +
    scale_x_continuous(
      name = "Window Sequence",
      breaks = scales::breaks_pretty(6),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(
      name = "Gene Count",
      sec.axis = sec_axis(~./ratio, name = "TE Count"),
      breaks = scales::breaks_pretty(5)
    ) +
    labs(
      title = paste("Species:", sp),
      subtitle = paste("Total Windows:", nrow(sp_data)),
      caption = paste("Data source:", paste(data_files, collapse = "; "))
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
      axis.title.y.left = element_text(color = "#d95f02", margin = margin(r = 10)),
      axis.title.y.right = element_text(color = "#1b9e77", margin = margin(l = 10)),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
      panel.border = element_rect(fill = NA, color = "grey70"),
      legend.position = "bottom",
      plot.margin = margin(15, 15, 15, 15)
    )
}

# --------------------------
# 批量输出模块
# --------------------------

# 获取有效物种列表
valid_species <- unique(all_data$species)
if(length(valid_species) == 0) stop("数据中未找到有效物种信息")

# 生成并保存图形
walk(valid_species, ~{
  sp <- .x
  message("正在处理物种：", sp)
  
  # 生成图形
  plot <- tryCatch({
    create_species_plot(sp)
  }, error = function(e) {
    message(sp, "图形生成失败：", e$message)
    NULL
  })
  
  # 安全保存输出
  if(!is.null(plot)) {
    filename <- file.path(output_dir, paste0(gsub("[^[:alnum:]]", "_", sp), ".pdf"))
    
    ggsave(
      filename = filename,
      plot = plot,
      device = cairo_pdf,
      width = 20, 
      height = 15,
      units = "cm",
      dpi = 300
    )
    message("成功保存：", filename)
  }
})

message("处理完成！输出目录：", normalizePath(output_dir))
