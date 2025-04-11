# 加载必要的包
library(tidyverse)
library(rtracklayer)
library(dplyr)

# 创建一个空数据框用于存储各个物种的基因组大小
genome_sizes = data.frame(Species = character(), Genome_Size = numeric())

# 定义物种列表
species_list1 <- c("Aedes_aegypti", "Drosophila_melanogaster", "Daphnia_pulex", "Tigriopus_copepod", 
                   "Lepeophtheirus_salmonis", "Hyalella_azteca", "Penaeus_chinensis", "Penaeus_monodon", 
                   "Penaeus_vannamei", "Portunus_trituberculatus", "Procambarus_clarkii", "Eriocheir_sinensis", 
                   "Homarus_americanus", "Argiope_bruennichi")


# 循环计算每个物种的基因组大小
for (species in species_list1) {
  # 从对应的 fai 文件中获取基因组大小
  genome_size = sum(read.table(file = paste0("12.ref/", species, ".genome.fa.fai"))$V2)
  
  # 将当前物种和基因组大小添加到数据框中
  genome_sizes = rbind(genome_sizes, data.frame(Species = species, Genome_Size = genome_size))
}


#替换下划线为空格
genome_sizes <- genome_sizes %>%
  mutate(Species = gsub("_", " ", Species))
# 创建一个空数据框用于存储所有物种的 TE 注释数据
allTE = data.frame()

#定义23中的路径物种名称

species_list23 <- c("Aed", "Dro", "Dap", "Tig", 
                    "Lep", "Hya", "Pen.c", "Pen.m", 
                    "Pen", "Por", "Pro", "Eri", 
                    "Hom", "Arg")


# 循环处理每个物种
for (species in species_list23) {
  # 构建文件路径
  te_file_path = paste0("23.EDTA/", species, "/", species, ".genome.fa.mod.EDTA.TEanno.gff3")
  
  # 检查文件是否存在
  if (file.exists(te_file_path)) {
    # 读取物种注释文件
    te_data = import(te_file_path) %>%
      as_tibble() %>%
      mutate(Species = species)
    
    # 将当前物种的 TE 注释数据添加到总体数据框中
    allTE = rbind(allTE, te_data)
  } else {
    message(paste("Warning: File not found for species", species))
  }
}


# 处理合并后的数据并添加基因组大小
allTE_processed = allTE %>%
  mutate(Classification0 = str_split(Classification, "/", simplify = TRUE)[, 1]) %>%
  left_join(genome_sizes, by = "Species")

# 从 allTE_processed 中删除 phase、TSD 和 TIR 等列
allTE_processed_filtered = allTE_processed %>%
  select(-phase, -TSD, -TIR, -motif, -tsd, -Parent)



#绘制图形
library(tidyverse)
library(ggplot2)
library(dplyr)
# 计算每个物种的转座子总数量
te_size <- allTE_summ %>%
  group_by(Species) %>%
  summarise(TE_size = sum(size))


# 合并基因组大小和转座子总数量
genome_te_data <- left_join(te_size, genome_sizes, by = "Species")

genome_te_data_sorted_size <-genome_te_data[order(genome_te_data$Genome_Size), ]


# 绘制基因组大小和转座子数量的关联图，添加物种名称到X轴并显示置信区间
library(ggrepel)
p1 <- ggplot(genome_te_data_sorted_size, aes(x = log2(Genome_Size), y = log2(TE_size))) +
  geom_point(color = "black", size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "#8074AC", size = 1, fill = "#0055A5") +
  geom_text_repel(aes(label = Species), size = 3, box.padding = 0.35, point.padding = 0.5,
                  segment.color = 'grey50') +
  theme_classic() +
  labs(x = "genome size (log2(bp))", y = "TE content (log2(bp))") +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))


p1
#计算相关系数
# 先计算两者的对数
log_genome_size <- log2(genome_te_data_sorted_size$Genome_Size)
log_te_size <- log2(genome_te_data_sorted_size$TE_size)

# 计算相关性系数
correlation_coefficient <- cor(log_genome_size, log_te_size, method = "pearson")

# 打印结果
print(correlation_coefficient)


#绘制superfamily分类


#数据统计
superfamily_data <- allTE_summ_Class %>%
  mutate(Classification = ifelse(str_detect(Classification, "/"), str_split(Classification, "/", simplify = TRUE)[, 2], Classification)) %>% # 保留/后面的内容，如果没有/则保留全部 
  group_by(Species, Classification, Classification0) %>%
  summarise(size = sum(size, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(Species) %>%
  mutate(proportion = size / sum(size)) %>%
  ungroup() 

superfamily_data$Species <- factor(superfamily_data$Species, levels = custom_species_order)

library(RColorBrewer)
#定义调色板数量
colourCount <-  length(unique(superfamily_data$Classification))
#绘图
p3 <- ggplot(superfamily_data, aes(x = Species, y = proportion, fill = as.factor(Classification), colour = Classification0)) +
  geom_bar(stat = "identity", position = "stack", size = 0.7, width = 0.9) +
  theme_bw() +
  labs(title = "The proportion of transposable elements(sumperfamily) in sixteen species", x = NULL ,y = NULL) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(colourCount)) +  # 使用自定义颜色为 Classification 提供填充颜色
  scale_colour_manual(values = class0_colors, name = "Classification") +  # 使用自定义颜色为 Classification0 提供边框颜色, 并设置图例标题
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
        legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"),  # 调整图例键的大小
        legend.text = element_text(size = 8)) +
  guides(fill = guide_legend(order = 1, title = NULL),  # 移除 Classification 图例标题并设置顺序为1
         colour = guide_legend(order = 2, title = "Classification"))  # 调整 Classification0 图例标题并设置顺序为2
p3



#不同家族转座子数目跟基因组大小的关系


#添加基因组xinx

# 去除 Species 列中为 "Tigriopus copepod" 的行
genome_sizes <- genome_sizes %>%
  filter(Species != "Capitulum mitella")




allTE_summ_Class <- allTE_summ_Class %>%
  filter(Species != "Tigriopus californicus")

allTE_summ_Class <- allTE_summ_Class %>%
  left_join(genome_sizes, by = "Species")

# 提取需要的列
te_count_data <- allTE_summ_Class %>%
  select(Species, Classification, count, Genome_Size) %>%
  mutate(log_count = log10(count + 1), log_genome_size = log10(Genome_Size))

correlation_results_count <- te_count_data %>%
  filter(!is.na(log_count) & !is.na(log_genome_size)) %>% # 移除缺失值
  group_by(Classification) %>%
  summarize(
    correlation = if(n() >= 3) cor(log_count, log_genome_size, method = "pearson") else NA,
    p_value = if(n() >= 3) cor.test(log_count, log_genome_size, method = "pearson")$p.value else NA
  ) %>%
  ungroup() %>%
  filter(Classification != "TIR/Sola2") %>%
  filter(Classification != "pararetrovirus") 


######size
te_size_data <- allTE_summ_Class %>%
  select(Species, Classification, size, Genome_Size) %>%
  mutate(log_size = log10(size + 1), log_genome_size = log10(Genome_Size))

correlation_results_size <- te_size_data %>%
  filter(!is.na(log_size) & !is.na(log_genome_size)) %>% # 移除缺失值
  group_by(Classification) %>%
  summarize(
    correlation = if(n() >= 3) cor(log_size, log_genome_size, method = "pearson") else NA,
    p_value = if(n() >= 3) cor.test(log_size, log_genome_size, method = "pearson")$p.value else NA
  ) %>%
  ungroup()

# 绘制气泡图
library(ggplot2)
library(ggrepel)
correlation_results_count <- correlation_results_count[correlation_results_count$correlation > 0.5, ]
# 绘制气泡图并添加标签
ggplot(correlation_results_count, aes(x = reorder(Classification, correlation), y = correlation, size = -log10(p_value), color = correlation)) +
  geom_text_repel(aes(label = sprintf("r = %.2f\np = %.2e", correlation, p_value)),size = 2.5, box.padding = 0.35, point.padding = 0.5 ,
                  fontface = "bold") +
  geom_point(alpha = 0.5) +
  scale_size_continuous(range = c(3, 15)) +
  scale_color_gradient2(midpoint = 0.5, low = "green", mid = "blue", high = "red") +
  scale_y_continuous(limits = c(0.5, 1)) +
  labs(x = NULL,
       size = "-log10(p-value)",
       color = "Correlation") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1 ,size = 6), 
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.key.size = unit(0.5, "cm"), # 调整图例键的大小
        legend.text = element_text(size = 6), # 调整图例文本的大小
        panel.grid.major = element_line(color = "gray", linetype = "dashed"),
        panel.grid.minor = element_line(color = "gray", linetype = "dashed")) 







#绘制转座子数量和基因大小的相关性

TE_count <- te_count_data %>%
  group_by(Species) %>%
  summarise(total_count = sum(count)) %>%
  left_join(genome_sizes, by = "Species")

ggplot(TE_count, aes(x = log10(Genome_Size), y = log10(total_count))) +
  geom_point(color = "black", size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "#8074AC", size = 1, fill = "#0055A5") +
  geom_text_repel(aes(label = Species), size = 3, box.padding = 0.2, point.padding = 0.1,
                  segment.color = 'grey50') +
  theme_classic() +
  labs(x = "genome size (log10(bp))", y = "TE content (log10(count))") +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
p4


# 数据准备
log_genome_size1 <- log10(TE_count$Genome_Size)
log_te_size1 <- log10(TE_count$total_count)

correlation_test1 <- cor.test(log_genome_size1, log_te_size1, method = "pearson")

# 输出相关性系数和P值
correlation_coefficient1 <- correlation_test1$estimate
p_value1<- correlation_test1$p.value
print(correlation_coefficient1) # 0.9052379
print(p_value1) #  2.083162e-05

########size
log_genome_size1 <- log10(TE_size$Genome_Size)
log_te_size1 <- log10(TE_size$total_size)

correlation_test1 <- cor.test(log_genome_size1, log_te_size1, method = "pearson")

# 输出相关性系数和P值
correlation_coefficient1 <- correlation_test1$estimate
p_value1<- correlation_test1$p.value
print(correlation_coefficient1) # 0.9052379
print(p_value1) #  2.083162e-05



###########################################################
#绘制转座子size各个家族与基因组大小的相关性图
te_size_data <- allTE_summ_Class %>%
  select(Species, Classification, size, Genome_Size) %>%
  mutate(log_size = log10(size + 1), log_genome_size = log10(Genome_Size))

TE_size <- te_size_data %>%
  group_by(Species) %>%
  summarise(total_size = sum(size)) %>%
  left_join(genome_sizes, by = "Species")

ggplot(TE_size, aes(x = log10(Genome_Size), y = log10(total_size))) +
  geom_point(color = "black", size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "#8074AC", size = 1, fill = "#0055A5") +
  geom_text_repel(aes(label = Species), size = 3, box.padding = 0.2, point.padding = 0.1,
                  segment.color = 'grey50') +
  theme_classic() +
  labs(x = "genome size (log10(bp))", y = "TE size (log10(bp))") +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))


te_size_data <- allTE_summ_Class %>%
  select(Species, Classification, size, Genome_Size) %>%
  mutate(log_size = log10(size + 1), log_genome_size = log10(Genome_Size)) %>%
  select(Species, Classification,log_size,log_genome_size )


correlation_results_size <- te_size_data %>%
  filter(!is.na(log_size) & !is.na(log_genome_size)) %>% # 移除缺失值
  group_by(Classification) %>%
  summarize(
    correlation = if(n() >= 3) cor(log_size, log_genome_size, method = "pearson") else NA,
    p_value = if(n() >= 3) cor.test(log_size, log_genome_size, method = "pearson")$p.value else NA
  ) %>%
  ungroup() %>%
  filter(Classification != "TIR/Sola2") %>%
  filter(Classification != "pararetrovirus")

# 打印结果
print(correlation_results_size)


# 绘制气泡图
library(ggplot2)
library(ggrepel)
correlation_results_size <- correlation_results_size[correlation_results_size$correlation > 0.5, ]
# 绘制气泡图并添加标签
ggplot(correlation_results_size, aes(x = reorder(Classification, correlation), y = correlation, size = -log10(p_value), color = correlation)) +
  geom_text_repel(aes(label = sprintf("r = %.2f\np = %.2e", correlation, p_value)),size = 2.5, box.padding = 0.1, point.padding = 0.1 ,
                  fontface = "bold") +
  geom_point(alpha = 0.5) +
  scale_size_continuous(range = c(3, 15), guide = guide_legend(override.aes = list(alpha = 0.5))) +
  scale_color_gradient2(midpoint = 0.5, low = "green", mid = "blue", high = "red", guide = guide_colorbar(barwidth = 1, barheight = 2.5)) +
  scale_y_continuous(limits = c(0.5, 1)) +
  labs(x = NULL,
       size = "-log10(p-value)",
       color = "Correlation") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6), 
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.key.size = unit(0.5, "cm"), # 调整图例键的大小
        legend.text = element_text(size = 6), # 调整图例文本的大小
        panel.grid.major = element_line(color = "gray", linetype = "dashed"),
        panel.grid.minor = element_line(color = "gray", linetype = "dashed"))



##########count
# 绘制气泡图
library(ggplot2)
library(ggrepel)

correlation_results_count <- correlation_results_count[correlation_results_count$correlation > 0.5, ] 
# 绘制气泡图并添加标签
ggplot(correlation_results_count, aes(x = reorder(Classification, correlation), y = correlation, size = -log10(p_value), color = correlation)) +
  geom_text_repel(aes(label = sprintf("r = %.2f\np = %.2e", correlation, p_value)),size = 2.5, box.padding = 0.1, point.padding = 0.1 ,
                  fontface = "bold") +
  geom_point(alpha = 0.5) +
  scale_size_continuous(range = c(3, 15), guide = guide_legend(override.aes = list(alpha = 0.5))) +
  scale_color_gradient2(midpoint = 0.5, low = "green", mid = "blue", high = "red", guide = guide_colorbar(barwidth = 1, barheight = 2.5)) +
  scale_y_continuous(limits = c(0.5, 1)) +
  labs(x = NULL,
       size = "-log10(p-value)",
       color = "Correlation") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6), 
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.key.size = unit(0.5, "cm"), # 调整图例键的大小
        legend.text = element_text(size = 6), # 调整图例文本的大小
        panel.grid.major = element_line(color = "gray", linetype = "dashed"),
        panel.grid.minor = element_line(color = "gray", linetype = "dashed"))
