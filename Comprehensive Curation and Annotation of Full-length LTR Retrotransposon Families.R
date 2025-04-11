#完整的LTR
#数据整理


# 加载必要的R包
library(dplyr)
library(readr)
library(rtracklayer)  # 用于导入GFF文件
library(stringr)
# 定义需要处理的物种列表
species_list <- c("Aed", "Arg", "Dap", "Eri","Dro", "Hom", "Hya", "Lep", "Pen", "Pen.m", "Pen.c", "Por", "Pro")

# 初始化一个列表来存储每个物种的数据框
species_data_frames <- list()

# 批量处理每个物种的数据并将结果保存在列表中
for (species in species_list) {
  # 导入并处理每个物种的 GFF 文件
  intactLTR_gff_file <- paste0("31.intactLTR/", species, "/intact.LTR.gff3")
  intactLTR_gff <- import(intactLTR_gff_file)
  intactLTR_gff <- as_tibble(intactLTR_gff) %>%
    select(seqnames, start, end, width, strand, ID, Name, ltr_identity, motif, tsd)
  
  # 读取并处理每个物种的 LTR_age 数据
  intactLTR_age_file <- paste0("23.EDTA/", species, "/", species, ".genome.fa.mod.EDTA.raw/LTR/", species, ".genome.fa.mod.pass.list")
  intactLTR_age <- read_table(intactLTR_age_file) %>%
    select(LTR_loc = `#LTR_loc`, Insertion_Time) %>%
    distinct()
  
  # 变异 intactLTR 并与 intactLTR_age 连接
  intactLTR <- mutate(intactLTR_gff, 
                      LTR_loc = str_c(seqnames, ":", start, "..", end)) %>%
    inner_join(intactLTR_age, by = "LTR_loc")
  
  # 读取每个物种的分类信息文件
  intactLTR_cls_file <- paste0("31.intactLTR/", species, "/intact.LTR.fa.rexdb.cls.tsv")
  intactLTR_cls <- read_delim(intactLTR_cls_file, 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)
  
  # 执行左连接并选择相关列
  intactLTR_final <- left_join(intactLTR, 
                               intactLTR_cls, 
                               by = c("ID" = "#TE")) %>%
    select(ID, Name, seqnames, start, end, width, ltr_identity, motif, tsd, Order, 
           Superfamily, Clade, Complete, Domains, Insertion_Time)
  
  # 将数据框保存到列表中，并命名为 物种名.species_data
  species_data_name <- paste0(species, ".species_data")
  species_data_frames[[species_data_name]] <- intactLTR_final
}


all.species_data <- bind_rows(species_data_frames, .id = "Species")


all.species_data.filter <- filter(all.species_data, !is.na(Superfamily) & !(Superfamily %in% c("unknown" ,"NA", "Tc1_Marine" ,"PiggyBac")))


selected_data <- select(all.species_data.filter, Species, Superfamily, Insertion_Time)


# 修改物种名称
selected_data <- selected_data %>%
  mutate(Species = case_when(
    grepl("Aed", Species) ~ "(A.aeg)",
    grepl("Arg", Species) ~ "(A.bru)",
    grepl("Dap", Species) ~ "(D.pul)",
    grepl("Eri", Species) ~ "(E.sin)",
    grepl("Hom", Species) ~ "(H.ame)",
    grepl("Hya", Species) ~ "(H.azt)",
    grepl("Lep", Species) ~ "(L.sal)",
    grepl("Pen\\.c", Species) ~ "(P.chi)",
    grepl("Pen\\.m", Species) ~ "(P.mon)",
    grepl("Pen", Species) ~ "(P.van)",
    grepl("Por", Species) ~ "(P.tri)",
    grepl("Pro", Species) ~ "(P.cla)",
    grepl("Dro", Species) ~ "(D.mel)",
    TRUE ~ Species
  )) %>%
mutate(Species = case_when(
  grepl("(A.aeg)", Species) ~ "A.aeg",
  grepl("(A.bru)", Species) ~ "A.bru",
  grepl("(D.pul)", Species) ~ "D.pul",
  grepl("(E.sin)", Species) ~ "E.sin",
  grepl("(H.ame)", Species) ~ "H.ame",
  grepl("(H.azt)", Species) ~ "H.azt",
  grepl("(L.sal)", Species) ~ "L.sal",
  grepl("(P.chi)", Species) ~ "P.chi",
  grepl("(P.mon)", Species) ~ "P.mon",
  grepl("(P.van)", Species) ~ "P.van",
  grepl("(P.tri)", Species) ~ "P.tri",
  grepl("(P.cla)", Species) ~ "P.cla",
  grepl("(D.mel)", Species) ~ "D.mel",
  TRUE ~ Species
))
############################
#绘图
library(ggplot2)
library(ggsci)
library(dplyr)

# 将插入时间转换为百万年前
selected_data <- selected_data %>%
  mutate(Insertion_Time_mya = Insertion_Time / 1e6) %>%
  filter(!is.na(Superfamily) & !(Superfamily %in% c("NA", "mixture")))

# 获取插入时间的实际范围
insertion_time_range <- range(selected_data$Insertion_Time_mya, na.rm = TRUE)

# 设置Y轴的上下限，适当增加一些边界
y_axis_limits <- c(0, 4)

custom_species_order <- c("E.sin", "P.tri", "P.cla", "H.ame", "P.van", "P.chi", "P.mon", "H.azt", "L.sal", "D.pul", "D.mel", "A.aeg", "A.bru")
custom_supmerfamily_order <- c("mixture","Bel-Pao", "Copia",  "Gypsy")
selected_data$Species <- factor(selected_data$Species, levels = custom_species_order)
selected_data$Superfamily <- factor(selected_data$Superfamily, levels = custom_supmerfamily_order)

violin_color <- c( "Gypsy" = "#D55740" , "Copia" = "#6CB9D2", "Bel-Pao" = "#479E88", "mixture"= "#f4a01a")
# 绘制提琴图，将 X 轴和 Y 轴颠倒，并统一 Y 轴刻度
selected_data<- selected_data %>%
  filter(Insertion_Time_mya <= 3)

ggplot(selected_data, aes(x = Species, y = Insertion_Time_mya, fill = Superfamily)) +
  geom_violin(trim = T , size = 0.1) +
  scale_fill_manual(values = violin_color) +
  facet_grid(Superfamily ~ ., scales = "fixed") +  # 每个Superfamily排列为一行，Y轴刻度固定
  labs(title = NULL,
       x = "Species",
       y = "Insertion Time (mya)") +  # 更新Y轴标签
  scale_y_reverse(limits = c(3, 0), breaks = c(0, 1,2,3)) +  # 设置Y轴刻度并颠倒
 # scale_y_continuous(limits = c(0, 4), breaks = c(0, 2, 4)) +  # 设置Y轴刻度
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # 旋转X轴标签以提高可读性
        panel.border = element_rect(color = "black", fill = NA, size = 0.1),# 为每个模块添加边框
        legend.position = "top",
        panel.spacing = unit(0.3, "lines")) + # 缩小两个物种之间的距离
  guides(fill = guide_legend(title = NULL ,ncol = 6))  # 移除图例的标题


##############################
library(ggplot2)
library(ggsci)
library(dplyr)

# 假设 selected_data 是你的数据框，并且插入时间单位是年
# 将插入时间转换为百万年前
selected_data <- selected_data %>%
  mutate(Insertion_Time_mya = Insertion_Time / 1e6)

# 计算每个物种的插入时间的中位数
ordered_species <- selected_data %>%
  group_by(Species) %>%
  summarise(median_insertion_time = median(Insertion_Time_mya, na.rm = TRUE)) %>%
  arrange(median_insertion_time) %>%
  pull(Species)

# 将 Species 转换为因子，并按中位数插入时间的顺序排列
selected_data$Species <- factor(selected_data$Species, levels = ordered_species)

# 计算每个超家族的插入时间的中位数
ordered_superfamily <- selected_data %>%
  group_by(Superfamily) %>%
  summarise(median_insertion_time = median(Insertion_Time_mya, na.rm = TRUE)) %>%
  arrange(median_insertion_time) %>%
  pull(Superfamily)

# 将 Superfamily 转换为因子，并按中位数插入时间的顺序排列
selected_data$Superfamily <- factor(selected_data$Superfamily, levels = ordered_superfamily)

# 绘制提琴图
ggplot(selected_data, aes(x = Species, y = Insertion_Time_mya, fill = Superfamily)) +
  geom_violin(trim = TRUE) +
  scale_fill_npg() +
  facet_wrap(~ Superfamily, scales = "free_y") +  # 按超家族分面显示，允许每个面的Y轴自由缩放
  labs(title = "Distribution of Insertion Times by Species and TE Superfamily",
       x = "Species",
       y = "Insertion Time (mya)") +  # 更新Y轴标签
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # 旋转X轴标签以提高可读性


###############################

# 现在可以通过物种名访问任何一个数据框，例如：
species_data_frames[["Arg.species_data"]]  # 访问物种 Arg 的数据框




#######################################################################################################

#绘图

#转化时间
selected_data <- selected_data %>%
  mutate(Insertion_Time_mya = Insertion_Time / 1e6)
# 计算Insertion_Time的全局最小和最大值
min_time <- min(sapply(species_data_frames, function(df) min(df$Insertion_Time, na.rm = TRUE)))
max_time <- max(sapply(species_data_frames, function(df) max(df$Insertion_Time, na.rm = TRUE)))

min_time
max_time
# 加载必要的R包
library(ggplot2)
library(dplyr)
library(ggsci)

# 定义颜色映射
superfamily_colors <- c(
  "Bel-Pao" = "#1D4787",  # 蓝色
  "Copia" = "#FF0000",    # 红色
  "Gypsy" = "#69B657"    # 绿色
  )




# 假设列表已经加载到R环境中，我们直接访问Arg的数据框
arg_data <- species_data_frames[["Arg.species_data"]]


#Arg
# 对数据进行必要的过滤
arg_filtered <- filter(arg_data, !is.na(Superfamily) & !(Superfamily %in% c("unknown" ,"NA" , "mixture", "Tc1_Mariner")))

# 绘制密度图 
Arg.p <- ggplot(arg_filtered, aes(x = Insertion_Time, color = Superfamily)) +
  geom_density(aes(y = ..density..), adjust=1) +  # adjust参数可以调整平滑程度
  scale_x_continuous(name = 'Insertion Time (Mys)', 
                     breaks = seq(0, 5000000, 1000000), labels = seq(0, 5, 1), limits = c(0, 5000000)) +  # 移除X轴刻度和标签
  #scale_y_continuous(labels = NULL) +  # 移除Y轴标签
  scale_color_manual(values = superfamily_colors) +  # 应用自定义颜色
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank()   # 去掉次网格线
         #axis.text.x = element_blank(),  # 隐藏X轴刻度文本
        #axis.title.x = element_blank(),  # 隐藏X轴标题
        #axis.text.y = element_blank(),  # 隐藏Y轴刻度文本
        #axis.title.y = element_blank())  # 隐藏Y轴标题
)
# 显示图形
print(Arg.p)



#Aed
aed_data <- species_data_frames[["Aed.species_data"]]

# 对数据进行必要的过滤
aed_filtered <- filter(aed_data, !is.na(Superfamily) & !(Superfamily %in% c("unknown" ,"NA" , "mixture", "Tc1_Mariner")))

# 绘制密度图 
Aed.p <- ggplot(aed_filtered, aes(x = Insertion_Time, color = Superfamily)) +
  geom_density(aes(y = ..density..), adjust=1) +  # adjust参数可以调整平滑程度
  scale_x_continuous(name = 'Insertion Time (Mys)', 
                     breaks = seq(0, 5000000, 1000000), labels = seq(0, 5, 1)) +  
  #scale_y_continuous(labels = NULL) +  # 移除Y轴标签
  scale_color_manual(values = superfamily_colors) +  # 应用自定义颜色
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank()   # 去掉次网格线
        #axis.text.x = element_blank(),  # 隐藏X轴刻度文本
        #axis.title.x = element_blank(),  # 隐藏X轴标题
        #axis.text.y = element_blank(),  # 隐藏Y轴刻度文本
        #axis.title.y = element_blank())  # 隐藏Y轴标题
  )
# 显示图形
print(Aed.p)





#Cap

cap_data <- species_data_frames[["Cap.species_data"]]
cap_filtered <- filter(cap_data, !is.na(Superfamily) & !(Superfamily %in% c("unknown", "NA")))
Cap_plot <- ggplot(cap_filtered, aes(x = Insertion_Time, color = Superfamily)) +
  geom_density(aes(y = ..density..), adjust=1) +
  scale_x_continuous(name = 'Insertion Time (Mys)', limits = c(min_time, max_time)) +
  scale_color_manual(values = superfamily_colors) +
  theme_bw() +
  theme(legend.position = 'none',
        #axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
print(Cap_plot)


#Dap
dap_data <- species_data_frames[["Dap.species_data"]]
# 对数据进行必要的过滤
dap_filtered <- filter(dap_data, !is.na(Superfamily) & !(Superfamily %in% c("unknown" ,"NA" , "mixture", "Tc1_Mariner", "PiggyBac")))

# 绘制密度图 
dap.p <- ggplot(dap_filtered, aes(x = Insertion_Time, color = Superfamily)) +
  geom_density(aes(y = ..density..), adjust=1) +  # adjust参数可以调整平滑程度
  scale_x_continuous(name = 'Insertion Time (Mys)', 
                     breaks = seq(0, 5000000, 1000000), labels = seq(0, 5, 1),limits = c(0, 5000000)) +  # 移除X轴刻度和标签
  #scale_y_continuous(labels = NULL) +  # 移除Y轴标签
  scale_color_manual(values = superfamily_colors) +  # 应用自定义颜色
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank()   # 去掉次网格线
        #axis.text.x = element_blank(),  # 隐藏X轴刻度文本
        #axis.title.x = element_blank(),  # 隐藏X轴标题
        #axis.text.y = element_blank(),  # 隐藏Y轴刻度文本
        #axis.title.y = element_blank())  # 隐藏Y轴标题
  )
# 显示图形
print(dap.p)


#Eri

eri_data <- species_data_frames[["Eri.species_data"]]
# 对数据进行必要的过滤
eri_filtered <- filter(eri_data, !is.na(Superfamily) & !(Superfamily %in% c("unknown" ,"NA" , "mixture", "Tc1_Mariner", "PiggyBac")))

# 绘制密度图 
eri.p <- ggplot(eri_filtered, aes(x = Insertion_Time, color = Superfamily)) +
  geom_density(aes(y = ..density..), adjust=1) +  # adjust参数可以调整平滑程度
  scale_x_continuous(name = 'Insertion Time (Mys)', 
                     breaks = seq(0, 5000000, 1000000), labels = seq(0, 5, 1), limits = c(0, 5000000)) +  # 移除X轴刻度和标签
  #scale_y_continuous(labels = NULL) +  # 移除Y轴标签
  scale_color_manual(values = superfamily_colors) +  # 应用自定义颜色
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank()   # 去掉次网格线
        #axis.text.x = element_blank(),  # 隐藏X轴刻度文本
        #axis.title.x = element_blank(),  # 隐藏X轴标题
        #axis.text.y = element_blank(),  # 隐藏Y轴刻度文本
        #axis.title.y = element_blank())  # 隐藏Y轴标题
  )
# 显示图形
print(eri.p)

#Dro
# 对数据进行必要的过滤
dro_data <- species_data_frames[["Dro.species_data"]]
dro_filtered <- filter(dro_data, !is.na(Superfamily) & !(Superfamily %in% c("unknown", "NA", "mixture", "Tc1_Mariner", "PiggyBac")))

# 绘制密度图 
dro.p <- ggplot(dro_filtered, aes(x = Insertion_Time, color = Superfamily)) +
  geom_density(aes(y = ..density..), adjust=1) +  # adjust参数可以调整平滑程度
  scale_x_continuous(name = 'Insertion Time (Mys)', 
                     breaks = seq(0, 5000000, 1000000), labels = seq(0, 5, 1), limits = c(0, 5000000)) +  # 设置X轴刻度和标签
  scale_color_manual(values = superfamily_colors) +  # 应用自定义颜色
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank()   # 去掉次网格线
  )
# 显示图形
print(dro.p)
 
#Hom
# 对数据进行必要的过滤
hom_data <- species_data_frames[["Hom.species_data"]]
hom_filtered <- filter(hom_data, !is.na(Superfamily) & !(Superfamily %in% c("unknown", "NA", "mixture", "Tc1_Mariner", "PiggyBac")))

# 绘制密度图 
hom.p <- ggplot(hom_filtered, aes(x = Insertion_Time, color = Superfamily)) +
  geom_density(aes(y = ..density..), adjust=1) +  # adjust参数可以调整平滑程度
  scale_x_continuous(name = 'Insertion Time (Mys)', 
                     breaks = seq(0, 5000000, 1000000), labels = seq(0, 5, 1), limits = c(0, 5000000)) +  # 设置X轴刻度和标签
  scale_color_manual(values = superfamily_colors) +  # 应用自定义颜色
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank()   # 去掉次网格线
  )
# 显示图形
print(hom.p)


#Hya
# 对数据进行必要的过滤
hya_data <- species_data_frames[["Hya.species_data"]]
hya_filtered <- filter(hya_data, !is.na(Superfamily) & !(Superfamily %in% c("unknown", "NA", "mixture", "Tc1_Mariner", "PiggyBac")))

# 绘制密度图 
hya.p <- ggplot(hya_filtered, aes(x = Insertion_Time, color = Superfamily)) +
  geom_density(aes(y = ..density..), adjust=1) +  # adjust参数可以调整平滑程度
  scale_x_continuous(name = 'Insertion Time (Mys)', 
                     breaks = seq(0, 5000000, 1000000), labels = seq(0, 5, 1), limits = c(0, 5000000)) +  # 设置X轴刻度和标签
  scale_color_manual(values = superfamily_colors) +  # 应用自定义颜色
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank()   # 去掉次网格线
  )
# 显示图形
print(hya.p)


#Lep
# 对数据进行必要的过滤
lep_data <- species_data_frames[["Lep.species_data"]]
lep_filtered <- filter(lep_data, !is.na(Superfamily) & !(Superfamily %in% c("unknown", "NA", "mixture", "Tc1_Mariner", "PiggyBac")))

# 绘制密度图 
lep.p <- ggplot(lep_filtered, aes(x = Insertion_Time, color = Superfamily)) +
  geom_density(aes(y = ..density..), adjust=1) +  # adjust参数可以调整平滑程度
  scale_x_continuous(name = 'Insertion Time (Mys)', 
                     breaks = seq(0, 5000000, 1000000), labels = seq(0, 5, 1), limits = c(0, 5000000)) +  # 设置X轴刻度和标签
  scale_color_manual(values = superfamily_colors) +  # 应用自定义颜色
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank()   # 去掉次网格线
  )
# 显示图形
print(lep.p)



#Pen

# 对数据进行必要的过滤
pen_data <- species_data_frames[["Pen.species_data"]]
pen_filtered <- filter(pen_data, !is.na(Superfamily) & !(Superfamily %in% c("unknown", "NA", "mixture", "Tc1_Mariner", "PiggyBac")))

# 绘制密度图 
pen.p <- ggplot(pen_filtered, aes(x = Insertion_Time, color = Superfamily)) +
  geom_density(aes(y = ..density..), adjust=1) +  # adjust参数可以调整平滑程度
  scale_x_continuous(name = 'Insertion Time (Mys)', 
                     breaks = seq(0, 5000000, 1000000), labels = seq(0, 5, 1), limits = c(0, 5000000)) +  # 设置X轴刻度和标签
  scale_color_manual(values = superfamily_colors) +  # 应用自定义颜色
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank()   # 去掉次网格线
  )
# 显示图形
print(pen.p)


#Pen.m
# 对数据进行必要的过滤
pen_m_data <- species_data_frames[["Pen.m.species_data"]]
pen_m_filtered <- filter(pen_m_data, !is.na(Superfamily) & !(Superfamily %in% c("unknown", "NA", "mixture", "Tc1_Mariner", "PiggyBac")))

# 绘制密度图 
pen_m.p <- ggplot(pen_m_filtered, aes(x = Insertion_Time, color = Superfamily)) +
  geom_density(aes(y = ..density..), adjust=1) +  # adjust参数可以调整平滑程度
  scale_x_continuous(name = 'Insertion Time (Mys)', 
                     breaks = seq(0, 5000000, 1000000), labels = seq(0, 5, 1), limits = c(0, 5000000)) +  # 设置X轴刻度和标签
  scale_color_manual(values = superfamily_colors) +  # 应用自定义颜色
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank()   # 去掉次网格线
  )
# 显示图形
print(pen_m.p)
 


#Pen.c
# 对数据进行必要的过滤
pen_c_data <- species_data_frames[["Pen.c.species_data"]]
pen_c_filtered <- filter(pen_c_data, !is.na(Superfamily) & !(Superfamily %in% c("unknown", "NA", "mixture", "Tc1_Mariner", "PiggyBac")))

# 绘制密度图 
pen_c.p <- ggplot(pen_c_filtered, aes(x = Insertion_Time, color = Superfamily)) +
  geom_density(aes(y = ..density..), adjust=1) +  # adjust参数可以调整平滑程度
  scale_x_continuous(name = 'Insertion Time (Mys)', 
                     breaks = seq(0, 5000000, 1000000), labels = seq(0, 5, 1), limits = c(0, 5000000)) +  # 设置X轴刻度和标签
  scale_color_manual(values = superfamily_colors) +  # 应用自定义颜色
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank()   # 去掉次网格线
  )
# 显示图形
print(pen_c.p)


#Por
# 对数据进行必要的过滤
por_data <- species_data_frames[["Por.species_data"]]
por_filtered <- filter(por_data, !is.na(Superfamily) & !(Superfamily %in% c("unknown", "NA", "mixture", "Tc1_Mariner", "PiggyBac")))

# 绘制密度图 
por.p <- ggplot(por_filtered, aes(x = Insertion_Time, color = Superfamily)) +
  geom_density(aes(y = ..density..), adjust=1) +  # adjust参数可以调整平滑程度
  scale_x_continuous(name = 'Insertion Time (Mys)', 
                     breaks = seq(0, 5000000, 1000000), labels = seq(0, 5, 1), limits = c(0, 5000000)) +  # 设置X轴刻度和标签
  scale_color_manual(values = superfamily_colors) +  # 应用自定义颜色
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank()   # 去掉次网格线
  )
# 显示图形
print(por.p)


#Pro
# 对数据进行必要的过滤
pro_data <- species_data_frames[["Pro.species_data"]]
pro_filtered <- filter(pro_data, !is.na(Superfamily) & !(Superfamily %in% c("unknown", "NA", "mixture", "Tc1_Mariner", "PiggyBac")))

# 绘制密度图 
pro.p <- ggplot(pro_filtered, aes(x = Insertion_Time, color = Superfamily)) +
  geom_density(aes(y = ..density..), adjust=1) +  # adjust参数可以调整平滑程度
  scale_x_continuous(name = 'Insertion Time (Mys)', 
                     breaks = seq(0, 5000000, 1000000), labels = seq(0, 5, 1), limits = c(0, 5000000)) +  # 设置X轴刻度和标签
  scale_color_manual(values = superfamily_colors) +  # 应用自定义颜色
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank()   # 去掉次网格线
  )
# 显示图形
print(pro.p)


#Tac

tac_data <- species_data_frames[["Tac.species_data"]]
tac_filtered <- filter(tac_data, !is.na(Superfamily) & !(Superfamily %in% c("unknown", "NA")))
Tac_plot <- ggplot(tac_filtered, aes(x = Insertion_Time, color = Superfamily)) +
  geom_density(aes(y = ..density..), adjust=1) +
  scale_x_continuous(name = 'Insertion Time (Mys)', limits = c(min_time, max_time)) +
  scale_color_manual(values = superfamily_colors) +
  theme_bw() +
  theme(legend.position = 'none',
       # axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
print(Tac_plot)

#Tig
tig_data <- species_data_frames[["Tig.species_data"]]
tig_filtered <- filter(tig_data, !is.na(Superfamily) & !(Superfamily %in% c("unknown", "NA")))
Tig_plot <- ggplot(tig_filtered, aes(x = Insertion_Time, color = Superfamily)) +
  geom_density(aes(y = ..density..), adjust=1) +
  scale_x_continuous(name = 'Insertion Time (Mys)', limits = c(min_time, max_time)) +
  scale_color_manual(values = superfamily_colors) +
  theme_bw() +
  theme(legend.position = 'top',
        #axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
print(Tig_plot)


####################################
# 对每个物种进行绘图
for (species_name in names(species_data_frames)) {
  species_data <- species_data_frames[[species_name]]
  species_filtered <- filter(species_data, !is.na(Superfamily) & !(Superfamily %in% c("unknown", "NA")))
  
  plot <- ggplot(species_filtered, aes(x = Insertion_Time, color = Superfamily)) +
    geom_density(aes(y = ..density..), adjust=1) +
    scale_x_continuous(name = 'Insertion Time (Mys)', limits = c(min_time, max_time)) +
    # scale_y_continuous(limits = c(min_density, max_density)) # 如果你决定设置Y轴范围
    scale_color_manual(values = superfamily_colors ) +
    theme_classic() +
    theme(legend.position = "right") +
    ggtitle(paste("Insertion Time Density for", species_name))
  
  print(plot)
}


