


#数据整理
# 物种名向量
species_names <- c("Aed","Arg", "Cap", "Dap", "Dro", "Eri", "Hom", "Hya", "Lep", "Pen", "Pen.m", "Pen.c", "Por", "Pro", "Tac", "Tig")

# 通过循环处理每个物种的数据
for (species in species_names) {
  # 动态构建原始数据的变量名
  data_var_name <- paste0(species, ".TE_divsum")
  
  # 动态构建处理后数据的变量名
  summarized_var_name <- paste0(species, "_summarized")
  
  # 从全局环境中获取对应的数据框
  species_data <- get(data_var_name, envir = .GlobalEnv)
  
  # 修改分类信息，将包含"SINE"的分类统一标记为"SINE"
  species_data$Classification <- ifelse(grepl("SINE", species_data$Classification), "SINE", species_data$Classification)
  
  # 按照Div和Classification进行分组，并对percent列进行求和
  summarized_data <- species_data %>% 
    group_by(Div, Classification) %>%
  summarise(percent = sum(percent), .groups = 'drop')  %>%
  filter(Classification != "Unknown")
  # 将处理后的数据框存储回全局环境，使用新的变量名
  assign(summarized_var_name, summarized_data, envir = .GlobalEnv)
}

library(ggsci)
library(scales)

#创建颜色分类
categories <- c("DNA.DTA", "DNA.DTC", "DNA.DTH", "DNA.DTT","DNA.DTM",
                "DNA.Helitron", "LINE.unknown", "LTR.Copia", "LTR.Gypsy",
                "LTR.unknown", "LTR.unknown_BEL", "MITE.DTA", "MITE.DTC",
                "MITE.DTH", "MITE.DTM", "MITE.DTT", "Penelope", "polinton",
                "SINE", "TIR.Kolobok", "TIR.PiggyBac", "TIR.Sola2", "TIR.Tc1_Mariner",
                "Unknown")
 
igv_colors_24 <- pal_igv("default")(24)
color_map <- setNames(igv_colors_24, categories)



#Aed1
library(tidyverse)
Aed.TE_divsum =  read.table("25.TEdiv/Aed/TE.divsum", header = TRUE) %>%
  pivot_longer(-Div, names_to = "Classification", 
               values_to = "size") %>%
  mutate(percent = (size/1278732104) * 100,
         Species = "Aed") %>% 
  filter(percent != 0)

# 修改分类信息，将包含"SINE"的分类统一标记为"SINE"
Aed.TE_divsum$Classification <- ifelse(grepl("SINE", Aed.TE_divsum$Classification), "SINE", Aed.TE_divsum$Classification)
# 按照新的Classification列进行分组，并对percent列进行求和
Aed_summarized <- Aed.TE_divsum %>% 
  group_by(Div, Classification) %>%
summarise(percent = sum(percent), .groups = 'drop')

#unknown
library(dplyr)

Aed_unknown_data <- Aed_summarized %>%
  filter(Classification == "Unknown")

ggplot(data = Aed_unknown_data, aes(x = Div, y = percent)) +
  geom_col(aes(fill = Classification), width = 0.8, color = "white") +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_fill_lancet() +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "Kimura 2‐Parameter divergence(%)", 
       y = "TE contents (%)", fill = NULL) + 
  theme_classic() +
  theme(legend.position = NULL) +
  guides(fill = guide_legend(nrow = 12 ,# 设置行数
                             ncol = 2))# 设置列数




#######
Aed_TE.DIV <- ggplot(data = Aed_summarized, aes(x = Div, y = percent)) +
  geom_col(aes(fill = Classification), width = 0.8, color = "white") +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_fill_manual(values = color_map) +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "Kimura 2‐Parameter divergence(%)", 
       y = "TE contents (%)", fill = NULL) + 
  theme_classic() +
  theme(legend.position = "right" ) +
  guides(fill = guide_legend(nrow = 12 ,# 设置行数
                             ncol = 2))# 设置列数



Aed_TE.DIV


#Arg2
library(tidyverse)
Arg.TE_divsum =  read.table("25.TEdiv/Arg/TE.divsum", header = TRUE) %>%
  pivot_longer(-Div, names_to = "Classification", 
               values_to = "size") %>%
  mutate(percent = (size/1778384155) * 100,
         Species = "Arg") %>% 
  filter(percent != 0)
# 修改分类信息，将包含"SINE"的分类统一标记为"SINE"
Arg.TE_divsum$Classification <- ifelse(grepl("SINE", Arg.TE_divsum$Classification), "SINE", Aed.TE_divsum$Classification)
# 按照新的Classification列进行分组，并对percent列进行求和
Arg_summarized <- Arg.TE_divsum %>% 
  group_by(Div, Classification) %>%
  summarise(percent = sum(percent), .groups = 'drop')


Arg_TE.DIV <- ggplot(data = Arg_summarized, aes(x = Div, y = percent)) +
  geom_col(aes(fill = Classification), width = 0.8, color = "white") +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_fill_manual(values = color_map) +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "Kimura 2‐Parameter divergence(%)", 
       y = "TE contents (%)", fill = NULL) + 
  theme_classic() +
  theme(legend.position = "right")

Arg_TE.DIV

#unknown
library(dplyr)

Arg_unknown_data <- Arg.TE_divsum %>%
  filter(Classification == "Unknown")

ggplot(data = Arg_unknown_data, aes(x = Div, y = percent)) +
  geom_col(aes(fill = Classification), width = 0.8, color = "white") +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_fill_lancet() +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "Kimura 2‐Parameter divergence(%)", 
       y = "TE contents (%)", fill = NULL) + 
  theme_classic() +
  theme(legend.position = "right" ) +
  guides(fill = guide_legend(nrow = 12 ,# 设置行数
                             ncol = 2))# 设置列数


#Cap3
library(tidyverse)
Cap.TE_divsum =  read.table("25.TEdiv/Cap/TE.divsum", header = TRUE) %>%
  pivot_longer(-Div, names_to = "Classification", 
               values_to = "size") %>%
  mutate(percent = (size/477053020) * 100,
         Species = "Cap") %>% 
  filter(percent != 0)

# 修改分类信息，将包含"SINE"的分类统一标记为"SINE"
Cap.TE_divsum$Classification <- ifelse(grepl("SINE", Cap.TE_divsum$Classification), "SINE", Aed.TE_divsum$Classification)
# 按照新的Classification列进行分组，并对percent列进行求和
Cap_summarized <- Cap.TE_divsum %>% 
  group_by(Div, Classification) %>%
  summarise(percent = sum(percent), .groups = 'drop')

Cap_TE.DIV <- ggplot(data = Cap_summarized, aes(x = Div, y = percent)) +
  geom_col(aes(fill = Classification), width = 0.8, color = "white") +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_fill_manual(values = color_map) +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "Kimura 2‐Parameter divergence(%)", 
       y = "TE contents (%)", fill = NULL) + 
  theme_classic() +
  theme(legend.position = "right")

Cap_TE.DIV

#unknown
library(dplyr)

Cap_unknown_data <- Cap.TE_divsum %>%
  filter(Classification == "Unknown")

ggplot(data = Cap_unknown_data, aes(x = Div, y = percent)) +
  geom_col(aes(fill = Classification), width = 0.8, color = "white") +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_fill_lancet() +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "Kimura 2‐Parameter divergence(%)", 
       y = "TE contents (%)", fill = NULL) + 
  theme_classic() +
  theme(legend.position = "right" ) +
  guides(fill = guide_legend(nrow = 12 ,# 设置行数
                             ncol = 2))# 设置列数
#Dap4
library(tidyverse)
Dap.TE_divsum =  read.table("25.TEdiv/Dap/TE.divsum", header = TRUE) %>%
  pivot_longer(-Div, names_to = "Classification", 
               values_to = "size") %>%
  mutate(percent = (size/133196385) * 100,
         Species = "Dap") %>% 
  filter(percent != 0)



Dap_TE.DIV <- ggplot(data = Dap_summarized, aes(x = Div, y = percent)) +
  geom_col(aes(fill = Classification), width = 0.8, color = "white") +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_fill_manual(values = color_map) +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "Kimura 2‐Parameter divergence(%)", 
       y = "TE contents (%)", fill = NULL) + 
  theme_classic() +
  theme(legend.position = "none")

Dap_TE.DIV
#unknown
library(dplyr)

Dap_unknown_data <- Dap.TE_divsum %>%
  filter(Classification == "Unknown")

ggplot(data = Dap_unknown_data, aes(x = Div, y = percent)) +
  geom_col(aes(fill = Classification), width = 0.8, color = "white") +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_fill_lancet() +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "Kimura 2‐Parameter divergence(%)", 
       y = "TE contents (%)", fill = NULL) + 
  theme_classic() +
  theme(legend.position = "right" ) +
  guides(fill = guide_legend(nrow = 12 ,# 设置行数
                             ncol = 2))# 设置列数
#Dro5
library(tidyverse)
Dro.TE_divsum =  read.table("25.TEdiv/Dro/TE.divsum", header = TRUE) %>%
  pivot_longer(-Div, names_to = "Classification", 
               values_to = "size") %>%
  mutate(percent = (size/143726002) * 100,
         Species = "Dro"
  )



Dro_TE.DIV <- ggplot(data = Dro_summarized, aes(x = Div, y = percent)) +
  geom_col(aes(fill = Classification), width = 0.8, color = "white") +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_fill_manual(values = color_map) +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "Kimura 2‐Parameter divergence(%)", 
       y = "TE contents (%)", fill = NULL) + 
  theme_classic() +
  theme(legend.position = "none")

Dro_TE.DIV
#unknown
library(dplyr)

Dro_unknown_data <- Dro.TE_divsum %>%
  filter(Classification == "Unknown")

ggplot(data = Dro_unknown_data, aes(x = Div, y = percent)) +
  geom_col(aes(fill = Classification), width = 0.8, color = "white") +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_fill_lancet() +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "Kimura 2‐Parameter divergence(%)", 
       y = "TE contents (%)", fill = NULL) + 
  theme_classic() +
  theme(legend.position = "right" ) +
  guides(fill = guide_legend(nrow = 12 ,# 设置行数
                             ncol = 2))# 设置列数
#Eri6
library(tidyverse)
Eri.TE_divsum =  read.table("25.TEdiv/Eri/TE.divsum", header = TRUE) %>%
  pivot_longer(-Div, names_to = "Classification", 
               values_to = "size") %>%
  mutate(percent = (size/1765827443) * 100,
         Species = "Eri"
  )



Eri_TE.DIV <- ggplot(data = Eri_summarized, aes(x = Div, y = percent)) +
  geom_col(aes(fill = Classification), width = 0.8, color = "white") +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_fill_manual(values = color_map) +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "Kimura 2‐Parameter divergence(%)", 
       y = "TE contents (%)", fill = NULL) + 
  theme_classic() +
  theme(legend.position = "none")

# 显示图形
print(Eri_TE.DIV)
#unknown
library(dplyr)

Eri_unknown_data <- Eri.TE_divsum %>%
  filter(Classification == "Unknown")

ggplot(data = Eri_unknown_data, aes(x = Div, y = percent)) +
  geom_col(aes(fill = Classification), width = 0.8, color = "white") +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_fill_lancet() +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "Kimura 2‐Parameter divergence(%)", 
       y = "TE contents (%)", fill = NULL) + 
  theme_classic() +
  theme(legend.position = "right" ) +
  guides(fill = guide_legend(nrow = 12 ,# 设置行数
                             ncol = 2))# 设置列数
#Hom7
library(tidyverse)
Hom.TE_divsum =  read.table("25.TEdiv/Hom/TE.divsum", header = TRUE) %>%
  pivot_longer(-Div, names_to = "Classification", 
               values_to = "size") %>%
  mutate(percent = (size/2292076018) * 100,
         Species = "Hom") %>% 
  filter(percent != 0)



Hom_TE.DIV <- ggplot(data = Hom_summarized, aes(x = Div, y = percent)) +
  geom_col(aes(fill = Classification), width = 0.8, color = "white") +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_fill_manual(values = color_map) +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "Kimura 2‐Parameter divergence(%)", 
       y = "TE contents (%)", fill = NULL) + 
  theme_classic() +
  theme(legend.position = "none")

Hom_TE.DIV

#unknown
library(dplyr)

Hom_unknown_data <- Hom.TE_divsum %>%
  filter(Classification == "Unknown")

ggplot(data = Hom_unknown_data, aes(x = Div, y = percent)) +
  geom_col(aes(fill = Classification), width = 0.8, color = "white") +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_fill_lancet() +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "Kimura 2‐Parameter divergence(%)", 
       y = "TE contents (%)", fill = NULL) + 
  theme_classic() +
  theme(legend.position = "right" ) +
  guides(fill = guide_legend(nrow = 12 ,# 设置行数
                             ncol = 2))# 设置列数
#Hya8

library(tidyverse)
Hya.TE_divsum =  read.table("25.TEdiv/Hya/TE.divsum", header = TRUE) %>%
  pivot_longer(-Div, names_to = "Classification", 
               values_to = "size") %>%
  mutate(percent = (size/542364171) * 100,
         Species = "Hya") %>% 
  filter(percent != 0)



Hya_TE.DIV <- ggplot(data = Hya_summarized, aes(x = Div, y = percent)) +
  geom_col(aes(fill = Classification), width = 0.8, color = "white") +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_fill_manual(values = color_map) +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "Kimura 2‐Parameter divergence(%)", 
       y = "TE contents (%)", fill = NULL) + 
  theme_classic() +
  theme(legend.position = "none")


Hya_TE.DIV
#unknown
library(dplyr)

Hya_unknown_data <- Hya.TE_divsum %>%
  filter(Classification == "Unknown")

ggplot(data = Hya_unknown_data, aes(x = Div, y = percent)) +
  geom_col(aes(fill = Classification), width = 0.8, color = "white") +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_fill_lancet() +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "Kimura 2‐Parameter divergence(%)", 
       y = "TE contents (%)", fill = NULL) + 
  theme_classic() +
  theme(legend.position = "right" ) +
  guides(fill = guide_legend(nrow = 12 ,# 设置行数
                             ncol = 2))# 设置列数 
#Lep
library(tidyverse)
Lep.TE_divsum =  read.table("25.TEdiv/Lep/TE.divsum", header = TRUE) %>%
  pivot_longer(-Div, names_to = "Classification", 
               values_to = "size") %>%
  mutate(percent = (size/647191672) * 100,
         Species = "Lep") %>% 
  filter(percent != 0)



Lep_TE.DIV <- ggplot(data = Lep_summarized, aes(x = Div, y = percent)) +
  geom_col(aes(fill = Classification), width = 0.8, color = "white") +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_fill_manual(values = color_map) +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "Kimura 2‐Parameter divergence(%)", 
       y = "TE contents (%)", fill = NULL) + 
  theme_classic() +
  theme(legend.position = "none")

Lep_TE.DIV
#unknown
library(dplyr)

Lep_unknown_data <- Lep.TE_divsum %>%
  filter(Classification == "Unknown")

ggplot(data = Lep_unknown_data, aes(x = Div, y = percent)) +
  geom_col(aes(fill = Classification), width = 0.8, color = "white") +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_fill_lancet() +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "Kimura 2‐Parameter divergence(%)", 
       y = "TE contents (%)", fill = NULL) + 
  theme_classic() +
  theme(legend.position = "right" ) +
  guides(fill = guide_legend(nrow = 12 ,# 设置行数
                             ncol = 2))# 设置列数
#Pen9
library(tidyverse)
Pen.TE_divsum =  read.table("25.TEdiv/Pen/TE.divsum", header = TRUE) %>%
  pivot_longer(-Div, names_to = "Classification", 
               values_to = "size") %>%
  mutate(percent = (size/1663581301) * 100,
         Species = "Pen") %>% 
  filter(percent != 0)



Pen_TE.DIV <- ggplot(data = Pen_summarized, aes(x = Div, y = percent)) +
  geom_col(aes(fill = Classification), width = 0.8, color = "white") +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_fill_manual(values = color_map) +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "Kimura 2‐Parameter divergence(%)", 
       y = "TE contents (%)", fill = NULL) + 
  theme_classic() +
  theme(legend.position = "none")

Pen_TE.DIV
#unknown
library(dplyr)

Pen_unknown_data <- Pen.TE_divsum %>%
  filter(Classification == "Unknown")

ggplot(data = Pen_unknown_data, aes(x = Div, y = percent)) +
  geom_col(aes(fill = Classification), width = 0.8, color = "white") +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_fill_lancet() +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "Kimura 2‐Parameter divergence(%)", 
       y = "TE contents (%)", fill = NULL) + 
  theme_classic() +
  theme(legend.position = "right" ) +
  guides(fill = guide_legend(nrow = 12 ,# 设置行数
                             ncol = 2))# 设置列数
#Pen.c10

library(tidyverse)
Pen.c.TE_divsum =  read.table("25.TEdiv/Pen.c/TE.divsum", header = TRUE) %>%
  pivot_longer(-Div, names_to = "Classification", 
               values_to = "size") %>%
  mutate(percent = (size/1466079641) * 100,
         Species = "Pen.c") %>% 
  filter(percent != 0)



Pen.c_TE.DIV <- ggplot(data = Pen.c_summarized, aes(x = Div, y = percent)) +
  geom_col(aes(fill = Classification), width = 0.8, color = "white") +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_fill_manual(values = color_map) +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "Kimura 2‐Parameter divergence(%)", 
       y = "TE contents (%)", fill = NULL) + 
  theme_classic() +
  theme(legend.position = "none")

Pen.c_TE.DIV


#unknown
library(dplyr)

Pen.c_unknown_data <- Pen.c.TE_divsum %>%
  filter(Classification == "Unknown")

ggplot(data = Pen.c_unknown_data, aes(x = Div, y = percent)) +
  geom_col(aes(fill = Classification), width = 0.8, color = "white") +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_fill_lancet() +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "Kimura 2‐Parameter divergence(%)", 
       y = "TE contents (%)", fill = NULL) + 
  theme_classic() +
  theme(legend.position = "right" ) +
  guides(fill = guide_legend(nrow = 12 ,# 设置行数
                             ncol = 2))# 设置列数
#Pen.m11
library(tidyverse)
Pen.m.TE_divsum =  read.table("25.TEdiv/Pen.m/TE.divsum", header = TRUE) %>%
  pivot_longer(-Div, names_to = "Classification", 
               values_to = "size") %>%
  mutate(percent = (size/2391056164) * 100,
         Species = "Pen.m") %>% 
  filter(percent != 0)



Pen.m_TE.DIV <- ggplot(data = Pen.m_summarized, aes(x = Div, y = percent)) +
  geom_col(aes(fill = Classification), width = 0.8, color = "white") +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_fill_manual(values = color_map) +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "Kimura 2‐Parameter divergence(%)", 
       y = "TE contents (%)", fill = NULL) + 
  theme_classic() +
  theme(legend.position = "none")

Pen.m_TE.DIV



#unknown
library(dplyr)

Pen.m_unknown_data <- Pen.m.TE_divsum %>%
  filter(Classification == "Unknown")

ggplot(data = Pen.m_unknown_data, aes(x = Div, y = percent)) +
  geom_col(aes(fill = Classification), width = 0.8, color = "white") +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_fill_lancet() +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "Kimura 2‐Parameter divergence(%)", 
       y = "TE contents (%)", fill = NULL) + 
  theme_classic() +
  theme(legend.position = "right" ) +
  guides(fill = guide_legend(nrow = 12 ,# 设置行数
                             ncol = 2))# 设置列数


#Por12

library(tidyverse)
Por.TE_divsum =  read.table("25.TEdiv/Por/TE.divsum", header = TRUE) %>%
  pivot_longer(-Div, names_to = "Classification", 
               values_to = "size") %>%
  mutate(percent = (size/1005062047) * 100,
         Species = "Por") %>% 
  filter(percent != 0)



Por_TE.DIV <- ggplot(data = Por_summarized, aes(x = Div, y = percent)) +
  geom_col(aes(fill = Classification), width = 0.8, color = "white") +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_fill_manual(values = color_map) +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "Kimura 2‐Parameter divergence(%)", 
       y = "TE contents (%)", fill = NULL) + 
  theme_classic() +
  theme(legend.position = "none")

Por_TE.DIV

#unknown
library(dplyr)

Por_unknown_data <- Por.TE_divsum %>%
  filter(Classification == "Unknown")

ggplot(data = Por_unknown_data, aes(x = Div, y = percent)) +
  geom_col(aes(fill = Classification), width = 0.8, color = "white") +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_fill_lancet() +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "Kimura 2‐Parameter divergence(%)", 
       y = "TE contents (%)", fill = NULL) + 
  theme_classic() +
  theme(legend.position = "right" ) +
  guides(fill = guide_legend(nrow = 12 ,# 设置行数
                             ncol = 2))# 设置列数

#Pro13
library(tidyverse)
Pro.TE_divsum =  read.table("25.TEdiv/Pro/TE.divsum", header = TRUE) %>%
  pivot_longer(-Div, names_to = "Classification", 
               values_to = "size") %>%
  mutate(percent = (size/2735361232) * 100,
         Species = "Pro") %>% 
  filter(percent != 0)



Pro_TE.DIV <- ggplot(data = Pro_summarized, aes(x = Div, y = percent)) +
  geom_col(aes(fill = Classification), width = 0.8, color = "white") +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_fill_manual(values = color_map) +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "Kimura 2‐Parameter divergence(%)", 
       y = "TE contents (%)", fill = NULL) + 
  theme_classic() +
  theme(legend.position = "none")

Pro_TE.DIV
#unknown
library(dplyr)

Pro_unknown_data <- Pro.TE_divsum %>%
  filter(Classification == "Unknown")

ggplot(data = Pro_unknown_data, aes(x = Div, y = percent)) +
  geom_col(aes(fill = Classification), width = 0.8, color = "white") +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_fill_lancet() +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "Kimura 2‐Parameter divergence(%)", 
       y = "TE contents (%)", fill = NULL) + 
  theme_classic() +
  theme(legend.position = "right" ) +
  guides(fill = guide_legend(nrow = 12 ,# 设置行数
                             ncol = 2))# 设置列数
#Tac14

library(tidyverse)
Tac.TE_divsum =  read.table("25.TEdiv/Tac/TE.divsum", header = TRUE) %>%
  pivot_longer(-Div, names_to = "Classification", 
               values_to = "size") %>%
  mutate(percent = (size/2167470406) * 100,
         Species = "Tac") %>% 
  filter(percent != 0)



Tac_TE.DIV <- ggplot(data = Tac_summarized, aes(x = Div, y = percent)) +
  geom_col(aes(fill = Classification), width = 0.8, color = "white") +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_fill_manual(values = color_map) +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "Kimura 2‐Parameter divergence(%)", 
       y = "TE contents (%)", fill = NULL) + 
  theme_classic() +
  theme(legend.position = "none")

Tac_TE.DIV 

#unknown
library(dplyr)

Tac_unknown_data <- Tac.TE_divsum %>%
  filter(Classification == "Unknown")

ggplot(data = Tac_unknown_data, aes(x = Div, y = percent)) +
  geom_col(aes(fill = Classification), width = 0.8, color = "white") +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_fill_lancet() +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "Kimura 2‐Parameter divergence(%)", 
       y = "TE contents (%)", fill = NULL) + 
  theme_classic() +
  theme(legend.position = "right" ) +
  guides(fill = guide_legend(nrow = 12 ,# 设置行数
                             ncol = 2))# 设置列数
#Tig15

library(tidyverse)
Tig.TE_divsum =  read.table("25.TEdiv/Tig/TE.divsum", header = TRUE) %>%
  pivot_longer(-Div, names_to = "Classification", 
               values_to = "size") %>%
  mutate(percent = (size/191142546) * 100,
         Species = "Tig") %>% 
  filter(percent != 0)



Tig_TE.DIV <- ggplot(data = Tig_summarized, aes(x = Div, y = percent)) +
  geom_col(aes(fill = Classification), width = 0.8, color = "white") +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_fill_manual(values = color_map) +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "Kimura 2‐Parameter divergence(%)", 
       y = "TE contents (%)", fill = NULL) + 
  theme_classic() +
  theme(legend.position = "none")

Tig_TE.DIV
#unknown
library(dplyr)

Tig_unknown_data <- Tig.TE_divsum %>%
  filter(Classification == "Unknown")

ggplot(data =Tig_unknown_data, aes(x = Div, y = percent)) +
  geom_col(aes(fill = Classification), width = 0.8, color = "white") +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_fill_lancet() +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "Kimura 2‐Parameter divergence(%)", 
       y = "TE contents (%)", fill = NULL) + 
  theme_classic() +
  theme(legend.position = "right" ) +
  guides(fill = guide_legend(nrow = 12 ,# 设置行数
                             ncol = 2))# 设置列数
