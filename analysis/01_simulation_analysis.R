setwd("/data/muriel/Projects/LMS/LiverModels")

# load packages
library(R.utils)
library(tidyverse)
library(plotly)

# Define cell types
#celltypes <- data.frame(celltype = seq(0,7), name = c("kupffer","hepatocyte","PV","CV","necrotic","macrophage","senescent","proliferated"))

# Create Figure folder
if (!dir.exists("Figures")) {
  dir.create("Figures")
}

### ### ### ### ### ### ###
### Default Mouse Liver ###
### ### ### ### ### ### ###

# get directory names
model_dir <- "08_Mouseliver_default_1-10_10x/"
simu_dirs <- dir(paste0(model_dir,"data/"))

# Define folder code
codes <- data.frame(code = c("0001", "0002","0003","0004", "0005","0006","0007", "0008","0009","0010"),
                    exposure = c(0,250,300,350,400,450,500,550,600,650))

### ### ### ### ### 
###   Figure 1C ### 
### ### ### ### ### 
# Extracellular concentration
file.name <- "log_blood_apap.csv"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")
if (!file.exists(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))) {gunzip(file.path, remove=FALSE)}

# Get first data frame
data <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
data$exposure <- codes %>% 
  filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
  pull(exposure)
data$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))

# Bind all other data frames
for (d in simu_dirs[2:length(simu_dirs)]) {
  file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
  if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
  data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
  data_tmp$exposure <- codes %>% 
    filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
    pull(exposure)
  data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
  data <- bind_rows(data,data_tmp)
}

# Calculate means and sd's
data_sum <- data %>% group_by(time, real_time, exposure) %>% 
  summarise(p_mean = mean(P_pv_local),
            p_sd = sd(P_pv_local)) %>% mutate(exposure = as.factor(exposure))

ggplot(data = data_sum %>% filter(real_time < 15 & real_time > -1 & !(exposure == 0)), aes(x = real_time, y = p_mean, color = exposure)) + 
  geom_point() +
  geom_line() + 
  xlab("Time (h)") + ylab("Blood APAP concentration (a.u.)") +
  scale_color_viridis_d(name = "APAP\ndose (a.u.)") +
  theme_classic()
ggsave(paste0("Figures/","Fig1C_blood_APAP.pdf"), width = 5.5, height = 3)

### ### ### ### ### 
###   Figure 1D ### 
### ### ### ### ### 
# Efflux and uptake concentrations
file.name <- "log_Elimination.csv"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")
if (!file.exists(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))) {gunzip(file.path, remove=FALSE)}

# Get first data frame
data <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
data$exposure <- codes %>% 
  filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
  pull(exposure)
data$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))

# Bind all other data frames
for (d in simu_dirs[2:length(simu_dirs)]) {
  file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
  if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
  data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
  data_tmp$exposure <- codes %>% 
    filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
    pull(exposure)
  data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
  data <- bind_rows(data,data_tmp)
}

# Get a single time point
data_t1 <- data %>% filter(real_time == 1, exposure %in% c(300,450,600))  %>%
  mutate(exposure = as.factor(exposure)) %>% 
  rename(efflux = efflux_perc,
         uptake = uptake_perc) %>%
  pivot_longer(cols = c(efflux,uptake), names_to = "Type", values_to = "Percentage") %>% 
  group_by(real_time,time,exposure, Type) %>% summarise(mean = mean(Percentage),
                                                        sd = sd(Percentage))

ggplot(data_t1, aes(x=Type, y = mean, fill = exposure)) +
  geom_bar(stat = "identity", position = "dodge2", colour = "black") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width=.5,
                position=position_dodge(.9)) +
  geom_text(aes(label = paste0(as.character(round(mean,1)),"%")),vjust=0.4,hjust=-0.2,
            position=position_dodge(0.95), size = 3, angle = 90) +
  scale_color_viridis_d(name = "APAP\ndose (a.u.)", aesthetics = "fill") +
  xlab("Elimination route") + ylab("Percentage (%)") +
  theme_classic() + 
  ggtitle("1 h") + 
  ylim(0,119) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("Figures/","Fig1D_APAP_elimination.pdf"), width = 3, height = 3)


### ### ### ### ### 
###   Figure 1E ### 
### ### ### ### ### 
# Intracellular concentrations
file.name <- "log_SingleCellConc.csv"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")
if (!file.exists(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))) {gunzip(file.path, remove=FALSE)}

# Bind all other data frames
if (interactive()) {
  if (askYesNo("Do you want to load ALL data files?")) {
    print("OK!")
    
    # Get first data frame
    data <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
    data$exposure <- codes %>% 
      filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
      pull(exposure)
    data$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))
    
    for (d in simu_dirs[2:length(simu_dirs)]) {
      file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
      if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
      data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
      data_tmp$exposure <- codes %>% 
        filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
        pull(exposure)
      data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
      data <- bind_rows(data,data_tmp)}
    
    if (interactive()) {
      if (askYesNo("Do you want to save these data? This may overwrite existing data.")) {
        print("OK, saving...")
        save(data, file = paste0(model_dir,"data_log_SingleCellConc.RData"))
        print("Saved!")}}
  } else {
    if (interactive()) {
      if (askYesNo("Do you want to load these data from file?")) {
        load(paste0(model_dir,"data_log_SingleCellConc.RData"))}}}
}

# summarise
data_sum <- data %>% filter(celltype == 1) %>% 
  mutate(exposure = as.factor(exposure)) %>%
  group_by(exposure,run, time, real_time) %>%
  summarize(mean = mean(Pin)) %>% ungroup %>%
  group_by(exposure, time, real_time) %>%
  summarize(Pin_mean = mean(mean),
            Pin_sd = sd(mean))

ggplot(data_sum %>% filter(real_time<50 & real_time>-1, !exposure %in% c(0)),aes(x = real_time, y= Pin_mean, color = exposure)) +
  geom_line() +
  geom_errorbar(aes(ymin = Pin_mean-Pin_sd, ymax = Pin_mean+Pin_sd)) +
  scale_color_viridis_d(name = "APAP\ndose (a.u.)") +
  xlab("Time (h)") + ylab("Intracellular APAP concentration (a.u.)") +
  theme_classic()
ggsave(paste0("Figures/","Fig1E_APAP_intracellular.pdf"), width = 5.5, height = 3)

### ### ### ### ### 
###Figure 1G,H,I### 
### ### ### ### ### 
# GSH graph
data_dist <- data %>% filter(celltype == 1) %>% 
  mutate(exposure = as.factor(exposure)) %>%
  select(-c("cell.id","celltype")) %>%
  group_by(exposure,run, time, real_time, distance) %>%
  summarize_all(list(mean = mean,sd = sd))

subset <- data_dist %>% filter(exposure == 450, run == 1,real_time < 80 & real_time >= -1)

plot_ly(subset, y = ~distance*3, x = ~real_time, z = ~G_cell_mean,
        type = "scatter3d", mode = 'markers',
        color = ~G_cell_mean, marker = list(size = 2)) %>%
  layout(
    scene = list(
      yaxis = list(title = "Distance to CV (\u03bcm)"),
      xaxis = list(title = "Time (h)"),
      zaxis = list(title = "GSH (a.u.)"),
      camera=list(
        eye = list(x=-1, y=-0.88, z=0.64))))

plot_ly(subset, y = ~distance*3, x = ~real_time, z = ~N_mean,
        type = "scatter3d", mode = 'markers',
        color = ~N_mean, marker = list(size = 2)) %>%
  layout(
    scene = list(
      yaxis = list(title = "Distance to CV (\u03bcm)"),
      xaxis = list(title = "Time (h)"),
      zaxis = list(title = "NAPQI (a.u.)"),
      camera=list(
        eye = list(x=-1, y=-0.88, z=0.64))))

plot_ly(subset, y = ~distance*3, x = ~real_time, z = ~C_mean,
        type = "scatter3d", mode = 'markers',
        color = ~C_mean, marker = list(size = 2)) %>%
  layout(
    scene = list(
      yaxis = list(title = "Distance to CV (\u03bcm)"),
      xaxis = list(title = "Time (h)"),
      zaxis = list(title = "NAPQI-cys (a.u.)"),
      camera=list(
        eye = list(x=-1, y=-0.88, z=0.64))))

# plot_ly(subset, x = ~distance*3, y = ~real_time, z = ~N_mean)  %>%
#   add_markers(color = ~N_mean) %>%
#   layout(
#     scene = list(
#       xaxis = list(title = "Distance to CV (\u03bcm)"),
#       yaxis = list(title = "Time (h)"),
#       zaxis = list(title = "NAPQI (a.u.)"),
#       camera=list(
#         eye = list(x=-2, y=-0.88, z=0.64))
#     ), showlegend=FALSE)
# 
# plot_ly(subset, x = ~distance*3, y = ~real_time, z = ~C_mean)  %>%
#   add_markers(color = ~C_mean) %>%
#   layout(
#     scene = list(
#       xaxis = list(title = "Distance to CV (\u03bcm)"),
#       yaxis = list(title = "Time (h)"),
#       zaxis = list(title = "NAPQI-cys (a.u.)"),
#       camera=list(
#         eye = list(x=-2, y=-0.88, z=0.64))
#     ), showlegend=FALSE)
# 

### ### ### ### ### 
###   Figure 1F ### 
### ### ### ### ### 
# Sulfation, glucuronidation and oxidation concentrations
file.name <- "log_oxi_gluc.csv"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")
if (!file.exists(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))) {gunzip(file.path, remove=FALSE)}

# Get first data frame
data <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
data$exposure <- codes %>% 
  filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
  pull(exposure)
data$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))

# Bind all other data frames
for (d in simu_dirs[2:length(simu_dirs)]) {
  file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
  if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
  data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
  data_tmp$exposure <- codes %>% 
    filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
    pull(exposure)
  data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
  data <- bind_rows(data,data_tmp)
}

# Get percentages
data_perc <- data %>% group_by(real_time, time, exposure, run) %>% 
  mutate(total = oxidation_total+gluc_total+sulf_total,
         oxid_perc = oxidation_total/total*100,
         gluc_perc = gluc_total/total*100,
         sulf_perc = sulf_total/total*100) %>% ungroup()

# Get a single time point
data_t6 <- data_perc %>% filter(real_time == 6, !exposure %in% c(0,650))  %>%
  mutate(exposure = as.factor(exposure)) %>% 
  rename(glucuronidation = gluc_total,
         sulfation = sulf_total,
         oxidation = oxidation_total) %>%
  pivot_longer(cols = c(glucuronidation,sulfation,oxidation), names_to = "Type", values_to = "Amount") %>% 
  group_by(real_time,time,exposure, Type) %>% summarise(mean = mean(Amount),
                                                        sd = sd(Amount))
data_t6 <- data_t6 %>% group_by(exposure, time, real_time) %>%
  mutate(total = sum(mean)) %>% mutate(perc = mean/total*100,
                                       perc_sd = sd/total*100) %>%
  mutate(Type = factor(Type, levels = c("glucuronidation","sulfation","oxidation")))
ggplot(data_t6, aes(x=Type, y = perc, fill = exposure)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  geom_text(aes(label = paste0(as.character(round(perc,1)),"%")),vjust=0.2,hjust=-0.5,
            position=position_dodge(0.95), size = 4, angle = 90) +
  geom_errorbar(aes(ymin = perc - perc_sd, ymax = perc + perc_sd), width=.5,
                position=position_dodge(.9)) +
  scale_color_viridis_d(name = "APAP\ndose (a.u.)", aesthetics = "fill") +
  xlab("Elimination route") + ylab("Percentage (%)") +
  theme_classic() + 
  ylim(0,119) +
  ggtitle("6 h") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
ggsave(paste0("Figures/","Fig1F_GlucSulf_Oxi.pdf"), width = 6.5, height = 4)


### ### ### ### ### 
###  Figure 2D  ### 
### ### ### ### ### 
# Cell counts
file.name <- "log_cellcounts.csv"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")
if (!file.exists(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))) {gunzip(file.path, remove=FALSE)}

# Get first data frame
data <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
data$exposure <- codes %>% 
  filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
  pull(exposure)
data$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))

# Bind all other data frames
for (d in simu_dirs[2:length(simu_dirs)]) {
  file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
  if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
  data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
  data_tmp$exposure <- codes %>% 
    filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
    pull(exposure)
  data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
  data <- bind_rows(data,data_tmp)
}

# Rename columns
data <- data %>% dplyr::rename(necrotic = "necrotic_total",
                               total = "hepatocytes_total",
                               senescent = "senescent_total",
                               Kupffer = "kupffer_total",
                               macrophages = "macrophages_total",
                               Time_h = "real_time") %>% 
  mutate(healthy = total - senescent)

# Pivot to long format
dl <- pivot_longer(data, c(necrotic,healthy,senescent,Kupffer,macrophages, total), names_to = "Cells", values_to = "Count") %>% 
  mutate(exposure = as.factor(exposure))

# Calculate mean and sd
dl_sum <- dl %>% filter(!exposure == 0) %>% group_by(time, Time_h, exposure, Cells) %>% summarise(mean_count = mean(Count),
                                                                                                  sd_count = sd(Count)) %>% 
  ungroup() 

# Make labeller
exp_names <- c("250" = "APAP = 250",
               "300" = "APAP = 300",
               "350" = "APAP = 350",
               "400" = "APAP = 400",
               "450" = "APAP = 450",
               "500" = "APAP = 500",
               "550" = "APAP = 550",
               "600" = "APAP = 600",
               "650" = "APAP = 650")

# Make plots
gplot <- ggplot(dl_sum %>% filter(!Cells %in% c("Kupffer","total")), aes(x= Time_h, y = mean_count, color = Cells)) + 
  geom_line(size = 0.5) +
  geom_errorbar(aes(ymin= mean_count-sd_count, ymax = mean_count+sd_count), width = 0.1, alpha = 0.1) +
  #geom_line(data = dl_sum %>% filter(Cells %in% c("total")), aes(x= Time_h, y = mean_count), color = 'grey') + 
  xlab("Time (h)") + ylab("Cell count") +
  theme_classic() +
  scale_color_viridis_d(name = "Cell type") +
  facet_grid(cols = vars(exposure),labeller = as_labeller(exp_names))
ggsave(paste0("Figures/","Fig2D_cell_counts_withoutTotal.pdf"), plot = gplot, width = 11.5, height = 1.5)

colors <- c("grey",viridisLite::viridis(4))
names(colors) <- c("total","healthy", "macrophages","necrotic", "senescent")

gplot <- ggplot(dl_sum %>% filter(!Cells == "Kupffer"), aes(x = Time_h, color = Cells, linetype = Cells)) + 
  geom_line(aes(y = mean_count), size = 0.5) +
  geom_errorbar(aes(ymin= mean_count-sd_count, ymax = mean_count+sd_count), width = 0.1, alpha = 0.1) +
  scale_linetype_manual(name = "Cell Type", values = c("total" = 2, "healthy" = 1, "macrophages" = 1,
                                                       "necrotic" = 1, "senescent" = 1)) +
  scale_color_manual(name = "Cell Type", values = colors) +
  labs(x = "Time (h)", y = "Cell count") + 
  theme_classic() +
  facet_grid(cols = vars(exposure),labeller = as_labeller(exp_names))
ggsave(paste0("Figures/","Fig2D_cell_counts.pdf"), plot = gplot, width = 11.5, height = 1.75)

gplot <- ggplot(dl_sum %>% filter(!Cells == "Kupffer"), aes(x = Time_h, color = Cells, linetype = Cells)) + 
  geom_line(aes(y = mean_count), size = 0.5) +
  geom_errorbar(aes(ymin= mean_count-sd_count, ymax = mean_count+sd_count), width = 0.1, alpha = 0.1) +
  scale_linetype_manual(name = "Cell Type", values = c("total" = 2, "healthy" = 1, "macrophages" = 1,
                                                       "necrotic" = 1, "senescent" = 1)) +
  scale_color_manual(name = "Cell Type", values = colors) +
  labs(x = "Time (h)", y = "Cell count") + 
  theme_classic() +
  theme(legend.position="top", 
        legend.title = element_text(size=5.75, face = "bold"),
        legend.text = element_text(size=5.75),
        legend.margin=margin(t=0, r=0, b=-0.3, l=-0.3, unit="cm")) +
  guides(color = guide_legend(nrow=1,byrow=TRUE)) +
  facet_wrap(~exposure, ncol = 3,labeller = as_labeller(exp_names))
ggsave(paste0("Figures/","Fig2D_cell_counts_long.pdf"), plot = gplot, width = 4, height = 4)

### ### ### ### ### 
###  Figure 2C  ### 
### ### ### ### ### 
# Calculate onset
onset <- dl %>% filter(Cells == "necrotic", Count > 0) %>% group_by(exposure,run) %>%
  summarise(onset = min(Time_h))

# Plot onset
ggplot(onset) +
  geom_line(aes(x = as.numeric(as.character(exposure)), y = onset), color = "grey") + 
  geom_point(aes(x = as.numeric(as.character(exposure)), y = onset, color = exposure), size = 3) + 
  xlab("APAP dose (a.u.)") + ylab("Onset necrosis (h)") +
  theme_classic() +
  scale_color_viridis_d() +
  scale_x_continuous(breaks = seq(250, 650, by = 50)) +
  theme(legend.position = "none")
ggsave(paste0("Figures/","Fig2C_necr_onset.pdf"), width = 3, height = 2)

### ### ### ### ### 
###  Figure 2E  ### 
### ### ### ### ### 
prop <- data %>% filter(Time_h == 168, !exposure %in% c(0,650)) %>% 
  mutate(exposure = as.factor(exposure),
         senescent = senescent/600,
         healthy = healthy/600) %>%
  pivot_longer(cols = c(healthy, senescent), names_to = "Type", values_to = "Proportion") %>%
  group_by(exposure, time, Type) %>% 
  summarise(mean = mean(Proportion),
            sd = sd(Proportion))

ggplot() +
  geom_bar(data = prop, aes(x = exposure, y = mean, fill = Type), stat="identity", position = "dodge") +
  geom_errorbar(data = prop, aes(x = exposure, ymin = mean-sd, ymax = mean+sd, color = Type), position = "dodge") +
  scale_fill_viridis_d() + 
  scale_color_viridis_d() +
  theme_classic() + 
  labs(x = "APAP dose (a.u.)", y = "Proportion of cells") +
  #ggtitle("168 hours") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position=c(0.85,0.85),
        legend.background = element_rect(linetype="solid", 
                                         colour ="lightgrey"),
        legend.title = element_text(size=6, face = "bold"),
        legend.text = element_text(size=6),
        legend.key.size = unit(0.2,"cm"))
ggsave(paste0("Figures/","Fig2E_Proportion.pdf"), width = 3, height = 2)

### ### ### ### ### 
### Figure 2E.2 ### 
### ### ### ### ### 
# get directory names
model_dir <- "08_Mouseliver_tippingpoint_1-9_10x/"
simu_dirs <- dir(paste0(model_dir,"data/"))

# Define folder code
codes <- data.frame(code = c("0001", "0002","0003","0004", "0005","0006","0007", "0008","0009"),
                    exposure = c(505,510,515,520,525,530,535,540,545))

# Cell counts
file.name <- "log_cellcounts.csv"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")
if (!file.exists(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))) {gunzip(file.path, remove=FALSE)}

# Get first data frame
data <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
data$exposure <- codes %>% 
  filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
  pull(exposure)
data$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))

# Bind all other data frames
for (d in simu_dirs[2:length(simu_dirs)]) {
  file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
  if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
  data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
  data_tmp$exposure <- codes %>% 
    filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
    pull(exposure)
  data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
  data <- bind_rows(data,data_tmp)
}

# Rename columns
data <- data %>% dplyr::rename(necrotic = "necrotic_total",
                               total = "hepatocytes_total",
                               senescent = "senescent_total",
                               Kupffer = "kupffer_total",
                               macrophages = "macrophages_total",
                               Time_h = "real_time") %>% 
  mutate(healthy = total - senescent)

# Pivot to long format
dl <- pivot_longer(data, c(necrotic,healthy,senescent,Kupffer,macrophages, total), names_to = "Cells", values_to = "Count") %>% 
  mutate(exposure = as.factor(exposure))

# Calculate mean and sd
dl_sum <- dl %>% filter(!exposure == 0) %>% group_by(time, Time_h, exposure, Cells) %>% summarise(mean_count = mean(Count),
                                                                                                  sd_count = sd(Count)) %>% 
  ungroup() 

# Make labeller
exp_names <- c("505" = "APAP = 505",
               "510" = "APAP = 510",
               "515" = "APAP = 515",
               "520" = "APAP = 520",
               "525" = "APAP = 525",
               "530" = "APAP = 530",
               "535" = "APAP = 535",
               "540" = "APAP = 540",
               "545" = "APAP = 545")

# Make plots
gplot <- ggplot(dl_sum %>% filter(!Cells %in% c("Kupffer","total")), aes(x= Time_h, y = mean_count, color = Cells)) + 
  geom_line(size = 0.5) +
  geom_errorbar(aes(ymin= mean_count-sd_count, ymax = mean_count+sd_count), width = 0.1, alpha = 0.1) +
  #geom_line(data = dl_sum %>% filter(Cells %in% c("total")), aes(x= Time_h, y = mean_count), color = 'grey') + 
  xlab("Time (h)") + ylab("Cell count") +
  theme_classic() +
  scale_color_viridis_d(name = "Cell type") +
  facet_grid(cols = vars(exposure),labeller = as_labeller(exp_names))
#ggsave(paste0("Figures/","FigSX_cell_counts_withoutTotal.pdf"), plot = gplot, width = 11.5, height = 1.5)

colors <- c("grey",viridisLite::viridis(4))
names(colors) <- c("total","healthy", "macrophages","necrotic", "senescent")

gplot <- ggplot(dl_sum %>% filter(!Cells == "Kupffer"), aes(x = Time_h, color = Cells, linetype = Cells)) + 
  geom_line(aes(y = mean_count), size = 0.5) +
  geom_errorbar(aes(ymin= mean_count-sd_count, ymax = mean_count+sd_count), width = 0.1, alpha = 0.1) +
  scale_linetype_manual(name = "Cell Type", values = c("total" = 2, "healthy" = 1, "macrophages" = 1,
                                                       "necrotic" = 1, "senescent" = 1)) +
  scale_color_manual(name = "Cell Type", values = colors) +
  labs(x = "Time (h)", y = "Cell count") + 
  theme_classic() +
  facet_grid(cols = vars(exposure),labeller = as_labeller(exp_names))
ggsave(paste0("Figures/","FigS1_cell_counts.pdf"), plot = gplot, width = 11.5, height = 1.75)

# colors <- c("grey",viridisLite::viridis(4))
# names(colors) <- c("total","healthy", "macrophages","necrotic", "senescent")
# 
# gplot <- ggplot() + 
#   geom_line(data = dl_sum %>% filter(Cells %in% c("total")), aes(x= Time_h, y = mean_count, color = 'total'), linetype = "dashed") + 
#   #geom_errorbar(data = dl_sum %>% filter(Cells %in% c("total")), aes(x = Time_h, ymin= mean_count-sd_count, ymax = mean_count+sd_count, color = "total"), width = 0.1, alpha = 0.1) +
#   geom_line(data = dl_sum %>% filter(Cells %in% c("macrophages")), aes(x = Time_h, y = mean_count, color = "macrophages"), size = 0.5) +
#   geom_errorbar(data = dl_sum %>% filter(Cells %in% c("macrophages")), aes(x = Time_h, ymin= mean_count-sd_count, ymax = mean_count+sd_count, color = "macrophages"), width = 0.1, alpha = 0.1) +
#   geom_line(data = dl_sum %>% filter(Cells %in% c("healthy")), aes(x = Time_h, y = mean_count, color = "healthy"), size = 0.5) +
#   geom_errorbar(data = dl_sum %>% filter(Cells %in% c("healthy")), aes(x = Time_h, ymin= mean_count-sd_count, ymax = mean_count+sd_count, color = "healthy"), width = 0.1, alpha = 0.1) +
#   geom_line(data = dl_sum %>% filter(Cells %in% c("necrotic")), aes(x = Time_h, y = mean_count, color = "necrotic"), size = 0.5) +
#   geom_errorbar(data = dl_sum %>% filter(Cells %in% c("necrotic")), aes(x = Time_h, ymin= mean_count-sd_count, ymax = mean_count+sd_count, color = "necrotic"), width = 0.1, alpha = 0.1) +
#   geom_line(data = dl_sum %>% filter(Cells %in% c("senescent")), aes(x = Time_h, y = mean_count, color = "senescent"), size = 0.5) +
#   geom_errorbar(data = dl_sum %>% filter(Cells %in% c("senescent")), aes(x = Time_h, ymin= mean_count-sd_count, ymax = mean_count+sd_count, color = "senescent"), width = 0.1, alpha = 0.1) +
#   labs(x = "Time (h)", y = "Cell count", color = "Cell Type") +
#   theme_classic() +
#   scale_color_manual(values = colors) +
#   facet_grid(cols = vars(exposure),labeller = as_labeller(exp_names))
# ggsave(paste0("Figures/","FigS1_cell_counts.pdf"), plot = gplot, width = 11.5, height = 1.75)

gplot <- ggplot() + 
  geom_line(data = dl_sum %>% filter(Cells %in% c("total")), aes(x= Time_h, y = mean_count, color = 'total'), linetype = "dashed") + 
  #geom_errorbar(data = dl_sum %>% filter(Cells %in% c("total")), aes(x = Time_h, ymin= mean_count-sd_count, ymax = mean_count+sd_count, color = "total"), width = 0.1, alpha = 0.1) +
  geom_line(data = dl_sum %>% filter(Cells %in% c("macrophages")), aes(x = Time_h, y = mean_count, color = "macrophages"), size = 0.5) +
  geom_errorbar(data = dl_sum %>% filter(Cells %in% c("macrophages")), aes(x = Time_h, ymin= mean_count-sd_count, ymax = mean_count+sd_count, color = "macrophages"), width = 0.1, alpha = 0.1) +
  geom_line(data = dl_sum %>% filter(Cells %in% c("healthy")), aes(x = Time_h, y = mean_count, color = "healthy"), size = 0.5) +
  geom_errorbar(data = dl_sum %>% filter(Cells %in% c("healthy")), aes(x = Time_h, ymin= mean_count-sd_count, ymax = mean_count+sd_count, color = "healthy"), width = 0.1, alpha = 0.1) +
  geom_line(data = dl_sum %>% filter(Cells %in% c("necrotic")), aes(x = Time_h, y = mean_count, color = "necrotic"), size = 0.5) +
  geom_errorbar(data = dl_sum %>% filter(Cells %in% c("necrotic")), aes(x = Time_h, ymin= mean_count-sd_count, ymax = mean_count+sd_count, color = "necrotic"), width = 0.1, alpha = 0.1) +
  geom_line(data = dl_sum %>% filter(Cells %in% c("senescent")), aes(x = Time_h, y = mean_count, color = "senescent"), size = 0.5) +
  geom_errorbar(data = dl_sum %>% filter(Cells %in% c("senescent")), aes(x = Time_h, ymin= mean_count-sd_count, ymax = mean_count+sd_count, color = "senescent"), width = 0.1, alpha = 0.1) +
  labs(x = "Time (h)", y = "Cell count", color = "Cell Type") +
  scale_color_manual(values = colors,guide = guide_legend()) +
  theme_classic() +
  theme(legend.position="top", 
        legend.title = element_text(size=5.75, face = "bold"),
        legend.text = element_text(size=5.75),
        legend.margin=margin(t=0, r=0, b=-0.3, l=-0.3, unit="cm")) +
  guides(color = guide_legend(nrow=1,byrow=TRUE)) +
  facet_wrap(~exposure, ncol = 3,labeller = as_labeller(exp_names))
ggsave(paste0("Figures/","FigS1_cell_counts_long.pdf"), plot = gplot, width = 4, height = 4)

### ### ### ### ### 
###  Figure 2C  ### 
### ### ### ### ### 
# Calculate onset
onset <- dl %>% filter(Cells == "necrotic", Count > 0) %>% group_by(exposure,run) %>%
  summarise(onset = min(Time_h))

# Plot onset
ggplot(onset) +
  geom_line(aes(x = as.numeric(as.character(exposure)), y = onset), color = "grey") + 
  geom_point(aes(x = as.numeric(as.character(exposure)), y = onset, color = exposure), size = 3) + 
  xlab("APAP dose (a.u.)") + ylab("Onset necrosis (h)") +
  theme_classic() +
  scale_color_viridis_d() +
  scale_x_continuous(breaks = seq(250, 650, by = 50)) +
  theme(legend.position = "none")
ggsave(paste0("Figures/","Fig2C_necr_onset.pdf"), width = 3, height = 2)

### ### ### ### ### 
### Figure 2E.2 ### 
### ### ### ### ### 
prop <- data %>% filter(Time_h == 168, !exposure %in% c(0,650)) %>% 
  mutate(exposure = as.factor(exposure),
         senescent = senescent/600,
         healthy = healthy/600) %>%
  pivot_longer(cols = c(healthy, senescent), names_to = "Type", values_to = "Proportion") %>%
  group_by(exposure, time, Type) %>% 
  summarise(mean = mean(Proportion),
            sd = sd(Proportion))

ggplot() +
  geom_bar(data = prop, aes(x = exposure, y = mean, fill = Type), stat="identity", position = "dodge") +
  geom_errorbar(data = prop, aes(x = exposure, ymin = mean-sd, ymax = mean+sd, color = Type), position = "dodge") +
  scale_fill_viridis_d() + 
  scale_color_viridis_d() +
  theme_classic() + 
  labs(x = "APAP dose (a.u.)", y = "Proportion of cells") +
  #ggtitle("168 hours") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position=c(0.5,0.85),
        legend.background = element_rect(linetype="solid", 
                                         colour ="lightgrey"),
        legend.title = element_text(size=6, face = "bold"),
        legend.text = element_text(size=6),
        legend.key.size = unit(0.2,"cm"))
ggsave(paste0("Figures/","Fig2E2_Proportion.pdf"), width = 3, height = 2)

# ------------------------------------------------------------------------------------------------------------------------------------------

### ### ### ### ### 
###  Figure SX  ### 
### ### ### ### ### 

# Plot Recovery
dl_recovery <- dl %>% filter(Cells %in% c("healthy","senescent")) %>% group_by(Time_h, time, exposure, run) %>% 
  summarise(hepatocytes = sum(Count)) %>% ungroup() %>% group_by(Time_h, time, exposure) %>% 
  summarise(mean_count = mean(hepatocytes),
            sd_count = sd(hepatocytes))
ggplot(dl_recovery, aes(x=Time_h, y = mean_count, color = exposure)) +
  geom_line() + 
  geom_point(size = 1) + 
  geom_errorbar(aes(ymin= mean_count-sd_count, ymax = mean_count+sd_count), width = 0.1, alpha = 0.8) +
  scale_color_viridis_d(name = "Exposure") +
  xlab("Time (h)") + ylab("Cell count") +
  theme_classic()
ggsave(paste0("Figures/","FigSX_cell_counts_default.pdf"), width = 4, height = 3)

ggplot() +
  geom_point(data = (dl %>% filter(Cells %in% c("healthy","senescent")) %>% group_by(Time_h, time, exposure, run) %>%
                       summarise(total = sum(Count))),
             aes(x=Time_h, y = total),size = 0.2, color = "grey", alpha = 0.25) +
  geom_point(data = (dl %>% filter(Cells %in% c("healthy","senescent"))), 
             aes(x=Time_h, y = Count, color = Cells), size = 0.2) + 
  scale_color_viridis_d(name = "Exposure") +
  xlab("Time (h)") + ylab("Cell count") +
  theme_classic() +
  facet_grid(~exposure)
ggsave(paste0("Figures/","FigSX_cell_counts_default_perRun.pdf"), width = 11.5, height = 1.5)



### ### ### ### ### ### ### ### ###
### No GSH zonation Mouse Liver ###
### ### ### ### ### ### ### ### ###
#rm(list = ls())
codes <- data.frame(code = c("0001", "0002","0003","0004", "0005","0006","0007", "0008","0009"), 
                    exposure = c(250,300,350,400,450,500,550,600,650))

# get directory names
model_dir <- "10_Mouseliver_nozon_1-9_10x/"
simu_dirs <- dir(paste0(model_dir,"data/"))
file.name <- "log_cellcounts.csv"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")

# Get first data frame
gunzip(file.path, remove=FALSE)
data_nozon <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
data_nozon$exposure <- codes %>% 
  filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
  pull(exposure)
data_nozon$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))

# Bind all other data frames
for (d in simu_dirs[2:length(simu_dirs)]) {
  file.name <- "log_cellcounts.csv"
  file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
  if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
  data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
  data_tmp$exposure <- codes %>% 
    filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
    pull(exposure)
  data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
  data_nozon <- bind_rows(data_nozon,data_tmp)
}

# Rename columns
data_nozon <- data_nozon %>% dplyr::rename(necrotic = "necrotic_total",
                                           healthy = "hepatocytes_total",
                                           senescent = "senescent_total",
                                           Kupffer = "kupffer_total",
                                           macrophages = "macrophages_total",
                                           Time_h = "real_time") %>% 
  mutate(healthy = healthy - senescent)

# Pivot to long format
dl_nozon <- pivot_longer(data_nozon, c(necrotic,healthy,senescent,Kupffer,macrophages), names_to = "Cells", values_to = "Count") %>% 
  mutate(exposure = as.factor(exposure))

# Calculate mean and sd
dl_nozon_sum <- dl_nozon %>% group_by(time, Time_h, exposure, Cells) %>% 
  summarise(mean_count = mean(Count),
            sd_count = sd(Count)) %>% 
  ungroup()

# Plot Recovery
dl_nozon_recovery <- dl_nozon %>% filter(Cells %in% c("healthy","senescent")) %>% group_by(Time_h, time, exposure, run) %>% 
  summarise(hepatocytes = sum(Count)) %>% ungroup() %>% group_by(Time_h, time, exposure) %>% 
  summarise(mean_count = mean(hepatocytes),
            sd_count = sd(hepatocytes))
ggplot(dl_nozon_recovery, aes(x=Time_h, y = mean_count, color = exposure)) +
  geom_line() + 
  geom_point(size = 1) + 
  geom_errorbar(aes(ymin= mean_count-sd_count, ymax = mean_count+sd_count), width = 0.1, alpha = 0.8) +
  scale_color_viridis_d(name = "Exposure") +
  xlab("Time (h)") + ylab("Cell count") +
  theme_classic()
ggsave(paste0("Figures/","FigSX_cell_counts_NoZonation.pdf"), width = 4, height = 3)


### ### ### ### ### ### ### ### ###
### No GSH zonation Mouse Liver ###
### ### ### ### ### ### ### ### ###
#rm(list = ls())

# get directory names
model_dir <- "09_Mouseliver_nozonGSH_1-9_10x/"
simu_dirs <- dir(paste0(model_dir,"data/"))
file.name <- "log_cellcounts.csv"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")

# Get first data frame
gunzip(file.path, remove=FALSE)
data_nozongsh <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
data_nozongsh$exposure <- codes %>% 
  filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
  pull(exposure)
data_nozongsh$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))

# Bind all other data frames
for (d in simu_dirs[2:length(simu_dirs)]) {
  file.name <- "log_cellcounts.csv"
  file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
  if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
  data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
  data_tmp$exposure <- codes %>% 
    filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
    pull(exposure)
  data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
  data_nozongsh <- bind_rows(data_nozongsh,data_tmp)
}

# Rename columns
data_nozongsh <- data_nozongsh %>% dplyr::rename(necrotic = "necrotic_total",
                                                 healthy = "hepatocytes_total",
                                                 senescent = "senescent_total",
                                                 Kupffer = "kupffer_total",
                                                 macrophages = "macrophages_total",
                                                 Time_h = "real_time") %>% 
  mutate(healthy = healthy - senescent)

# Pivot to long format
dl_nozongsh <- pivot_longer(data_nozongsh, c(necrotic,healthy,senescent,Kupffer,macrophages), names_to = "Cells", values_to = "Count") %>% 
  mutate(exposure = as.factor(exposure))

# Calculate mean and sd
dl_nozongsh_sum <- dl_nozongsh %>% group_by(time, Time_h, exposure, Cells) %>% summarise(mean_count = mean(Count),
                                                                                         sd_count = sd(Count)) %>% 
  ungroup()
# Plot
# ggplot(dl_sum %>% filter(!Cells %in% c("Kupffer", "macrophages")), aes(x= Time_h, y = mean_count, color = Cells)) + 
#   geom_line(size = 1) +
#   geom_errorbar(aes(ymin= mean_count-sd_count, ymax = mean_count+sd_count), width = 0.1, alpha = 0.1) +
#   xlab("Time (h)") + ylab("Cell count") +
#   theme_classic() +
#   scale_color_viridis_d(name = "Hepatocytes\nstate") +
#   facet_grid(rows = vars(exposure))
# ggsave(paste0("Figures/","Fig2X_cell_counts.pdf"), width = 3.5, height = 5.5)

# # Plot Cell Death
# ggplot(dl_sum %>% filter(!Cells %in% c("Kupffer", "macrophages"),
#                          Time_h %in% c(48))) +
#   geom_bar(aes(x=Cells, y = mean_count, fill = exposure), stat = "identity", position = "dodge2", colour = "black") + 
#   scale_color_viridis_d(name = "Exposure", aesthetics = "fill") +
#   xlab("Hepatocyte state") + ylab("Cell count") +
#   ggtitle("48 h") +
#   theme_classic() + 
#   theme(plot.title = element_text(hjust = 0.5))
# ggsave(paste0("Figures/","Fig2C_cell_counts_NoZonationP450_48h.pdf"), width = 4, height = 2.5)


# Plot Recovery
dl_nozongsh_recovery <- dl_nozongsh %>% filter(Cells %in% c("healthy","senescent")) %>% group_by(Time_h, time, exposure, run) %>% 
  summarise(hepatocytes = sum(Count)) %>% ungroup() %>% group_by(Time_h, time, exposure) %>% 
  summarise(mean_count = mean(hepatocytes),
            sd_count = sd(hepatocytes))
ggplot(dl_nozongsh_recovery, aes(x=Time_h, y = mean_count, color = exposure)) +
  geom_line() + 
  geom_point(size = 1) + 
  geom_errorbar(aes(ymin= mean_count-sd_count, ymax = mean_count+sd_count), width = 0.1, alpha = 0.8) +
  scale_color_viridis_d(name = "Exposure") +
  xlab("Time (h)") + ylab("Cell count") +
  theme_classic()
ggsave(paste0("Figures/","FigSX_cell_counts_NoZonationGSH.pdf"), width = 4, height = 3)

# # Plot Recovery
# dl_recovery <- dl_sum %>% filter(Cells %in% c("healthy","senescent")) %>% group_by(Time_h, time, exposure) %>% 
#   summarise(hepatocytes = sum(mean_count))
# ggplot(dl_recovery %>% filter(Time_h > -1 & Time_h < 100)) +
#   geom_line(aes(x=Time_h, y = hepatocytes, color = exposure)) + 
#   geom_point(aes(x=Time_h, y = hepatocytes, color = exposure), size = 1) + 
#   scale_color_viridis_d(name = "Exposure") +
#   xlab("Time (h)") + ylab("Cell count") +
#   theme_classic()
# ggsave(paste0("Figures/","Fig3C_cell_counts_NoZonationGSH.pdf"), width = 4, height = 3)

### ### ### ### ### ### ### ### ###
###No P450 zonation Mouse Liver ###
### ### ### ### ### ### ### ### ###
#rm(list = ls())

# get directory names
model_dir <- "11_Mouseliver_nozonP450_1-9_10x/"
simu_dirs <- dir(paste0(model_dir,"data/"))
file.name <- "log_cellcounts.csv"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")

# Get first data frame
gunzip(file.path, remove=FALSE)
data_nozonp450 <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
data_nozonp450$exposure <- codes %>% 
  filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
  pull(exposure)
data_nozonp450$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))

# Bind all other data frames
for (d in simu_dirs[2:length(simu_dirs)]) {
  file.name <- "log_cellcounts.csv"
  file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
  if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
  data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
  data_tmp$exposure <- codes %>% 
    filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
    pull(exposure)
  data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
  data_nozonp450 <- bind_rows(data_nozonp450,data_tmp)
}

# Rename columns
data_nozonp450 <- data_nozonp450 %>% dplyr::rename(necrotic = "necrotic_total",
                                                   healthy = "hepatocytes_total",
                                                   senescent = "senescent_total",
                                                   Kupffer = "kupffer_total",
                                                   macrophages = "macrophages_total",
                                                   Time_h = "real_time") %>% 
  mutate(healthy = healthy - senescent)

# Pivot to long format
dl_nozonp450 <- pivot_longer(data_nozonp450, c(necrotic,healthy,senescent,Kupffer,macrophages), names_to = "Cells", values_to = "Count") %>% 
  mutate(exposure = as.factor(exposure))

# Calculate mean and sd
dl_nozonp450_sum <- dl_nozonp450 %>% group_by(time, Time_h, exposure, Cells) %>% summarise(mean_count = mean(Count),
                                                                                           sd_count = sd(Count)) %>% 
  ungroup()
# Plot
# ggplot(dl_sum %>% filter(!Cells %in% c("Kupffer", "macrophages")), aes(x= Time_h, y = mean_count, color = Cells)) + 
#   geom_line(size = 1) +
#   geom_errorbar(aes(ymin= mean_count-sd_count, ymax = mean_count+sd_count), width = 0.1, alpha = 0.1) +
#   xlab("Time (h)") + ylab("Cell count") +
#   theme_classic() +
#   scale_color_viridis_d(name = "Hepatocytes\nstate") +
#   facet_grid(rows = vars(exposure))
# ggsave(paste0("Figures/","Fig2X_cell_counts.pdf"), width = 3.5, height = 5.5)

# # Plot Cell Death
# ggplot(dl_sum %>% filter(!Cells %in% c("Kupffer", "macrophages"),
#                          Time_h %in% c(48))) +
#   geom_bar(aes(x=Cells, y = mean_count, fill = exposure), stat = "identity", position = "dodge2", colour = "black") + 
#   scale_color_viridis_d(name = "Exposure", aesthetics = "fill") +
#   xlab("Hepatocyte state") + ylab("Cell count") +
#   ggtitle("48 h") +
#   theme_classic() + 
#   theme(plot.title = element_text(hjust = 0.5))
# ggsave(paste0("Figures/","Fig2C_cell_counts_NoZonationP450_48h.pdf"), width = 4, height = 2.5)

# Plot Recovery
dl_nozonp450_recovery <- dl_nozonp450 %>% filter(Cells %in% c("healthy","senescent")) %>% group_by(Time_h, time, exposure, run) %>% 
  summarise(hepatocytes = sum(Count)) %>% ungroup() %>% group_by(Time_h, time, exposure) %>% 
  summarise(mean_count = mean(hepatocytes),
            sd_count = sd(hepatocytes))
ggplot(dl_nozonp450_recovery, aes(x=Time_h, y = mean_count, color = exposure)) +
  geom_line() + 
  geom_point(size = 1) + 
  geom_errorbar(aes(ymin= mean_count-sd_count, ymax = mean_count+sd_count), width = 0.1, alpha = 0.8) +
  scale_color_viridis_d(name = "Exposure") +
  xlab("Time (h)") + ylab("Cell count") +
  theme_classic()
ggsave(paste0("Figures/","FigSX_cell_counts_NoZonationP450.pdf"), width = 4, height = 3)

# # Plot Recovery
# dl_recovery <- dl_sum %>% filter(Cells %in% c("healthy","senescent")) %>% group_by(Time_h, time, exposure) %>% 
#   summarise(hepatocytes = sum(mean_count))
# ggplot(dl_recovery %>% filter(Time_h > -1 & Time_h < 100)) +
#   geom_line(aes(x=Time_h, y = hepatocytes, color = exposure)) + 
#   geom_point(aes(x=Time_h, y = hepatocytes, color = exposure), size = 1) + 
#   scale_color_viridis_d(name = "Exposure") +
#   xlab("Time (h)") + ylab("Cell count") +
#   theme_classic()
# ggsave(paste0("Figures/","Fig3D_cell_counts_NoZonationP450.pdf"), width = 4, height = 3)
# 


### ### ### ### ### 
###  Figure 3B  ### 
### ### ### ### ### 
colors <- c(viridisLite::viridis(4))
names(colors) <- c("Both","None", "No GSH","No P450")

pchs <- c(15,16,17,18)
names(pchs) <- c("Both","None", "No GSH","No P450")

ggplot() +
  geom_point(data = dl_sum %>% filter(Cells == c("necrotic"), Time_h == 14), 
             aes(x= as.numeric(as.character(exposure)), y = mean_count, color = "Both", shape = "Both"), size = 2) +
  geom_errorbar(data = dl_sum %>% filter(Cells == c("necrotic"), Time_h == 14), 
                aes(x= as.numeric(as.character(exposure)), ymin = mean_count-sd_count, ymax =  mean_count+sd_count, color = "Both"), width = 10) +
  geom_line(data = dl_sum %>% filter(Cells == c("necrotic"), Time_h == 14), 
            aes(x= as.numeric(as.character(exposure)), y = mean_count, color = "Both")) +
  
  geom_point(data = dl_nozon_sum %>% filter(Cells == c("necrotic"), Time_h == 14), 
             aes(x= as.numeric(as.character(exposure)), y = mean_count, color = "None", shape = "None"), size = 2) +
  geom_errorbar(data = dl_nozon_sum %>% filter(Cells == c("necrotic"), Time_h == 14), 
                aes(x= as.numeric(as.character(exposure)), ymin = mean_count-sd_count, ymax =  mean_count+sd_count, color = "None"), width = 10) +
  geom_line(data = dl_nozon_sum %>% filter(Cells == c("necrotic"), Time_h == 14), 
            aes(x= as.numeric(as.character(exposure)), y = mean_count, color = "None")) +
  
  geom_point(data = dl_nozongsh_sum %>% filter(Cells == c("necrotic"), Time_h == 14), 
             aes(x= as.numeric(as.character(exposure)), y = mean_count, color = "No GSH", shape = "No GSH"), size = 2) +
  geom_errorbar(data = dl_nozongsh_sum %>% filter(Cells == c("necrotic"), Time_h == 14), 
                aes(x= as.numeric(as.character(exposure)), ymin = mean_count-sd_count, ymax =  mean_count+sd_count, color = "No GSH"), width = 10) +
  geom_line(data = dl_nozongsh_sum %>% filter(Cells == c("necrotic"), Time_h == 14), 
            aes(x= as.numeric(as.character(exposure)), y = mean_count, color = "No GSH")) +
  
  geom_point(data = dl_nozonp450_sum %>% filter(Cells == c("necrotic"), Time_h == 14),
             aes(x= as.numeric(as.character(exposure)), y = mean_count, color = "No P450", shape = "No P450"), size = 2) +
  geom_errorbar(data = dl_nozonp450_sum %>% filter(Cells == c("necrotic"), Time_h == 14),
                aes(x= as.numeric(as.character(exposure)), ymin = mean_count-sd_count, ymax =  mean_count+sd_count, color = "No P450"), width = 10) +
  geom_line(data = dl_nozonp450_sum %>% filter(Cells == c("necrotic"), Time_h == 14),
            aes(x= as.numeric(as.character(exposure)), y = mean_count, color = "No P450")) +
  
  scale_color_manual(values = colors) +
  scale_shape_manual(values = pchs) +
  theme_classic() + 
  labs(x = "APAP dose (a.u.)", y = "Number of cells", color = "Zonation", shape = "Zonation") +
  ggtitle("Necrotic cells at 14 hours") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("Figures/","Fig3B_Zonation_sensitivity_damage.pdf"), width = 3, height = 2)


### ### ### ### ### 
###  Figure 3C  ### 
### ### ### ### ### 
colors <- c(viridisLite::viridis(4))
names(colors) <- c("Both","None", "No GSH","No P450")

pchs <- c(15,16,17,18)
names(pchs) <- c("Both","None", "No GSH","No P450")

ggplot() +
  geom_point(data = dl_sum %>% filter(Cells == c("senescent"), Time_h == 48), 
             aes(x= as.numeric(as.character(exposure)), y = mean_count, color = "Both", shape = "Both"), size = 2) +
  geom_errorbar(data = dl_sum %>% filter(Cells == c("senescent"), Time_h == 48), 
                aes(x= as.numeric(as.character(exposure)), ymin = mean_count-sd_count, ymax =  mean_count+sd_count, color = "Both"), width = 10) +
  geom_line(data = dl_sum %>% filter(Cells == c("senescent"), Time_h == 48), 
            aes(x= as.numeric(as.character(exposure)), y = mean_count, color = "Both")) +
  
  geom_point(data = dl_nozon_sum %>% filter(Cells == c("senescent"), Time_h == 48), 
             aes(x= as.numeric(as.character(exposure)), y = mean_count, color = "None", shape = "None"), size = 2) +
  geom_errorbar(data = dl_nozon_sum %>% filter(Cells == c("senescent"), Time_h == 48), 
                aes(x= as.numeric(as.character(exposure)), ymin = mean_count-sd_count, ymax =  mean_count+sd_count, color = "None"), width = 10) +
  geom_line(data = dl_nozon_sum %>% filter(Cells == c("senescent"), Time_h == 48), 
            aes(x= as.numeric(as.character(exposure)), y = mean_count, color = "None")) +
  
  geom_point(data = dl_nozongsh_sum %>% filter(Cells == c("senescent"), Time_h == 48), 
             aes(x= as.numeric(as.character(exposure)), y = mean_count, color = "No GSH", shape = "No GSH"), size = 2) +
  geom_errorbar(data = dl_nozongsh_sum %>% filter(Cells == c("senescent"), Time_h == 48), 
                aes(x= as.numeric(as.character(exposure)), ymin = mean_count-sd_count, ymax =  mean_count+sd_count, color = "No GSH"), width = 10) +
  geom_line(data = dl_nozongsh_sum %>% filter(Cells == c("senescent"), Time_h == 48), 
            aes(x= as.numeric(as.character(exposure)), y = mean_count, color = "No GSH")) +
  
  geom_point(data = dl_nozonp450_sum %>% filter(Cells == c("senescent"), Time_h == 48),
             aes(x= as.numeric(as.character(exposure)), y = mean_count, color = "No P450", shape = "No P450"), size = 2) +
  geom_errorbar(data = dl_nozonp450_sum %>% filter(Cells == c("senescent"), Time_h == 48),
                aes(x= as.numeric(as.character(exposure)), ymin = mean_count-sd_count, ymax =  mean_count+sd_count, color = "No P450"), width = 10) +
  geom_line(data = dl_nozonp450_sum %>% filter(Cells == c("senescent"), Time_h == 48),
            aes(x= as.numeric(as.character(exposure)), y = mean_count, color = "No P450")) +
  
  scale_color_manual(values = colors) +
  scale_shape_manual(values = pchs) +
  theme_classic() + 
  labs(x = "APAP dose (a.u.)", y = "Number of cells", color = "Zonation", shape = "Zonation") +
  ggtitle("Senescent cells at 48 hours") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("Figures/","Fig3C_Zonation_sensitivity_senescence.pdf"), width = 3, height = 2)

### ### ### ### ### 
###  Figure 3D  ### 
### ### ### ### ### 
colors <- c(viridisLite::viridis(4))
names(colors) <- c("Both","None", "No GSH","No P450")

pchs <- c(15,16,17,18)
names(pchs) <- c("Both","None", "No GSH","No P450")

ggplot() +
  geom_point(data = dl_sum %>% filter(Cells == c("healthy"), Time_h == 168), 
             aes(x= as.numeric(as.character(exposure)), y = mean_count, color = "Both", shape = "Both"), size = 2) +
  geom_errorbar(data = dl_sum %>% filter(Cells == c("healthy"), Time_h == 168), 
                aes(x= as.numeric(as.character(exposure)), ymin = mean_count-sd_count, ymax =  mean_count+sd_count, color = "Both"), width = 10) +
  geom_line(data = dl_sum %>% filter(Cells == c("healthy"), Time_h == 168), 
            aes(x= as.numeric(as.character(exposure)), y = mean_count, color = "Both")) +
  
  geom_point(data = dl_nozon_sum %>% filter(Cells == c("healthy"), Time_h == 168), 
             aes(x= as.numeric(as.character(exposure)), y = mean_count, color = "None", shape = "None"), size = 2) +
  geom_errorbar(data = dl_nozon_sum %>% filter(Cells == c("healthy"), Time_h == 168), 
                aes(x= as.numeric(as.character(exposure)), ymin = mean_count-sd_count, ymax =  mean_count+sd_count, color = "None"), width = 10) +
  geom_line(data = dl_nozon_sum %>% filter(Cells == c("healthy"), Time_h == 168), 
            aes(x= as.numeric(as.character(exposure)), y = mean_count, color = "None")) +
  
  geom_point(data = dl_nozongsh_sum %>% filter(Cells == c("healthy"), Time_h == 168), 
             aes(x= as.numeric(as.character(exposure)), y = mean_count, color = "No GSH", shape = "No GSH"), size = 2) +
  geom_errorbar(data = dl_nozongsh_sum %>% filter(Cells == c("healthy"), Time_h == 168), 
                aes(x= as.numeric(as.character(exposure)), ymin = mean_count-sd_count, ymax =  mean_count+sd_count, color = "No GSH"), width = 10) +
  geom_line(data = dl_nozongsh_sum %>% filter(Cells == c("healthy"), Time_h == 168), 
            aes(x= as.numeric(as.character(exposure)), y = mean_count, color = "No GSH")) +
  
  geom_point(data = dl_nozonp450_sum %>% filter(Cells == c("healthy"), Time_h == 168),
             aes(x= as.numeric(as.character(exposure)), y = mean_count, color = "No P450", shape = "No P450"), size = 2) +
  geom_errorbar(data = dl_nozonp450_sum %>% filter(Cells == c("healthy"), Time_h == 168),
                aes(x= as.numeric(as.character(exposure)), ymin = mean_count-sd_count, ymax =  mean_count+sd_count, color = "No P450"), width = 10) +
  geom_line(data = dl_nozonp450_sum %>% filter(Cells == c("healthy"), Time_h == 168),
            aes(x= as.numeric(as.character(exposure)), y = mean_count, color = "No P450")) +
  
  scale_color_manual(values = colors) +
  scale_shape_manual(values = pchs) +
  theme_classic() + 
  labs(x = "APAP dose (a.u.)", y = "Number of cells", color = "Zonation", shape = "Zonation") +
  ggtitle("Healthy cells at 168 hours") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("Figures/","Fig3D_Zonation_sensitivity_recovery.pdf"), width = 3, height = 2)


# ------------------------------------------------------------------------------------------------------------------------------------------
### ### ### ### ### 
###  Figure 4A  ### 
### ### ### ### ### 
model_dir <- "08_Mouseliver_default_1-10_10x/"
simu_dirs <- dir(paste0(model_dir,"data/"))

# Define folder code
codes <- data.frame(code = c("0001", "0002","0003","0004", "0005","0006","0007", "0008","0009","0010"),
                    exposure = c(0,250,300,350,400,450,500,550,600,650))

# Intracellular concentrations
file.name <- "log_SingleCellConc.csv"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")
if (!file.exists(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))) {gunzip(file.path, remove=FALSE)}

# Bind all other data frames
if (interactive()) {
  if (askYesNo("Do you want to load ALL data files?")) {
    print("OK!")
    
    # Get first data frame
    data <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
    data$exposure <- codes %>% 
      filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
      pull(exposure)
    data$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))
    
    for (d in simu_dirs[2:length(simu_dirs)]) {
      file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
      if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
      data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
      data_tmp$exposure <- codes %>% 
        filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
        pull(exposure)
      data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
      data <- bind_rows(data,data_tmp)}
    
    if (interactive()) {
      if (askYesNo("Do you want to save these data? This may overwrite existing data.")) {
        print("OK, saving...")
        save(data, file = paste0(model_dir,"data_log_SingleCellConc.RData"))
        print("Saved!")}}
  } else {
    if (interactive()) {
      if (askYesNo("Do you want to load these data from file?")) {
        load(paste0(model_dir,"data_log_SingleCellConc.RData"))}}}
}

# summarise
data_sum <- data %>% filter(celltype %in% c(1,6)) %>% 
  mutate(exposure = as.factor(exposure)) %>%
  group_by(exposure,run, time, real_time) %>% 
  mutate(p53_total = p53 + p53p) %>%
  summarize(m_p53_total = mean(p53_total),
            m_p53p = mean(p53p),
            m_p53 = mean(p53),
            m_p21 = mean(p21),
            m_mdm2 = mean(MDM2)) %>% ungroup %>%
  group_by(exposure, time, real_time) %>%
  summarize(mean_p53_total = mean(m_p53_total),
            sd_p53_total = sd(m_p53_total),
            mean_p53p = mean(m_p53p),
            sd_p53p = sd(m_p53p),
            mean_p53 = mean(m_p53),
            sd_p53 = sd(m_p53),
            mean_p21 = mean(m_p21),
            sd_p21 = sd(m_p21),
            mean_mdm2 = mean(m_mdm2),
            sd_mdm2 = sd(m_mdm2))


ggplot(data_sum %>% filter(real_time>-1, !exposure %in% c(0, 650))) +
  geom_line(aes(x = real_time, y= mean_p53p, color = exposure)) +
  geom_errorbar(aes(x = real_time, ymin = mean_p53p-sd_p53p, ymax = mean_p53p+sd_p53p, color = exposure)) +
  scale_color_viridis_d(name = "APAP\ndose (a.u.)") +
  xlab("Time (h)") + ylab("Concentration (a.u.)") +
  theme_classic() +
  ggtitle("phospho-p53") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("Figures/","Fig4A_p53p_default.pdf"), width = 4, height = 3)

ggplot(data_sum %>% filter(real_time>-1, !exposure %in% c(0, 650))) +
  geom_line(aes(x = real_time, y= mean_p53, color = exposure)) +
  geom_errorbar(aes(x = real_time, ymin = mean_p53-sd_p53, ymax = mean_p53+sd_p53, color = exposure)) +
  scale_color_viridis_d(name = "APAP\ndose (a.u.)") +
  xlab("Time (h)") + ylab("Concentration (a.u.)") +
  theme_classic() +
  ggtitle("p53") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("Figures/","FigS2A_p53_default.pdf"), width = 4, height = 3)

ggplot(data_sum %>% filter(real_time>-1, !exposure %in% c(0, 650))) +
  geom_line(aes(x = real_time, y= mean_p53_total, color = exposure)) +
  geom_errorbar(aes(x = real_time, ymin = mean_p53_total-sd_p53_total, ymax = mean_p53_total+sd_p53_total, color = exposure)) +
  scale_color_viridis_d(name = "APAP\ndose (a.u.)") +
  xlab("Time (h)") + ylab("Concentration (a.u.)") +
  theme_classic() +
  ggtitle("Total p53") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("Figures/","FigS2B_p53_total_default.pdf"), width = 4, height = 3)

ggplot(data_sum %>% filter(real_time>-1, !exposure %in% c(0, 650))) +
  geom_line(aes(x = real_time, y= mean_p21, color = exposure)) +
  geom_errorbar(aes(x = real_time, ymin = mean_p21-sd_p21, ymax = mean_p21+sd_p21, color = exposure)) +
  scale_color_viridis_d(name = "APAP\ndose (a.u.)") +
  xlab("Time (h)") + ylab("Concentration (a.u.)") +
  theme_classic() +
  ggtitle("p21") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("Figures/","Fig4A_p21_default.pdf"), width = 4, height = 3)

ggplot(data_sum %>% filter(real_time>-1, !exposure %in% c(0, 650))) +
  geom_line(aes(x = real_time, y= mean_mdm2, color = exposure)) +
  geom_errorbar(aes(x = real_time, ymin = mean_mdm2-sd_mdm2, ymax = mean_mdm2+sd_mdm2, color = exposure)) +
  scale_color_viridis_d(name = "APAP\ndose (a.u.)") +
  xlab("Time (h)") + ylab("Concentration (a.u.)") +
  theme_classic() +
  ggtitle("MDM2") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("Figures/","Fig4A_mdm2_default.pdf"), width = 4, height = 3)

### ### ### ### ### 
###  Figure 4B  ### 
### ### ### ### ### 

model_dir <- "13_Mouseliver_Mdm2_1-5_10x/"
simu_dirs <- dir(paste0(model_dir,"data/"))

# Define folder code
codes <- data.frame(code = c("0001", "0002","0003","0004", "0005"),
                    r = c(0.01, 0.1, 1, 10, 100))

# Intracellular concentrations
file.name <- "log_SingleCellConc.csv"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")
if (!file.exists(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))) {gunzip(file.path, remove=FALSE)}

# Bind all other data frames
if (interactive()) {
  if (askYesNo("Do you want to load ALL data files?")) {
    print("OK!")
    
    # Get first data frame
    data <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
    data$r <- codes %>% 
      filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
      pull(r)
    data$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))
    
    for (d in simu_dirs[2:length(simu_dirs)]) {
      file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
      if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
      data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
      data_tmp$r <- codes %>% 
        filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
        pull(r)
      data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
      data <- bind_rows(data,data_tmp)}
    
    if (interactive()) {
      if (askYesNo("Do you want to save these data? This may overwrite existing data.")) {
        print("OK, saving...")
        save(data, file = paste0(model_dir,"data_log_SingleCellConc.RData"))
        print("Saved!")}}
  } else {
    if (interactive()) {
      if (askYesNo("Do you want to load these data from file?")) {
        load(paste0(model_dir,"data_log_SingleCellConc.RData"))}}}
}

# summarise
data_sum <- data %>% filter(celltype %in% c(1,6)) %>% 
  mutate(r = as.factor(r)) %>%
  group_by(r,run, time, real_time) %>% 
  mutate(p53_total = p53 + p53p) %>%
  summarize(m_p53_total = mean(p53_total),
            m_p53p = mean(p53p),
            m_p53 = mean(p53),
            m_p21 = mean(p21),
            m_mdm2 = mean(MDM2)) %>% ungroup %>%
  group_by(r, time, real_time) %>%
  summarize(mean_p53_total = mean(m_p53_total),
            sd_p53_total = sd(m_p53_total),
            mean_p53p = mean(m_p53p),
            sd_p53p = sd(m_p53p),
            mean_p53 = mean(m_p53),
            sd_p53 = sd(m_p53),
            mean_p21 = mean(m_p21),
            sd_p21 = sd(m_p21),
            mean_mdm2 = mean(m_mdm2),
            sd_mdm2 = sd(m_mdm2))

# # Intracellular concentrations
# file.name <- "DDR.csv"
# file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")
# if (!file.exists(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))) {gunzip(file.path, remove=FALSE)}
# 
# # Get first data frame
# data <- read_csv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
# data$r <- codes %>% 
#   filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
#   pull(r)
# data$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))
# 
# # Bind all other data frames
# for (d in simu_dirs[2:length(simu_dirs)]) {
#   file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
#   if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
#   data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
#   data_tmp$r <- codes %>% 
#     filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
#     pull(r)
#   data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
#   data <- bind_rows(data,data_tmp)
# }
# 
# data_sum <- data %>% mutate(r = as.factor(r)) %>%
#   group_by(time, real_time, r) %>%
#   summarise(mean_p21 = mean(p21_mean),
#             mean_p53 = mean(p53_mean),
#             mean_mdm2 = mean(mdm2_mean),
#             sd_p21 = sd(p21_mean),
#             sd_p53 = sd(p53_mean),
#             sd_mdm2 = sd(mdm2_mean)) 

ggplot(data_sum %>% filter(real_time > -1)) +
  geom_point(aes(x = real_time, y = mean_p21, color = r), pch = 20) +
  geom_errorbar(aes(x = real_time, ymin = mean_p21-sd_p21, ymax = mean_p21+sd_p21, color = r), alpha = 0.4) +
  scale_color_viridis_d(name = "Factor r") +
  xlab("Time (h)") + ylab("Concentration (a.u.)") +
  theme_classic() +
  ggtitle("p21") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("Figures/","Fig4B_p21_feedback_mdm2.pdf"), width = 4, height = 3)

ggplot(data_sum %>% filter(real_time > -1)) +
  geom_point(aes(x = real_time, y = mean_p53p, color = r), pch = 20) +
  geom_errorbar(aes(x = real_time, ymin = mean_p53p-sd_p53p, ymax = mean_p53p+sd_p53p, color = r), alpha = 0.4) +
  scale_color_viridis_d(name = "Factor r") +
  xlab("Time (h)") + ylab("Concentration (a.u.)") +
  theme_classic() +
  ggtitle("phospho-p53") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("Figures/","Fig4B_p53p_feedback_mdm2.pdf"), width = 4, height = 3)

ggplot(data_sum %>% filter(real_time > -1)) +
  geom_point(aes(x = real_time, y = mean_p53, color = r), pch = 20) +
  geom_errorbar(aes(x = real_time, ymin = mean_p53-sd_p53, ymax = mean_p53+sd_p53, color = r), alpha = 0.4) +
  scale_color_viridis_d(name = "Factor r") +
  xlab("Time (h)") + ylab("Concentration (a.u.)") +
  theme_classic() +
  ggtitle("p53") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("Figures/","FigS2C_p53_feedback_mdm2.pdf"), width = 4, height = 3)

ggplot(data_sum %>% filter(real_time > -1)) +
  geom_point(aes(x = real_time, y = mean_p53_total, color = r), pch = 20) +
  geom_errorbar(aes(x = real_time, ymin = mean_p53_total-sd_p53_total, ymax = mean_p53_total+sd_p53_total, color = r), alpha = 0.4) +
  scale_color_viridis_d(name = "Factor r") +
  xlab("Time (h)") + ylab("Concentration (a.u.)") +
  theme_classic() +
  ggtitle("Total p53") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("Figures/","FigS2D_p53_total_feedback_mdm2.pdf"), width = 4, height = 3)

### ### ### ### ### 
###  Figure 4C  ### 
### ### ### ### ### 

file.name <- "log_cellcounts.csv"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")
if (!file.exists(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))) {gunzip(file.path, remove=FALSE)}

# Get first data frame
data <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
data$r <- codes %>% 
  filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
  pull(r)
data$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))

# Bind all other data frames
for (d in simu_dirs[2:length(simu_dirs)]) {
  file.name <- "log_cellcounts.csv"
  file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
  if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
  data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
  data_tmp$r <- codes %>% 
    filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
    pull(r)
  data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
  data <- bind_rows(data,data_tmp)
}

# Rename columns
data <- data %>% dplyr::rename(necrotic = "necrotic_total",
                               healthy = "hepatocytes_total",
                               senescent = "senescent_total",
                               Kupffer = "kupffer_total",
                               macrophages = "macrophages_total",
                               Time_h = "real_time") %>% 
  mutate(healthy = healthy - senescent)

data_sum <- data %>% select(-c(run)) %>% group_by(Time_h, time, r) %>%
  summarise_all(mean)
ggplot() +
  geom_point(data = data_sum, aes(x = Time_h, y = senescent, color = r))

# Pivot to long format
dl <- pivot_longer(data, c(necrotic,healthy,senescent,Kupffer,macrophages), names_to = "Cells", values_to = "Count") %>% 
  mutate(r = as.factor(r))

# Calculate mean and sd
dl_sum <- dl %>% group_by(time, Time_h, r, Cells) %>% summarise(mean_count = mean(Count),
                                                                sd_count = sd(Count)) %>% 
  ungroup()


ggplot() +
  geom_line(data = dl_sum %>% filter(Cells == c("senescent"), Time_h == 48), 
            aes(x= as.numeric(as.character(r)), y = mean_count), color = "grey") +
  geom_point(data = dl_sum %>% filter(Cells == c("senescent"), Time_h == 48), 
             aes(x= as.numeric(as.character(r)), y = mean_count, color = r), pch = 15, size = 2) +
  geom_errorbar(data = dl_sum %>% filter(Cells == c("senescent"), Time_h == 48), 
                aes(x= as.numeric(as.character(r)), ymin = mean_count-sd_count, ymax =  mean_count+sd_count, color = r), width = 1) +
  scale_x_continuous(trans='log2',breaks = c(0.01,0.1,1,10,100)) +
  scale_color_viridis_d(name = "factor r") +
  theme_classic() + 
  labs(x = "Factor r", y = "Number of cells") +
  ggtitle("Senescent cells at 48 hours") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")
ggsave(paste0("Figures/","Fig4C_MDM2_sensitivity_senescence.pdf"), width = 3, height = 2)

### ### ### ### ### 
###  Figure 4D  ### 
### ### ### ### ### 

model_dir <- "13_Mouseliver_Mdm2_1-5_10x/"
simu_dirs <- dir(paste0(model_dir,"data/"))

# Define folder code
codes <- data.frame(code = c("0001", "0002","0003","0004", "0005"),
                    r = c(0.01, 0.1, 1, 10, 100))

file.name <- "cell_division_hepatocytes.txt"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")
if (!file.exists(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))) {gunzip(file.path, remove=FALSE)}

# Get first data frame
data <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
data$r <- codes %>% 
  filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
  pull(r)
data <- data %>% mutate(run = sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1]))),
                        step = seq(1,nrow(data)),
                        Time_h = (Time - 96)/12)

# Bind all other data frames
for (d in simu_dirs[2:length(simu_dirs)]) {
  file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
  if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
  data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
  data_tmp$r <- codes %>% 
    filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
    pull(r)
  data_tmp <- data_tmp %>% mutate(run = sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d))),
                                  step = seq(1,nrow(data_tmp)),
                                  Time_h = (Time - 96)/12)
  data <- bind_rows(data,data_tmp)
}

labels <- c(
  "0.01" = "r = 0.01",
  "0.1"="r = 0.1",
  "1"="r = 1",
  "10"="r = 10",
  "100" = "r = 100")
ggplot() +
  geom_step(data = data, aes(x = Time_h, y = step, color = run)) +
  theme_classic() +
  labs(x = "Time (h)", y = "Proliferation") +
  scale_color_viridis_d() +
  theme(legend.position = "none") +
  facet_grid(~r, labeller = as_labeller(labels))
ggsave(paste0("Figures/","Fig4D_MDM2_sensitivity_proliferation.pdf"), width = 10, height = 2)

### ### ### ### ### ### 
###  Figure 5B left ### 
### ### ### ### ### ###  

# get directory names
model_dir <- "08_Mouseliver_default_1-10_10x/"
simu_dirs <- dir(paste0(model_dir,"data/"))

# Define folder code
codes <- data.frame(code = c("0001", "0002","0003","0004", "0005","0006","0007", "0008","0009","0010"),
                    exposure = c(0,250,300,350,400,450,500,550,600,650))


# Cell counts
file.name <- "log_cellcounts.csv"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")
if (!file.exists(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))) {gunzip(file.path, remove=FALSE)}

# Get first data frame
data <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
data$exposure <- codes %>% 
  filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
  pull(exposure)
data$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))

# Bind all other data frames
for (d in simu_dirs[2:length(simu_dirs)]) {
  file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
  if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
  data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
  data_tmp$exposure <- codes %>% 
    filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
    pull(exposure)
  data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
  data <- bind_rows(data,data_tmp)
}

# Rename columns
data <- data %>% dplyr::rename(necrotic = "necrotic_total",
                               total = "hepatocytes_total",
                               senescent = "senescent_total",
                               Kupffer = "kupffer_total",
                               macrophages = "macrophages_total",
                               Time_h = "real_time") %>% 
  mutate(healthy = total - senescent)

# Pivot to long format
dl <- pivot_longer(data, c(necrotic,healthy,senescent,Kupffer,macrophages, total), names_to = "Cells", values_to = "Count") %>% 
  mutate(exposure = as.factor(exposure))

# Calculate mean and sd
dl_sum <- dl %>% filter(!exposure == 0) %>% group_by(time, Time_h, exposure, Cells) %>% summarise(mean_count = mean(Count),
                                                                                                  sd_count = sd(Count)) %>% 
  ungroup() 

# Make labeller
exp_names <- c("250" = "APAP = 250",
               "300" = "APAP = 300",
               "350" = "APAP = 350",
               "400" = "APAP = 400",
               "450" = "APAP = 450",
               "500" = "APAP = 500",
               "550" = "APAP = 550",
               "600" = "APAP = 600",
               "650" = "APAP = 650")

colors <- c("grey",viridisLite::viridis(4))
names(colors) <- c("total","healthy", "macrophages","necrotic", "senescent")

# Option 1
gplot <- ggplot() + 
  geom_line(data = dl_sum %>% filter(Cells %in% c("total"), exposure == 450), aes(x= Time_h, y = mean_count, color = 'total'), linetype = "dashed") + 
  #geom_errorbar(data = dl_sum %>% filter(Cells %in% c("total")), aes(x = Time_h, ymin= mean_count-sd_count, ymax = mean_count+sd_count, color = "total"), width = 0.1, alpha = 0.1) +
  geom_line(data = dl_sum %>% filter(Cells %in% c("macrophages"), exposure == 450), aes(x = Time_h, y = mean_count, color = "macrophages"), size = 0.5) +
  geom_errorbar(data = dl_sum %>% filter(Cells %in% c("macrophages"), exposure == 450), aes(x = Time_h, ymin= mean_count-sd_count, ymax = mean_count+sd_count, color = "macrophages"), width = 0.1, alpha = 0.1) +
  geom_line(data = dl_sum %>% filter(Cells %in% c("healthy"), exposure == 450), aes(x = Time_h, y = mean_count, color = "healthy"), size = 0.5) +
  geom_errorbar(data = dl_sum %>% filter(Cells %in% c("healthy"), exposure == 450), aes(x = Time_h, ymin= mean_count-sd_count, ymax = mean_count+sd_count, color = "healthy"), width = 0.1, alpha = 0.1) +
  geom_line(data = dl_sum %>% filter(Cells %in% c("necrotic"), exposure == 450), aes(x = Time_h, y = mean_count, color = "necrotic"), size = 0.5) +
  geom_errorbar(data = dl_sum %>% filter(Cells %in% c("necrotic"), exposure == 450), aes(x = Time_h, ymin= mean_count-sd_count, ymax = mean_count+sd_count, color = "necrotic"), width = 0.1, alpha = 0.1) +
  geom_line(data = dl_sum %>% filter(Cells %in% c("senescent"), exposure == 450), aes(x = Time_h, y = mean_count, color = "senescent"), size = 0.5) +
  geom_errorbar(data = dl_sum %>% filter(Cells %in% c("senescent"), exposure == 450), aes(x = Time_h, ymin= mean_count-sd_count, ymax = mean_count+sd_count, color = "senescent"), width = 0.1, alpha = 0.1) +
  labs(x = "Time (h)", y = "Cell count", color = "Cell Type") +
  theme_classic() +
  scale_color_manual(values = colors)
ggsave(paste0("Figures/","Fig5B_cell_counts_withMacrophages.pdf"), plot = gplot, width = 4, height = 2)

# Option 2
ggplot() + 
  geom_line(data = dl_sum %>% filter(Cells %in% c("senescent"), exposure %in% c(400,450,500)), aes(x = Time_h, y = mean_count, color = exposure), size = 0.5) +
  geom_errorbar(data = dl_sum %>% filter(Cells %in% c("senescent"), exposure %in% c(400,450,500)), aes(x = Time_h, ymin= mean_count-sd_count, ymax = mean_count+sd_count, color = exposure), width = 0.1, alpha = 0.1) +
  labs(x = "Time (h)", y = "Senescent cells") +
  theme_classic() +
  scale_color_viridis_d(name = "APAP\ndose (a.u.)") +
  ylim(0,200)
ggsave(paste0("Figures/","Fig5B_senescent_counts_withMacrophages.pdf"), width = 3, height = 2)

# Number of senescent cells after 168 hours
withmacro <- dl_sum %>% filter(Cells == "senescent", Time_h == 168) %>% select(-c(time, Time_h,Cells))

### ### ### ### ### ### 
### Figure 5B right ### 
### ### ### ### ### ###  

# get directory names
model_dir <- "15_Mouseliver_noMf_1-5_10x/"
simu_dirs <- dir(paste0(model_dir,"data/"))

# Define folder code
codes <- data.frame(code = c("0001", "0002","0003","0004", "0005"),
                    exposure = c(400,450,500,550,600))


# Cell counts
file.name <- "log_cellcounts.csv"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")
if (!file.exists(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))) {gunzip(file.path, remove=FALSE)}

# Get first data frame
data <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
data$exposure <- codes %>% 
  filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
  pull(exposure)
data$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))

# Bind all other data frames
for (d in simu_dirs[2:length(simu_dirs)]) {
  file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
  if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
  data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
  data_tmp$exposure <- codes %>% 
    filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
    pull(exposure)
  data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
  data <- bind_rows(data,data_tmp)
}

# Rename columns
data <- data %>% dplyr::rename(necrotic = "necrotic_total",
                               total = "hepatocytes_total",
                               senescent = "senescent_total",
                               Kupffer = "kupffer_total",
                               macrophages = "macrophages_total",
                               Time_h = "real_time") %>% 
  mutate(healthy = total - senescent)

# Pivot to long format
dl <- pivot_longer(data, c(necrotic,healthy,senescent,Kupffer,macrophages, total), names_to = "Cells", values_to = "Count") %>% 
  mutate(exposure = as.factor(exposure))

# Calculate mean and sd
dl_sum <- dl %>% filter(!exposure == 0) %>% group_by(time, Time_h, exposure, Cells) %>% summarise(mean_count = mean(Count),
                                                                                                  sd_count = sd(Count)) %>% 
  ungroup() 

# Make labeller
exp_names <- c("250" = "APAP = 250",
               "300" = "APAP = 300",
               "350" = "APAP = 350",
               "400" = "APAP = 400",
               "450" = "APAP = 450",
               "500" = "APAP = 500",
               "550" = "APAP = 550",
               "600" = "APAP = 600",
               "650" = "APAP = 650")

colors <- c("grey",viridisLite::viridis(4))
names(colors) <- c("total","healthy", "macrophages","necrotic", "senescent")

# Option 1
gplot <- ggplot() + 
  geom_line(data = dl_sum %>% filter(Cells %in% c("total"), exposure == 450), aes(x= Time_h, y = mean_count, color = 'total'), linetype = "dashed") + 
  #geom_errorbar(data = dl_sum %>% filter(Cells %in% c("total")), aes(x = Time_h, ymin= mean_count-sd_count, ymax = mean_count+sd_count, color = "total"), width = 0.1, alpha = 0.1) +
  geom_line(data = dl_sum %>% filter(Cells %in% c("macrophages"), exposure == 450), aes(x = Time_h, y = mean_count, color = "macrophages"), size = 0.5) +
  geom_errorbar(data = dl_sum %>% filter(Cells %in% c("macrophages"), exposure == 450), aes(x = Time_h, ymin= mean_count-sd_count, ymax = mean_count+sd_count, color = "macrophages"), width = 0.1, alpha = 0.1) +
  geom_line(data = dl_sum %>% filter(Cells %in% c("healthy"), exposure == 450), aes(x = Time_h, y = mean_count, color = "healthy"), size = 0.5) +
  geom_errorbar(data = dl_sum %>% filter(Cells %in% c("healthy"), exposure == 450), aes(x = Time_h, ymin= mean_count-sd_count, ymax = mean_count+sd_count, color = "healthy"), width = 0.1, alpha = 0.1) +
  geom_line(data = dl_sum %>% filter(Cells %in% c("necrotic"), exposure == 450), aes(x = Time_h, y = mean_count, color = "necrotic"), size = 0.5) +
  geom_errorbar(data = dl_sum %>% filter(Cells %in% c("necrotic"), exposure == 450), aes(x = Time_h, ymin= mean_count-sd_count, ymax = mean_count+sd_count, color = "necrotic"), width = 0.1, alpha = 0.1) +
  geom_line(data = dl_sum %>% filter(Cells %in% c("senescent"), exposure == 450), aes(x = Time_h, y = mean_count, color = "senescent"), size = 0.5) +
  geom_errorbar(data = dl_sum %>% filter(Cells %in% c("senescent"), exposure == 450), aes(x = Time_h, ymin= mean_count-sd_count, ymax = mean_count+sd_count, color = "senescent"), width = 0.1, alpha = 0.1) +
  labs(x = "Time (h)", y = "Cell count", color = "Cell Type") +
  theme_classic() +
  scale_color_manual(values = colors)
ggsave(paste0("Figures/","Fig5B_cell_counts_withoutMacrophages.pdf"), plot = gplot, width = 4, height = 2)

# Option 2
ggplot() + 
  geom_line(data = dl_sum %>% filter(Cells %in% c("senescent"), exposure %in% c(400,450,500)), aes(x = Time_h, y = mean_count, color = exposure), size = 0.5) +
  geom_errorbar(data = dl_sum %>% filter(Cells %in% c("senescent"), exposure %in% c(400,450,500)), aes(x = Time_h, ymin= mean_count-sd_count, ymax = mean_count+sd_count, color = exposure), width = 0.1, alpha = 0.1) +
  labs(x = "Time (h)", y = "Senescent cells") +
  theme_classic() +
  scale_color_viridis_d(name = "APAP\ndose (a.u.)") +
  ylim(0,200)
ggsave(paste0("Figures/","Fig5B_senescent_counts_withoutMacrophages.pdf"), width = 3, height = 2)

# Number of senescent cells after 168 hours
withoutmacro <- dl_sum %>% filter(Cells == "senescent", Time_h == 168) %>% select(-c(time, Time_h,Cells))

# Percentage increase in senescent cells
senesence <- left_join(withmacro, withoutmacro, by = c("exposure"), suffix = c("_with","_without"))
senesence <- senesence %>% group_by(exposure) %>% mutate(perc_increase = (mean_count_with/mean_count_without) * 100)

### ### ### ### ### ### ###
###  Figure 5C and 5D   ###
### ### ### ### ### ### ###

### part 1 ###
# k_inhibitor
# get directory names
model_dir <- "14_Mouseliver_k_inhib_1-3_10x/"
simu_dirs <- dir(paste0(model_dir,"data/"))

# Define folder code
codes <- data.frame(code = c("0001", "0002","0003"),
                    r = c(0.1,1,10))


# Cell counts
file.name <- "log_cellcounts.csv"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")
if (!file.exists(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))) {gunzip(file.path, remove=FALSE)}

# Get first data frame
data_k_inhib <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
data_k_inhib$r <- codes %>% 
  filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
  pull(r)
data_k_inhib$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))

# Bind all other data frames
for (d in simu_dirs[2:length(simu_dirs)]) {
  file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
  if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
  data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
  data_tmp$r <- codes %>% 
    filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
    pull(r)
  data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
  data_k_inhib <- bind_rows(data_k_inhib,data_tmp)
}

# Rename columns
data_k_inhib <- data_k_inhib %>% dplyr::rename(necrotic = "necrotic_total",
                                               total = "hepatocytes_total",
                                               senescent = "senescent_total",
                                               Kupffer = "kupffer_total",
                                               macrophages = "macrophages_total",
                                               Time_h = "real_time") %>% 
  mutate(healthy = total - senescent)

# Pivot to long format
dl_k_inhib <- pivot_longer(data_k_inhib, c(necrotic,healthy,senescent,Kupffer,macrophages, total), names_to = "Cells", values_to = "Count")# %>% 
#mutate(r = log10(r))

dl_k_inhib <- dl_k_inhib %>% mutate(Parameter = "k_inhibitor")

# Calculate mean and sd
dl_sum_k_inhib <- dl_k_inhib %>% group_by(time, Time_h, r, Cells, Parameter) %>% 
  summarise(mean_count = mean(Count),
            sd_count = sd(Count)) %>% 
  ungroup()

# k_stim
# get directory names
model_dir <- "14_Mouseliver_k_stim_1-4_10x/"
simu_dirs <- dir(paste0(model_dir,"data/"))

# Define folder code
codes <- data.frame(code = c("0001", "0002","0003","0004"),
                    r = c(0.01, 0.1,1,10))


# Cell counts
file.name <- "log_cellcounts.csv"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")
if (!file.exists(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))) {gunzip(file.path, remove=FALSE)}

# Get first data frame
data_k_stim <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
data_k_stim$r <- codes %>% 
  filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
  pull(r)
data_k_stim$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))

# Bind all other data frames
for (d in simu_dirs[2:length(simu_dirs)]) {
  file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
  if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
  data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
  data_tmp$r <- codes %>% 
    filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
    pull(r)
  data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
  data_k_stim <- bind_rows(data_k_stim,data_tmp)
}

# Rename columns
data_k_stim <- data_k_stim %>% dplyr::rename(necrotic = "necrotic_total",
                                             total = "hepatocytes_total",
                                             senescent = "senescent_total",
                                             Kupffer = "kupffer_total",
                                             macrophages = "macrophages_total",
                                             Time_h = "real_time") %>% 
  mutate(healthy = total - senescent)

# Pivot to long format
dl_k_stim <- pivot_longer(data_k_stim, c(necrotic,healthy,senescent,Kupffer,macrophages, total), names_to = "Cells", values_to = "Count")# %>% 
#mutate(r = log10(r))

dl_k_stim <- dl_k_stim %>% mutate(Parameter = "k_stim")

# Calculate mean and sd
dl_sum_k_stim <- dl_k_stim %>% group_by(time, Time_h, r, Cells, Parameter) %>% 
  summarise(mean_count = mean(Count),
            sd_count = sd(Count)) %>% 
  ungroup() 

# k_mitogen
# get directory names
model_dir <- "14_Mouseliver_k_mitogen_1-3_10x/"
simu_dirs <- dir(paste0(model_dir,"data/"))

# Define folder code
codes <- data.frame(code = c("0001", "0002","0003"),
                    r = c(0.1,1,10))


# Cell counts
file.name <- "log_cellcounts.csv"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")
if (!file.exists(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))) {gunzip(file.path, remove=FALSE)}

# Get first data frame
data_k_mitogen <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
data_k_mitogen$r <- codes %>% 
  filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
  pull(r)
data_k_mitogen$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))

# Bind all other data frames
for (d in simu_dirs[2:length(simu_dirs)]) {
  file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
  if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
  data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
  data_tmp$r <- codes %>% 
    filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
    pull(r)
  data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
  data_k_mitogen <- bind_rows(data_k_mitogen,data_tmp)
}

# Rename columns
data_k_mitogen <- data_k_mitogen %>% dplyr::rename(necrotic = "necrotic_total",
                                                   total = "hepatocytes_total",
                                                   senescent = "senescent_total",
                                                   Kupffer = "kupffer_total",
                                                   macrophages = "macrophages_total",
                                                   Time_h = "real_time") %>% 
  mutate(healthy = total - senescent)

# Pivot to long format
dl_k_mitogen <- pivot_longer(data_k_mitogen, c(necrotic,healthy,senescent,Kupffer,macrophages, total), names_to = "Cells", values_to = "Count")# %>% 
#mutate(r = log10(r))

dl_k_mitogen <- dl_k_mitogen %>% mutate(Parameter = "k_mitogen")

# Calculate mean and sd
dl_sum_k_mitogen <- dl_k_mitogen %>% group_by(time, Time_h, r, Cells,Parameter) %>% 
  summarise(mean_count = mean(Count),
            sd_count = sd(Count)) %>% 
  ungroup()

### part 2 ###
# k_inhibitor
# get directory names
model_dir <- "14_Mouseliver_k_inhib_1-8_10x/"
simu_dirs <- dir(paste0(model_dir,"data/"))

# Define folder code
codes <- data.frame(code = c("0001", "0002","0003","0004", "0005","0006","0007", "0008"),
                    r = c(0.125,0.167,0.25,0.5,2,4,6,8))

# Cell counts
file.name <- "log_cellcounts.csv"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")
if (!file.exists(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))) {gunzip(file.path, remove=FALSE)}

# Get first data frame
data_k_inhib_2 <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
data_k_inhib_2$r <- codes %>% 
  filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
  pull(r)
data_k_inhib_2$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))

# Bind all other data frames
for (d in simu_dirs[2:length(simu_dirs)]) {
  file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
  if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
  data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
  data_tmp$r <- codes %>% 
    filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
    pull(r)
  data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
  data_k_inhib_2 <- bind_rows(data_k_inhib_2,data_tmp)
}

# Rename columns
data_k_inhib_2 <- data_k_inhib_2 %>% dplyr::rename(necrotic = "necrotic_total",
                                                   total = "hepatocytes_total",
                                                   senescent = "senescent_total",
                                                   Kupffer = "kupffer_total",
                                                   macrophages = "macrophages_total",
                                                   Time_h = "real_time") %>% 
  mutate(healthy = total - senescent)

# Pivot to long format
dl_k_inhib_2 <- pivot_longer(data_k_inhib_2, c(necrotic,healthy,senescent,Kupffer,macrophages, total), names_to = "Cells", values_to = "Count")# %>% 
#mutate(r = log10(r))

dl_k_inhib_2 <- dl_k_inhib_2 %>% mutate(Parameter = "k_inhibitor")

# Calculate mean and sd
dl_sum_k_inhib_2 <- dl_k_inhib_2 %>% group_by(time, Time_h, r, Cells, Parameter) %>% 
  summarise(mean_count = mean(Count),
            sd_count = sd(Count)) %>% 
  ungroup()

# k_stim
# get directory names
model_dir <- "14_Mouseliver_k_stim_1-8_10x/"
simu_dirs <- dir(paste0(model_dir,"data/"))

# Define folder code
codes <- data.frame(code = c("0001", "0002","0003","0004", "0005","0006","0007", "0008"),
                    r = c(0.125,0.167,0.25,0.5,2,4,6,8))

# Cell counts
file.name <- "log_cellcounts.csv"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")
if (!file.exists(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))) {gunzip(file.path, remove=FALSE)}

# Get first data frame
data_k_stim_2 <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
data_k_stim_2$r <- codes %>% 
  filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
  pull(r)
data_k_stim_2$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))

# Bind all other data frames
for (d in simu_dirs[2:length(simu_dirs)]) {
  file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
  if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
  data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
  data_tmp$r <- codes %>% 
    filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
    pull(r)
  data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
  data_k_stim_2 <- bind_rows(data_k_stim_2,data_tmp)
}

# Rename columns
data_k_stim_2 <- data_k_stim_2 %>% dplyr::rename(necrotic = "necrotic_total",
                                                 total = "hepatocytes_total",
                                                 senescent = "senescent_total",
                                                 Kupffer = "kupffer_total",
                                                 macrophages = "macrophages_total",
                                                 Time_h = "real_time") %>% 
  mutate(healthy = total - senescent)

# Pivot to long format
dl_k_stim_2 <- pivot_longer(data_k_stim_2, c(necrotic,healthy,senescent,Kupffer,macrophages, total), names_to = "Cells", values_to = "Count")# %>% 
#mutate(r = log10(r))

dl_k_stim_2 <- dl_k_stim_2 %>% mutate(Parameter = "k_stim")

# Calculate mean and sd
dl_sum_k_stim_2 <- dl_k_stim_2 %>% group_by(time, Time_h, r, Cells, Parameter) %>% 
  summarise(mean_count = mean(Count),
            sd_count = sd(Count)) %>% 
  ungroup()

# k_mitogen
# get directory names
model_dir <- "14_Mouseliver_k_mitogen_1-8_10x/"
simu_dirs <- dir(paste0(model_dir,"data/"))

# Define folder code
codes <- data.frame(code = c("0001", "0002","0003","0004", "0005","0006","0007", "0008"),
                    r = c(0.125,0.167,0.25,0.5,2,4,6,8))

# Cell counts
file.name <- "log_cellcounts.csv"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")
if (!file.exists(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))) {gunzip(file.path, remove=FALSE)}

# Get first data frame
data_k_mitogen_2 <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
data_k_mitogen_2$r <- codes %>% 
  filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
  pull(r)
data_k_mitogen_2$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))

# Bind all other data frames
for (d in simu_dirs[2:length(simu_dirs)]) {
  file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
  if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
  data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
  data_tmp$r <- codes %>% 
    filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
    pull(r)
  data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
  data_k_mitogen_2 <- bind_rows(data_k_mitogen_2,data_tmp)
}

# Rename columns
data_k_mitogen_2 <- data_k_mitogen_2 %>% dplyr::rename(necrotic = "necrotic_total",
                                                       total = "hepatocytes_total",
                                                       senescent = "senescent_total",
                                                       Kupffer = "kupffer_total",
                                                       macrophages = "macrophages_total",
                                                       Time_h = "real_time") %>% 
  mutate(healthy = total - senescent)

# Pivot to long format
dl_k_mitogen_2 <- pivot_longer(data_k_mitogen_2, c(necrotic,healthy,senescent,Kupffer,macrophages, total), names_to = "Cells", values_to = "Count")# %>% 
#mutate(r = log10(r))

dl_k_mitogen_2 <- dl_k_mitogen_2 %>% mutate(Parameter = "k_mitogen")

# Calculate mean and sd
dl_sum_k_mitogen_2 <- dl_k_mitogen_2 %>% group_by(time, Time_h, r, Cells, Parameter) %>% 
  summarise(mean_count = mean(Count),
            sd_count = sd(Count)) %>% 
  ungroup() 

### ### ### ###

# Bind data frames
dl_all <- bind_rows(dl_k_inhib, dl_k_stim, dl_k_mitogen,
                    dl_k_inhib_2, dl_k_stim_2, dl_k_mitogen_2)

# Bind data frames
dl_sum <- bind_rows(dl_sum_k_inhib, dl_sum_k_stim, dl_sum_k_mitogen,
                    dl_sum_k_inhib_2, dl_sum_k_stim_2, dl_sum_k_mitogen_2)

dl_sum_tmp <- dl_sum %>% filter(Time_h == 48, Cells == "senescent", !r == 0.01)
# 4-fold increase in k_stim
dl_sum_tmp %>% filter(r == 4, Parameter == "k_stim") %>% pull(mean_count) /
  dl_sum_tmp %>% filter(r == 1, Parameter == "k_stim") %>% pull(mean_count)
# 10-fold increase in k_inhib
dl_sum_tmp %>% filter(r == 10, Parameter == "k_inhibitor") %>% pull(mean_count) /
  dl_sum_tmp %>% filter(r == 1, Parameter == "k_inhibitor") %>% pull(mean_count)

ggplot(dl_sum %>% filter(Time_h == 48, Cells == "senescent", !r == 0.01)) + 
  geom_point(aes(x = r, y = mean_count, color = Parameter)) +
  geom_errorbar(aes(x = r, ymin = mean_count-sd_count, ymax = mean_count+sd_count, color = Parameter), width = 0.05) +
  geom_line(aes(x = r, y = mean_count, color = Parameter)) +
  labs(x = "Multiplication factor", y = "Senescent cells") +
  theme_classic() +
  #theme(axis.text.x = element_text(angle = 90)) +
  scale_x_continuous(trans='log2', breaks = c(0.1,0.125,0.167,0.25,0.5,1,2,4,6,8,10),
                     labels = c("0.1","","","0.25","0.5","1","2","4","","","10")) +
  scale_color_viridis_d(labels = c(expression(k[inhib]), 
                                   expression(k[mit]*"  "), 
                                   expression(k[stim]*" ")))
ggsave(paste0("Figures/","Fig5C_senescent_counts_q.pdf"), width = 3.5, height = 2)

d1 <- dl_sum %>% filter(Cells %in% c("senescent"), !r == 0.01) %>% 
  mutate(r = as.factor(round(r,3)),
         Parameter_2 = sapply(Parameter,function (p){
           out <- ifelse(p=="k_inhibitor", "k[inhib]", 
                         ifelse(p=="k_mitogen","k[mit]","k[stim]"))
           out
         }))
ggplot() + 
  geom_line(data = d1, aes(x = Time_h, y = mean_count, color = r)) +
  geom_errorbar(data = d1, aes(x = Time_h, ymin = mean_count-sd_count, ymax = mean_count+sd_count, color = r), width = 0.05) +
  labs(x = "Time (h)", y = "Senescent hepatocytes") +
  theme_classic() +
  scale_color_viridis_d(name = "Multiplication\nfactor") +
  facet_grid(~Parameter_2, labeller = label_parsed) + 
  theme(strip.text.x = element_text(size = 14),
        axis.title = element_text(size = 16))
ggsave(paste0("Figures/","Fig5D_senes_prol.pdf"), width = 10, height = 3.25)

data_senescence <- dl_all %>% filter(Time_h == 48, !r == 0.01) %>% 
  pivot_wider(values_from = c(Count), id_cols = c(r, Parameter, Time_h,run), names_from = Cells)

### ### ### ### ### 
###  Figure 5E  ###
### ### ### ### ###

### part 1 ###
# k_inhibitor
# get directory names
model_dir <- "14_Mouseliver_k_inhib_1-3_10x/"
simu_dirs <- dir(paste0(model_dir,"data/"))

# Define folder code
codes <- data.frame(code = c("0001", "0002","0003"),
                    r = c(0.1,1,10))


# Cell counts
file.name <- "cell_division_hepatocytes.txt"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")
if (!file.exists(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))) {gunzip(file.path, remove=FALSE)}

# Get first data frame
prol_k_inhib <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
if (!nrow(prol_k_inhib) == 0) {
  prol_k_inhib <- prol_k_inhib %>% filter(Time <=672)
  prol_k_inhib$r <- codes %>% 
    filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
    pull(r)
  prol_k_inhib$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))
  prol_k_inhib <- prol_k_inhib %>% group_by(r,run) %>% summarise(proliferated = n())
  if (nrow(prol_k_inhib) == 0) {
    prol_k_inhib <- tibble(proliferated = 0, 
                           r = codes %>% 
                             filter(code == regmatches(simu_dirs[51],regexpr("[0-9]{4}",simu_dirs[51]))) %>% 
                             pull(r),
                           run = sub('.', '', regmatches(simu_dirs[51],regexpr("-[0-9]{1,3}",simu_dirs[51]))))
  }
  
} else {
  prol_k_inhib <- tibble(proliferated = 0, 
                         r = codes %>% 
                           filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
                           pull(r),
                         run = sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1]))))
}

# Bind all other data frames
for (d in simu_dirs[2:length(simu_dirs)]) {
  file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
  if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
  data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
  if (!nrow(data_tmp) == 0) {
    data_tmp <- data_tmp %>% filter(Time <=672)
    data_tmp$r <- codes %>% 
      filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
      pull(r)
    data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
    data_tmp <- data_tmp %>% group_by(r,run) %>% summarise(proliferated = n())
    if (nrow(data_tmp) == 0) {
      data_tmp <- tibble(proliferated = 0, 
                         r = codes %>% 
                           filter(code == regmatches(simu_dirs[51],regexpr("[0-9]{4}",simu_dirs[51]))) %>% 
                           pull(r),
                         run = sub('.', '', regmatches(simu_dirs[51],regexpr("-[0-9]{1,3}",simu_dirs[51]))))
    }
    
  } else {
    data_tmp <- tibble(proliferated = 0, 
                       r = codes %>% 
                         filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
                         pull(r),
                       run = sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d))))
  }
  prol_k_inhib <- bind_rows(prol_k_inhib,data_tmp)
}

sum_prol_k_inhib <- prol_k_inhib %>% group_by(r) %>% 
  # summarise(mean_count_proliferated = mean(proliferated),
  #           sd_count_proliferated = sd(proliferated)) %>% 
  ungroup() %>% mutate(Parameter = "k_inhibitor", Time_h = 48)

# k_stim
# get directory names
model_dir <- "14_Mouseliver_k_stim_1-4_10x/"
simu_dirs <- dir(paste0(model_dir,"data/"))

# Define folder code
codes <- data.frame(code = c("0001", "0002","0003","0004"),
                    r = c(0.01, 0.1,1,10))


# Cell counts
file.name <- "cell_division_hepatocytes.txt"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")
if (!file.exists(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))) {gunzip(file.path, remove=FALSE)}

# Get first data frame
prol_k_stim <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
if (!nrow(prol_k_stim) == 0) {
  prol_k_stim <- prol_k_stim %>% filter(Time <=672)
  prol_k_stim$r <- codes %>% 
    filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
    pull(r)
  prol_k_stim$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))
  prol_k_stim <- prol_k_stim %>% group_by(r,run) %>% summarise(proliferated = n())
  if (nrow(prol_k_stim) == 0) {
    prol_k_stim <- tibble(proliferated = 0, 
                          r = codes %>% 
                            filter(code == regmatches(simu_dirs[51],regexpr("[0-9]{4}",simu_dirs[51]))) %>% 
                            pull(r),
                          run = sub('.', '', regmatches(simu_dirs[51],regexpr("-[0-9]{1,3}",simu_dirs[51]))))
  }
  
} else {
  prol_k_stim <- tibble(proliferated = 0, 
                        r = codes %>% 
                          filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
                          pull(r),
                        run = sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1]))))
}

# Bind all other data frames
for (d in simu_dirs[2:length(simu_dirs)]) {
  file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
  if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
  data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
  if (!nrow(data_tmp) == 0) {
    data_tmp <- data_tmp %>% filter(Time <=672)
    data_tmp$r <- codes %>% 
      filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
      pull(r)
    data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
    data_tmp <- data_tmp %>% group_by(r,run) %>% summarise(proliferated = n())
    if (nrow(data_tmp) == 0) {
      data_tmp <- tibble(proliferated = 0, 
                         r = codes %>% 
                           filter(code == regmatches(simu_dirs[51],regexpr("[0-9]{4}",simu_dirs[51]))) %>% 
                           pull(r),
                         run = sub('.', '', regmatches(simu_dirs[51],regexpr("-[0-9]{1,3}",simu_dirs[51]))))
    }
    
  } else {
    data_tmp <- tibble(proliferated = 0, 
                       r = codes %>% 
                         filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
                         pull(r),
                       run = sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d))))
  }
  prol_k_stim <- bind_rows(prol_k_stim,data_tmp)
}

sum_prol_k_stim <- prol_k_stim %>% group_by(r) %>% 
  # summarise(mean_count_proliferated = mean(proliferated),
  #           sd_count_proliferated = sd(proliferated)) %>% 
  ungroup() %>% mutate(Parameter = "k_stim", Time_h = 48)

# k_mitogen
# get directory names
model_dir <- "14_Mouseliver_k_mitogen_1-3_10x/"
simu_dirs <- dir(paste0(model_dir,"data/"))

# Define folder code
codes <- data.frame(code = c("0001", "0002","0003"),
                    r = c(0.1,1,10))


# Cell counts
file.name <- "cell_division_hepatocytes.txt"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")
if (!file.exists(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))) {gunzip(file.path, remove=FALSE)}

# Get first data frame
prol_k_mitogen <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
if (!nrow(prol_k_mitogen) == 0) {
  prol_k_mitogen <- prol_k_mitogen %>% filter(Time <=672)
  prol_k_mitogen$r <- codes %>% 
    filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
    pull(r)
  prol_k_mitogen$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))
  prol_k_mitogen <- prol_k_mitogen %>% group_by(r,run) %>% summarise(proliferated = n())
  if (nrow(prol_k_mitogen) == 0) {
    prol_k_mitogen <- tibble(proliferated = 0, 
                             r = codes %>% 
                               filter(code == regmatches(simu_dirs[51],regexpr("[0-9]{4}",simu_dirs[51]))) %>% 
                               pull(r),
                             run = sub('.', '', regmatches(simu_dirs[51],regexpr("-[0-9]{1,3}",simu_dirs[51]))))
  }
  
} else {
  prol_k_mitogen <- tibble(proliferated = 0, 
                           r = codes %>% 
                             filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
                             pull(r),
                           run = sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1]))))
}

# Bind all other data frames
for (d in simu_dirs[2:length(simu_dirs)]) {
  file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
  if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
  data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
  if (!nrow(data_tmp) == 0) {
    data_tmp <- data_tmp %>% filter(Time <=672)
    data_tmp$r <- codes %>% 
      filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
      pull(r)
    data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
    data_tmp <- data_tmp %>% group_by(r,run) %>% summarise(proliferated = n())
    if (nrow(data_tmp) == 0) {
      data_tmp <- tibble(proliferated = 0, 
                         r = codes %>% 
                           filter(code == regmatches(simu_dirs[51],regexpr("[0-9]{4}",simu_dirs[51]))) %>% 
                           pull(r),
                         run = sub('.', '', regmatches(simu_dirs[51],regexpr("-[0-9]{1,3}",simu_dirs[51]))))
    }
    
  } else {
    data_tmp <- tibble(proliferated = 0, 
                       r = codes %>% 
                         filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
                         pull(r),
                       run = sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d))))
  }
  prol_k_mitogen <- bind_rows(prol_k_mitogen,data_tmp)
}

sum_prol_k_mitogen <- prol_k_mitogen %>% group_by(r) %>% 
  # summarise(mean_count_proliferated = mean(proliferated),
  #           sd_count_proliferated = sd(proliferated)) %>% 
  ungroup() %>% mutate(Parameter = "k_mitogen", Time_h = 48)

### part 2 ###
# k_inhibitor
# get directory names
model_dir <- "14_Mouseliver_k_inhib_1-8_10x/"
simu_dirs <- dir(paste0(model_dir,"data/"))

# Define folder code
codes <- data.frame(code = c("0001", "0002","0003","0004", "0005","0006","0007", "0008"),
                    r = c(0.125,0.167,0.25,0.5,2,4,6,8))


# Cell counts
file.name <- "cell_division_hepatocytes.txt"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")
if (!file.exists(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))) {gunzip(file.path, remove=FALSE)}

# Get first data frame
prol_k_inhib_2 <- read_tsv(paste0(model_dir,"data/",simu_dirs[51],"/",file.name))
if (!nrow(prol_k_inhib_2) == 0) {
  prol_k_inhib_2 <- prol_k_inhib_2 %>% filter(Time <=672)
  prol_k_inhib_2$r <- codes %>% 
    filter(code == regmatches(simu_dirs[51],regexpr("[0-9]{4}",simu_dirs[51]))) %>% 
    pull(r)
  prol_k_inhib_2$run <- sub('.', '', regmatches(simu_dirs[51],regexpr("-[0-9]{1,3}",simu_dirs[51])))
  prol_k_inhib_2 <- prol_k_inhib_2 %>% group_by(r,run) %>% summarise(proliferated = n())
  if (nrow(prol_k_inhib_2) == 0) {
    prol_k_inhib_2 <- tibble(proliferated = 0, 
                             r = codes %>% 
                               filter(code == regmatches(simu_dirs[51],regexpr("[0-9]{4}",simu_dirs[51]))) %>% 
                               pull(r),
                             run = sub('.', '', regmatches(simu_dirs[51],regexpr("-[0-9]{1,3}",simu_dirs[51]))))
  }
} else {
  prol_k_inhib_2 <- tibble(proliferated = 0, 
                           r = codes %>% 
                             filter(code == regmatches(simu_dirs[51],regexpr("[0-9]{4}",simu_dirs[51]))) %>% 
                             pull(r),
                           run = sub('.', '', regmatches(simu_dirs[51],regexpr("-[0-9]{1,3}",simu_dirs[51]))))
}

# Bind all other data frames
for (d in simu_dirs[2:length(simu_dirs)]) {
  file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
  if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
  data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
  if (!nrow(data_tmp) == 0) {
    data_tmp <- data_tmp %>% filter(Time <=672)
    data_tmp$r <- codes %>% 
      filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
      pull(r)
    data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
    data_tmp <- data_tmp %>% group_by(r,run) %>% summarise(proliferated = n())
    if (nrow(data_tmp) == 0) {
      data_tmp <- tibble(proliferated = 0, 
                         r = codes %>% 
                           filter(code == regmatches(simu_dirs[51],regexpr("[0-9]{4}",simu_dirs[51]))) %>% 
                           pull(r),
                         run = sub('.', '', regmatches(simu_dirs[51],regexpr("-[0-9]{1,3}",simu_dirs[51]))))
    }
    
  } else {
    data_tmp <- tibble(proliferated = 0, 
                       r = codes %>% 
                         filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
                         pull(r),
                       run = sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d))))
  }
  prol_k_inhib_2 <- bind_rows(prol_k_inhib_2,data_tmp)
}

sum_prol_k_inhib_2 <- prol_k_inhib_2 %>% group_by(r) %>% 
  # summarise(mean_count_proliferated = mean(proliferated),
  #           sd_count_proliferated = sd(proliferated)) %>% 
  ungroup() %>% mutate(Parameter = "k_inhibitor", Time_h = 48)

# k_stim
# get directory names
model_dir <- "14_Mouseliver_k_stim_1-8_10x/"
simu_dirs <- dir(paste0(model_dir,"data/"))

# Define folder code
codes <- data.frame(code = c("0001", "0002","0003","0004", "0005","0006","0007", "0008"),
                    r = c(0.125,0.167,0.25,0.5,2,4,6,8))


# Cell counts
file.name <- "cell_division_hepatocytes.txt"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")
if (!file.exists(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))) {gunzip(file.path, remove=FALSE)}

# Get first data frame
prol_k_stim_2 <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
if (!nrow(prol_k_stim_2) == 0) {
  prol_k_stim_2 <- prol_k_stim_2 %>% filter(Time <=672)
  prol_k_stim_2$r <- codes %>% 
    filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
    pull(r)
  prol_k_stim_2$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))
  prol_k_stim_2 <- prol_k_stim_2 %>% group_by(r,run) %>% summarise(proliferated = n())
  if (nrow(prol_k_stim_2) == 0) {
    prol_k_stim_2 <- tibble(proliferated = 0, 
                            r = codes %>% 
                              filter(code == regmatches(simu_dirs[51],regexpr("[0-9]{4}",simu_dirs[51]))) %>% 
                              pull(r),
                            run = sub('.', '', regmatches(simu_dirs[51],regexpr("-[0-9]{1,3}",simu_dirs[51]))))
  }
  
} else {
  prol_k_stim_2 <- tibble(proliferated = 0, 
                          r = codes %>% 
                            filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
                            pull(r),
                          run = sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1]))))
}

# Bind all other data frames
for (d in simu_dirs[2:length(simu_dirs)]) {
  file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
  if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
  data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
  if (!nrow(data_tmp) == 0) {
    data_tmp <- data_tmp %>% filter(Time <=672)
    data_tmp$r <- codes %>% 
      filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
      pull(r)
    data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
    data_tmp <- data_tmp %>% group_by(r,run) %>% summarise(proliferated = n())
    if (nrow(data_tmp) == 0) {
      data_tmp <- tibble(proliferated = 0, 
                         r = codes %>% 
                           filter(code == regmatches(simu_dirs[51],regexpr("[0-9]{4}",simu_dirs[51]))) %>% 
                           pull(r),
                         run = sub('.', '', regmatches(simu_dirs[51],regexpr("-[0-9]{1,3}",simu_dirs[51]))))
    }
    
  } else {
    data_tmp <- tibble(proliferated = 0, 
                       r = codes %>% 
                         filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
                         pull(r),
                       run = sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d))))
  }
  prol_k_stim_2 <- bind_rows(prol_k_stim_2,data_tmp)
}

sum_prol_k_stim_2 <- prol_k_stim_2 %>% group_by(r) %>% 
  # summarise(mean_count_proliferated = mean(proliferated),
  #           sd_count_proliferated = sd(proliferated)) %>% 
  ungroup() %>% mutate(Parameter = "k_stim", Time_h = 48)

# k_mitogen
# get directory names
model_dir <- "14_Mouseliver_k_mitogen_1-8_10x/"
simu_dirs <- dir(paste0(model_dir,"data/"))

# Define folder code
codes <- data.frame(code = c("0001", "0002","0003","0004", "0005","0006","0007", "0008"),
                    r = c(0.125,0.167,0.25,0.5,2,4,6,8))

# Cell counts
file.name <- "cell_division_hepatocytes.txt"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")
if (!file.exists(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))) {gunzip(file.path, remove=FALSE)}

# Get first data frame
prol_k_mitogen_2 <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
if (!nrow(prol_k_mitogen_2) == 0) {
  prol_k_mitogen_2 <- prol_k_mitogen_2 %>% filter(Time <=672)
  prol_k_mitogen_2$r <- codes %>% 
    filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
    pull(r)
  prol_k_mitogen_2$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))
  prol_k_mitogen_2 <- prol_k_mitogen_2 %>% group_by(r,run) %>% summarise(proliferated = n())
  if (nrow(prol_k_mitogen_2) == 0) {
    prol_k_mitogen_2 <- tibble(proliferated = 0, 
                               r = codes %>% 
                                 filter(code == regmatches(simu_dirs[51],regexpr("[0-9]{4}",simu_dirs[51]))) %>% 
                                 pull(r),
                               run = sub('.', '', regmatches(simu_dirs[51],regexpr("-[0-9]{1,3}",simu_dirs[51]))))
  }
  
} else {
  prol_k_mitogen_2 <- tibble(proliferated = 0, 
                             r = codes %>% 
                               filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
                               pull(r),
                             run = sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1]))))
}

# Bind all other data frames
for (d in simu_dirs[2:length(simu_dirs)]) {
  file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
  if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
  data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
  if (!nrow(data_tmp) == 0) {
    data_tmp <- data_tmp %>% filter(Time <=672)
    data_tmp$r <- codes %>% 
      filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
      pull(r)
    data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
    data_tmp <- data_tmp %>% group_by(r,run) %>% summarise(proliferated = n())
    if (nrow(data_tmp) == 0) {
      data_tmp <- tibble(proliferated = 0, 
                         r = codes %>% 
                           filter(code == regmatches(simu_dirs[51],regexpr("[0-9]{4}",simu_dirs[51]))) %>% 
                           pull(r),
                         run = sub('.', '', regmatches(simu_dirs[51],regexpr("-[0-9]{1,3}",simu_dirs[51]))))
    }
    
  } else {
    data_tmp <- tibble(proliferated = 0, 
                       r = codes %>% 
                         filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
                         pull(r),
                       run = sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d))))
  }
  prol_k_mitogen_2 <- bind_rows(prol_k_mitogen_2,data_tmp)
}

sum_prol_k_mitogen_2 <- prol_k_mitogen_2 %>% group_by(r) %>% 
  # summarise(mean_count_proliferated = mean(proliferated),
  #           sd_count_proliferated = sd(proliferated)) %>% 
  ungroup() %>% mutate(Parameter = "k_mitogen", Time_h = 48)

# Bind all data

data_proliferation <- bind_rows(sum_prol_k_inhib,sum_prol_k_mitogen,sum_prol_k_stim,
                                sum_prol_k_inhib_2,sum_prol_k_mitogen_2,sum_prol_k_stim_2)

data_sen_prol <- left_join(data_senescence, data_proliferation, by = c("Parameter", "Time_h", "r", "run"))

# Make figure
data_sen_prol <- data_sen_prol %>% mutate(Parameter_2 = sapply(Parameter,function (p){
  out <- ifelse(p=="k_inhibitor", "k[inhib]", 
                ifelse(p=="k_mitogen","k[mit]","k[stim]"))
  out
}))
ggplot(data_sen_prol %>% mutate(r = as.factor(r))) +
  geom_point(aes(x = senescent , y = proliferated, color = r)) +
  theme_classic() +
  scale_color_viridis_d(name = "Multiplication\nfactor") +
  labs(x = "Senescent hepatocytes", y = "Proliferated hepatocytes") +
  facet_grid(~Parameter_2, labeller = label_parsed) + 
  theme(strip.text.x = element_text(size = 14),
        axis.title = element_text(size = 16))
ggsave(paste0("Figures/","Fig5E_senesc_prolif_tradeoff_48h.pdf"), width = 10, height = 3.25)

### ### ### ### ### ### ###
### ### Treatments  ### ###
### ### ### ### ### ### ###

# get directory names
model_dir <- "08_Mouseliver_default_1-10_10x/"
simu_dirs <- dir(paste0(model_dir,"data/"))

# Define folder code
codes <- data.frame(code = c("0001", "0002","0003","0004", "0005","0006","0007", "0008","0009","0010"),
                    exposure = c(0,250,300,350,400,450,500,550,600,650))

file.name <- "log_cellcounts.csv"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")

# Get first data frame
gunzip(file.path, remove=FALSE)
data_control <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
data_control$exposure <- codes %>% 
  filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
  pull(exposure)
data_control$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))

# Bind all other data frames
for (d in simu_dirs[2:length(simu_dirs)]) {
  file.name <- "log_cellcounts.csv"
  file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
  if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
  data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
  data_tmp$exposure <- codes %>% 
    filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
    pull(exposure)
  data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
  data_control <- bind_rows(data_control,data_tmp)
}

# Rename columns
data_control <- data_control %>% dplyr::rename(necrotic = "necrotic_total",
                                               healthy = "hepatocytes_total",
                                               senescent = "senescent_total",
                                               Kupffer = "kupffer_total",
                                               macrophages = "macrophages_total",
                                               Time_h = "real_time") %>% 
  mutate(healthy = healthy - senescent)

# Pivot to long format
dl_control <- pivot_longer(data_control, c(necrotic,healthy,senescent,Kupffer,macrophages), names_to = "Cells", values_to = "Count") %>% 
  mutate(exposure = as.factor(exposure))

# Calculate mean and sd
dl_control_sum <- dl_control %>% group_by(time, Time_h, exposure, Cells) %>% summarise(mean_count = mean(Count),
                                                                                       sd_count = sd(Count)) %>% 
  ungroup()

dl_control_filter <- dl_control_sum %>% filter(exposure == 500)

### ### ### ### ### 
###   Figure 6E ### 
### ### ### ### ### 

# NAC part 1 #

# get directory names
model_dir <- "16_Mouseliver_NAC_1-5_10x/"
simu_dirs <- dir(paste0(model_dir,"data/"))

# Define folder code
codes <- data.frame(code = c("0001", "0002","0003","0004", "0005"),
                    trtmt = c(3,6,12,24,48))

file.name <- "log_cellcounts.csv"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")

# Get first data frame
gunzip(file.path, remove=FALSE)
data_NAC <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
data_NAC$trtmt <- codes %>% 
  filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
  pull(trtmt)
data_NAC$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))

# Bind all other data frames
for (d in simu_dirs[2:length(simu_dirs)]) {
  file.name <- "log_cellcounts.csv"
  file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
  if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
  data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
  data_tmp$trtmt <- codes %>% 
    filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
    pull(trtmt)
  data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
  data_NAC <- bind_rows(data_NAC,data_tmp)
}

# NAC part 2 #

# get directory names
model_dir <- "16_Mouseliver_NAC_1-1_10x/"
simu_dirs <- dir(paste0(model_dir,"data/"))

# Define folder code
codes <- data.frame(code = c("0001"),
                    trtmt = c(18))

file.name <- "log_cellcounts.csv"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")

# Get first data frame
gunzip(file.path, remove=FALSE)
data_NAC_p2 <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
data_NAC_p2$trtmt <- codes %>% 
  filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
  pull(trtmt)
data_NAC_p2$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))

# Bind all other data frames
for (d in simu_dirs[2:length(simu_dirs)]) {
  file.name <- "log_cellcounts.csv"
  file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
  if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
  data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
  data_tmp$trtmt <- codes %>% 
    filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
    pull(trtmt)
  data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
  data_NAC_p2 <- bind_rows(data_NAC_p2,data_tmp)
}

data_NAC <- bind_rows(data_NAC, data_NAC_p2)

# Rename columns
data_NAC <- data_NAC %>% dplyr::rename(necrotic = "necrotic_total",
                                       healthy = "hepatocytes_total",
                                       senescent = "senescent_total",
                                       Kupffer = "kupffer_total",
                                       macrophages = "macrophages_total",
                                       Time_h = "real_time") %>% 
  mutate(healthy = healthy - senescent)

# Pivot to long format
dl_NAC <- pivot_longer(data_NAC, c(necrotic,healthy,senescent,Kupffer,macrophages), names_to = "Cells", values_to = "Count") %>% 
  mutate(trtmt = as.factor(trtmt))

# Calculate mean and sd
dl_NAC_sum <- dl_NAC %>% group_by(time, Time_h, trtmt, Cells) %>% summarise(mean_count = mean(Count),
                                                                            sd_count = sd(Count)) %>% 
  ungroup()

# colors <- c(viridisLite::viridis(3))
# names(colors) <- c("necrotic","senescent", "healthy")
# 
# pchs <- c(15,16,17)
# names(pchs) <- c("necrotic","senescent", "healthy")
# 
ggplot() +
  geom_point(data = dl_NAC_sum %>% filter(Cells %in% c("necrotic","senescent","healthy"), Time_h == 24), 
             aes(x= trtmt, y = mean_count, color = Cells, shape = Cells), size = 1) +
  geom_errorbar(data = dl_NAC_sum %>% filter(Cells %in% c("necrotic","senescent","healthy"), Time_h == 24),
                aes(x= trtmt, ymin = mean_count-sd_count, ymax =  mean_count+sd_count, color = Cells), width = 0.3,show.legend=FALSE) +
  geom_hline(data = dl_control_filter %>% filter(Cells %in% c("necrotic","senescent","healthy"), Time_h == 24),
             aes(yintercept = mean_count, color = Cells), linetype = 'dotted', alpha = 0.5) +
  scale_color_viridis_d() +
  labs(x = "Time of treatment (h)", y = "Number of cells", color = "Cell type", shape = 'Cell type') +
  ggtitle("NAC") +
  theme_classic() + 
  ylim(0,600) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("Figures/","Fig6E_NAC_24h.pdf"), width = 3, height = 2)

ggplot() +
  geom_point(data = dl_NAC_sum %>% filter(Cells %in% c("necrotic","senescent","healthy"), Time_h == 168), 
             aes(x= trtmt, y = mean_count, color = Cells, shape = Cells), size = 1) +
  geom_errorbar(data = dl_NAC_sum %>% filter(Cells %in% c("necrotic","senescent","healthy"), Time_h == 168),
                aes(x= trtmt, ymin = mean_count-sd_count, ymax =  mean_count+sd_count, color = Cells), width = 0.3,show.legend=FALSE) +
  geom_hline(data = dl_control_filter %>% filter(Cells %in% c("necrotic","senescent","healthy"), Time_h == 168),
             aes(yintercept = mean_count, color = Cells), linetype = 'dotted', alpha = 0.5) +
  scale_color_viridis_d() +
  labs(x = "Time of treatment (h)", y = "Number of cells", color = "Cell type", shape = 'Cell type') +
  ggtitle("NAC") +
  theme_classic() + 
  ylim(0,600) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("Figures/","Fig6E_NAC_168h.pdf"), width = 3, height = 2)

### ### ### ### ### 
###  Figure 6A  ### 
### ### ### ### ### 

# get directory names
model_dir <- "16_Mouseliver_NAC_1-5_10x/"
simu_dirs <- dir(paste0(model_dir,"data/"))

# Define folder code
codes <- data.frame(code = c("0001", "0002","0003","0004", "0005"),
                    trtmt = c(3,6,12,24,48))

# Intracellular concentrations
file.name <- "log_SingleCellConc.csv"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")
if (!file.exists(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))) {gunzip(file.path, remove=FALSE)}

# Bind all other data frames
if (interactive()) {
  if (askYesNo("Do you want to load ALL data files?")) {
    print("OK!")
    
    # Get first data frame
    data <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
    data$trtmt <- codes %>% 
      filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
      pull(trtmt)
    data$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))
    
    for (d in simu_dirs[2:length(simu_dirs)]) {
      file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
      if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
      data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
      data_tmp$trtmt <- codes %>% 
        filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
        pull(trtmt)
      data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
      data <- bind_rows(data,data_tmp)}
    
    if (interactive()) {
      if (askYesNo("Do you want to save these data? This may overwrite existing data.")) {
        print("OK, saving...")
        save(data, file = paste0(model_dir,"data_log_SingleCellConc.RData"))
        print("Saved!")}}
  } else {
    if (interactive()) {
      if (askYesNo("Do you want to load these data from file?")) {
        load(paste0(model_dir,"data_log_SingleCellConc.RData"))}}}
}

data_p1 <- data

# get directory names
model_dir <- "16_Mouseliver_NAC_1-1_10x/"
simu_dirs <- dir(paste0(model_dir,"data/"))

# Define folder code
codes <- data.frame(code = c("0001"),
                    trtmt = c(18))

# Intracellular concentrations
file.name <- "log_SingleCellConc.csv"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")
if (!file.exists(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))) {gunzip(file.path, remove=FALSE)}

# Bind all other data frames
if (interactive()) {
  if (askYesNo("Do you want to load ALL data files?")) {
    print("OK!")
    
    # Get first data frame
    data <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
    data$trtmt <- codes %>% 
      filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
      pull(trtmt)
    data$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))
    
    for (d in simu_dirs[2:length(simu_dirs)]) {
      file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
      if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
      data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
      data_tmp$trtmt <- codes %>% 
        filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
        pull(trtmt)
      data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
      data <- bind_rows(data,data_tmp)}
    
    if (interactive()) {
      if (askYesNo("Do you want to save these data? This may overwrite existing data.")) {
        print("OK, saving...")
        save(data, file = paste0(model_dir,"data_log_SingleCellConc.RData"))
        print("Saved!")}}
  } else {
    if (interactive()) {
      if (askYesNo("Do you want to load these data from file?")) {
        load(paste0(model_dir,"data_log_SingleCellConc.RData"))}}}
}

data <- bind_rows(data_p1, data)

# summarise
data_sum <- data %>% filter(celltype == 1) %>% 
  mutate(trtmt = as.factor(trtmt)) %>%
  group_by(trtmt,run, time, real_time) %>%
  summarize(mean = mean(G_cell)) %>% ungroup %>%
  group_by(trtmt, time, real_time) %>%
  summarize(G_mean = mean(mean),
            G_sd = sd(mean))

ggplot(data_sum %>% filter(real_time >= 0 & real_time <= 100),
       aes(x = real_time, y= G_mean, color = trtmt)) +
  geom_line() +
  #geom_errorbar(aes(ymin = G_mean-G_sd, ymax = G_mean+G_sd)) +
  scale_color_viridis_d(name = "Time of\ntreatment (h)") +
  xlab("Time (h)") + ylab("GSH (a.u.)") +
  ggtitle("NAC") +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("Figures/","Fig6A_NAC_GSH.pdf"), width = 3, height = 2)

### ### ### ### ### 
###   Figure 6F ### 
### ### ### ### ### 

# get directory names
#model_dir <- "17_Mouseliver_4MP_1-6_10x/"
model_dir <- "17_Mouseliver_4MP_v2revision_1-6_10x/"
simu_dirs <- dir(paste0(model_dir,"data/"))

# Define folder code
codes <- data.frame(code = c("0001", "0002","0003","0004", "0005", "0006"),
                    trtmt = c(3,6,12,18,24,48))

file.name <- "log_cellcounts.csv"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")

# Get first data frame
gunzip(file.path, remove=FALSE)
data_4MP <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
data_4MP$trtmt <- codes %>% 
  filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
  pull(trtmt)
data_4MP$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))

# Bind all other data frames
for (d in simu_dirs[2:length(simu_dirs)]) {
  file.name <- "log_cellcounts.csv"
  file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
  if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
  data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
  data_tmp$trtmt <- codes %>% 
    filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
    pull(trtmt)
  data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
  data_4MP <- bind_rows(data_4MP,data_tmp)
}

# Rename columns
data_4MP <- data_4MP %>% dplyr::rename(necrotic = "necrotic_total",
                                       healthy = "hepatocytes_total",
                                       senescent = "senescent_total",
                                       Kupffer = "kupffer_total",
                                       macrophages = "macrophages_total",
                                       Time_h = "real_time") %>% 
  mutate(healthy = healthy - senescent)

# Pivot to long format
dl_4MP <- pivot_longer(data_4MP, c(necrotic,healthy,senescent,Kupffer,macrophages), names_to = "Cells", values_to = "Count") %>% 
  mutate(trtmt = as.factor(trtmt))

# Calculate mean and sd
dl_4MP_sum <- dl_4MP %>% group_by(time, Time_h, trtmt, Cells) %>% summarise(mean_count = mean(Count),
                                                                            sd_count = sd(Count)) %>% 
  ungroup()

ggplot() +
  geom_point(data = dl_4MP_sum %>% filter(Cells %in% c("necrotic","senescent","healthy"), Time_h == 24), 
             aes(x= trtmt, y = mean_count, color = Cells, shape = Cells), size = 1) +
  geom_errorbar(data = dl_4MP_sum %>% filter(Cells %in% c("necrotic","senescent","healthy"), Time_h == 24),
                aes(x= trtmt, ymin = mean_count-sd_count, ymax =  mean_count+sd_count, color = Cells), width = 0.3,show.legend=FALSE) +
  geom_hline(data = dl_control_filter %>% filter(Cells %in% c("necrotic","senescent","healthy"), Time_h == 24),
             aes(yintercept = mean_count, color = Cells), linetype = 'dotted', alpha = 0.5) +
  scale_color_viridis_d() +
  labs(x = "Time of treatment (h)", y = "Number of cells", color = "Cell type", shape = 'Cell type') +
  ggtitle("4MP") +
  theme_classic() + 
  ylim(0,600) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("Figures/","Fig6F_4MP_24h.pdf"), width = 3, height = 2)

ggplot() +
  geom_point(data = dl_4MP_sum %>% filter(Cells %in% c("necrotic","senescent","healthy"), Time_h == 168), 
             aes(x= trtmt, y = mean_count, color = Cells, shape = Cells), size = 1) +
  geom_errorbar(data = dl_4MP_sum %>% filter(Cells %in% c("necrotic","senescent","healthy"), Time_h == 168),
                aes(x= trtmt, ymin = mean_count-sd_count, ymax =  mean_count+sd_count, color = Cells), width = 0.3,show.legend=FALSE) +
  geom_hline(data = dl_control_filter %>% filter(Cells %in% c("necrotic","senescent","healthy"), Time_h == 168),
             aes(yintercept = mean_count, color = Cells), linetype = 'dotted', alpha = 0.5) +
  scale_color_viridis_d() +
  labs(x = "Time of treatment (h)", y = "Number of cells", color = "Cell type", shape = 'Cell type') +
  ggtitle("4MP") +
  ylim(0,600) +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("Figures/","Fig6F_4MP_168h.pdf"), width = 3, height = 2)

### ### ### ### ### 
###  Figure 6B  ### 
### ### ### ### ### 

# get directory names
model_dir <- "17_Mouseliver_4MP_v2revision_1-6_10x/"
simu_dirs <- dir(paste0(model_dir,"data/"))

# Define folder code
codes <- data.frame(code = c("0001", "0002","0003","0004", "0005", "0006"),
                    trtmt = c(3,6,12,18,24,48))

# Intracellular concentrations
file.name <- "log_SingleCellConc.csv"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")
if (!file.exists(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))) {gunzip(file.path, remove=FALSE)}

# Bind all other data frames
if (interactive()) {
  if (askYesNo("Do you want to load ALL data files?")) {
    print("OK!")
    
    # Get first data frame
    data <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
    data$trtmt <- codes %>% 
      filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
      pull(trtmt)
    data$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))
    
    for (d in simu_dirs[2:length(simu_dirs)]) {
      file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
      if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
      data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
      data_tmp$trtmt <- codes %>% 
        filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
        pull(trtmt)
      data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
      data <- bind_rows(data,data_tmp)}
    
    if (interactive()) {
      if (askYesNo("Do you want to save these data? This may overwrite existing data.")) {
        print("OK, saving...")
        save(data, file = paste0(model_dir,"data_log_SingleCellConc.RData"))
        print("Saved!")}}
  } else {
    if (interactive()) {
      if (askYesNo("Do you want to load these data from file?")) {
        load(paste0(model_dir,"data_log_SingleCellConc.RData"))}}}
}

# summarise
data_sum <- data %>% filter(celltype == 1) %>% 
  mutate(trtmt = as.factor(trtmt)) %>%
  group_by(trtmt,run, time, real_time) %>%
  summarize(mean = mean(C)) %>% ungroup %>%
  group_by(trtmt, time, real_time) %>%
  summarize(C_mean = mean(mean),
            C_sd = sd(mean))

ggplot(data_sum %>% filter(real_time >= 0 & real_time <= 100),
       aes(x = real_time, y= C_mean, color = trtmt)) +
  geom_line() +
  #geom_errorbar(aes(ymin = C_mean-C_sd, ymax = C_mean+C_sd)) +
  scale_color_viridis_d(name = "Time of\ntreatment (h)") +
  xlab("Time (h)") + ylab("NAPQI-cys (a.u.)") +
  ggtitle("4MP") +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("Figures/","Fig6B_4MP_C.pdf"), width = 3, height = 2)

### ### ### ### ###
###   Figure 6G ###
### ### ### ### ###

# get directory names
model_dir <- "18_Mouseliver_PFTa_1-6_10x/"
simu_dirs <- dir(paste0(model_dir,"data/"))

# Define folder code
codes <- data.frame(code = c("0001", "0002","0003","0004", "0005", "0006"),
                    trtmt = c(3,6,12,18,24,48))

file.name <- "log_cellcounts.csv"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")

# Get first data frame
gunzip(file.path, remove=FALSE)
data_p53inh <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
data_p53inh$trtmt <- codes %>%
  filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>%
  pull(trtmt)
data_p53inh$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))

# Bind all other data frames
for (d in simu_dirs[2:length(simu_dirs)]) {
  file.name <- "log_cellcounts.csv"
  file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
  if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
  data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
  data_tmp$trtmt <- codes %>%
    filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>%
    pull(trtmt)
  data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
  data_p53inh <- bind_rows(data_p53inh,data_tmp)
}

# Rename columns
data_p53inh <- data_p53inh %>% dplyr::rename(necrotic = "necrotic_total",
                                             healthy = "hepatocytes_total",
                                             senescent = "senescent_total",
                                             Kupffer = "kupffer_total",
                                             macrophages = "macrophages_total",
                                             Time_h = "real_time") %>%
  mutate(healthy = healthy - senescent)

# Pivot to long format
dl_p53inh <- pivot_longer(data_p53inh, c(necrotic,healthy,senescent,Kupffer,macrophages), names_to = "Cells", values_to = "Count") %>%
  mutate(trtmt = as.factor(trtmt))

# Calculate mean and sd
dl_p53inh_sum <- dl_p53inh %>% group_by(time, Time_h, trtmt, Cells) %>% summarise(mean_count = mean(Count),
                                                                                  sd_count = sd(Count)) %>%
  ungroup()

ggplot() +
  geom_point(data = dl_p53inh_sum %>% filter(Cells %in% c("necrotic","senescent","healthy"), Time_h == 24),
             aes(x= trtmt, y = mean_count, color = Cells, shape = Cells), size = 1) +
  geom_errorbar(data = dl_p53inh_sum %>% filter(Cells %in% c("necrotic","senescent","healthy"), Time_h == 24),
                aes(x= trtmt, ymin = mean_count-sd_count, ymax =  mean_count+sd_count, color = Cells), width = 0.3,show.legend=FALSE) +
  geom_hline(data = dl_control_filter %>% filter(Cells %in% c("necrotic","senescent","healthy"), Time_h == 24),
             aes(yintercept = mean_count, color = Cells), linetype = 'dotted', alpha = 0.5) +
  scale_color_viridis_d() +
  labs(x = "Time of treatment (h)", y = "Number of cells", color = "Cell type", shape = 'Cell type') +
  ggtitle(expression("Pifithrin-"*alpha)) +
  theme_classic() +
  ylim(0,600) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("Figures/","Fig6G_PFT_24h.pdf"), width = 3, height = 2)

ggplot() +
  geom_point(data = dl_p53inh_sum %>% filter(Cells %in% c("necrotic","senescent","healthy"), Time_h == 168),
             aes(x= trtmt, y = mean_count, color = Cells, shape = Cells), size = 1) +
  geom_errorbar(data = dl_p53inh_sum %>% filter(Cells %in% c("necrotic","senescent","healthy"), Time_h == 168),
                aes(x= trtmt, ymin = mean_count-sd_count, ymax =  mean_count+sd_count, color = Cells), width = 0.3,show.legend=FALSE) +
  geom_hline(data = dl_control_filter %>% filter(Cells %in% c("necrotic","senescent","healthy"), Time_h == 168),
             aes(yintercept = mean_count, color = Cells), linetype = 'dotted', alpha = 0.5) +
  scale_color_viridis_d() +
  labs(x = "Time of treatment (h)", y = "Number of cells", color = "Cell type", shape = 'Cell type') +
  ggtitle(expression("Pifithrin-"*alpha)) +
  theme_classic() +
  ylim(0,600) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("Figures/","Fig6G_PFT_168h.pdf"), width = 3, height = 2)

### ### ### ### ###
###  Figure 6C  ###
### ### ### ### ###

# get directory names
model_dir <- "18_Mouseliver_PFTa_1-6_10x/"
simu_dirs <- dir(paste0(model_dir,"data/"))

# Define folder code
codes <- data.frame(code = c("0001", "0002","0003","0004", "0005", "0006"),
                    trtmt = c(3,6,12,18,24,48))

# Intracellular concentrations
file.name <- "log_SingleCellConc.csv"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")
if (!file.exists(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))) {gunzip(file.path, remove=FALSE)}

# Bind all other data frames
if (interactive()) {
  if (askYesNo("Do you want to load ALL data files?")) {
    print("OK!")
    
    # Get first data frame
    data <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
    data$trtmt <- codes %>%
      filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>%
      pull(trtmt)
    data$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))
    
    for (d in simu_dirs[2:length(simu_dirs)]) {
      file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
      if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
      data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
      data_tmp$trtmt <- codes %>%
        filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>%
        pull(trtmt)
      data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
      data <- bind_rows(data,data_tmp)}
    
    if (interactive()) {
      if (askYesNo("Do you want to save these data? This may overwrite existing data.")) {
        print("OK, saving...")
        save(data, file = paste0(model_dir,"data_log_SingleCellConc.RData"))
        print("Saved!")}}
  } else {
    if (interactive()) {
      if (askYesNo("Do you want to load these data from file?")) {
        load(paste0(model_dir,"data_log_SingleCellConc.RData"))}}}
}

# summarise
data_sum1 <- data %>% filter(celltype == 1) %>%
  mutate(trtmt = as.factor(trtmt),
         p53_total = p53 + p53p) %>%
  group_by(trtmt,run, time, real_time) %>%
  summarize(mean_p53 = mean(p53),
            mean_p53p = mean(p53p),
            mean_p53_total = mean(p53_total),
            mean_p21 = mean(p21),
            mean_MDM2 = mean(MDM2),
            sd_p53 = sd(p53),
            sd_p53p = sd(p53p),
            sd_p53_total = sd(p53_total),
            sd_p21 = sd(p21),
            sd_MDM2 = sd(MDM2)) %>% ungroup
data_sum <-  data_sum1  %>% group_by(trtmt, time, real_time) %>%
  summarize(p53_mean = mean(mean_p53),
            p53_sd = sd(mean_p53),
            p53p_mean = mean(mean_p53p),
            p53p_sd = sd(mean_p53p),
            p53_total_mean = mean(mean_p53_total),
            p53_total_sd = sd(mean_p53_total),
            p21_mean = mean(mean_p21),
            p21_sd = sd(mean_p21),
            MDM2_mean = mean(mean_MDM2),
            MDM2_sd = sd(mean_MDM2))

# Check out variability between single cells
tmp <- data %>% filter(run == 1, trtmt == 48, celltype == 1, cell.id %in% sample(0:600,20))
ggplot(tmp,
       aes(x = real_time, y= p21, color = cell.id)) +
  #geom_line() +
  geom_point(size = 0.2) +
  #geom_errorbar(aes(ymin = mean_p21-sd_p21, ymax = mean_p21+sd_p21)) +
  scale_color_viridis_c(name = "Time of\ntreatment (h)") +
  xlab("Time (h)") + ylab("p21 (a.u.)") +
  ggtitle(expression("Pifithrin-"*alpha)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'none')

ggplot(data_sum1 %>% filter(real_time >= 0 & real_time <= 48),
       aes(x = real_time, y= mean_p21, color = trtmt)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean_p21-sd_p21, ymax = mean_p21+sd_p21)) +
  scale_color_viridis_d(name = "Time of\ntreatment (h)") +
  xlab("Time (h)") + ylab("p21 (a.u.)") +
  ggtitle(expression("Pifithrin-"*alpha)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_grid(trtmt ~ run)

ggplot(data_sum %>% filter(real_time >= 0 & real_time <= 100),
       aes(x = real_time, y= p53p_mean, color = trtmt)) +
  geom_line() +
  #geom_point() +
  #geom_errorbar(aes(ymin = p53p_mean-p53p_sd, ymax = p53p_mean+p53p_sd)) +
  scale_color_viridis_d(name = "Time of\ntreatment (h)") +
  xlab("Time (h)") + ylab("phospho-p53 (a.u.)") +
  ggtitle(expression("Pifithrin-"*alpha)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("Figures/","Fig6C_PFT_p53p.pdf"), width = 3, height = 2)

ggplot(data_sum %>% filter(real_time >= 0 & real_time <= 100),
       aes(x = real_time, y= p53_mean, color = trtmt)) +
  geom_line() +
  #geom_errorbar(aes(ymin = p53p_mean-p53p_sd, ymax = p53p_mean+p53p_sd)) +
  scale_color_viridis_d(name = "Time of\ntreatment (h)") +
  xlab("Time (h)") + ylab("p53 (a.u.)") +
  ggtitle(expression("Pifithrin-"*alpha)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("Figures/","FigS3A_PFT_p53.pdf"), width = 3, height = 2)

ggplot(data_sum %>% filter(real_time >= 0 & real_time <= 100),
       aes(x = real_time, y= p53_total_mean, color = trtmt)) +
  geom_line() +
  #geom_errorbar(aes(ymin = p53_total_mean-p53_total_sd, ymax = p53_total_mean+p53_total_sd)) +
  scale_color_viridis_d(name = "Time of\ntreatment (h)") +
  xlab("Time (h)") + ylab("Total p53 (a.u.)") +
  ggtitle(expression("Pifithrin-"*alpha)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("Figures/","FigS3B_PFT_p53_total.pdf"), width = 3, height = 2)

ggplot(data_sum %>% filter(real_time >= 0 & real_time <= 100),
       aes(x = real_time, y= MDM2_mean, color = trtmt)) +
  geom_line() +
  geom_errorbar(aes(ymin = MDM2_mean-MDM2_sd, ymax = MDM2_mean+MDM2_sd), size = 0.2) +
  scale_color_viridis_d(name = "Time of\ntreatment (h)") +
  xlab("Time (h)") + ylab("MDM2 (a.u.)") +
  ggtitle(expression("Pifithrin-"*alpha)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("Figures/","FigS3C_PFT_MDM2.pdf"), width = 4, height = 2)

ggplot(data_sum %>% filter(real_time >= 0 & real_time <= 100),
       aes(x = real_time, y= p21_mean, color = trtmt)) +
  geom_line() +
  geom_errorbar(aes(ymin = p21_mean-p21_sd, ymax = p21_mean+p21_sd), size = 0.2) +
  scale_color_viridis_d(name = "Time of\ntreatment (h)") +
  xlab("Time (h)") + ylab("p21 (a.u.)") +
  ggtitle(expression("Pifithrin-"*alpha)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("Figures/","FigS3D_PFT_p21.pdf"), width = 4, height = 2)

### ### ### ### ### 
###   Figure 6H ### 
### ### ### ### ### 

# get directory names
model_dir <- "19_Mouseliver_CDN_1-6_10x/"
simu_dirs <- dir(paste0(model_dir,"data/"))

# Define folder code
codes <- data.frame(code = c("0001", "0002","0003","0004", "0005", "0006"),
                    trtmt = c(3,6,12,18,24,48))

file.name <- "log_cellcounts.csv"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")

# Get first data frame
gunzip(file.path, remove=FALSE)
data_cdn <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
data_cdn$trtmt <- codes %>% 
  filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
  pull(trtmt)
data_cdn$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))

# Bind all other data frames
for (d in simu_dirs[2:length(simu_dirs)]) {
  file.name <- "log_cellcounts.csv"
  file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
  if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
  data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
  data_tmp$trtmt <- codes %>% 
    filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
    pull(trtmt)
  data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
  data_cdn <- bind_rows(data_cdn,data_tmp)
}

# Rename columns
data_cdn <- data_cdn %>% dplyr::rename(necrotic = "necrotic_total",
                                       healthy = "hepatocytes_total",
                                       senescent = "senescent_total",
                                       Kupffer = "kupffer_total",
                                       macrophages = "macrophages_total",
                                       Time_h = "real_time") %>% 
  mutate(healthy = healthy - senescent)

# Pivot to long format
dl_cdn <- pivot_longer(data_cdn, c(necrotic,healthy,senescent,Kupffer,macrophages), names_to = "Cells", values_to = "Count") %>% 
  mutate(trtmt = as.factor(trtmt))

# Calculate mean and sd
dl_cdn_sum <- dl_cdn %>% group_by(time, Time_h, trtmt, Cells) %>% summarise(mean_count = mean(Count),
                                                                            sd_count = sd(Count)) %>% 
  ungroup()

ggplot() +
  geom_point(data = dl_cdn_sum %>% filter(Cells %in% c("necrotic","senescent","healthy"), Time_h == 24), 
             aes(x= trtmt, y = mean_count, color = Cells, shape = Cells), size = 1) +
  geom_errorbar(data = dl_cdn_sum %>% filter(Cells %in% c("necrotic","senescent","healthy"), Time_h == 24),
                aes(x= trtmt, ymin = mean_count-sd_count, ymax =  mean_count+sd_count, color = Cells), width = 0.3,show.legend=FALSE) +
  geom_hline(data = dl_control_filter %>% filter(Cells %in% c("necrotic","senescent","healthy"), Time_h == 24),
             aes(yintercept = mean_count, color = Cells), linetype = 'dotted', alpha = 0.5) +
  scale_color_viridis_d() +
  labs(x = "Time of treatment (h)", y = "Number of cells", color = "Cell type", shape = 'Cell type') +
  ggtitle("Clodronate") +
  theme_classic() + 
  ylim(0,600) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("Figures/","Fig6H_CDN_24h.pdf"), width = 3, height = 2)

ggplot() +
  geom_point(data = dl_cdn_sum %>% filter(Cells %in% c("necrotic","senescent","healthy"), Time_h == 168), 
             aes(x= trtmt, y = mean_count, color = Cells, shape = Cells), size = 1) +
  geom_errorbar(data = dl_cdn_sum %>% filter(Cells %in% c("necrotic","senescent","healthy"), Time_h == 168),
                aes(x= trtmt, ymin = mean_count-sd_count, ymax =  mean_count+sd_count, color = Cells), width = 0.3,show.legend=FALSE) +
  geom_hline(data = dl_control_filter %>% filter(Cells %in% c("necrotic","senescent","healthy"), Time_h == 168),
             aes(yintercept = mean_count, color = Cells), linetype = 'dotted', alpha = 0.5) +
  scale_color_viridis_d() +
  labs(x = "Time of treatment (h)", y = "Number of cells", color = "Cell type", shape = 'Cell type') +
  ggtitle("Clodronate") +
  theme_classic() + 
  ylim(0,600) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("Figures/","Fig6H_CDN_168h.pdf"), width = 3, height = 2)

ggplot() +
  geom_line(data = dl_cdn_sum %>% filter(Cells %in% c("macrophages"),Time_h >= 0 & Time_h <= 100), 
            aes(x= Time_h, y = mean_count, color = trtmt), size = 1) +
  # geom_errorbar(data = dl_cdn_sum %>% filter(Cells %in% c("macrophages")),
  #               aes(x= Time_h, ymin = mean_count-sd_count, ymax =  mean_count+sd_count, color = trtmt), width = 0.3,show.legend=FALSE) +
  scale_color_viridis_d() +
  labs(x = "Time (h)", y = "Number of macrophages", color = "Time of\ntreatment (h)") +
  ggtitle("Clodronate") +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("Figures/","Fig6D_CDN_macrophages.pdf"), width = 3, height = 2)

### ### ### ### ### 
###   Figure p21 expr ### 
### ### ### ### ### 

# get directory names
model_dir <- "19_Mouseliver_CDN_1-6_10x/"
simu_dirs <- dir(paste0(model_dir,"data/"))

# Define folder code
codes <- data.frame(code = c("0001", "0002","0003","0004", "0005", "0006"),
                    trtmt = c(3,6,12,18,24,48))

# Intracellular concentrations
file.name <- "log_SingleCellConc.csv"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")
if (!file.exists(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))) {gunzip(file.path, remove=FALSE)}

# Bind all other data frames
if (interactive()) {
  if (askYesNo("Do you want to load ALL data files?")) {
    print("OK!")
    
    # Get first data frame
    data <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
    data$trtmt <- codes %>%
      filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>%
      pull(trtmt)
    data$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))
    
    for (d in simu_dirs[2:length(simu_dirs)]) {
      file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
      if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
      data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
      data_tmp$trtmt <- codes %>%
        filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>%
        pull(trtmt)
      data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
      data <- bind_rows(data,data_tmp)}
    
    if (interactive()) {
      if (askYesNo("Do you want to save these data? This may overwrite existing data.")) {
        print("OK, saving...")
        save(data, file = paste0(model_dir,"data_log_SingleCellConc.RData"))
        print("Saved!")}}
  } else {
    if (interactive()) {
      if (askYesNo("Do you want to load these data from file?")) {
        load(paste0(model_dir,"data_log_SingleCellConc.RData"))}}}
}

# summarise
data_sum <- data %>% filter(celltype == 1) %>%
  mutate(trtmt = as.factor(trtmt),
         p53_total = p53 + p53p) %>%
  group_by(trtmt,run, time, real_time) %>%
  summarize(mean_p53 = mean(p53),
            mean_p53p = mean(p53p),
            mean_p53_total = mean(p53_total),
            mean_p21 = mean(p21),
            mean_MDM2 = mean(MDM2)) %>% ungroup %>%
  group_by(trtmt, time, real_time) %>%
  summarize(p53_mean = mean(mean_p53),
            p53_sd = sd(mean_p53),
            p53p_mean = mean(mean_p53p),
            p53p_sd = sd(mean_p53p),
            p53_total_mean = mean(mean_p53_total),
            p53_total_sd = sd(mean_p53_total),
            p21_mean = mean(mean_p21),
            p21_sd = sd(mean_p21),
            MDM2_mean = mean(mean_MDM2),
            MDM2_sd = sd(mean_MDM2))

ggplot(data_sum %>% filter(real_time >= 0 & real_time <= 100),
       aes(x = real_time, y= p21_mean, color = trtmt)) +
  geom_line() +
  geom_errorbar(aes(ymin = p21_mean-p21_sd, ymax = p21_mean+p21_sd)) +
  scale_color_viridis_d(name = "Time of\ntreatment (h)") +
  xlab("Time (h)") + ylab("p21 (a.u.)") +
  ggtitle(expression("Clodronate")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("Figures/","FigSX_CDN_p21.pdf"), width = 3.5, height = 2)


### ### ### ### ### 
###  Figure SX  ### 
### ### ### ### ### 

# get directory names
model_dir <- "20_Mouseliver_NecrProb_1-5_10x/"
simu_dirs <- dir(paste0(model_dir,"data/"))

# Cell counts
file.name <- "log_cellcounts.csv"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")
if (!file.exists(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))) {gunzip(file.path, remove=FALSE)}

# Define folder code
codes <- data.frame(code = c("0001", "0002","0003","0004", "0005"),
                    factor = c(0.00001,0.0001,0.001,0.01,0.1))

# Get first data frame
data_NecProb <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
data_NecProb$factor <- codes %>% 
  filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
  pull(factor)
data_NecProb$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))

# Bind all other data frames
for (d in simu_dirs[2:length(simu_dirs)]) {
  file.name <- "log_cellcounts.csv"
  file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
  if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
  data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
  data_tmp$factor <- codes %>% 
    filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
    pull(factor)
  data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
  data_NecProb <- bind_rows(data_NecProb,data_tmp)
}

# Rename columns
data_NecProb <- data_NecProb %>% dplyr::rename(necrotic = "necrotic_total",
                                               healthy = "hepatocytes_total",
                                               senescent = "senescent_total",
                                               Kupffer = "kupffer_total",
                                               macrophages = "macrophages_total",
                                               Time_h = "real_time") %>% 
  mutate(healthy = healthy - senescent)

# Pivot to long format
# Pivot to long format
dl_NecProb <- pivot_longer(data_NecProb, c(necrotic,healthy,senescent,Kupffer,macrophages), names_to = "Cells", values_to = "Count") %>% 
  mutate(factor = as.factor(factor))

# Calculate mean and sd
dl_NecProb_sum <- dl_NecProb %>% group_by(time, Time_h, factor, Cells) %>% summarise(mean_count = mean(Count),
                                                                     sd_count = sd(Count)) %>% 
  ungroup()

# Make plots
exp_names <- c("1e-05" = "p[necrosis]==1e-5",
               "1e-04" = "p[necrosis]==1e-4",
               "0.001" = "p[necrosis]==0.001",
               "0.01" =  "p[necrosis] == 0.01",
               "0.1" =  "p[necrosis] == 0.1")

colors <- c("grey",viridisLite::viridis(4))
names(colors) <- c("total","healthy", "macrophages","necrotic", "senescent")

ggplot(dl_NecProb_sum %>% filter(!Cells == "Kupffer", !factor == 10), aes(x = Time_h, color = Cells, linetype = Cells)) + 
  geom_line(aes(y = mean_count), size = 0.5) +
  geom_errorbar(aes(ymin= mean_count-sd_count, ymax = mean_count+sd_count), width = 0.1, alpha = 0.1) +
  scale_linetype_manual(name = "Cell Type", values = c("total" = 2, "healthy" = 1, "macrophages" = 1,
                                                       "necrotic" = 1, "senescent" = 1)) +
  scale_color_manual(name = "Cell Type", values = colors) +
  labs(x = "Time (h)", y = "Cell count") + 
  theme_classic() +
  facet_grid(cols = vars(factor),labeller = as_labeller(exp_names,default = label_parsed))
ggsave(paste0("Figures/","FigS1B_cell_counts_NecrosisProb.pdf"), width = 8.5, height = 1.5)

### ### ### ### ### 
###  Figure SX  ### 
### ### ### ### ### 

# get directory names
model_dir <- "21_Mouseliver_DAMPdeg_1-5_10x/"
simu_dirs <- dir(paste0(model_dir,"data/"))

# Cell counts
file.name <- "log_cellcounts.csv"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")
if (!file.exists(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))) {gunzip(file.path, remove=FALSE)}

# Define folder code
codes <- data.frame(code = c("0001", "0002","0003","0004", "0005"),
                    factor = c(0.001,0.01,0.1,1,10))

# Get first data frame
data_DAMPdeg <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
data_DAMPdeg$factor <- codes %>% 
  filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
  pull(factor)
data_DAMPdeg$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))

# Bind all other data frames
for (d in simu_dirs[2:length(simu_dirs)]) {
  file.name <- "log_cellcounts.csv"
  file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
  if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
  data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
  data_tmp$factor <- codes %>% 
    filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
    pull(factor)
  data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
  data_DAMPdeg <- bind_rows(data_DAMPdeg,data_tmp)
}

# Rename columns
data_DAMPdeg <- data_DAMPdeg %>% dplyr::rename(necrotic = "necrotic_total",
                                               healthy = "hepatocytes_total",
                                               senescent = "senescent_total",
                                               Kupffer = "kupffer_total",
                                               macrophages = "macrophages_total",
                                               Time_h = "real_time") %>% 
  mutate(healthy = healthy - senescent)

# Pivot to long format
# Pivot to long format
dl_DAMPdeg <- pivot_longer(data_DAMPdeg, c(necrotic,healthy,senescent,Kupffer,macrophages), names_to = "Cells", values_to = "Count") %>% 
  mutate(factor = as.factor(factor))

# Calculate mean and sd
dl_DAMPdeg_sum <- dl_DAMPdeg %>% group_by(time, Time_h, factor, Cells) %>% summarise(mean_count = mean(Count),
                                                                     sd_count = sd(Count)) %>% 
  ungroup()

# Make plots
exp_names <- c("0.001" = "d[D]==0.012*~h^{-1}",
               "0.01" =  "d[D] == 0.12*~h^{-1}",
               "0.1" =  "d[D] == 1.2*~h^{-1}",
               "1" =  "d[D] == 12*~h^{-1}")

colors <- c("grey",viridisLite::viridis(4))
names(colors) <- c("total","healthy", "macrophages","necrotic", "senescent")

ggplot(dl_DAMPdeg_sum %>% filter(!Cells == "Kupffer", !factor == 10), aes(x = Time_h, color = Cells, linetype = Cells)) + 
  geom_line(aes(y = mean_count), size = 0.5) +
  geom_errorbar(aes(ymin= mean_count-sd_count, ymax = mean_count+sd_count), width = 0.1, alpha = 0.1) +
  scale_linetype_manual(name = "Cell Type", values = c("total" = 2, "healthy" = 1, "macrophages" = 1,
                                                       "necrotic" = 1, "senescent" = 1)) +
  scale_color_manual(name = "Cell Type", values = colors) +
  labs(x = "Time (h)", y = "Cell count") + 
  theme_classic() +
  facet_grid(cols = vars(factor),labeller = as_labeller(exp_names,default = label_parsed))
ggsave(paste0("Figures/","FigS1C_cell_counts_DAMPdeg.pdf"), width = 6.5, height = 1.5)

### ### ### ### ### 
###  Figure SX  ### 
### ### ### ### ### 

# get directory names
model_dir <- "22_Mouseliver_MCP1deg_1-5_10x/"
simu_dirs <- dir(paste0(model_dir,"data/"))

# Cell counts
file.name <- "log_cellcounts.csv"
file.path <- paste0(model_dir,"data/",simu_dirs[1],"/",file.name,".gz")
if (!file.exists(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))) {gunzip(file.path, remove=FALSE)}

# Define folder code
codes <- data.frame(code = c("0001", "0002","0003","0004", "0005"),
                    factor = c(0.002,0.02,0.2,2,20))

# Get first data frame
data_MCP1deg <- read_tsv(paste0(model_dir,"data/",simu_dirs[1],"/",file.name))
data_MCP1deg$factor <- codes %>% 
  filter(code == regmatches(simu_dirs[1],regexpr("[0-9]{4}",simu_dirs[1]))) %>% 
  pull(factor)
data_MCP1deg$run <- sub('.', '', regmatches(simu_dirs[1],regexpr("-[0-9]{1,3}",simu_dirs[1])))

# Bind all other data frames
for (d in simu_dirs[2:length(simu_dirs)]) {
  file.name <- "log_cellcounts.csv"
  file.path <- paste0(model_dir,"data/",d,"/",file.name,".gz")
  if (!file.exists(paste0(model_dir,"data/",d,"/",file.name))) {gunzip(file.path, remove=FALSE)}
  data_tmp <- read_tsv(paste0(model_dir,"data/",d,"/",file.name))
  data_tmp$factor <- codes %>% 
    filter(code == regmatches(d,regexpr("[0-9]{4}",d))) %>% 
    pull(factor)
  data_tmp$run <- sub('.', '', regmatches(d,regexpr("-[0-9]{1,3}",d)))
  data_MCP1deg <- bind_rows(data_MCP1deg,data_tmp)
}

# Rename columns
data_MCP1deg <- data_MCP1deg %>% dplyr::rename(necrotic = "necrotic_total",
                                               healthy = "hepatocytes_total",
                                               senescent = "senescent_total",
                                               Kupffer = "kupffer_total",
                                               macrophages = "macrophages_total",
                                               Time_h = "real_time") %>% 
  mutate(healthy = healthy - senescent)

# Pivot to long format
# Pivot to long format
dl_MCP1deg <- pivot_longer(data_MCP1deg, c(necrotic,healthy,senescent,Kupffer,macrophages), names_to = "Cells", values_to = "Count") %>% 
  mutate(factor = as.factor(factor))

# Calculate mean and sd
dl_MCP1deg_sum <- dl_MCP1deg %>% group_by(time, Time_h, factor, Cells) %>% summarise(mean_count = mean(Count),
                                                                                     sd_count = sd(Count)) %>% 
  ungroup()

# Make plots
# Make labeller
exp_names <- c("0.002" = "d[M]==0.024*~h^{-1}",
               "0.02" =  "d[M] == 0.24*~h^{-1}",
               "0.2" =  "d[M] == 2.4*~h^{-1}",
               "2" =  "d[M] == 24*~h^{-1}")

colors <- c("grey",viridisLite::viridis(4))
names(colors) <- c("total","healthy", "macrophages","necrotic", "senescent")

ggplot(dl_MCP1deg_sum %>% filter(!Cells == "Kupffer", !factor == 20), aes(x = Time_h, color = Cells, linetype = Cells)) + 
  geom_line(aes(y = mean_count), size = 0.5) +
  geom_errorbar(aes(ymin= mean_count-sd_count, ymax = mean_count+sd_count), width = 0.1, alpha = 0.1) +
  scale_linetype_manual(name = "Cell Type", values = c("total" = 2, "healthy" = 1, "macrophages" = 1,
                                                       "necrotic" = 1, "senescent" = 1)) +
  scale_color_manual(name = "Cell Type", values = colors) +
  labs(x = "Time (h)", y = "Cell count") + 
  theme_classic() +
  facet_grid(cols = vars(factor),labeller = as_labeller(exp_names,default = label_parsed))
ggsave(paste0("Figures/","FigS1D_cell_counts_MCP1deg.pdf"), width = 6.5, height = 1.5)
