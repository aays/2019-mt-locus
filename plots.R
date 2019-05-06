# mt locus paper plots

library(readr)
library(dplyr)
library(magrittr)
library(tidyr)
library(ggplot2)
library(purrr)
library(patchwork)
library(wesanderson)
library(stringr)
library(fs)
library(depmixS4)
library(viridis)
library(ggalt) # for relative annotations
library(broom)
library(car)
select <- dplyr::select

setwd('~/Desktop/Coding/GitHub/2019-mt-locus/')

# figure 1
# ZnS across the mt locus

fnames <- fs::dir_ls('data/ld-windowed', glob = '*zns*txt')
zns_files <- map_dfr(fnames, read_delim, delim = ' ', col_types = cols(), .id = 'name') %>% 
  mutate(region = str_extract(name, '(mt|plus|chromosome_6)')) %>% 
  mutate(midpoint = start + 500) %>% 
  filter(site_count != 0) %>% 
  dplyr::select(-name)

zns_plot_theme <- function(font_size = 16) {
  theme(axis.title = element_text(family = "Helvetica", size = font_size),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text = element_text(family = "Helvetica", size = font_size, color = 'black'),
        axis.line = element_line(colour = 'black', linetype = 'solid'),
        axis.line.x = element_line(size = 0.9),
        axis.line.y = element_line(size = 0.9),
        panel.background = element_blank())
}

mt_only <- zns_files %>% 
  filter(region == 'mt') %>% 
  filter(site_count > 30) %>% 
  mutate(
    domain = case_when(
      start <= 420000 ~ 'T',
      start > 420000 & start <= 826000 ~ 'R',
      start > 826000 ~ 'C'
    )
  )

chr6_only <- zns_files %>% 
  filter(site_count > 30) %>% 
  filter(region == 'chromosome_6') %>% 
  filter(start < 298000 | end > 943000)

plot_data <- bind_rows(mt_only, chr6_only) %>% 
  arrange(start)

# HMM fit
hmm_fit <- function(d, nstates) {
  model <- d %>% 
    depmix(
      list(zns ~ 1), data = .,
      family = list(gaussian()), nstates = nstates)
  model_fit <- fit(model)
  posterior_states <- posterior(model_fit)
  return(posterior_states)
}

posterior_3 <- hmm_fit(plot_data, 3)

# add HMM fit to plot_data
plot_data <- bind_cols(plot_data, dplyr::select(posterior_3, state))
plot_data$state <- as.factor(plot_data$state)

# zns across first ~3 Mb of chr6
hmm_colors <- viridis(3)
zns_plot_all <- plot_data %>% # only gametologs and chr6
  ggplot(aes(x = midpoint, y = zns, color = state)) +
  geom_rect(aes(xmin = 298000, xmax = 420000, ymin = -0.5, ymax = 1.5),  # 0-175k in bed - previously 473000
            color = '#eaecef', fill = '#eaecef') +
  geom_rect(aes(xmin = 421000, xmax = 826000, ymin = -0.5, ymax = 1.5),
            color = '#d9dce0', fill = '#d9dce0') +
  geom_rect(aes(xmin = 827000, xmax = 943000, ymin = -0.5, ymax = 1.5), 
            color = '#eaecef', fill = '#eaecef') +
  geom_point(size = 1.2) +
  labs(
    x = 'position on chromosome 6 (bp)',
    y = expression(paste(Z[nS]))
  ) +
  scale_color_manual(
    values = c(hmm_colors[3], hmm_colors[2], hmm_colors[1])) +
  scale_x_continuous(
    breaks = seq(0, 40, by = 5)*1e5, 
    labels = scales::comma,
    limits = c(0, 31.5)*1e5) +
  coord_cartesian(
    x = c(0, 30.5)*1e5, y = c(0, 1)) +
  zns_plot_theme(16) +
  annotate('text', label = 'T', x = 360000, y = 0.05, size = 8) +
  annotate('text', label = 'R', x = 638000, y = 0.05, size = 8) +
  annotate('text', label = 'C', x = 884000, y = 0.05, size = 8) +
  guides(color = FALSE) +
  NULL

zns_plot_all

ggsave('plots/fig_1.pdf', plot = zns_plot_all,
       width = par('din')[1] * 1.75, height = par('din')[1])

# figure 2 -
# Fst and ZnS over shared regions

polymorphism_shared <- read_csv('data/cds-popgen/polymorphism_all.shared.20190225.zns.csv', col_types = cols()) %>% 
  mutate(piN_piS = theta_pi_0fold_ALL / theta_pi_4fold_ALL) %>% 
  mutate(piN_piS_PLUS = theta_pi_0fold_MTPLUS / theta_pi_4fold_MTPLUS) %>% 
  mutate(piN_piS_MINUS = theta_pi_0fold_MTMINUS / theta_pi_4fold_MTMINUS) %>% 
  mutate(chrom = 'shared') # for boxplot later
polymorphism_limited <- read_csv('data/cds-popgen/polymorphism.mtLimited.zns.csv', col_types = cols()) %>% 
  mutate(piN_piS = theta_pi_0fold_ALL / theta_pi_4fold_ALL) %>% 
  mutate(chrom = 'limited')
polymorphism_autosomal <- read_csv('data/cds-popgen/polymorphism.nonMT.csv', col_types = cols()) %>% 
  separate(mtPlus_coords, into = c('chrom', 'pos'), sep = ':') %>% 
  mutate(piN_piS = theta_pi_0fold_ALL / theta_pi_4fold_ALL)
polymorphism_chr6 <- polymorphism_autosomal %>% 
  filter(chrom == 'chromosome_6') %>% 
  sample_n(100)
polymorphism_other <- polymorphism_autosomal %>% 
  filter(chrom != 'chromosome_6') %>% 
  mutate(chrom = 'other') %>% 
  sample_n(100)
windowed <- read_csv(
  'data/cds-popgen/polymorphism_all_windowed_script.gametolog.csv',
  col_types = cols()
  ) %>% 
  rename(start = start_coords)

windowed$start[1] <- 298000 # to match zns file

corr_plot_theme <- function(font_size = 16) {
  theme(axis.title = element_text(family = "Helvetica", size = font_size),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text = element_text(family = "Helvetica", size = font_size, color = 'black'),
        axis.line = element_line(colour = 'black', linetype = 'solid'),
        axis.line.x = element_line(size = 0.9),
        axis.line.y = element_line(size = 0.9),
        axis.ticks = element_line(colour = 'black'),
        plot.tag = element_text(colour = 'black', size = font_size, face = 'bold'),
        panel.background = element_blank())
}
  
windowed_all <- left_join(
  windowed, filter(zns_files, region == 'mt'), by = 'start'
  ) %>% 
  mutate(
    Fst = as.numeric(Fst), 
    GCsilent = as.numeric(GCsilent)
  )

# Fst plot
windowed_fst_plot <- windowed_all %>% 
  ggplot(aes(x = zns, y = Fst)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE) +
  labs(
    x = expression(Z[nS]),
    y = expression(F[ST])
  ) +
  coord_cartesian(
    y = c(0, 1)
  ) +
  ylim(0, 1) +
  corr_plot_theme(16) +
  annotate('text', x = 0.2, y = 0.9, size = 6,
           label = 'italic(R) ^ 2 == 0.534',
           parse = TRUE) +
  annotate('text', x = 0.2, y = 0.8, size = 6,
           label = 'p < 2.2 %*% 10 ^ -16', parse = TRUE)

windowed_fst_plot

ggsave('plots/fig_2.pdf', plot = windowed_fst_plot,
       width = par('din')[1] * 0.75, height = par('din')[1] * 0.75)

lm(Fst ~ zns, data = windowed_all) %>% 
  summary() # r2 = 0.5338, p < 2.2e-16

# figure 3 -
# a) pi and ZnS over both mating types
# b) piN piS and ZnS over shared genes

# pi ~ ZnS - all samples pooled
ggplot(polymorphism_shared, aes(x = zns_all, y = theta_pi_4fold_ALL)) +
  geom_point() +
  geom_smooth(method = 'lm')
lm(theta_pi_4fold_ALL ~ zns_all, data = polymorphism_shared) %>% 
  summary() # R2 = 0.000176, p = 0.934

# 3A - mt-separated
pi_mt_plot <- ggplot(polymorphism_shared, aes(x = zns_all)) +
  geom_point(
    data = polymorphism_shared,
    aes(y = theta_pi_4fold_MTPLUS),
    size = 2,
    color = wes_palette('Darjeeling1')[2]
  ) + # plus
  geom_point(
    data = polymorphism_shared,
    aes(y = theta_pi_4fold_MTMINUS), 
    size = 2,
    color = wes_palette('Darjeeling1')[3]) +
  geom_smooth(
    data = polymorphism_shared,
    aes(y = theta_pi_4fold_MTPLUS), 
    method = 'lm',
    se = FALSE,
    color = wes_palette('Darjeeling1')[2]) + # plus
  geom_smooth(
    data = polymorphism_shared,
    aes(y = theta_pi_4fold_MTMINUS), method = 'lm',
    se = FALSE,
    color = wes_palette('Darjeeling1')[3]) +
  corr_plot_theme(16) +
  labs(
    x = expression(paste(Z[nS])),
    y = expression(paste(pi[silent])),
    tags = 'a'
  ) +
  coord_cartesian(
    x = c(0, 1),
    y = c(0, 0.045)
  ) +
  annotate(
    'text', x = 0.9, y = 0.04, 
    size = 6, label = 'MT+',
    color = wes_palette('Darjeeling1')[2]) +
  annotate(
    'text', x = 0.9, y = 0.035, 
    size = 6, label = 'MT-',
    color = wes_palette('Darjeeling1')[3]) +
  NULL
  
pi_mt_plot

# fits
lm(theta_pi_4fold_MTPLUS ~ zns_all, data = polymorphism_shared) %>% 
  summary() # R2 = 0.1395, p = 0.016
  
lm(theta_pi_4fold_MTMINUS ~ zns_all, data = polymorphism_shared) %>% 
  summary() # R2 = 0.187, p = 0.0047

# 3B - piN/piS
pi_all <- list(polymorphism_shared, polymorphism_limited) %>% 
  map_dfr(~ select(., gene, contains('piN_piS'), sites_ALL, zns_all), .id = 'name') %>% 
  mutate(
    is_gametolog = case_when(
      name == 1 ~ 1,
      name == 2 ~ 0
    )
  ) %>% 
  select(-name)

pi_all %>% 
  filter(is_gametolog == 1) %>% 
  filter(gene != '522915', gene != 'MT0796') %>% 
  select(piN_piS, zns_all, sites_ALL) %>% 
  lm(piN_piS ~ zns_all, data = ., weights = sites_ALL) %>% 
  summary() # R2 = 0.106, p = 0.043

# incl. weird genes
pi_all %>% 
  filter(is_gametolog == 1) %>%
  select(piN_piS, zns_all, sites_ALL) %>% 
  lm(piN_piS ~ zns_all, data = ., weights = sites_ALL) %>% 
  summary() # R2 = 0.094, p = 0.051

pi_zns_plot <- polymorphism_shared %>% 
  filter(gene != 'MT0796', gene != '522915') %>% 
  ggplot(aes(x = zns_all, y = piN_piS)) +
  geom_point(size = 2) +
  geom_smooth(method = 'lm', se = FALSE) +
  corr_plot_theme(16) +
  labs(
    x = expression(paste(Z[nS])),
    y = expression(paste(pi[N]/pi[S])),
    tag = 'b'
  ) +
  coord_cartesian(
    x = c(0, 1),
    y = c(0, 1)
  ) +
  annotate(
    'text', x = 0.2, y = 0.9, size = 6, 
     label = 'italic(R) ^ 2 == 0.106',
     parse = TRUE) +
  annotate(
    'text', x = 0.2, y = 0.8, size = 6,
    label = 'p == 0.043', parse = TRUE)

pi_zns_plot

ggsave('plots/piNpiS_zns_genes.pdf', plot = pi_zns_plot)

# combined plot
fig_3 <- pi_mt_plot + pi_zns_plot
fig_3

ggsave('plots/fig_3.pdf', plot = fig_3,
       width = par('din')[1] * 1.5, height = par('din')[1] * 1.5 * 0.5)

# diversity stats

### silent pi in all, MT+, and MT- in shared regions

polymorphism_other %>% 
  summarise_at(vars(starts_with('theta_pi_4fold')), mean, na.rm = TRUE) %>% 
  gather(measure, value)

polymorphism_chr6 %>% 
  summarise_at(vars(starts_with('theta_pi_4fold')), mean, na.rm = TRUE) %>% 
  gather(measure, value)

polymorphism_shared %>% 
  summarise_at(vars(starts_with('theta_pi_4fold')), mean, na.rm = TRUE) %>% 
  gather(measure, value) # all three are lower than autosomal levels, but MT+ especially low

### silent pi in MT+ and MT- limited regions

polymorphism_limited %>% 
  mutate(allele = case_when(
    !is.na(mtPlus_coords) ~ 'plus',
    !is.na(mtMinus_coords) ~ 'minus')
  ) %>% 
  group_by(allele) %>% 
  summarise(mean_theta = mean(theta_pi_4fold_ALL))

# GC4
windowed_all %>% 
  lm(GCsilent ~ zns, data = .) %>% 
  summary() # r2 = 0.023, p = 0.00028, m = -0.0477

# figure 4 - organelle linkage
# A - cp + mta/fus
# B - mt + mtd1/mid
# C - mt with itself
# D - randomly selected pairs

linkage_plot_theme <- function(font_size = 16) {
  theme(axis.title = element_text(family = "Helvetica", size = font_size),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text = element_text(family = "Helvetica", size = font_size, color = 'black'),
        axis.line = element_line(colour = 'black', linetype = 'solid'),
        axis.line.x = element_line(size = 0.9),
        axis.line.y = element_line(size = 0.9),
        axis.ticks = element_blank(),
        plot.tag = element_text(colour = 'black', size = font_size, face = 'bold'),
        panel.background = element_blank())
}

ld_dflist <- fs::dir_ls('data/organelle-linkage/', regexp = '\\.txt') %>%
  map_dfr(read_delim, delim = ' ', col_types = cols(), .id = 'name') %>% 
  mutate(name = str_replace(name, 'data/organelle-linkage/', '')) %>% 
  mutate(name = str_replace(name, '.txt', '')) %>% 
  split(.$name)

ld_dflist %>% 
  map(~ summarise(., mean_dprime = mean(dprime))) %>% 
  unlist()
# allpairs = 0.781
# minus = 1.0
# mt only = 1.0
# plus = 0.611

cp_plot <- ggplot(ld_dflist$plus, aes(x = dprime)) +
  geom_histogram(binwidth = 0.1, fill = wes_palette('Darjeeling1')[2]) +
  linkage_plot_theme(16) +
  coord_cartesian(x = c(-0.05, 1.05)) +
  labs(x = "D'", y = 'count', tag = 'a') +
  annotate_textp(
    'text', x = 0.05, y = 0.9, size = 16, 
     label = "cpDNA x MT+")
cp_plot

mt_plot <- ggplot(ld_dflist$minus, aes(x = dprime)) +
  geom_histogram(binwidth = 0.1, fill = wes_palette('Darjeeling1')[3]) +
  linkage_plot_theme(16) +
  coord_cartesian(x = c(-0.05, 1.05)) +
  labs(x = "D'", y = 'count', tag = 'b') +
  annotate_textp(
    'text', x = 0.05, y = 0.9, size = 16,
    label = 'mtDNA x MT-'
  )
mt_plot

mt_only_plot <- ggplot(ld_dflist$mt_only, aes(x = dprime)) +
  geom_histogram(binwidth = 0.1, fill = wes_palette('Darjeeling1')[4]) +
  linkage_plot_theme(16) +
  coord_cartesian(x = c(-0.05, 1.05)) +
  labs(x = "D'", y = 'count', tag = 'c') +
  annotate_textp(
    'text', x = 0.05, y = 0.9, size = 16,
    label = 'mtDNA'
  )
mt_only_plot

allpairs_plot <- ggplot(ld_dflist$allpairs, aes(x = dprime)) +
  geom_histogram(binwidth = 0.1, fill = wes_palette('Darjeeling1')[5]) +
  linkage_plot_theme(16) +
  labs(x = "D'", y = expression(paste('count (x', 10^3, ')')), tag = 'd') +
  coord_cartesian(x = c(-0.05, 1.05)) +
  scale_y_continuous(labels = seq(0, 400, by = 100)) +
  annotate_textp(
    'text', x = 0.05, y = 0.9, size = 16,
    label = 'autosomal'
  )
allpairs_plot

fig_4 <- cp_plot + mt_plot + mt_only_plot + allpairs_plot
fig_4

ggsave('plots/fig_4.pdf', plot = fig_4,
       width = par('din')[1] * 1.25, height = par('din')[1])

