# Script microbiome data

#---- Install packages-----

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager") 
a
BiocManager::install("phyloseq")
install.packages("ggplot2")
install.packages("vegan")
install.packages("dplyr")
install.packages("tidyverse")
install.packages("dada2")
install.packages("hillR")
install.packages("ranacapa")
install.packages("ggpubr")
install.packages("scales")
install.packages("pheatmap")
install.packages("ggvenn")
install.packages("ggVennDiagramm")
install.packages("randomcoloR")
install.packages("shiny")
install.packages("ranacapa", dependencies = TRUE)
install.packages("ggrepel")
install.packages("agricolae")

library(phyloseq)
library(ggplot2)
library(vegan)
library(dplyr)
library(tidyverse)
library(dada2)
library(hillR)
library(ranacapa) #cannot get this to work yet
library(ggpubr)
library(scales)
library(pheatmap)
library(ggvenn)
library(ggVennDiagram) 
library(randomcoloR)
library(ggrepel)
library(agricolae)

#-----Set pathway and colour package----

setwd("C:/Users/raabf/Nextcloud/Cloud/Documents/Pretests/PT4/Sequencing/Sequencing_PT4_test2_phyloseq_files")

#Colours
#Colour package with lot's of contrast. Enough colours for phyla.

Mycolors1 <- c(
  "#CBD588", "#5F7FC7","powderblue", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", 
  "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#66bd63",
  "#f46d43", "#4292c6", "#fa9fb5", "#8c510a", "#80cdc1",
  "#c2a5cf", "#bababa", "#fdbf6f", "#fb9a99", "#ffffb3",
  "#5ab4ac", "#a1d76a", "#af8dc3", "#d95f02", "#fbb4ae",
  "#b3e2cd", "#fdcdac", "#cbd5e8", "#e41a1c", "#377eb8",
  "#8dd3c7", "#ffffb3", "#fee090", "#313695", "#a50026",
  "#1E90FF", "#228B22", "#9400D3", "#D2691E", "#4682B4", 
  "#48D1CC", "#FFA500", "#FF4500")




Mycolors2 <- distinctColorPalette(80)

#---Load data + metadata------

#read files
OTU <- read.csv2("OTU_test2_v3.csv", row.names = 1) #row.names = 1
TAX <- read.csv2("taxonomy_test2_v2.csv", row.names = 1) #row.names = 1

#change long headers in taxa table to "Family", "Genus" etc.
colnames(TAX) <- sub("\\..*", "", colnames(TAX))

#
str(OTU)


#create physeq
OTU_mat <- otu_table(as.matrix(OTU), taxa_are_rows = TRUE) 
TAX_mat <- tax_table(as.matrix(TAX))

physeq_dada <- phyloseq(OTU_mat, TAX_mat)
physeq_dada


#test row names
colnames(tax_table(physeq_dada))

# make metadata
meta <- read.csv2("metadata_PT4_v3.csv", row.names = 1, head = TRUE)

meta$treatment_timepoint <- paste(meta$treatment, meta$timepoint, sep = "_")

meta$treatment_timepoint_replicate <- paste(
  meta$treatment, meta$timepoint, meta$replicate, sep = "_"
)

#test metadata
colnames(meta)

# Convert metadata data.frame to sample_data first
sample_metadata <- sample_data(meta)

# Then merge
physeq <- merge_phyloseq(physeq_dada, sample_metadata)



##-----Plotting all data (does not work well, clean first)-----

#plot family
#plot_bar(physeq, fill = "Phylum") 


#save plot
#fig1 <- plot_bar(physeq, fill = "Phylum") 
#fig1 + scale_fill_manual(values = Mycolors2)

#ggsave("Phylum.png", plot = fig1, width = 8, height = 6)

#dim(OTU)
#head(OTU)

#all(rownames(TAX) == rownames(OTU))
#colnames(TAX_mat)

#-----Plot OTU read distribution----

reads_per_OTU = data.frame(nreads = sort (taxa_sums(physeq), TRUE), SVs = 1:ntaxa(physeq_dada))

ggplot(data=reads_per_OTU, aes(x = SVs, y = nreads))+
  geom_area(stat="identity", fill= "#69BE28")+ scale_x_log10()



#---Plot sample read distribution----
reads_per_sample = data.frame(nreads = sort(sample_sums(physeq), TRUE),
                              samples = 1:nsamples(physeq), 
                              type = "samples")
ggplot(data = reads_per_sample, 
       aes(x=samples, y=nreads))+
  geom_bar(stat= "identity", fill = "#69BE28")



#-----Rarefy-----

rarefied <- rarefy_even_depth(physeq, sample.size = 30000, rngseed=1)
rarefied
#180OTUs were removed because they are no longer present in any sample after random subsampling

saveRDS(rarefied, "20260413_rarefied.RDS")   #change date when necessary 
rarecurve <- ggrare(rarefied, step= 1000, se = FALSE) #does not work yet



## Check for empty ASVs - prune

sum(taxa_sums(rarefied) == 0)

#remove empty ASVs when present
rarefied = prune_taxa(taxa_sums(rarefied) > 0,
                      rarefied)


## Relative abundances

rarefied.rel = transform_sample_counts(rarefied, function(x)x/sum(x))


## Filter and transform 

physeq.rel_fil = filter_taxa(rarefied.rel, function(x) mean(x) > 1e-5, TRUE)                              
physeq.rel_fil 

physeq.rel1 = transform_sample_counts(physeq.rel_fil, function(x) x/ sum(x))  
physeq.rel1

sample_sums(physeq.rel_fil)
sample_sums(rarefied)
saveRDS (physeq.rel1, "20260422_Rarefied_test2.RDS")


## Agglomeration + smelting

#Agglomerate Phylum
tax_glom = tax_glom(physeq.rel1, taxrank = "Phylum")
#Psmelt phylum
melted.p = psmelt(tax_glom)
View(melted.p)

# Agglomerate Class
tax_glom_class = tax_glom(physeq.rel1, taxrank = "Class")
#Psmelt class
melted.p_class = psmelt(tax_glom_class)
View(melted.p_class)

#Agglomerate Order
tax_glom_order = tax_glom(physeq.rel1, taxrank = "Order")
#Psmelt Order
melted.p_order = psmelt(tax_glom_order)
View(melted.p_order)

#Agglomerate Family
tax_glom_family = tax_glom(physeq.rel1, taxrank = "Family")
#Psmelt genus
melted.p_family = psmelt(tax_glom_family)
View(melted.p_family)

# Agglomerate Genus
tax_glom_genus = tax_glom(physeq.rel1, taxrank = "Genus")
#Psmelt genus
melted.p_genus = psmelt(tax_glom_genus)
View(melted.p_genus)


## Figure abundance per phyla

ggplot(melted.p, aes(treatment, y= Abundance))+
  geom_point()+
  facet_wrap(~Phylum)+
  scale_x_discrete(guide = guide_axis(angle = 50))+
  theme_classic ()+
  xlab("treatment")+
  ylab("Abundance") 

ggsave("20260409_abundance per phyla.png", width=6, height=5, dpi=300)



#--Abundance plots on phylum level----

plot_abundance_phylum_grouped_timepoint <-ggplot(melted.p, aes(x=as.factor(treatment), 
                                                      y = Abundance/3, 
                                                      fill = Phylum)) + 
  geom_bar(stat = "identity")+
  scale_x_discrete(guide = guide_axis(angle = 50))+
  facet_wrap(~timepoint)+
  theme_classic()+
  xlab("Treatment")+
  ylab("Abundance")+
  scale_y_continuous(labels = percent, limits = c(0, 1))+
  scale_fill_manual(values= Mycolors1)

plot_abundance_phylum_grouped_timepoint
ggsave("abundance_phylum_level_grouped_timepoints.png", width=6, height=5, dpi=300)

plot_abundance_phylum_grouped_treatment <-ggplot(melted.p, aes(x=as.factor(timepoint), 
                                      y = Abundance/3, 
                                      fill = Phylum)) + 
  geom_bar(stat = "identity")+
  scale_x_discrete(guide = guide_axis(angle = 50))+
  facet_wrap(~treatment)+
  theme_classic()+
  xlab("Timepoints")+
  ylab("Abundance")+
  scale_y_continuous(labels = percent, limits = c(0, 1))+
  scale_fill_manual(values= Mycolors1)

plot_abundance_phylum_grouped_treatment
ggsave("abundance_phylum_level_grouped_treatments.png", width=6, height=5, dpi=300)


#--Abundance plots on class level----

plot_abundance_class_grouped_timepoint <-ggplot(melted.p_class, aes(x=as.factor(treatment), 
                                                               y = Abundance/3, 
                                                               fill = Class)) + 
  geom_bar(stat = "identity")+
  scale_x_discrete(guide = guide_axis(angle = 50))+
  facet_wrap(~timepoint)+
  theme_classic()+
  xlab("Treatment")+
  ylab("Abundance")+
  scale_y_continuous(labels = percent, limits = c(0, 1))+
  scale_fill_manual(values= Mycolors2)

plot_abundance_class_grouped_timepoint
ggsave("abundance_class_level_grouped_timepoints.png", width=12, height=10, dpi=300)

plot_abundance_class_grouped_treatment <-ggplot(melted.p_class, aes(x=as.factor(timepoint), 
                                                               y = Abundance/3, 
                                                               fill = Class)) + 
  geom_bar(stat = "identity")+
  scale_x_discrete(guide = guide_axis(angle = 50))+
  facet_wrap(~treatment)+
  theme_classic()+
  xlab("Timepoints")+
  ylab("Abundance")+
  scale_y_continuous(labels = percent, limits = c(0, 1))+
  scale_fill_manual(values= Mycolors2)

plot_abundance_class_grouped_treatment
ggsave("abundance_class_level_grouped_treatments.png", width=12, height=10, dpi=300)


#--Abundance plots on order level----
## too many taxa to show in a graph

plot_abundance_order_grouped_timepoint <-ggplot(melted.p_order, aes(x=as.factor(treatment), 
                                                                    y = Abundance/3, 
                                                                    fill = Order)) + 
  geom_bar(stat = "identity")+
  scale_x_discrete(guide = guide_axis(angle = 50))+
  facet_wrap(~timepoint)+
  theme_classic()+
  xlab("Treatment")+
  ylab("Abundance")+
  scale_y_continuous(labels = percent, limits = c(0, 1))+
  scale_fill_manual(values= Mycolors2)

plot_abundance_order_grouped_timepoint
ggsave("abundance_order_level_grouped_timepoints.png", width=6, height=5, dpi=300)

plot_abundance_order_grouped_treatment <-ggplot(melted.p_order, aes(x=as.factor(timepoint), 
                                                                    y = Abundance/3, 
                                                                    fill = Order)) + 
  geom_bar(stat = "identity")+
  scale_x_discrete(guide = guide_axis(angle = 50))+
  facet_wrap(~treatment)+
  theme_classic()+
  xlab("Treatment")+
  ylab("Abundance")+
  scale_y_continuous(labels = percent, limits = c(0, 1))+
  scale_fill_manual(values= Mycolors2)

plot_abundance_order_grouped_treatment
ggsave("abundance_order_level_grouped_treatments.png", width=6, height=5, dpi=300)


#--old Alpha and beta diversity----


### Alpha diversity measurements

plot_rich <- plot_richness(rarefied, x = "treatment_timepoint", 
                               measures = c("Observed", "Chao1", "Shannon", "InvSimpson"), 
                               title = "Shannon")+
  #facet_wrap(~timepoint)+
  geom_boxplot(aes(fill = factor(treatment), group = treatment))+
  theme(aspect.ratio = 1) 
#ylim(1800, 3300)

plot_rich
ggsave("plot_rich_alpha.png", width=12, height=10, dpi=300)

#Shannon richness

plot_rich_new <- plot_richness(rarefied, x = "treatment", 
                                    measures = c("Shannon"), 
                                    title = "Shannon")+
  facet_wrap(~timepoint)+
  geom_boxplot(aes(fill = factor(treatment), group = treatment))+
  theme(aspect.ratio = 1) 
  #ylim(1800, 3300)

plot_rich_new
ggsave("plot_rich_alpha_shannon.png", width=12, height=10, dpi=300)


## Beta diversity
#Difference between two communities
#Ordination plot

rarefied.o = ordinate (rarefied, method = "PCoA", distance = "bray")
plot_ordination(rarefied, rarefied.o)

#ordinate
rarefied@sam_data$timepoint = as.factor(rarefied@sam_data$timepoint)
plot_ordination(rarefied, rarefied.o, color = "treatment_timepoint_replicate", label= "treatment_timepoint_replicate")+
  geom_point(size=4)+
  theme_classic()

ggsave("20260409_ordinate_beta_diversity.png", width = 6, height=5, dpi = 300)

plot_scree(rarefied.o)
str(rarefied.o)

df = plot_ordination(rarefied, rarefied.o, justDF = TRUE)
df

#-----remove outlier-----

# create subset

ps <- prune_samples(sample_names(physeq) != "F11", physeq)

ps <- prune_samples(!(sample_names(physeq) %in% c("F11", "F7")), physeq)

#-----Plot OTU read distribution----

reads_per_OTU = data.frame(nreads = sort (taxa_sums(ps), TRUE), SVs = 1:ntaxa(physeq_dada))

ggplot(data=reads_per_OTU, aes(x = SVs, y = nreads))+
  geom_area(stat="identity", fill= "#69BE28")+ scale_x_log10()



#---Plot sample read distribution----
reads_per_sample = data.frame(nreads = sort(sample_sums(ps), TRUE),
                              samples = 1:nsamples(ps), 
                              type = "samples")
ggplot(data = reads_per_sample, 
       aes(x=samples, y=nreads))+
  geom_bar(stat= "identity", fill = "#69BE28")



#-----Rarefy-----

rarefied <- rarefy_even_depth(ps, sample.size = 30000, rngseed=1)
rarefied
#259 OTUs were removed because they are no longer present in any sample after random subsampling
#264 OTUs removed with F11 and F7

saveRDS(rarefied, "20260422_rarefied_wo_F11_F7.RDS")   #change date when necessary 
rarecurve <- ggrare(rarefied, step= 1000, se = FALSE) #does not work yet



## Check for empty ASVs - prune

sum(taxa_sums(rarefied) == 0)

#remove empty ASVs when present
rarefied = prune_taxa(taxa_sums(rarefied) > 0,
                      rarefied)


## Relative abundances

rarefied.rel = transform_sample_counts(rarefied, function(x)x/sum(x))


## Filter and transform 

ps.rel_fil = filter_taxa(rarefied.rel, function(x) mean(x) > 1e-5, TRUE)                              
ps.rel_fil 

ps.rel1 = transform_sample_counts(ps.rel_fil, function(x) x/ sum(x))  
ps.rel1

sample_sums(ps.rel_fil)
sample_sums(rarefied)
saveRDS (ps.rel1, "20260422_Rarefied_wo_F11_F11.RDS")


## Agglomeration + smelting

#Agglomerate Phylum
tax_glom = tax_glom(ps.rel1, taxrank = "Phylum")
#Psmelt phylum
melted.p = psmelt(tax_glom)
View(melted.p)

# Agglomerate Class
tax_glom_class = tax_glom(ps.rel1, taxrank = "Class")
#Psmelt class
melted.p_class = psmelt(tax_glom_class)
View(melted.p_class)

#Agglomerate Order
tax_glom_order = tax_glom(ps.rel1, taxrank = "Order")
#Psmelt Order
melted.p_order = psmelt(tax_glom_order)
View(melted.p_order)

#Agglomerate Family
tax_glom_family = tax_glom(ps.rel1, taxrank = "Family")
#Psmelt genus
melted.p_family = psmelt(tax_glom_family)
View(melted.p_family)

# Agglomerate Genus
tax_glom_genus = tax_glom(ps.rel1, taxrank = "Genus")
#Psmelt genus
melted.p_genus = psmelt(tax_glom_genus)
View(melted.p_genus)


## Figure abundance per phyla

ggplot(melted.p, aes(treatment, y= Abundance))+
  geom_point()+
  facet_wrap(~Phylum)+
  scale_x_discrete(guide = guide_axis(angle = 50))+
  theme_classic ()+
  xlab("treatment")+
  ylab("Abundance") 

ggsave("20260409_abundance per phyla_wo_F11_F7.png", width=6, height=5, dpi=300)



#--Abundance plots on phylum level----

melted_summary_Phylum_abund <- melted.p %>% 
  group_by(timepoint, treatment, Phylum) %>%
  summarize_at("Abundance", sum)#


ggplot(melted_summary_Phylum_abund, aes(x = as.factor(treatment),
                           y = Abundance,
                           fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~timepoint) +
  theme_classic() +
  xlab("Treatment") +
  ylab("Abundance") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_x_discrete(guide = guide_axis(angle = 50)) +
  scale_fill_manual(values = Mycolors1)+
  theme(axis.text = element_text(size = 22))+
  theme(axis.title = element_text(size = 22))+
  theme(strip.text = element_text(size = 22))

ggsave("20260427_EGU_abundance_phylum_grouped_timepoints_wo_F11_F7.png", width=12, height=10, dpi=300)



#--Abundance plots on class level----

melted_summary_class_abund <- melted.p_class %>% 
  group_by(timepoint, treatment, Class) %>%
  summarize_at("Abundance", sum)#


ggplot(melted_summary_class_abund, aes(x = as.factor(treatment),
                                        y = Abundance,
                                        fill = Class)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~timepoint) +
  theme_classic() +
  xlab("Treatment") +
  ylab("Abundance") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_x_discrete(guide = guide_axis(angle = 50)) +
  scale_fill_manual(values = Mycolors2)+
  theme(axis.text = element_text(size = 14))+
  theme(axis.title = element_text(size = 20))+
  theme(strip.text = element_text(size = 16))

ggsave("new__v2_abundance_class_grouped_timepoints_wo_F11_F7.png", width=12, height=10, dpi=300)



#--Alpha and statistics----


### Alpha diversity measurements Luis


physeq_general_diversity <- estimate_richness(rarefied, measures = c("Shannon","Observed"))


physeq_general_diversity_meta <- cbind(sample_data(rarefied), physeq_general_diversity)


#boxplot
ggplot(physeq_general_diversity_meta, 
       aes(x = as.factor(treatment), 
           y = Shannon, 
           fill = treatment_qualitative)) + 
  geom_boxplot(aes(colour = treatment_qualitative)) +
  facet_wrap(~timepoint) +
  theme_classic() +
  xlab("Treatment") +
  ylab("Shannon") +
  labs(fill = "Treatment", colour = "Treatment")
  #scale_y_continuous(labels = percent, limits = c(0, 1))+
  #scale_fill_manual(values= Mycolors2)

ggsave("shannon_wo_F11_F7.png", width=12, height=10, dpi=300)

 

#pointbar old
ggplot(physeq_general_diversity_meta,
       aes(x = as.factor(treatment),
           y = Shannon,
           color = treatment_qualitative)) +
  geom_point(size = 3, position = position_jitter(width = 0.1)) +
  stat_summary(fun = mean, geom = "point", size = 5, color = "black") +
  facet_wrap(~timepoint) +
  theme_classic() +
  xlab("Treatment") +
  ylab("Shannon") +
  labs(color = "Treatment")

ggsave("shannon_wo_F11_F7_pointbar.png", width=6, height=5, dpi=300)

#pointbar new
ggplot(physeq_general_diversity_meta,
       aes(x = as.factor(treatment),
           y = Shannon,
           color = treatment_qualitative,
           group = treatment_qualitative)) +
  geom_point(size = 3, position = position_jitter(width = 0.1)) +
  
  # mean line
  stat_summary(fun = mean, geom = "point", aes(group = 1), size = 6, shape = 17) +
  
  # mean points (optional but usually helpful)
 # stat_summary(fun = mean, geom = "point", size = 3, color = "black") +
  
  # standard deviation error bars
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "errorbar", width = 0.2) +
  
  facet_wrap(~timepoint) +
  theme_classic() +
  xlab("Treatment") +
  ylab("Shannon") #+
  #labs(color = "Treatment")

ggsave("shannon_wo_F11_F7_pointbar_new.png", width=6, height=5, dpi=300)


# modification for axis size and legend title

ggplot(physeq_general_diversity_meta,
       aes(x = as.factor(treatment),
           y = Shannon,
           color = treatment_qualitative,
           group = treatment_qualitative)) +
  geom_point(size = 3, position = position_jitter(width = 0.1)) +
  
  # mean line
  stat_summary(fun = mean, geom = "point", aes(group = 1), size = 6, shape = 17) +
  
  # mean points (optional but usually helpful)
  # stat_summary(fun = mean, geom = "point", size = 3, color = "black") +
  
  # standard deviation error bars
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "errorbar", width = 0.2) +
  
  facet_wrap(~timepoint) +
  theme_classic() +
# 1. Achsenbeschriftungen und Legendentitel festlegen
labs(x = "Treatment", 
     y = "Shannon Index", 
     color = expression(""^2*H[2]*O ~ "concentration")) + 
  # Hier die Farbanpassung:
  scale_color_viridis_d(option = "plasma") +# Optionen: "viridis", "magma", "plasma", "inferno", "cividis"
  
  # 2. Schriftgrößen im Theme anpassen
  theme(
    axis.title.x = element_text(size = 14), # Größe X-Achsenbeschriftung
    axis.title.y = element_text(size = 14), # Größe Y-Achsenbeschriftung
    axis.text = element_text(size = 12),    # Größe der Zahlen an den Achsen
    legend.title = element_text(size = 13), # Größe Legendentitel
    legend.text = element_text(size = 11),  # Größe Legendentext
    strip.text = element_text(size = 12) )   # Größe der Facet-Überschriften

ggsave("20260428_EGU_shannon_wo_F11_F7_pointbar.png", width=6, height=5, dpi=300)

##----anova----

anova_test <- aov(Shannon ~ treatment_timepoint, physeq_general_diversity_meta)

summary(anova_test)
#Df  Sum Sq  Mean Sq F value Pr(>F)
#treatment_timepoint 13 0.02528 0.001945   0.612  0.822
#Residuals           26 0.08266 0.003179

qqnorm(residuals(anova_test))
hist(residuals(anova_test))
shapiro.test(physeq_general_diversity_meta$Shannon)
#data:  physeq_general_diversity_meta$Shannon
#W = 0.96987, p-value = 0.3564 -> normalverteilt

tukey_result <- TukeyHSD(anova_test)

print(tukey_result)


##---observed----
ggplot(physeq_general_diversity_meta,
       aes(x = as.factor(treatment),
           y = Observed,
           color = treatment_qualitative,
           group = treatment_qualitative)) +
  geom_point(size = 3, position = position_jitter(width = 0.1)) +
  
  # mean line
  stat_summary(fun = mean, geom = "point", aes(group = 1), size = 6, shape = 17) +
  
  # mean points (optional but usually helpful)
  # stat_summary(fun = mean, geom = "point", size = 3, color = "black") +
  
  # standard deviation error bars
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "errorbar", width = 0.2) +
  
  facet_wrap(~timepoint) +
  theme_classic() +
  # 1. Achsenbeschriftungen und Legendentitel festlegen
  labs(x = "Treatment", 
       y = "Observed Richness", 
       color = expression(""^2*H[2]*O ~ "concentration")) + 
  scale_color_viridis_d(option = "plasma") + # Optionen: "viridis", "magma", "plasma", "inferno", "cividis"
  
  # 2. Schriftgrößen im Theme anpassen
  theme(
    axis.title.x = element_text(size = 14), # Größe X-Achsenbeschriftung
    axis.title.y = element_text(size = 14), # Größe Y-Achsenbeschriftung
    axis.text = element_text(size = 12),    # Größe der Zahlen an den Achsen
    legend.title = element_text(size = 13), # Größe Legendentitel
    legend.text = element_text(size = 11),  # Größe Legendentext
    strip.text = element_text(size = 12)    # Größe der Facet-Überschriften
  )

ggsave("20260429_EGU_observed_wo_F11_F7_pointbar.png", width=6, height=5, dpi=300)


ggplot(physeq_general_diversity_meta, 
       aes(x = as.factor(treatment), 
           y = Observed, 
           fill = treatment_qualitative)) + 
  geom_boxplot(aes(colour = treatment_qualitative)) +
  facet_wrap(~timepoint) +
  theme_classic() +
  xlab("Treatment") +
  ylab("Observed") +
  labs(fill = "Treatment", colour = "Treatment")

ggsave("observed_wo_F11_F7.png", width=12, height=10, dpi=300)

ggplot(physeq_general_diversity_meta,
       aes(x = as.factor(treatment),
           y = Observed,
           color = treatment_qualitative)) +
  geom_point(size = 3, position = position_jitter(width = 0.1)) +
  stat_summary(fun = mean, geom = "point", size = 5, color = "black") +
  facet_wrap(~timepoint) +
  theme_classic() +
  xlab("Treatment") +
  ylab("Observed") +
  labs(color = "Treatment")

ggsave("observed_wo_F11_F7_pointbar.png", width=6, height=5, dpi=300)

##----anova-----
anova_test_observed <- aov(Observed ~ treatment_timepoint, physeq_general_diversity_meta)


summary(anova_test_observed)
#Df Sum Sq Mean Sq F value Pr(>F)
#treatment_timepoint 13  60137    4626   0.746  0.705
#Residuals           26 161245    6202  

qqnorm(residuals(anova_test_observed))
hist(residuals(anova_test_observed))
shapiro.test(physeq_general_diversity_meta$Observed)
#data:  physeq_general_diversity_meta$Observed
#W = 0.98326, p-value = 0.8076 -> normalverteilt

tukey_result <- TukeyHSD(anova_test_observed)

print(tukey_result)

#-----Beta diversity----
#Difference between two communities
#Ordination plot PCOA

rarefied.o = ordinate (rarefied, method = "PCoA", distance = "bray")
plot_ordination(rarefied, rarefied.o)

#ordinate
rarefied@sam_data$timepoint = as.factor(rarefied@sam_data$timepoint)
plot_ordination(rarefied, rarefied.o, color = "treatment_qualitative", label= "treatment_timepoint")+
  labs(color = "Treatment")+
  geom_point(size=4)+
  theme_classic()+
  scale_color_viridis_d(option = "plasma") + # Optionen: "viridis", "magma", "plasma", "inferno", "cividis"
  guides(color = guide_legend(override.aes = list(label = "")))

ggsave("20260413_ordinate_beta_diversity.png", width = 6, height=5, dpi = 300)

colnames(sample_data(rarefied))

plot_scree(rarefied.o)
str(rarefied.o)

df = plot_ordination(rarefied, rarefied.o, justDF = TRUE)
df

#Ordination plot NMDS

rarefied.o = ordinate (rarefied, method = "NMDS", distance = "bray")
plot_ordination(rarefied, rarefied.o)

#ordinate
rarefied@sam_data$timepoint = as.factor(rarefied@sam_data$timepoint)

plot_ordination(rarefied, rarefied.o, color = "treatment_qualitative")+
  labs(color = "Treatment")+
  geom_text_repel(aes(label = timepoint)) +
  geom_point(size=4)+
  theme_classic()+
  
  guides(color = guide_legend(override.aes = list(label = "")))

ggsave("20260413_ordinate_timepoint_wo_F11_F7.png", width = 6, height=5, dpi = 300)

#new plos
# 1. Sicherstellen, dass timepoint ein Factor ist
rarefied@sam_data$timepoint = as.factor(rarefied@sam_data$timepoint)

# 2. Palette definieren
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

# 3. Plot erstellen
plot_ordination(rarefied, rarefied.o, 
                color = "treatment_qualitative", 
                shape = "timepoint") + # Form nach Zeitpunkt
  geom_point(size = 4) +
  scale_color_manual(values = cbPalette) + # Barrierefreie Farben
  labs(color = "Treatment", 
       shape = "Zeitpunkt") +
  theme_classic() +
  # Entfernt das "a" oder die Zahlen aus der Legende der Farben
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 4)))


# 1. Zeitpunkte als Faktoren sicherstellen
rarefied@sam_data$timepoint = as.factor(rarefied@sam_data$timepoint)

# 2. Plot mit Viridis-Palette
plot_ordination(rarefied, rarefied.o, 
                color = "treatment_qualitative", 
                shape = "timepoint") + 
  geom_point(size = 4) +
  
  # Hier die Viridis-Farbskala für diskrete Gruppen
  # 'option' kann sein: "viridis", "magma", "plasma", "inferno" oder "cividis"
  scale_color_viridis_d(option = "plasma") + 
  
  labs(color = "Treatment", 
       shape = "Timepoint") +
  theme_classic() +
  
  # Legende säubern
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 4)))


ggsave("20260428_ordinate_timepoint_wo_F11_F7.png", width = 6, height=5, dpi = 300)


#größere titel, achsen und legende

# 1. Zeitpunkte als Faktoren sicherstellen
rarefied@sam_data$timepoint = as.factor(rarefied@sam_data$timepoint)

# 2. Plot mit vergrößerten Elementen
plot_ordination(rarefied, rarefied.o, 
                color = "treatment_qualitative", 
                shape = "timepoint") + 
  geom_point(size = 4) +
  scale_color_viridis_d(option = "plasma") + 
  labs(color = "Treatment", 
       shape = "Zeitpunkt") +
  theme_classic() +
  
  # Hier erfolgt die Anpassung der Schriftgrößen
  theme(
    # Achsentitel (z.B. "Axis 1 [15.4%]")
    axis.title = element_text(size = 16, face = "bold"), 
    
    # Achsenbeschriftung (die Zahlen an den Achsen)
    axis.text = element_text(size = 14),
    
    # Titel der Legende ("Treatment" & "Timepoint")
    legend.title = element_text(size = 15, face = "bold"),
    
    # Text innerhalb der Legende (die Gruppen-Namen)
    legend.text = element_text(size = 13),
    
    # Optional: Abstand der Legende vergrößern
    legend.key.size = unit(1.2, "cm")
  ) +
  
  # Symbole in der Legende vergrößern (damit man die Formen gut erkennt)
  guides(color = guide_legend(override.aes = list(size = 5)),
         shape = guide_legend(override.aes = list(size = 5)))


ggsave("20260428_beta_timepoint_wo_F11_F7.png", width = 6, height=5, dpi = 300)




##----Transform phyloseq to vegan----

veganotu = function(physeq) {
  require("vegan")
  OTU = otu_table(physeq)
  if (taxa_are_rows(OTU)) {OTU = t(OTU)}
  return(as(OTU, "matrix"))
}

#FUN
#?vegdist()
  physeq_FUN_vegan <- veganotu(rarefied) #physeq after rarefied 

#The distance matrix
  physeq_FUN_ASV_dist <- vegdist(physeq_FUN_vegan, method = "bray") #gower
  physeq_FUN_dataframe <- as(sample_data(rarefied), "data.frame")
  
  names(physeq_FUN_dataframe)

###-----permanova from luis----
  
  
  mod1_PERMA_FUN1 <- adonis2(physeq_FUN_ASV_dist ~treatment_timepoint,
                             data = physeq_FUN_dataframe, 
                             by = "onedf")
  mod1_PERMA_FUN1
  
  mod1_PERMA_FUN2 <- adonis2(physeq_FUN_ASV_dist ~treatment*timepoint,
                             data = physeq_FUN_dataframe, 
                             by = "onedf")
  mod1_PERMA_FUN2
  
  mod1_PERMA_FUN3 <- adonis2(physeq_FUN_ASV_dist ~treatment+timepoint,
                             data = physeq_FUN_dataframe, 
                             by = "onedf")
  mod1_PERMA_FUN3
  
  
  
  

#---old analyses Lieke-----

plot_rich <- plot_richness(rarefied, x = "treatment_timepoint", 
                           measures = c("Observed", "Chao1", "Shannon", "InvSimpson"), 
                           title = "Shannon")+
  #facet_wrap(~timepoint)+
  #facet_grid(~timepoint)+
  geom_boxplot(aes(fill = factor(treatment), group = treatment))+
  theme(aspect.ratio = 1) 
#ylim(1800, 3300)

plot_rich
ggsave("plot_rich_alpha.png", width=12, height=10, dpi=300)

#Shannon richness

plot_rich_new <- plot_richness(rarefied, x = "treatment", 
                               measures = c("Shannon", "Observed"), 
                               title = "Shannon")+
  facet_wrap(~timepoint)+
  geom_boxplot(aes(fill = factor(treatment), group = treatment))+
  theme(aspect.ratio = 1) 
#ygeom_point()#ylim(1800, 3300)

plot_rich_new
ggsave("plot_rich_alpha_shannon.png", width=12, height=10, dpi=300)

rarefied@sam_data

