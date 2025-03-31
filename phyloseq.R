library("phyloseq"); packageVersion("phyloseq")  ##‘1.34.0’
library("ggplot2"); packageVersion("ggplot2") ##‘3.3.3’
theme_set(theme_bw())
library("plyr")
library("cluster")
library("readxl")    
library("dplyr")       
library("knitr")
library("BiocStyle")


sample_names(set1)
rank_names(set1)
sample_variables(set1)

#filtering
#rank_names(set1)
table(tax_table(set1)[, "Phylum"], exclude = NULL)
prevdf = apply(X = otu_table(set1),
               MARGIN = ifelse(taxa_are_rows(set1), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(set1),
                    tax_table(set1))
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
filterPhyla = c("p_Campylobacterota", "p_Deferribacterota")
# Filter entries with unidentified Phylum.
ps1 = subset_taxa(set1, !Phylum %in% filterPhyla)
ps1

prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))
#pdf("barplot/abund_prev.pdf")
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps1),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, linetype = 2) +  geom_point(size = 2) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
#dev.off()

prevalenceThreshold = 0.05 * nsamples(ps1)
prevalenceThreshold #1.05
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps1)
ps2

#require("genefilter")
#flist <- filterfun(kOverA(3, 1e-06))
#ps3=filter_taxa(ps2, flist, TRUE)
#ps3

set1n  = transform_sample_counts(ps2, function(x) x / sum(x) )
set1n

set2n  = transform_sample_counts(ps2, function(x) x / sum(x) * 1000000 + 1)

write.csv(file="set1n.csv",otu_table(set1n))
write.csv(file="set2n.csv",otu_table(set2n))

#saveRDS(set1n, "set1n.rds")
#saveRDS(set2n, "set2n.rds")
#saveRDS(ps2, "ps2.rds")
#saveRDS(set1, "set1.rds")
set1n <- readRDS("set1n.rds")
ps2 <- readRDS("ps2.rds")
set1 <- readRDS("set1.rds")

postscript("barplot/Family_filter.eps",hor=T)
plot_bar(set1n, "Order", fill="Phylum", facet_grid=group~day)
dev.off()

postscript("barplot/Order_filter.eps",hor=T)
plot_bar(set1n, "Order", fill="Phylum", facet_grid=group~day)
dev.off()

postscript("barplot/Phylum.eps",hor=T)
pdf("barplot/Phylum.pdf")
plot_bar(set1n, "Phylum", fill="Phylum", facet_grid=batch~group)
dev.off()


##p_Actinobacteriota
set1n_Actinobacteriota <- subset_taxa(set1n, Phylum %in% c("Actinobacteriota"))
pdf("barplot/group_Actinobacteriota.pdf")
plot_bar(set1n_Actinobacteriota, x="Species", fill = "Genus", facet_grid = group~day, title="Actinobacteriota") 
plot_bar(set1n_Actinobacteriota, x="Genus", fill = "Family", facet_grid = group~day, title="Actinobacteriota") 
#  geom_bar(aes(color=Species, fill=Species), stat="identity", position="stack")
dev.off()

set1n_Bacteroidota <- subset_taxa(set1n, Phylum %in% c("Bacteroidota"))
#plot_bar(set1n_Bacteroidota, x="Genus", fill = "Genus", facet_grid = .~group) +geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
pdf("barplot/group_Bacteroidota_genus.pdf")
plot_bar(set1n_Bacteroidota, x="Species", fill = "Genus", facet_grid = group~day, title="Bacteroidetes") 
plot_bar(set1n_Bacteroidota, x="Genus", fill = "Family", facet_grid = group~day, title="Bacteroidota") 
plot_bar(set1n_Bacteroidota, x="Family", fill = "Family", facet_grid = group~day, title="Bacteroidota") 
#  geom_bar(aes(color=Species, fill=Species), stat="identity", position="stack")
dev.off()

set1n_Cyanobacteria <- subset_taxa(set1n, Phylum %in% c("p_Cyanobacteria"))
#plot_bar(set1n_Bacteroidota, x="Genus", fill = "Genus", facet_grid = .~group) +geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
pdf("barplot/group_p_Cyanobacteria_family.pdf")
plot_bar(set1n_Cyanobacteria, x="Species", fill = "Genus", facet_grid = group~day, title="Cyanobacteria") 
plot_bar(set1n_Cyanobacteria, x="Genus", fill = "Family", facet_grid = group~day2, title="Cyanobacteria") 
plot_bar(set1n_Cyanobacteria, x="Family", fill = "Family", facet_grid = group~day, title="Cyanobacteria") 
#  geom_bar(aes(color=Species, fill=Species), stat="identity", position="stack")
dev.off()

set1n_Desulfobacterota <- subset_taxa(set1n, Phylum %in% c("p_Desulfobacterota"))
#plot_bar(set1n_Bacteroidota, x="Genus", fill = "Genus", facet_grid = .~group) +geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
pdf("barplot/group_p_Desulfobacterota.pdf")
plot_bar(set1n_Desulfobacterota, x="Species", fill = "Genus", facet_grid = group~day, title="p_Desulfobacterota") 
plot_bar(set1n_Desulfobacterota, x="Genus", fill = "Family", facet_grid = group~day, title="p_Desulfobacterota") 
plot_bar(set1n_Desulfobacterota, x="Family", fill = "Family", facet_grid = group~day, title="p_Desulfobacterota") 
#  geom_bar(aes(color=Species, fill=Species), stat="identity", position="stack")
dev.off()

##p_Firmicutes
set1n_Firmicutes <- subset_taxa(set1n, Phylum %in% c("Firmicutes"))
pdf("barplot/group_Firmicutes_order.pdf")
plot_bar(set1n_Firmicutes, x="Order", fill = "Order", facet_grid = group~day, title="Firmicutes") 
#  geom_bar(aes(color=Species, fill=Species), stat="identity", position="stack")
plot_bar(set1n_Firmicutes, x="Class", fill = "Class", facet_grid = group~day, title="Firmicutes")
plot_bar(set1n_Firmicutes, x="Genus", fill = "Genus", facet_grid = group~., title="Firmicutes") #  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
dev.off()
#Clostridia
#Lachnospiraceae
set1n_Clostridia <- subset_taxa(set1n, Class %in% c("Clostridia"))
pdf("barplot/group_Clostridia.eps")
plot_bar(set1n_Clostridia, x="Family", fill = "Family", facet_grid = group~day, title="Clostridia") 
dev.off()

#Lachnospiraceae
set1n_lactobacillales <- subset_taxa(set1n, Order %in% c("o_Lactobacillales"))
postscript("barplot/group_Lachnospiraceae_species.eps",hor=T)
plot_bar(set1n_lactobacillales, x="Species", fill = "Genus", facet_grid = group~day, title="Lachnospiraceae") 
dev.off()

#Lachnospiraceae
set1n_Lachnospiraceae <- subset_taxa(set1n, Family %in% c("Lachnospiraceae"))
postscript("barplot/group_Lachnospiraceae_species.eps",hor=T)
plot_bar(set1n_Lachnospiraceae, x="Genus", fill = "Family", facet_grid = group~day, title="Lachnospiraceae") 
dev.off()

#Lactobacillaceae
set1n_Lactobacillaceae <- subset_taxa(set1n, Family %in% c("f_Lactobacillaceae"))
postscript("barplot/group_Lactobacillales_species.eps",hor=T)
plot_bar(set1n_Lactobacillaceae, x="Genus", fill = "Genus", facet_grid = group~day, title="Lactobacillaceae") 
dev.off()

#p_Fusobacteriota
set1n_Fusobacteriota <- subset_taxa(set1n, Phylum %in% c("p_Fusobacteriota"))
postscript("barplot/group_p_Fusobacteriota.eps",hor=T)
plot_bar(set1n_Fusobacteriota, "Species", fill="Species", facet_grid=group~day, title="p_Fusobacteriota")
#plot_bar(set1n_proteobacteria, x="Genus", fill = "Species", facet_grid = group~., title="Proteobacteria") +
#  geom_bar(aes(color=Species, fill=Species), stat="identity", position="stack")
dev.off()

#p_Proteobacteria
set1n_proteobacteria <- subset_taxa(set1n, Phylum %in% c("p_Proteobacteria"))
postscript("barplot/group_proteobacteria.eps",hor=T)
plot_bar(set1n_proteobacteria, "Genus", fill="Species", facet_grid=group~day, title="Proteobacteria")
#plot_bar(set1n_proteobacteria, x="Genus", fill = "Species", facet_grid = group~., title="Proteobacteria") +
#  geom_bar(aes(color=Species, fill=Species), stat="identity", position="stack")
dev.off()

#p_Patescibacteria
set1n_Patescibacteria <- subset_taxa(set1n, Phylum %in% c("p_Patescibacteria"))
postscript("barplot/group_p_Patescibacteria.eps",hor=T)
plot_bar(set1n_Patescibacteria, "Species", fill="Species", facet_grid=group~day, title="p_Patescibacteria")
plot_bar(set1n_Patescibacteria, x="Genus", fill = "Family", facet_grid=group~day, title="Proteobacteria") 
#  geom_bar(aes(color=Species, fill=Species), stat="identity", position="stack")
dev.off()

##p_Verrucomicrobiota
set1n_Verrucomicrobiota <- subset_taxa(set1n, Phylum %in% c("p_Verrucomicrobiota"))
postscript("barplot/group_Verrucomicrobiota.eps",hor=T)
plot_bar(set1n_Verrucomicrobiota, x="Species", fill = "Species", facet_grid = group~day, title="Verrucomicrobiota") 
#  geom_bar(aes(color=Species, fill=Species), stat="identity", position="stack")
#plot_bar(set1n_Firmicutes, x="Genus", fill = "Genus", facet_grid = group~., title="Firmicutes") +
#  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
dev.off()

###alpha
##https://micca.readthedocs.io/en/latest/phyloseq.html
pdf("barplot/phyla_experiment.pdf")
ps.rarefied = rarefy_even_depth(ps2, rngseed = 1, replace=F)
plot_bar(ps.rarefied, fill="Phylum") + facet_wrap(group~day, scales= "free_x",nrow=2)
plot_bar(ps.rarefied, fill="Phylum") + facet_wrap(group~experiment, scales= "free_x",nrow=2)
dev.off()

level_order <- c('3d', '7d', '30d') 
pdf("alpha/diversity_shannon.pdf")
plot_richness(ps.rarefied, x="day2", color="group", measures=c("Shannon")) + geom_boxplot(aes(x = factor(day2, level = level_order)))
dev.off()

pdf("alpha/diversity_observed.pdf")
plot_richness(ps.rarefied, x="day2", color="group", measures=c("Observed")) + geom_boxplot(aes(x = factor(day2, level = level_order)))
dev.off()

pdf("alpha/diversity_alpha.pdf")
plot_richness(ps.rarefied, x="day2", color="group", measures=c("Observed","Shannon")) + geom_boxplot(aes(x = factor(day2, level = level_order)))  
dev.off()

# Store as a new data variable
pAlpha=plot_richness(ps.rarefied, x="day2", color="group", measures=c("Observed","Shannon"))
alphadt = data.table(pAlpha$data)
# Subset to just the Shannon index
alphadt <- alphadt[(variable == "Shannon")]
# Order by Days
alphadt <- alphadt[order(day)][(is.na(se))]
# Define the plot
level_order <- c('3d', '7d', '30d') 
pdf("alpha/diversity_shannon_day.pdf")
ggplot(data = alphadt, 
       mapping = aes(day2, value,
                     color = group)) +
  # shape = ReportedAntibioticUsage)) + 
  geom_point(size = 5) + 
  geom_path() +
  geom_point(data = alphadt, 
             size = 8, alpha = 0.35) +
  #facet_wrap(~group, ncol = 2, scales = "free_y") +
  ylab("Shannon Index") +
  geom_boxplot(aes(x = factor(day2, level = level_order))) +
  ggtitle("Shannon Index after Rapamycin treatment")






sample_data(ps.rarefied) <- data.frame(sample_data(ps.rarefied), alpha_div)

alpha_div <- estimate_richness(ps.rarefied, measures = c("Shannon", "Observed"))
shannon_diversity <- sample_data(ps.rarefied)$Shannon
observed_diversity <- sample_data(ps.rarefied)$Observed
groups <- sample_data(ps.rarefied)$group
# Performing Kruskal-Wallis test
test_result <- kruskal.test(shannon_diversity ~ groups)
# Get the p-value
p_value <- test_result$p.value

### 
ps.sub3 <- subset_samples(ps.rarefied, day == 3)
ps.day3 <- prune_taxa(taxa_sums(ps.sub3)>0, ps.sub3)
#day7n = transform_sample_counts(ps.day7, function(x) x / sum(x) )
ps.day3.rarefied = rarefy_even_depth(ps.day3, rngseed = 1, replace=F)
rich.day3 = estimate_richness(ps.day3)
pairwise.wilcox.test(rich.day3$Shannon, sample_data(ps.day3.rarefied)$group)

alpha_div.3d <- estimate_richness(ps.day3.rarefied, measures = c("Shannon", "Observed"))
shannon_diversity.3d <- sample_data(ps.day3.rarefied)$Shannon
observed_diversity.3d <- sample_data(ps.day3.rarefied)$Observed
groups <- sample_data(ps.day3.rarefied)$group
# Performing Kruskal-Wallis test
shannon_result.3d <- kruskal.test(shannon_diversity.3d ~ groups)
observed_result.3d <- kruskal.test(observed_diversity.3d ~ groups)
# Get the p-value
p_value.shannon.3d <- shannon_result.3d$p.value
p_value.shannon.3d
p_value.observed.3d <- observed_result.3d$p.value
p_value.observed.3d

###
ps.sub7 <- subset_samples(ps.rarefied, day == 7)
ps.day7 <- prune_taxa(taxa_sums(ps.sub7)>0, ps.sub7)
#day7n = transform_sample_counts(ps.day7, function(x) x / sum(x) )
ps.day7.rarefied = rarefy_even_depth(ps.day7, rngseed = 1, replace=F)
rich.day7 = estimate_richness(ps.day7)
pairwise.wilcox.test(rich.day7$Shannon, sample_data(ps.day7.rarefied)$group)

alpha_div.7d <- estimate_richness(ps.day7.rarefied, measures = c("Shannon", "Observed"))
shannon_diversity.7d <- sample_data(ps.day7.rarefied)$Shannon
observed_diversity.7d <- sample_data(ps.day7.rarefied)$Observed
groups <- sample_data(ps.day7.rarefied)$group
# Performing Kruskal-Wallis test
shannon_result.7d <- kruskal.test(shannon_diversity.7d ~ groups)
observed_result.7d <- kruskal.test(observed_diversity.7d ~ groups)
# Get the p-value
p_value.shannon.7d <- shannon_result.7d$p.value
p_value.shannon.7d
p_value.observed.7d <- observed_result.7d$p.value
p_value.observed.7d

##30d
ps.sub30 <- subset_samples(ps.rarefied, day == 30)
ps.day30 <- prune_taxa(taxa_sums(ps.sub30)>0, ps.sub30)
#day7n = transform_sample_counts(ps.day7, function(x) x / sum(x) )
ps.day30.rarefied = rarefy_even_depth(ps.day30, rngseed = 1, replace=F)
rich.day30 = estimate_richness(ps.day30.rarefied)
pairwise.wilcox.test(rich.day30$Shannon, sample_data(ps.day30)$group)

alpha_div.30d <- estimate_richness(ps.day30.rarefied, measures = c("Shannon", "Observed"))
shannon_diversity.30d <- sample_data(ps.day30.rarefied)$Shannon
observed_diversity.30d <- sample_data(ps.day30.rarefied)$Observed
groups <- sample_data(ps.day30.rarefied)$group
# Performing Kruskal-Wallis test
shannon_result.30d <- kruskal.test(shannon_diversity.30d ~ groups)
observed_result.30d <- kruskal.test(observed_diversity.30d ~ groups)
# Get the p-value
p_value.shannon.30d <- shannon_result.30d$p.value
p_value.shannon.30d
p_value.observed.30d <- observed_result.30d$p.value
p_value.observed.30d

##https://github.com/ying14/yingtools2?tab=readme-ov-file
library("yingtools2")
library("microbiomeMarker")
##day 3
ps.sub3 <- subset_samples(ps2, day == 3)
ps.day3 <- prune_taxa(taxa_sums(ps.sub3)>10, ps.sub3)
ps.day3
day3n = transform_sample_counts(ps.day3, function(x) x / sum(x) )
ps.day3.rarefied = rarefy_even_depth(ps.day3, rngseed = 1, replace=F)
rich.day3 = estimate_richness(ps.day3.rarefied)
pairwise.wilcox.test(rich.day3$Shannon, sample_data(ps.day3.rarefied)$group)
#phyloseq2lefse(ps = ps.day3, covars = c("group"), taxa.levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), file.name = "day3.txt")
saveRDS(ps.day3, "ps.day3.rds")
saveRDS(day3n, "day3n.rds")
#set1n <- readRDS("set1n.rds")
phyloseq2lefse(ps = ps.day3, covars = c("group"), taxa.levels = c("Species"), file.name = "day3.txt")

phyloseq2lefse(ps = day3n, covars = c("group"), taxa.levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), file.name = "day3n.txt")
ps.day3.norm=normalize(ps.day3, method = "TSS") ##relative log expression, RLE uses a pseudo-reference calculated using the geometric mean of the gene-specific abundances over all samples. The scaling factors are then calculated as the median of the gene counts ratios between the samples and the reference.
#run_deseq2(day3n, group = "group")
mm_edger <- run_edger(
  ps.day3.norm,
  wilcoxon_cutoff = 0.05,
  group = "group",
  kw_cutoff = 0.05,
  multigrp_strat = TRUE,
  lda_cutoff = 4
)
mm_edger
 
mm_lefse3 <- run_lefse(
  day3n,
  wilcoxon_cutoff = 0.05,
  group = "group",
  kw_cutoff = 0.05,
  multigrp_strat = FALSE,
  lda_cutoff = 4
)
mm_lefse3
head(marker_table(mm_lefse))
write.csv(file="lefse/lefse_marker_3d.csv", marker_table(mm_lefse3))
p_abd3 <- plot_abundance(mm_lefse3, group = "group")
p_abd3

##7day
ps.sub7 <- subset_samples(ps2, day == 7)
ps.day7 <- prune_taxa(taxa_sums(ps.sub7)>10, ps.sub7)
ps.day7
day7n = transform_sample_counts(ps.day7, function(x) x / sum(x) )
ps.day7.rarefied = rarefy_even_depth(ps.day7, rngseed = 1, replace=F)
rich.day7 = estimate_richness(ps.day7.rarefied)
pairwise.wilcox.test(rich.day7$Shannon, sample_data(ps.day7.rarefied)$group)
#phyloseq2lefse(ps = ps.day3, covars = c("group"), taxa.levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), file.name = "day3.txt")
saveRDS(ps.day7, "ps.day7.rds")
saveRDS(day7n, "day7n.rds")
#set1n <- readRDS("set1n.rds")
phyloseq2lefse(ps = ps.day3, covars = c("group"), taxa.levels = c("Species"), file.name = "day3.txt")

phyloseq2lefse(ps = day3n, covars = c("group"), taxa.levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), file.name = "day3n.txt")
ps.day7.norm=normalize(ps.day7, method = "TSS") ##relative log expression, RLE uses a pseudo-reference calculated using the geometric mean of the gene-specific abundances over all samples. The scaling factors are then calculated as the median of the gene counts ratios between the samples and the reference.
#run_deseq2(day3n, group = "group")

mm_lefse7 <- run_lefse(
  day7n,
  wilcoxon_cutoff = 0.05,
  group = "group",
  kw_cutoff = 0.05,
  multigrp_strat = FALSE,
  lda_cutoff = 4
)
mm_lefse7
head(marker_table(mm_lefse7))
write.csv(file="lefse/lefse_marker_7d.csv", marker_table(mm_lefse7))
p_abd7 <- plot_abundance(mm_lefse7, group = "group")
pdf("lefse/day7.pdf")
p_abd7
dev.off()

##30day
ps.sub30 <- subset_samples(ps2, day == 30)
ps.day30 <- prune_taxa(taxa_sums(ps.sub30)>10, ps.sub30)
ps.day30
day30n = transform_sample_counts(ps.day30, function(x) x / sum(x) )
ps.day30.rarefied = rarefy_even_depth(ps.day30, rngseed = 1, replace=F)
rich.day30 = estimate_richness(ps.day30.rarefied)
pairwise.wilcox.test(rich.day30$Shannon, sample_data(ps.day30.rarefied)$group)
#phyloseq2lefse(ps = ps.day3, covars = c("group"), taxa.levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), file.name = "day3.txt")
saveRDS(ps.day30, "ps.day30.rds")
saveRDS(day30n, "day30n.rds")
#set1n <- readRDS("set1n.rds")
phyloseq2lefse(ps = ps.day30, covars = c("group"), taxa.levels = c("Species"), file.name = "day3.txt")

phyloseq2lefse(ps = day3n, covars = c("group"), taxa.levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), file.name = "day3n.txt")
ps.day30.norm=normalize(ps.day30, method = "TSS") ##relative log expression, RLE uses a pseudo-reference calculated using the geometric mean of the gene-specific abundances over all samples. The scaling factors are then calculated as the median of the gene counts ratios between the samples and the reference.
#run_deseq2(day3n, group = "group")

mm_lefse30 <- run_lefse(
  day30n,
  wilcoxon_cutoff = 0.05,
  group = "group",
  kw_cutoff = 0.05,
  multigrp_strat = FALSE,
  lda_cutoff = 4
)
mm_lefse30
head(marker_table(mm_lefse30))
write.csv(file="lefse/lefse_marker_30d.csv", marker_table(mm_lefse30))
p_abd30 <- plot_abundance(mm_lefse30, group = "group")
pdf("lefse/day30.pdf")
p_abd30
dev.off()



#lda <- lda.effect(day7n,class="group3",n_boots=NULL)
lda <- lda.effect(day3n,class="group",n_boots=NULL)
lda.plot(lda)
pdf(file="lefse/day3-circular.pdf")
#lda.clado(lda)
lda.clado(lda, font.size = 3, clade.label.height = 0.8)
dev.off()
pdf(file="lefse/day3-rectangular.pdf")
lda.clado(lda, layout="rectangular", font.size = 2, clade.label.height = 1)
dev.off()

lefse = run_lefse(
  day3n,
  wilcoxon_cutoff = 0.01,
  group = "group",
  kw_cutoff = 0.05,
  multigrp_strat = TRUE,
  lda_cutoff = 2
)
marker_table(lefse)
write.csv(file="lefse/markers_3day.csv",marker_table(lefse))



##day 7
ps.sub7 <- subset_samples(ps2, day == 7)
ps.day7 <- prune_taxa(taxa_sums(ps.sub7)>0, ps.sub7)
day7n = transform_sample_counts(ps.day7, function(x) x / sum(x) )
ps.day7.rarefied = rarefy_even_depth(ps.day7, rngseed = 1, replace=F)

phyloseq2lefse(ps = day7n, covars = c("group"), taxa.levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), file.name = "day7n.txt")
normalize(day7n, method = "RLE") ##relative log expression, RLE uses a pseudo-reference calculated using the geometric mean of the gene-specific abundances over all samples. The scaling factors are then calculated as the median of the gene counts ratios between the samples and the reference.
run_deseq2(day7n, group = "group")
#lda <- lda.effect(day7n,class="group3",n_boots=NULL)
lda <- lda.effect(day7n,class="group",n_boots=NULL)
pdf("lefse/lda_7d.pdf")
lda.plot(lda)
dev.off()
pdf(file="lefse/day7-circular.pdf")
#lda.clado(lda)
lda.clado(lda, font.size = 3, clade.label.height = 0.8)
dev.off()
pdf(file="lefse/day7-rectangular.pdf")
lda.clado(lda, layout="rectangular", font.size = 2, clade.label.height = 1)
dev.off()

lefse = run_lefse(
  day7n,
  wilcoxon_cutoff = 0.01,
  group = "group",
  kw_cutoff = 0.05,
  multigrp_strat = TRUE,
  lda_cutoff = 2
)
marker_table(lefse)
write.csv(file="lefse/markers_7day.csv",marker_table(lefse))

#####30d
ps.sub30 <- subset_samples(ps2, day == 30)
ps.day30 <- prune_taxa(taxa_sums(ps.sub30)>0, ps.sub30)
day30n = transform_sample_counts(ps.day30, function(x) x / sum(x) )

phyloseq2lefse(ps = day30n, covars = c("group"), taxa.levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), file.name = "day30n.txt")
normalize(day30n, method = "RLE") ##relative log expression, RLE uses a pseudo-reference calculated using the geometric mean of the gene-specific abundances over all samples. The scaling factors are then calculated as the median of the gene counts ratios between the samples and the reference.
run_deseq2(day30n, group = "group")
#lda <- lda.effect(day7n,class="group3",n_boots=NULL)
lda <- lda.effect(day30n,class="group",n_boots=NULL)
pdf("lefse/lda_30d.pdf")
lda.plot(lda)
dev.off()
pdf(file="lefse/day30-circular.pdf")
#lda.clado(lda)
lda.clado(lda, font.size = 3, clade.label.height = 0.8)
dev.off()
pdf(file="lefse/day30-rectangular.pdf")
lda.clado(lda, layout="rectangular", font.size = 2, clade.label.height = 1)
dev.off()

lefse = run_lefse(
  day30n,
  wilcoxon_cutoff = 0.05,
  group = "group",
  kw_cutoff = 0.05,
  multigrp_strat = TRUE,
  lda_cutoff = 2
)
marker_table(lefse)
write.csv(file="lefse/markers_30day.csv",marker_table(lefse))

#https://yiluheihei.github.io/microbiomeMarker/articles/microbiomeMarker-vignette.html

plot_cladogram(lefse, color = c("darkgreen", "red"))




ps.sub7 <- subset_samples(ps2, day == 7)
ps.day7 <- prune_taxa(taxa_sums(ps.sub7)>0, ps.sub7)
#day7n = transform_sample_counts(ps.day7, function(x) x / sum(x) )
ps.day7.rarefied = rarefy_even_depth(ps.day7, rngseed = 1, replace=F)
rich.day7 = estimate_richness(ps.day7.rarefied)
pairwise.wilcox.test(rich.day7$Shannon, sample_data(ps.day7.rarefied)$group)

ps.sub30 <- subset_samples(ps2, day == 30)
ps.day30 <- prune_taxa(taxa_sums(ps.sub30)>0, ps.sub30)
#day7n = transform_sample_counts(ps.day7, function(x) x / sum(x) )
ps.day30.rarefied = rarefy_even_depth(ps.day30, rngseed = 1, replace=F)
rich.day30 = estimate_richness(ps.day30.rarefied)
pairwise.wilcox.test(rich.day30$Shannon, sample_data(ps.day30.rarefied)$group, p.adjust.method = "bonf")

ps.subRapa <- subset_samples(ps2, group == "Rapamycin")
ps.rapa <- prune_taxa(taxa_sums(ps.subRapa)>0, ps.subRapa)
#day7n = transform_sample_counts(ps.day7, function(x) x / sum(x) )
ps.rapa.rarefied = rarefy_even_depth(ps.rapa, rngseed = 1, replace=F)
rich.rapa = estimate_richness(ps.rapa.rarefied)
pairwise.wilcox.test(rich.rapa$Shannon, sample_data(ps.rapa.rarefied)$day, p.adjust.method = "bonf")

##beta
library(dplyr)
library(ggplot2)
library(patchwork)
library(phyloseq)
library(microbiome)
library(microViz)
set.seed(1)

ps2 %>%
  tax_fix() %>%
  tax_agg(rank = "Family")

ps2 %>%
  tax_fix() %>%
  tax_agg(rank = "Genus")

ps3=ps2 %>% tax_fix(unknowns = c("g_unassigned"))

shao4d_psX <- ps3 %>%
  # keep only taxa belonging to genera that have over 100 counts in at least 5% of samples
  tax_filter(min_prevalence = 0.05, undetected = 100, tax_level = "Genus") %>%
  # aggregate counts at genus-level & transform with robust CLR transformation
  tax_transform(trans = "rclr", rank = "Genus") %>%
  # no distances are needed for PCA: so skip dist_calc and go straight to ord_calc
  ord_calc(method = "PCA")
#> Proportional min_prevalence given: 0.05 --> min 16/306 samples.

pdf("beta/PCA.pdf")
PCA_plot <- shao4d_psX %>%
  ord_plot(
    colour = "day2", shape = "group", size = 5,
    plot_taxa = 10:1,
    tax_vec_length = 0.3,
    tax_lab_style = tax_lab_style(
      type = "label", max_angle = 90, aspect_ratio = 1,
      size = 4, alpha = 0.8, fontface = "bold", # style the labels
      label.r = unit(0, "mm") # square corners of labels - see ?geom_label
    )
  ) +
  coord_fixed(ratio = 1, clip = "off", xlim = c(-4, 4))
# match coord_fixed() ratio to tax_lab_style() aspect_ratio for correct text angles
PCA_plot
dev.off()

PCA_plot2 <- shao4d_psX %>%
  ord_plot(
    colour = "group", shape = "day2", size = 5,
    plot_taxa = 10:1,
    tax_vec_length = 0.3,
    tax_lab_style = tax_lab_style(
      type = "label", max_angle = 90, aspect_ratio = 1,
      size = 4, alpha = 0.8, fontface = "bold", # style the labels
      label.r = unit(0, "mm") # square corners of labels - see ?geom_label
    )
  ) +
  coord_fixed(ratio = 1, clip = "off", xlim = c(-4, 4))
# match coord_fixed() ratio to tax_lab_style() aspect_ratio for correct text angles
PCA_plot2

pdf("beta/PCA_frame.pdf")
PCA_plot_custom2 <- PCA_plot2 +
  # add a convex hull around the points for each group, to aid the eye
  stat_chull(
    mapping = aes(colour = group, linetype = day2),
    linewidth = 0.5, alpha = 0.5, show.legend = FALSE
  ) +
  scale_shape_manual(
    name = "Day",
    values = c("circle", "circle open", "circle plus"), labels = c("Day 3", "Day 7", "Day 30")
  ) +
  # set a custom colour scale and customise the legend order and appearance
  scale_color_manual(
    name = "Day",
    values = c("black", "grey45","grey86"), labels = c("Day 3", "Day 7", "Day 30"),
    guide = guide_legend(override.aes = list(size = 4))
  ) +
  # add a title and delete the automatic caption
  labs(title = "4-day-old gut microbiota and birth mode", caption = NULL) +
  # put the legend at the bottom and draw a border around it
  theme(legend.position = "bottom", legend.background = element_rect())
PCA_plot_custom2

PCA_plot_custom <- PCA_plot +
  # add a convex hull around the points for each group, to aid the eye
  stat_chull(
    mapping = aes(colour = day2, linetype = day2),
    linewidth = 0.5, alpha = 0.5, show.legend = FALSE
  ) +
  scale_shape_manual(
    name = "Group",
    values = c("circle", "circle open"), labels = c("Rapamycin","Ctl")
  ) +
  # set a custom colour scale and customise the legend order and appearance
  scale_color_manual(
    name = "Day",
    values = c("black", "grey45","grey86"), labels = c("Day 30", "Day 3", "Day 7"),
    guide = guide_legend(override.aes = list(size = 4))
  ) +
  # add a title and delete the automatic caption
  labs(title = "Gut microbiota after Rapamycin treatment", caption = NULL) +
  # put the legend at the bottom and draw a border around it
  theme(legend.position = "bottom", legend.background = element_rect())
PCA_plot_custom
dev.off()

ps.sub30 <- subset_samples(ps2, day == 30)
ps.day30 <- prune_taxa(taxa_sums(ps.sub30)>0, ps.sub30)
ps2.day30=ps.day30 %>% tax_fix(unknowns = c("g_unassigned"))
ps2.day30 %>%
  tax_filter(min_prevalence = 5 / 100, tax_level = "Genus") %>%
  tax_agg("Genus") %>%
  dist_calc("aitchison") %>%
  dist_permanova(
    variables = c("group"),
    n_perms = 99, # this is a low number of permutations for speed, you should set more e.g. 9999
    seed = 12345, complete_cases = TRUE, verbose = "max"
  )
phylo <- ps2.day30 %>%
  ps_filter(group %in% c("Rapamycin", "Ctl")) %>%
  tax_mutate(Species = NULL) %>%
  tax_fix()
phylo %>%
  samdat_tbl() %>%
  dplyr::mutate(across(where(is.character), as.factor)) %>%
  skimr::skim()

interestingGenera <- ps3 %>%
tax_select("p_Bacteroidota") %>%
  tax_top(n = 10, rank = "Genus")
interestingGenera
ps3 %>%
  ps_filter(day == 30) %>%
  tax_sort(by = sum, at = "Genus") %>% # put other taxa in a reasonable order
  comp_barplot(
    tax_level = "Genus", n_taxa = 15, merge_other = FALSE, other_name = "Other",
    palette = distinct_palette(15, pal = "kelly", add = "grey90"),
    tax_order = interestingGenera # this is the reordering magic
  ) +
  coord_flip()
ps3 %>%
  ps_select(group) %>% 
  ps_filter(day == 30) %>%
  phyloseq::merge_samples(group = "Rapamycin") %>%
  comp_barplot(tax_level = "Genus", n_taxa = 12, bar_width = 0.8) +
  coord_flip() + labs(x = NULL, y = NULL)

pdf("barplot/day3.pdf")
ps3 %>%
  ps_filter(day == 3) %>% # only one sample in this group
  # convert DiseaseState into ordered factor to control order of facets
  ps_mutate(
    groupState = factor(group, levels = c("Rapamycin", "Ctl"))
  ) %>%
  comp_barplot(
    tax_level = "Genus", n_taxa = 15,
    bar_outline_colour = NA, facet_by = "group"
  ) +
  coord_flip()
dev.off()

pdf("barplot/day7.pdf")
ps3 %>%
  ps_filter(day == 7) %>% # only one sample in this group
  # convert DiseaseState into ordered factor to control order of facets
  ps_mutate(
    groupState = factor(group, levels = c("Rapamycin", "Ctl"))
  ) %>%
  comp_barplot(
    tax_level = "Genus", n_taxa = 15,
    bar_outline_colour = NA, facet_by = "group"
  ) +
  coord_flip()
dev.off()

pdf("barplot/day30.pdf")
ps3 %>%
  ps_filter(day == 30) %>% # only one sample in this group
  # convert DiseaseState into ordered factor to control order of facets
  ps_mutate(
    groupState = factor(group, levels = c("Rapamycin", "Ctl"))
  ) %>%
  comp_barplot(
    tax_level = "Genus", n_taxa = 15,
    bar_outline_colour = NA, facet_by = "group"
  ) +
  coord_flip()
dev.off()
