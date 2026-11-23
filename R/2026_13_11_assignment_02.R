##***************************
## Assignment 1, Trichinella nativa comparing with Trichinella Spiralis
##
## Jordyn Anene  
## Luis Humberto Villasana Avila  
## Shariq Mumtaz Hashmi
## Sumayyah Khan
##
## 2026-11-02
##
##***************************

#comments made by me (mikayla) start #mk#
#downloads are from current wd, (working in a project)

## Cleaning the environment to make sure that the code can run on its own ----
rm(list = ls())

## _ Packages used -------
library(tidyverse)
conflicted::conflicts_prefer(dplyr::filter())
library(viridis)
library(Biostrings)
# + scale_color/fill_viridis_c/d()
theme_set(theme_light())

## Preparing CSV files ----
# Download T. nativa into script and create object 
df_te_nat = read.csv("../data/GCA_001447565.2.csv")
df_te_nat = read_csv("Data/GCA_001447565.2.csv")
# Download T. spiralis into script and create object 
df_trich_te = read.csv("../Data/trichinella_spiralis_te.csv")
df_trich_te = read_csv("Data/trichinella_spiralis_te.csv")
#mk# downloading Trichinella T9 (s9)
df_s9 <- read_csv("data/GCA_001447505.1(1).csv")

#fixing the issue with qry_left and in_repeat_end columns for Nativa:
df_te_nat = df_te_nat |>
  mutate(qry_left2 = as.numeric(str_replace_all(qry_left, "\\(|\\)", "")), .after = qry_left) |>
  mutate(in_repeat_end2 = as.numeric(str_replace_all(in_repeat_end, "\\(|\\)", "")), .after = in_repeat_end) |> 
  mutate(no_bp_in_complement2 = as.numeric(str_replace_all(no_bp_in_complement, "\\(|\\)", "")), .after = no_bp_in_complement)

#fixing the issue with qry_left and in_repeat_end columns for Spiralis:
df_trich_te = df_trich_te |>
  mutate(qry_left2 = as.numeric(str_replace_all(qry_left, "\\(|\\)", "")), .after = qry_left) |>
  mutate(in_repeat_end2 = as.numeric(str_replace_all(in_repeat_end, "\\(|\\)", "")), .after = in_repeat_end) |> 
  mutate(no_bp_in_complement2 = as.numeric(str_replace_all(no_bp_in_complement, "\\(|\\)", "")), .after = no_bp_in_complement)
#mk#added the same thing to no_bp_in_compliment

#mk#now doing the same for s9
df_s9 = df_s9 |>
  mutate(qry_left2 = as.numeric(str_replace_all(qry_left, "\\(|\\)", "")), .after = qry_left) |>
  mutate(in_repeat_end2 = as.numeric(str_replace_all(in_repeat_end, "\\(|\\)", "")), .after = in_repeat_end) |> 
  mutate(no_bp_in_complement2 = as.numeric(str_replace_all(no_bp_in_complement, "\\(|\\)", "")), .after = no_bp_in_complement)

## Removing non-TE elements and classifying matching classes into either NA transposons or retrotransposons----

# List of non-TE elements
ls_no_te = c("Simple_repeat", "simple_repeat", "Low_complexity", "Satellite", "rRNA", "snRNA", "ARTEFACT", "Unknown")
ls_no_te
#mk# these two vectors were the same, so i removed one for simplicity

#Creating a copy of the original file for possible further use before updating the original function for Nativa
df_te_nat_org = df_te_nat

#Creating a copy of the original file for possible further use before updating the original function for Spiralis
df_trich_te_org = df_trich_te

#mk#doing the same for s9
df_s9_org <- df_s9

# Updating dataset with only TEs for Nativa
df_te_nat = df_te_nat |>
  filter(matching_class %notin% ls_no_te)
df_te_nat

# Updating dataset with only TEs for Spiralis
df_trich_te = df_trich_te |>
  filter(matching_class %notin% ls_no_te)
df_trich_te

#mk#same for s9
df_s9 <- df_s9 |> 
  filter(matching_class %notin% ls_no_te)

# Making sure there are not any non-TE left in Nativa
df_te_nat$matching_class %in% ls_no_te |>
  sum()
df_te_nat |>
  filter(matching_class %notin% ls_no_te) |> View()

# Making sure there are not any non-TE left in Spiralis
df_trich_te$matching_class %in% ls_no_te |>
  sum()

# To observe what are the TE elements present in Nativa
df_te_nat$matching_class |> unique()

# To observe what are the TE elements present in Spiralis
df_trich_te$matching_class |> unique()

# To identify how many unique qry_ids we have in Nativa
df_te_nat$qry_id |> unique() |> length()
#=============> we have 438 unique contigs

# To identify how many unique qry_ids we have in Nativa
df_trich_te$qry_id |> unique() |> length()
#=============> we have 315 unique contigs 
#mk# my answer is 243, but i'm not sure why they're different

# To obtain a list of the dna transposons for T. nativa 
ls_dna_qry = df_te_nat |> 
  filter(str_starts(matching_class, "DNA/")) |>
  pull(qry_id) |> unique()
ls_dna_qry

# To obtain a list of the dna transposons for T. spiralis 
ls_ts_dna_qry = df_trich_te |>
  filter(str_starts(matching_class, "DNA/|RC")) |>
  pull(qry_id) |> unique()
ls_ts_dna_qry

# To obtain a list of a all the retrotransposons in Nativa
ls_retro_qry = df_te_nat |>
  filter(str_starts(matching_class, "PLE|PLE/|LTR/|LTR|LINE|LINE/|SINE/")) |>
  pull(qry_id) |> unique()
ls_retro_qry

# To obtain a list of a all the retrotransposons in Spiralis
ls_ts_retro_qry = df_trich_te |>
  filter(str_starts(matching_class, "PLE|PLE/|LTR/|LTR|LINE|LINE/|SINE/")) |>
  pull(qry_id) |> unique()
ls_ts_retro_qry

# Checking that the number of DNA transposons and retrotransposons match in Nativa
Biostrings::setdiff(ls_dna_qry, ls_retro_qry) |> length()
Biostrings::setdiff(ls_retro_qry, ls_dna_qry) |> length()

# Checking that the number of DNA transposons and retrotransposons match in Spiralis
Biostrings::setdiff(ls_ts_dna_qry, ls_ts_retro_qry) |> length()
Biostrings::setdiff(ls_ts_retro_qry, ls_ts_dna_qry) |> length()


## Creating different data frames according to the type of TE ----
#Data frame with only DNA transponsos in Nativa
df_dna_only = df_te_nat |>
  filter(qry_id %in% Biostrings::setdiff(ls_dna_qry, ls_retro_qry))
df_dna_only

#Data frame with only DNA transponsos in Spiralis
df_ts_dnaonly = df_trich_te |>
  filter(qry_id %in% Biostrings::setdiff(ls_ts_dna_qry, ls_ts_retro_qry))
df_ts_dnaonly

#Data frame with only retrotransponsos in Nativa
df_retro_only = df_te_nat |>
  filter(qry_id %in% Biostrings::setdiff(ls_retro_qry, ls_dna_qry))
df_retro_only

#Data frame with only retrotransponsos in spiralis
df_ts_retroonly = df_trich_te |>
  filter(qry_id %in% Biostrings::setdiff(ls_ts_retro_qry, ls_ts_dna_qry))
df_ts_retroonly

# Data frame with both DNA transponsos and retrotransponsos in Nativa
df_bothTE_nat = df_te_nat |>
  filter(qry_id %in% Biostrings::intersect(ls_dna_qry, ls_retro_qry))
df_bothTE_nat

# Data frame with both DNA transponsos and retrotransponsos in Spiralis
df_ts_both = df_trich_te |>
  filter(qry_id %in% Biostrings::intersect(ls_ts_dna_qry, ls_ts_retro_qry))
df_ts_both


## Finding the inflection point for large contigs ----

# For T.Spiralis
# Creating a plot to determine the inflection point to identify the size of large contigs in spiralis
df_trich_te_org |> 
  group_by(qry_id) |> 
  summarise(count = n()) |> 
  arrange(desc(count)) |> 
  pull(count) |> 
  plot()

# Using log10 to make more clear where the inflection point is located in spiralis
df_trich_te_org |> 
  group_by(qry_id) |> 
  summarise(count = n()) |> 
  arrange(desc(count)) |> 
  pull(count) |> 
  log10() |>
  plot()
#================> We determine that large contigs are above 50 in width


# Determine the inflection point to identify the size of large contigs in Nativa
df_te_nat_org |>
  group_by(qry_id) |>
  summarise(count = n()) |>
  arrange(desc(count)) |>
  pull(count) |>
  plot()

# Using log10 to make more clear where the inflection point is located in Nativa
df_te_nat_org |>
  group_by(qry_id) |>
  summarise(count = n()) |>
  arrange(desc(count)) |>
  pull(count) |>
  log10() |>
  plot()
#================> We determine that large contigs are above 50 in width

#================> Cutoff for repeat-dense contigs: >50 repeats per contig (both species).
#================> the curve bends from a gentle to a steep decline at ~50
#================> Repeats per contig (all repeat classes >>>> org tables) were sorted and plotted on a log10 scale;

#inflexion point in s9
df_s9_org |>
  group_by(qry_id) |>
  summarise(count = n()) |>
  arrange(desc(count)) |>
  pull(count) |>
  log10() |>
  plot()
#================> about 300

# To obtain from the original cvs file the contigs larger than 50 in nativa
ls_dense_nat = df_te_nat_org |>
  group_by(qry_id) |>
  summarise(count = n()) |>
  arrange(desc(count)) |>
  filter(count > 50) |>
  pull(qry_id)
ls_dense_nat

# To obtain number of large contigs present in nativa
ls_dense_nat |> length()
#=================> In T.nativa there are 268 contigs

# To obtain from the original cvs file the contigs larger than 50 in Spiralis 
ls_dense_ts = df_trich_te_org |>
  group_by(qry_id) |>
  summarise(count = n()) |>
  arrange(desc(count)) |>
  filter(count > 50) |>
  pull(qry_id)

# To obtain number of large contigs present in Spiralis
ls_dense_ts |> length()
#=================> In T.spiralis there are 220 contigs 
#mk# i only have 14

#=================>  The same cutoff is used for both species so the comparison is consistent.
#=================>  Note: we haven't used FASTA files, contigs are ranked by repeat count, not by length, so "repeat dense" rather than "largest".

#mk#same for s9
ls_dense_s9 = df_s9_org |>
  group_by(qry_id) |>
  summarise(count = n()) |>
  arrange(desc(count)) |>
  filter(count > 300) |>
  pull(qry_id)

# Creating a new table with TE rows that are in the ls_dense for which species, so we can plot them against perc_del and perc_div

df_dense_nat = df_te_nat |>
  filter(qry_id %in% ls_dense_nat) |> 
  select(!qry_left) |> 
  select(!no_bp_in_complement) |> 
  select(!in_repeat_end)
#mk#added the select to get rid of the columns

df_dense_ts = df_trich_te |>
  filter(qry_id %in% ls_dense_ts) |> 
  select(!qry_left) |> 
  select(!no_bp_in_complement) |> 
  select(!in_repeat_end)
#mk#same here

#mk#making the same for s9
df_dense_s9 = df_s9 |>
  filter(qry_id %in% ls_dense_s9) |> 
  select(!qry_left) |> 
  select(!no_bp_in_complement) |> 
  select(!in_repeat_end)

nrow(df_te_nat)
nrow(df_dense_nat)

nrow(df_trich_te)
nrow(df_dense_ts)

#====> Dense contigs hold 823 of 1319 TEs in T.Nativa and 838 of 1037 TEs in T.Spiralis.
#mk#946 of 1268 in spiralis

## Creating new columns classifying contigs as either DNA Transposons or Retrotransposons and Nativa or spiralis to later create a single table to make comparisons and plotting easier ----

#Adding a new column classifying contigs as either DNA Transposons or Retrotransposons in Nativa
df_dense_nat = df_dense_nat |>
  mutate(TE_type = if_else(str_starts(matching_class, "DNA/|RC"), "DNA Transposons", "Retrotransposons"))

#Adding a new column classifying contigs as either DNA Transposons or Retrotransposons in Spiralis
df_dense_ts = df_dense_ts |>
  mutate(TE_type = if_else(str_starts(matching_class, "DNA/|RC"), "DNA Transposons", "Retrotransposons"))

#mk# same for s9
df_dense_s9 = df_dense_s9 |>
  mutate(TE_type = if_else(str_starts(matching_class, "DNA/|RC"), "DNA Transposons", "Retrotransposons"))

# Identifying the amount of DNA Transposons and Retrotransposons in both Nativa and Spiralis to check values
df_dense_nat |>
  count(TE_type)

df_dense_ts |>
  count(TE_type)

df_dense_nat$qry_id |> unique()
df_dense_ts$qry_id |> unique()

df_dense_nat |>
  filter(TE_type == "Retrotransposons") |>
  count(matching_class)  

df_dense_ts |>
  filter(TE_type == "Retrotransposons") |>
  count(matching_class)

#creating a new table with both species labelled in the new column
df_dense_nat = df_dense_nat |>
  mutate(Species = "T.nativa")

df_dense_ts = df_dense_ts |>
  mutate(Species = "T.spiralis")

df_both_species = bind_rows(df_dense_nat, df_dense_ts)
#mk#changed some things upstream so I can run this

nrow(df_both_species)
df_both_species |> count(Species, TE_type)

#mk#adding s9
df_dense_s9 = df_dense_s9 |>
  mutate(Species = "s9")

df_all3 <- bind_rows(df_dense_nat, df_dense_s9, df_dense_ts)

df_all3 |> count(Species, TE_type)

## Determine the timings of the TEs ----
# Do DNA transposons and retrotransposons differ in age, and is that pattern the same in T.nativa and T.spiralis?
# Do older TEs appear more decayed? Use the relationship between perc_div and perc_del

# Exploratory: Checking for TE divergence
df_both_species |> 
  ggplot(aes(x = perc_div, fill = TE_type)) +
  geom_histogram(binwidth = 1, color = "white") +
  labs(x = "Divergence from reference sequence %", y = "Number of TE copies", title = "Distribution of TE divergence by type and species") +
  facet_grid(TE_type ~ Species, scales = 'free_y') +
  theme_light() +
  scale_fill_brewer(palette = "Pastel2") +
  theme(legend.position = "none")

ggsave("Output/Fig4_TE_Species_dist.png", width = 10, height = 6)


#==> T. nativa: no impossible values (all between ~2% and ~37% divergence); data passes the check.
#==> Retrotransposons: one main peak at ~20%; a few very old copies (~37%) separated from the main distribution,
#    and a thin tail of younger copies (2–8%) too small to indicate a recent burst. Outliers kept as genuine values.
#==> DNA transposons: patchy distribution (n = 85), a few younger copies at ~9–10%; gaps reflect small sample size.

#mk# part 1 for s9
df_all3 |> 
  ggplot(aes(x = perc_div, fill = TE_type)) +
  geom_histogram(binwidth = 1, color = "white") +
  labs(x = "Divergence from reference sequence %", y = "Number of TE copies", title = "Distribution of TE divergence by type and species") +
  facet_grid(TE_type ~ Species, scales = 'free_y') +
  theme_light() +
  scale_fill_brewer(palette = "Pastel2") +
  theme(legend.position = "none")

ggsave("Output/fig_spec_dist_s9.png", width = 10, height = 6)
#==> s9 seems to have fewer overall copies than either of the other two species


# Plotting perc_div (substitution, as a proxy for age) against perc_del (deletion) for both species.
# To examine whether older TE copies carry more deletions, and whether this relationship differs between DNA transposons and Retrotransposon, and T.spiralis species and T.nativa species.

df_both_species |>
  ggplot(aes(x = perc_div, y = perc_del, colour = TE_type)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = lm) +
  facet_wrap(~ Species) +
  labs(x = "Divergence from the reference sequence in %", y = "Deletion in %", color = "TE type")

ggsave("Output/fig_div_vs_del.png", width = 10, height = 6)

#mk# plotting so I can look at dna and retro against each other by species type
df_both_species |>
  ggplot(aes(x = perc_div, y = perc_del, colour = Species)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = lm) +
  facet_wrap(~ TE_type) +
  labs(x = "Divergence from the reference sequence in %", y = "Deletion in %", color = "TE type")

ggsave("Output/fig_divdel_spec.png", width = 10, height = 6)

#==> now we see the differences in the regression lines
#==> the line in retrotransposons is very similar between the two species. This graph shows a weak positive link between deletion and divergence
#==> dna transposons show opposing trends, which may indicate each species is evolving their TEs differently, but the confidence intervals are wide, so there is a lot of variability

# In both species, older retrotransposon copies show slightly more deletion, whereas DNA transposons show no clear relationship; the overlapping confidence bands indicate the differences between TE types and species are weak.

#mk# now trying it with s9 added
df_all3 |>
  ggplot(aes(x = perc_div, y = perc_del, colour = TE_type)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = lm) +
  facet_wrap(~ Species) +
  labs(x = "Divergence from the reference sequence in %", y = "Deletion in %", color = "TE type")

ggsave("Output/fig_div_vs_del_s9.png", width = 10, height = 6)

#mk# and with comparison by species
df_all3 |>
  ggplot(aes(x = perc_div, y = perc_del, colour = Species)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = lm) +
  facet_wrap(~ TE_type) +
  labs(x = "Divergence from the reference sequence in %", y = "Deletion in %", color = "TE type")

ggsave("Output/fig_div_vs_del_spec_s9.png", width = 10, height = 6)
#==> s9 has a larger trend in DNA transposons, but the retrotransposon trend is very similar to the other two species

# Part 2
# Are DNA transposons and retrotransposons different ages, and is that difference the same in T.native and T.spiralis?

#Statistic analysis: Using wilcoxon ranked sum tested, based on the histogram we generated we noticed that the distibution of DNA TEs is regular. With wilcoxon ranked sum we will check whether the divergence differ between the two classes of TE (DNA vs Retrotransposon).

help("wilcox.test")

wilcox.test(perc_div ~ TE_type, data = df_both_species |> filter(Species == "T.nativa"))
#In T.nativa, divergence did not differ between DNA transposons and retrotransposons (Wilcoxon rank-sum test, W = 29899, p = 0.48; medians ≈ 21% for both).

wilcox.test(perc_div ~ TE_type, data = df_both_species |> filter(Species == "T.spiralis"))

# DNA transposons and retrotransposons did not differ in divergence in T.nativa, but DNA transposons were slightly more diverged in T.spiralis (W = 41606, p = 0.038). The T. spiralis difference is small (~1.4 percentage points), so the TE-type age pattern is broadly similar between species.      #=======> The small T.spiralis difference (p = 0.038) does not survive correction for running two tests, so it should be treated as a tentative trend, not a firm result.
#mk#(w = 45006, p = 0.4518)

#mk# now calculating for s9
wilcox.test(perc_div ~ TE_type, data = df_all3 |> filter(Species == "s9"))
#==> (w = 2904, p = 0.4932)
#==> medians ~20.5% for both

#  Comparing divergence between species within each TE type, to test directly whether TE timing differs between T.nativa and T.spiralis.

wilcox.test(perc_div ~ Species, data = df_both_species |> filter(TE_type == "Retrotransposons"))

wilcox.test(perc_div ~ Species, data = df_both_species |> filter(TE_type == "DNA Transposons"))

df_median_div = df_both_species |>
  group_by(Species, TE_type) |>
  summarise(med = median(perc_div), n = n())

#mk# can't add s9 into the data, wilcox tests is a 2-sample test

df_median_div2 = df_all3 |>
  group_by(Species, TE_type) |>
  summarise(med = median(perc_div), n = n())

#=====> Medians were calculated per species and TE type to report the size of the age difference alongside the Wilcoxon p-values, since a p-value alone shows whether groups differ, not by how much. ANd to add captions to the Boxplot.

df_both_species |> filter(Species == "T.nativa") |>
  ggplot(aes(x = TE_type, y = perc_div, fill = TE_type)) +
  scale_fill_brewer(palette = "Pastel2") +
  geom_boxplot(width = 0.4, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.4, size = 0.5, colour = "grey30") +
  theme_light()

df_both_species |> filter(Species == "T.spiralis") |>
  ggplot(aes(x = TE_type, y = perc_div, fill = TE_type)) +
  scale_fill_brewer(palette = "Pastel2") +
  geom_boxplot(width = 0.4, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.4, size = 0.5, colour = "grey30") +
  theme_light()

df_both_species |>
  ggplot(aes(x = TE_type, y = perc_div, fill = TE_type)) +
  geom_violin(alpha = 0.4) +
  scale_fill_brewer(palette = "Pastel2") +
  geom_boxplot(width = 0.2, outlier.shape = NA, fill = "white") +
  geom_jitter(width = 0.1, alpha = 0.2, size = 0.8, colour = "grey30") +
  facet_wrap(~ Species) +
  labs(x = "TE type", y = "Divergence from reference sequence (%)", title = "TE divergence by type and species") +
  labs(caption = "T. nativa: DNA 21.2% (n = 85) vs Retro 21.3% (n = 738), Wilcoxon p = 0.48\nT. spiralis: DNA 22.8% (n = 100) vs Retro 21.4% (n = 738), Wilcoxon p = 0.038\nBetween species: Retro p = 0.76, DNA p = 0.030") +
  theme(legend.position = "none")

ggsave("Output/fig_div_TEtype.png", width = 10, height = 6)

#Histograms showed irregular distributions for DNA transposons and unequal group sizes, so a non-parametric Wilcoxon test was used instead of a t-test.

#==> Median divergence is around 21% for both TE types in T.nativa, so DNA transposons and retrotransposons are of similar age. In T.spiralis, DNA transposons sit slightly higher (~23% vs ~21%), suggesting slightly older DNA TE copies.
#==> Retrotransposons look almost identical in both species (median ~21%, range ~3–38%), and span a wider range of ages than DNA TEs, i.e. activity over a longer period.
#==> All violins have a single bulge at ~20–25% and no peak at low divergence: no evidence of a recent burst of TE activity in either species; most copies are old.
#==> DNA groups are small (n = 85, 100) vs retros (n = 738 each).

#mk# adding s9
df_all3 |>
  ggplot(aes(x = TE_type, y = perc_div, fill = TE_type)) +
  geom_violin(alpha = 0.4) +
  scale_fill_brewer(palette = "Pastel2") +
  geom_boxplot(width = 0.2, outlier.shape = NA, fill = "white") +
  geom_jitter(width = 0.1, alpha = 0.2, size = 0.8, colour = "grey30") +
  facet_wrap(~ Species) +
  labs(x = "TE type", y = "Divergence from reference sequence (%)", title = "TE divergence by type and species") +
  labs(caption = "T. nativa: DNA 21.2% (n = 85) vs Retro 21.3% (n = 738), Wilcoxon p = 0.48\nT. spiralis: DNA 22.8% (n = 100) vs Retro 21.4% (n = 738), Wilcoxon p = 0.038\ns9: DNA 20.5% (n = 27) vs Retro 20.85% (n = 234)") +
  theme(legend.position = "none")

ggsave("Output/fig_div_TEtype_s9.png", width = 10, height = 6)
#==> s9 has a lower divergence average in DNA transposons than the other species
#==> retrotransposon trends look similar

# Part 3

#Do older TE copies tend to be shorter, and does that differ between DNA transposons and retrotransposons, and between T. nativa and T. spiralis?

df_both_species |>
  ggplot(aes(x = perc_div, y = qry_width, colour = TE_type)) +
  scale_y_log10() +
  geom_point(size = 0.8, alpha = 0.3) +
  geom_smooth(method = lm) +
  facet_grid(TE_type ~ Species, scales = "free_y") +
  labs(x = "Divergence from the reference sequence (%)",
       y = "TE copy length (bp, log scale)",
       title = "Older TE copies appear longer",
       caption = "Each point is one TE copy; lines are linear fits with 95% confidence bands.\nn: T. nativa DNA = 85, Retro = 738; T. spiralis DNA = 100, Retro = 738") +
  theme_light() +
  theme(legend.position = "none")

ggsave("Output/fig_div_TEwidth_4.png", width = 10, height = 8)

#==> Most TE copies are short fragments (<200 bp), far below full length TE size; consistent with old, decayed TEs. Contrary to expectation, copy length increases with divergence, clearly for retrotransposons (both species, n = 738), weakly for DNA transposons (wide confidence bands, n = 85, 100).
#==> Likely bias: highly diverged copies must be longer to pass RepeatMasker's score threshold, so old short fragments are under-detected. Length vs divergence therefore reflects detectability as well as decay.

#mk#seeing if the trend changes in certain repeat families
#mk# taking the most commo one so the graph is still readable
top_families <- df_both_species |>
  count(matching_class, sort = TRUE) |>
  slice_head(n = 9) |>
  pull(matching_class)

df_both_species |>
  filter(matching_class %in% top_families) |>
  ggplot(aes(x = perc_div, y = qry_width, colour = Species)) +
  scale_y_log10() +
  geom_point(size = 0.8, alpha = 0.4) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ matching_class, scales = "free") +
  labs(x = "Divergence from consensus (%)",
       y = "TE copy length (bp, log scale)",
       title = "Length vs. Divergence Within Individual TE Families") +
  theme_light()
#==> most TE famlies are still showing a positive coorelation between deletion and divergence

ggsave("Output/fig_TE_fam_div.png", width = 10, height = 6)

#mk# adding s9 to the comparison
df_all3 |>
  ggplot(aes(x = perc_div, y = qry_width, colour = TE_type)) +
  scale_y_log10() +
  geom_point(size = 0.8, alpha = 0.3) +
  geom_smooth(method = lm) +
  facet_grid(TE_type ~ Species, scales = "free_y") +
  labs(x = "Divergence from the reference sequence (%)",
       y = "TE copy length (bp, log scale)",
       title = "Older TE copies appear longer",
       caption = "Each point is one TE copy; lines are linear fits with 95% confidence bands.\nn: T. nativa DNA = 85, Retro = 738; T. spiralis DNA = 100, Retro = 738; s9 DNA = 27, Retro = 234") +
  theme_light() +
  theme(legend.position = "none")

ggsave("Output/fig_div_TEwidth_s9.png", width = 10, height = 6)
#==> s9 seems to follow the same trend

#mk# now lookng at the TE families wih s9 included
top_families2 <- df_all3 |>
  count(matching_class, sort = TRUE) |>
  slice_head(n = 9) |>
  pull(matching_class)

df_all3 |>
  filter(matching_class %in% top_families2) |>
  ggplot(aes(x = perc_div, y = qry_width, colour = Species)) +
  scale_y_log10() +
  geom_point(size = 0.8, alpha = 0.4) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ matching_class, scales = "free") +
  labs(x = "Divergence from consensus (%)",
       y = "TE copy length (bp, log scale)",
       title = "Length vs. Divergence Within Individual TE Families") +
  theme_light()
#==> the relationship is only really different in LINE/CR1, s9 has a negative realtionship

ggsave("Output/fig_TE_fam_div_s9.png", width = 10, height = 6)

# Part 4: Do contigs in T.nativa and T.spiralis differ in their composition of TE families?
# PCA on TE family copy counts per contig (one row per contig, one column per family),
# coloured by species, to see whether contigs from the two species carry similar or different TE family mixes.

# Creating a wide table to run PCA

df_pca_wide = df_both_species |>
  group_by(qry_id, Species, matching_class) |>
  summarise(count = n()) |>
  pivot_wider(names_from = matching_class,
              values_from = count,
              values_fill = 0)

dim(df_pca_wide)

help("prcomp")

pca_te = df_pca_wide |>
  ungroup() |>
  select(-qry_id, -Species) |>
  prcomp(scale. = TRUE)

summary(pca_te)

# # Table of PC scores for each contig (from pca_te$x), with species added, for plotting PC1 vs PC2.

help("as_tibble")

df_pca_score = pca_te$x |>
  as_tibble() |>
  mutate(Species = df_pca_wide$Species)

df_pca_score

#Plotting
df_pca_score |>
  ggplot(aes(x = PC1, y = PC2, colour = Species)) +
  geom_point(alpha = 0.5) +
  stat_ellipse(level = 0.95, segments = 51, type = "t") +
  labs(x = "PC1 = 8.6% of variance", y = "PC2 = 7.1% of variance", title = "Contigs from both species have similar TE composition") +
  labs( caption = "Each point is one contig; ellipses cover ~95% of each species contigs.\nPCA on scaled TE family counts (23 families).") +
  theme_light()

pca_te$rotation[ ,1:2]

ggsave("Output/PCA_species.png", width = 8, height = 6)

#mk# doing pca for s9
df_pc_wide2 = df_all3 |>
  group_by(qry_id, Species, matching_class) |>
  summarise(count = n()) |>
  pivot_wider(names_from = matching_class,
              values_from = count,
              values_fill = 0)

pc_te = df_pc_wide2 |>
  ungroup() |>
  select(-qry_id, -Species) |>
  prcomp(scale. = TRUE)

df_pc_score2 = pc_te$x |>
  as_tibble() |>
  mutate(Species = df_pc_wide2$Species)

df_pc_score2 |>
  ggplot(aes(x = PC1, y = PC2, colour = Species)) +
  geom_point(alpha = 0.5) +
  stat_ellipse(level = 0.95, segments = 51, type = "t") +
  labs(x = "PC1", y = "PC2",) +
  labs( caption = "Each point is one contig; ellipses cover ~95% of each species contigs.\nPCA on scaled TE family counts (23 families).") +
  theme_light()

pc_te$rotation[ ,1:2]

ggsave("Output/PCA_species_s9.png", width = 8, height = 6)
