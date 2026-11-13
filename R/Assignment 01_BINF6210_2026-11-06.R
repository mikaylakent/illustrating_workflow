##***************************
## Assignment 1
##
## Group07_Merged_Final_File
##
## 2026-11-06
##
##***************************

## _ Packages used -------
library(tidyverse)
conflicted::conflicts_prefer(dplyr::filter())
library(viridis)
# + scale_color/fill_viridis_c/d()
theme_set(theme_light())

# Startup ends here

## _ Comment codes ------
# Coding explanations (#, often after the code, but not exclusively)
# Code organization (## XXXXX -----)
# Justification for a section of code ## XXX
# Dead end analyses because it did not work, or not pursuing this line of inquiry (but leave it in as a trace of it, to potentially solve this issue, or avoid making the same mistake in the future # (>_<) 
# Solutions/results/interpretations (#==> XXX)
# Reference to manuscript pieces, figures, results, tables, ... # (*_*)
# TODO items #TODO
# names for data frames (df_name), lists (ls_name), ... (Thanks Jacqueline May)

## code for assignment part 1 ----

df_trich_sp <- read.csv("../data/trichinella_spiralis_te.csv")
df_trich_sp

plot(sort(table(df_trich_sp$matching_class), decreasing = T)[1:10])

df_trich_s9 <- read.csv("../data/GCA_001447505.1.csv")
df_trich_s9

plot(sort(table(df_trich_s9$matching_class), decreasing = T)[1:10])

# filtering for only TEs
df_trich_sp$matching_class |> unique()

non_te_sp = c("simple_repeat", "Simple_repeat",
              "Low_complexity", "Satellite",
              "rRNA", "snRNA", "ARTEFACT",
              "Unknown")

df_trich_s9$matching_class |> unique()
non_te_s9 = c("Simple_repeat", "Satellite", "Low_complexity", "rRNA", "snRNA",
              "ARTEFACT", "Unknown")

# visualizing the number of TEs per contig
df_trich_sp |> 
  group_by(qry_id) |> 
  summarise(count = n()) |> # this is the number of TEs per contig
  arrange(desc(count)) |> 
  filter(count < 1260) |> # this may not be correct, s9 has max count of 1259
  pull(count) |> 
  hist()

df_trich_s9 |> 
  group_by(qry_id) |> 
  summarise(count = n()) |> # this is the number of TEs per contig
  arrange(desc(count)) |> 
  pull(count) |> 
  hist()

# histogram of the percentage deletion
df_trich_sp |> 
  filter(matching_class %notin% non_te_sp) |>
  ggplot(aes(x = perc_del)) +
  geom_histogram()

df_trich_s9 |> 
  filter(matching_class %notin% non_te_s9) |>
  ggplot(aes(x = perc_del)) +
  geom_histogram()

# getting the TEs
te_sp <- df_trich_sp |> 
  filter(matching_class %notin% non_te_sp)

te_s9 <- df_trich_s9 |>
  filter(matching_class %notin% non_te_s9)

# getting DNA- and retro- transposons
te_sp$matching_class |> unique()

te_s9$matching_class |> unique()

# defining retrotransposons
retro_sp <- c("PLE/Naiad", "PLE", "LINE/CR1", "LTR/Pao", "LTR/Gypsy", "LINE/RTE-BovB",
              "SINE/tRNA-RTE", "LTR/Ngaro", "SINE/tRNA-Deu-L2", "LTR/DIRS")
retro_s9 <- c("LTR/Gypsy", "PLE", "PLE/Naiad", "LINE/CR1", "LTR/Pao", "LTR/Ngaro",
              "SINE/tRNA-Deu-L2", "SINE/tRNA-RTE", "LINE/RTE-BovB", "LINE")

new_sp <- df_trich_sp |> 
  filter(matching_class %notin% non_te_sp) |>
  mutate(class_te = if_else(matching_class %in% retro_sp, "retrotrans", "DNAtrans")) 
new_sp |>
  ggplot(aes( x = perc_del, y = perc_div, colour = class_te)) +
  geom_point() +
  geom_smooth()

new_s9 <- df_trich_s9 |> 
  filter(matching_class %notin% non_te_s9) |>
  mutate(class_te = if_else(matching_class %in% retro_s9, "retrotrans", "DNAtrans")) 
new_s9 |>
  ggplot(aes( x = perc_del, y = perc_div, colour = class_te)) +
  geom_point() +
  geom_smooth()

# do TE families evolve differently
# compare  perc del vs perc div of DNA trans between species
# need to have a column with species
df_all_te <- bind_rows(
  new_sp |>
    select(perc_div, perc_del, matching_class, class_te) |>
    mutate(species = "spiralis"),
  new_s9 |>
    select(perc_div, perc_del, matching_class, class_te) |>
    mutate(species = "T9"))

names(df_all_te)

ggplot(df_all_te, aes(x = perc_del, y = perc_div, colour = species)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~class_te)

# trying to test it
dna_transposons <- df_all_te |>
  filter(class_te == "DNAtrans")

model <- lm(perc_div ~ 1 + perc_del + species + perc_del*species, data = dna_transposons)
summary(model)

# the species does have a significant effect on the divergence
model2 <- lm(perc_div ~ 1 + species, data = dna_transposons)
summary(model2)

# for retrotransposons
retro_transposons <- df_all_te |>
  filter(class_te == "retrotrans")

model3 <- lm(perc_div ~ 1 + perc_del + species + perc_del*species, data = retro_transposons)
summary(model3)
# the species does not have a significant effect on the divergence

model4 <- lm(perc_div ~ 1 + perc_del, data = retro_transposons)
summary(model4)
# the percentage deletion has a significant effect

df_all_te |>
  ggplot(aes(x = perc_div, fill = species)) +
  geom_density(bins = 30, alpha = 0.6, position = "identity") +
  scale_fill_manual(values = c("spiralis" = "steelblue", "T9" = "tomato")) +
  labs(x = "Percent Divergence", y = "Count", fill = "Species") 

common_te <- c("DNA/MULE-MuDR", "DNA/PIF-Harbinger", "DNA/TcMar-Tc1", "LINE/CR1", "LTR/Gypsy", "LTR/Ngaro", "LTR/Pao", "PLE", "PLE/Naiad")

df_common_te <- df_all_te |>
  filter(matching_class %in% common_te)
df_common_te

df_common_te |>
  ggplot(aes(x = matching_class, y = perc_div, fill = species)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "TE Class", y = "Pecrent Divergence", fill = "Species") 


## code for assignment part 2 ----

#TE Composition

# Loading datasets

df_trichinella_sp <- read.csv("../data/trichinella_spiralis_te.csv")

df_trichinella_t9 <- read.csv("../data/GCA_001447505.1.csv")


# Getting to know about the datasets

dim(df_trichinella_sp)
dim(df_trichinella_t9)

names(df_trichinella_sp)
names(df_trichinella_t9)

head(df_trichinella_sp)
head(df_trichinella_t9)

unique(df_trichinella_sp$matching_class)
unique(df_trichinella_t9$matching_class)

# Filtering for TE annotations

non_te <- c("Simple_repeat",
            "simple_repeat",
            "Low_complexity",
            "Satellite",
            "Unknown",
            "rRNA",
            "snRNA",
            "ARTEFACT")


df_sp_te <- df_trichinella_sp |>
  filter(!(matching_class %in% non_te))


df_t9_te <- df_trichinella_t9 |>
  filter(!(matching_class %in% non_te))

# Adding dataset names

df_sp_te <- df_sp_te |>
  mutate(dataset = "T. spiralis")


df_t9_te <- df_t9_te |>
  mutate(dataset = "T9")


# Combining datasets

df_all_tes <- bind_rows(df_sp_te, df_t9_te)


# Checking TE annotation counts

table(df_all_tes$dataset)

# Counting TE families in each dataset

te_composition <- df_all_tes |>
  count(dataset, matching_class)

te_composition

# Calculating TE family proportions

te_composition <- te_composition |>
  group_by(dataset) |>
  mutate(proportion = n / sum(n),
         percentage = proportion * 100) |>
  ungroup()

te_composition

# Visualizing TE family composition

te_composition |>
  ggplot(aes(x = matching_class,
             y = proportion,
             colour = dataset)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "TE Family/Class",
       y = "Proportion of TE Annotations",
       colour = "Dataset")

# Ranking TE families by abundance

te_composition |>
  arrange(dataset, n)
# Gives us from the lowest abundance to the highest

# Top five TE families in T. spiralis

te_composition |>
  filter(dataset == "T. spiralis") |>
  arrange(n) |>
  tail(5) |>
  summarise(top_five_percentage = sum(percentage))
# Gives us the combined percentage of top five, %88.5

# Top five TE families in T9

te_composition |>
  filter(dataset == "T9") |>
  arrange(n) |>
  tail(5) |>
  summarise(top_five_percentage = sum(percentage))
# Gives us the combined percentage of top five, %91.7