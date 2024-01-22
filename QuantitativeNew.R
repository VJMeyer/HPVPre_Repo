#Worldwide HPV Prevalence Among Women 50 Years and Older with Predominantly Normal Cytology: A Meta-Analysis
#Master Thesis - Mina Deol

# code QUANTITATIVE ANALYSIS ---------------------------------------------------

#packages required for 
install.packages("tidyverse")
install.packages("metafor")
install.packages("meta")
install.packages("devtools")
devtools::install_github("MathiasHarrer/dmetar") #installation manager might ask for direction select 1 to update all, and then n to not install from sources the compilation
install.packages("tidyr")


#load packages
library(tidyverse)
library(metafor)
library(meta)
library(dmetar)
library (dplyr)
library(tidyr)
library(stringr)
library(readr)


setwd("~/Library/CloudStorage/GoogleDrive-minadeol@gmail.com (2024-01-14 14:41)/My Drive/Thesis/R Code/Version1")
any_type <- read_csv("Anytype_OLDI.csv")  
hr_type <- read_csv("HR_OLDI.csv") 
gtype <- read_csv("gtype.csv") 


# DATA CLEAN UP----  
# Ensure there are no MISSING values in the important columns of ANYTYPE
class(any_type)
head(any_type)
str(any_type)
summary(any_type)
view(any_type)

# Filter rows where num_older_wom is NA and select study_id
na_in_num_older_wom <- any_type %>%
  filter(is.na(num_older_wom)) %>%
  select(study_id)

# Filter rows where num_hpv_pos is NA and select study_id
na_in_num_hpv_pos <-any_type %>%
  filter(is.na(num_hpv_pos)) %>%
  select(study_id)

# Filter rows where num_hpv_neg is NA and select study_id
na_in_num_hpv_neg <- any_type %>%
  filter(is.na(num_hpv_neg)) %>%
  select(study_id)

# View the results
na_in_num_older_wom
na_in_num_hpv_pos
na_in_num_hpv_neg

# there was missing data in one row (63O) of num_hpv_neg which was corrected 

# Filter out rows where the identifier ends with 'P'
any_type_filtered <- any_type %>%
  filter(!str_detect(study_id, "P$"))

summary(any_type_filtered)

nrow(any_type_filtered)

# Extract and view unique study ids in hr_type_filtered to make sure it is correct
unique_study_ids_any <- any_type_filtered %>%
  select(study_id) %>%
  distinct()

print(unique_study_ids_any, n = 83)

#Ensure there are no MISSING values in the columns  of  HR TYPE  where there should be data
class(hr_type)
head(hr_type)
str(hr_type)
summary(hr_type)
view(hr_type)
nrow(hr_type)

# Filter rows where num_older_wom is NA and select study_id
na_in_num_older_wom_hr <- hr_type %>%
  filter(is.na(num_older_wom)) %>%
  select(study_id)

# Filter rows where num_hr_hpv_pos is NA and select study_id
na_in_num_hr_hpv_pos <- hr_type %>%
  filter(is.na(num_hr_hpv_pos)) %>%
  select(study_id)

# Filter rows where hr_prev is NA and select study_id
na_in_hr_prev <- hr_type %>%
  filter(is.na(hr_prev)) %>%
  select(study_id)

# View the results
na_in_num_older_wom_hr
na_in_num_hr_hpv_pos
na_in_hr_prev

# we must exclude study id 23HG from the HR analysis because it has crude numbers for all gtypes but not for HR type
# Filter out rows where the identifier ends with 'P' or is '23HG'
hr_type_filtered <- hr_type %>%
  filter(!str_detect(study_id, "P$"), study_id != "23HG")

summary(hr_type_filtered)
nrow(hr_type_filtered)


# Extract and view unique study ids in hr_type_filtered to make sure it is correct
unique_study_ids_hr <- hr_type_filtered %>%
  select(study_id) %>%
  distinct()

print(unique_study_ids_hr, n = 112)

# Filter out rows where the identifier ends with 'P' for prevalance only for gtype_filtered 

#trim whitespaces
gtype <- gtype %>%
  mutate(study_id = trimws(study_id))

gtype_filtered <- gtype %>%
  filter(!str_detect(study_id, "P$"))

class(gtype_filtered)
head(gtype_filtered)
str(gtype_filtered)
summary(gtype_filtered)
view(gtype_filtered)

summary(gtype_filtered)
nrow(gtype_filtered)

# Extract and view unique study ids in gtype_filtered to make sure it is correct
unique_study_ids_gtype <- gtype_filtered %>%
  select(study_id) %>%
  distinct()

nrow(unique_study_ids_gtype)
print(unique_study_ids_gtype, n = 43)



# Data frames for quantitative analysis ----
# any hpv type analysis ---> any_type_filtered
#hr type analysis ---> hr_type_filtered
#type-specific analysis --->  gtype_filtered



# Meta-analysis model for any_type_filtered ----

summary(any_type_filtered) # there are zero event values in the num_hpv_pos

# Identifying zero-event cases
zero_event_cases <- any_type_filtered[any_type_filtered$num_hpv_pos == 0, ]

# Counting the number of zero-event cases
num_zero_event_cases <- nrow(zero_event_cases)

# Printing the number of zero-event cases
print(paste("Number of zero-event cases:", num_zero_event_cases))

# Printing the study IDs of zero-event cases
if(num_zero_event_cases > 0) {
  print("Study IDs with zero-event cases:")
  print(zero_event_cases$author_year)
  print(zero_event_cases$study_id)
} else {
  print("There are no zero-event cases.")
}

# there are three studies with zero event cases 
# [1] "Becker et al., 1991" "Eren et al., 2010"   "Rahmat et al., 2021"
# [1] "23O" "44O" "71O"

# Option 1: random effects MA, Freeman-Tukey double arcsine transformation, restricted maximum likelihood estimator, Knapp-Hartung adjustment
any_type_ma <- metaprop(event = num_hpv_pos, 
                        n = num_older_wom,
                        studlab = author_year,
                        data = any_type_filtered, 
                        sm = "PFT", 
                        method.tau = "REML",
                        method.ci = "NAsm",  # Specify the method for confidence intervals
                        add = 0,
                        fixed = FALSE,
                        random = TRUE,
                        hakn = TRUE,
                        title = "Anytype HPV Prevalence in Women (50+) with Predominantly Normal Cytology")

summary(any_type_ma)

# # Option 2: rma with the double arcsine transformation (and then converted back to the proportions for the output)
# #Calculate the effect sizes and their variances using the Freeman-Tukey double arcsine transformation.
# any_type_filtered$yi <- escalc(measure = "PFT", xi = num_hpv_pos, ni = num_older_wom, add = 0, data = any_type_filtered)$yi
# any_type_filtered$vi <- escalc(measure = "PFT", xi = num_hpv_pos, ni = num_older_wom, add = 0, data = any_type_filtered)$vi
# 
# #Perform the meta-analysis using REML for the estimation of between-study variance and applying the Knapp-Hartung adjustment for the confidence intervals
# any_type_ma3 <- rma(yi, vi, data = any_type_filtered, method = "REML", test = "knha")
# 
# #Display the summary of the meta-analysis, including the back-transformed pooled effect size and its confidence intervals.
# summary(any_type_ma3)
# 
# # Back-transformation to original scale
# any_type_ma3$pred <- predict(any_type_ma3, transf = transf.ipft.hm, targ = list(ni = any_type_filtered$num_older_wom))
# print(any_type_ma3$pred)

# ANYTYPE FOREST PLOT ----

# Create a forest plot based on the meta-analysis results

# Set up a larger plot window size in R
pdf("forestplot_anytype.pdf", width=10, height=20)  # Width and height in inches, adjust as needed

forest(any_type_ma, 
       common = FALSE, 
       print.tau2 = TRUE, 
       print.Q = TRUE, 
       print.pval.Q = TRUE, 
       print.I2 = TRUE, 
       rightcols = FALSE, 
       pooled.totals = FALSE, 
       weight.study = "random", 
       leftcols = c("studlab", "num_hpv_pos", "num_older_wom", "effect", "ci"), 
       leftlabs = c("Study (Year)", "HPV+ Cases", "Total", "Prevalence", "95% C.I."), 
       xlab = "Pooled Prevalence Rate", 
       smlab = "", 
       xlim = c(0,1), 
       pscale = 1, 
       squaresize = 0.5, 
       fs.hetstat = 10, 
       digits = 2, 
       col.square = "navy", 
       col.square.lines = "navy", 
       col.diamond = "maroon", 
       col.diamond.lines = "maroon")

# Close the PDF device
dev.off()


# Note: The argument 'comb.random' is not valid for the 'forest()' function.
# If you need to specify additional arguments for visual aspects like square sizes or colors,
# you will have to use the arguments specific to the 'forest()' function in the 'meta' package.

# Remember to replace 'any_type_ma$Author' and 'any_type_ma$Year' with the actual names of the columns in your data.


# Additional customization may be required to perfectly match the provided plot

# forest(any_type_ma,
#        comb.fixed = FALSE, # Do not combine studies using fixed effect model
#        comb.random = TRUE, # Combine studies using random effects model
#        xlim = c(-4, 4), # Set the limits of the x-axis
#        xlab = "Effect Sizes", # Label for the x-axis
#        alim = c(-4, 4), # Set the axis limits for the plot
#        at = seq(-4, 4, 0.5), # Tick marks on the axis
#        cex = 0.75, # Size of the text in the plot
#        col = "blue", # Color for the points and lines
#        pch = 19, # Type of point used to display the effect sizes
#        psize = 1, # Size of points to represent studies
#        lwd = 1, # Line width for the confidence interval lines
#        refline = 0, # Reference line, typically at 0 for no effect
#        digits = 2, # Digits for rounding the effect sizes and CIs
#        mlab = "RE Model", # Label for the summary effect in the model
#        backtransf = TRUE, # Back-transform the effect sizes and CIs
#        sortvar = any_type_ma1$TE, # Sort studies by effect size
#        smlab = "Prevalence" # Label for the size of the studies
# )


# ANYTYPE Sensitvity analysis ----

## create a vector of study identifiers that had zero events
studies_with_zero_events <- c("Becker et al., 1991", "Eren et al., 2010", "Rahmat et al., 2021")

# Filter out the studies with zero events from the data frame
any_type_filtered_sensitivity <- any_type_filtered[!any_type_filtered$author_year %in% studies_with_zero_events, ]

any_type_ma_sens <- metaprop(event = num_hpv_pos, 
                             n = num_older_wom,
                             studlab = author_year,
                             data = any_type_filtered_sensitivity, 
                             sm = "PFT", 
                             method.tau = "REML",
                             add = 0,
                             fixed = FALSE,
                             random = TRUE,
                             hakn = TRUE,
                             title = "Anytype HPV Prevalence in Women (50+) with Predominantly Normal Cytology")
summary(any_type_ma_sens)


# Meta-analysis model for hr_type_filtered ------------------------------------------------

summary(hr_type_filtered) # there are no zero event values in num_hr_hpv_pos 

# Option 1: random effects MA, Freeman-Tukey double arcsine transformation, restricted maximum likelihood estimator,, Freeman-Tukey double arcsine transformation 

hr_type_filtered$num_hr_hpv_pos <- as.integer(hr_type_filtered$num_hr_hpv_pos)

hr_type_ma <- metaprop(event = num_hr_hpv_pos, 
                       n = num_older_wom,
                       studlab = author_year,
                       data = hr_type_filtered, 
                       sm = "PFT", 
                       method.tau = "REML",
                       method.ci = "NAsm",  # Specify the method for confidence intervals
                       add = 0,
                       fixed = FALSE,
                       random = TRUE,
                       hakn = TRUE,
                       title = "HR-HPV Prevalence in Women (50+) with Predominantly Normal Cytology")

summary(hr_type_ma)


# # Option 2: rma with the double arcsine transformation (and then converted back to the proportions for the output)
# #Calculate the effect sizes and their variances using the Freeman-Tukey double arcsine transformation.
# hr_type_filtered$yi <- escalc(measure = "PFT", xi = num_hr_hpv_pos, ni = num_older_wom, add = 0, data = hr_type_filtered)$yi
# hr_type_filtered$vi <- escalc(measure = "PFT", xi = num_hr_hpv_pos, ni = num_older_wom, add = 0, data = hr_type_filtered)$vi
# 
# #Perform the meta-analysis using REML for the estimation of between-study variance and applying the Knapp-Hartung adjustment for the confidence intervals
# hr_type_ma2 <- rma(yi, vi, data = hr_type_filtered, method = "REML", test = "knha")
# 
# #Display the summary of the meta-analysis, including the back-transformed pooled effect size and its confidence intervals.
# summary(hr_type_ma2)
# 
# # Back-transformation to original scale
# hr_type_ma2$pred <- predict(hr_type_ma2, transf = transf.ipft.hm, targ = list(ni = hr_type_filtered$num_older_wom))
# print(hr_type_ma2$pred)

# HR-TYPE FOREST PLOT ----

# Create a forest plot based on the meta-analysis results

# Set up a larger plot window size in R
pdf("forestplot_hr_type.pdf", width=10, height=30)  # Width and height in inches, adjust as needed

forest(hr_type_ma, 
       common = FALSE, 
       print.tau2 = TRUE, 
       print.Q = TRUE, 
       print.pval.Q = TRUE, 
       print.I2 = TRUE, 
       rightcols = FALSE, 
       pooled.totals = FALSE, 
       weight.study = "random", 
       leftcols = c("studlab", "num_hr_hpv_pos", "num_older_wom", "effect", "ci"), 
       leftlabs = c("Study (Year)", "HR HPV+ Cases", "Total", "Prevalence", "95% C.I."), 
       xlab = "Pooled Prevalence Rate", 
       smlab = "", 
       xlim = c(0,1), 
       pscale = 1, 
       squaresize = 0.5, 
       fs.hetstat = 10, 
       digits = 2, 
       col.square = "navy", 
       col.square.lines = "navy", 
       col.diamond = "maroon", 
       col.diamond.lines = "maroon")

# Close the PDF device
dev.off()
# Meta-analysis model for gtype_filtered -----------------------------------------------

summary(gtype_filtered) # there are zero event values in num_type

# Identify unique study_id's in the gtype spreadsheet
unique_values_studyid <- unique(gtype_filtered$study_id)
# Print the unique values
print(unique_values_studyid) #43 studies 

# check the number of genotypes and the instances of them in the data frame "gtype_filtered" column "hpv_type" 
# Identify unique values
unique_values_type <- unique(gtype_filtered$hpv_type)
# Print the unique values
print(unique_values_type)


library(dplyr)
library(tidyr)
library(ggplot2)

# If "hpv_type" contains concatenated genotypes, separate them into individual rows.
gtype_figure <- gtype_filtered %>%
  group_by(hpv_type) %>%
  summarise(total_cases = sum(num_type, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(total_cases))

# Get the top 5 HPV types
top_hpv_types <- head(gtype_figure, 5)

top_hpv_types

# Create a bar plot for visual representation
ggplot(top_hpv_types, aes(x = hpv_type, y = total_cases, fill = hpv_type)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "HPV Type", y = "Total Number of Cases", 
       title = "Top 5 HPV Types in Women Over 50",
       fill = "HPV Type")

# Alternatively, create a pie chart for visual representation
ggplot(top_hpv_types, aes(x = "", y = total_cases, fill = hpv_type)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() +
  labs(fill = "HPV Type", 
       title = "Top 5 HPV Types in Women Over 50, N = 43")


# Identifying zero-event cases
zero_event_cases_gtype <- gtype_filtered[gtype_filtered$num_type == 0, ]

# Counting the number of zero-event cases
num_zero_event_cases_gtype <- nrow(zero_event_cases_gtype)

# Printing the number of zero-event cases
print(paste("Number of zero-event cases gtype:", num_zero_event_cases_gtype))

# Printing the study IDs of zero-event cases
if(num_zero_event_cases_gtype > 0) {
  print("Study IDs with zero-event cases gtype:")
  print(zero_event_cases_gtype$author_year)
  print(zero_event_cases_gtype$study_id)
} else {
  print("There are no zero-event cases in gtype.")
}

# [1] "Study IDs with zero-event cases gtype:"
# [1] "Cathro et al., 2009" "Debrah et al., 2021" "Klug et al., 2007"   "Kobetz et al., 2012"
# [1] "30OG"  "39OG"  "57OHG" "58OHG"

# need to make sure that all values for the analysis are in integer format






g_type_ma <- metaprop(event = num_hr_hpv_pos, 
                      n = num_older_wom,
                      studlab = author_year,
                      data = hr_type_filtered, 
                      s