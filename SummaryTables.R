#Worldwide HPV Prevalence Among Women 50 Years and Older with Predominantly Normal Cytology: A Meta-Analysis
#Master Thesis - Mina Deol

# code for TABLE 1 summary ---------------------------------------------------

#packages required for 
#install.packages("tidyverse")
#install.packages("metafor")
#install.packages("meta")
#install.packages("devtools")
#devtools::install_github("MathiasHarrer/dmetar")
#install.packages("tidyr")

#installation manager might ask for direction select 1 to update all, and then n to not install from sources the compilation
#load packages
library(tidyverse)
library(metafor)
library(meta)
library (dplyr)
library(tidyr)
library(stringr)

#confirm working directory
getwd()
list.files() #list files within it

#read files
file_list <- list.files(pattern = "*.csv")
data_list <- lapply(file_list, read.csv)

data_list[[1]]  # Access the first data frame etc (1= anytype, 2= HR)
#name the elements of data_list based on the file names
names(data_list) <- file_list
# check names
names(data_list) #[1] "Anytype_OLDI.csv"  "HR_OLDI.csv"       "HR_years.csv"      "overall_years.csv" "Type_OLDI.csv"  

#check the class of the data frame 
class(data_list[[1]])

#remove columns not need for analysis
data_list[[1]] <- data_list[[1]] %>% select(-Notes, -test_details)
data_list[[2]] <- data_list[[2]] %>% select(-Comments, -test_details)

#select the columns that you want to keep for data_list[[3]] 
# data_list[[3]] <- data_list[[3]] %>% select(study_id, author_year, title, study_design, ending_year, 
#                                             world_region, world_subregion, country, city_state, focus, 
#                                             overall_cytology, pap_method, recruitment_setting, 
#                                             hpv_types_reported, hpv_test, num_older_wom, num_hpv_pos, 
#                                             num_hr_hpv_pos, hpv_type, num_type, type_prev)

#check column names
colnames(data_list[[1]])
colnames(data_list[[2]])
colnames(data_list[[5]]) #for the "Type_OLDI.csv" 

#because there will be new NA values introduced when the data frames are merged,  
#we can replace the orginal NA values in the data frames  with something that can be identified 

# Replace NA in all columns of the data lists but without losing the class of data frame
# data_list[[1]] <- data_list[[1]] %>%
#   mutate_all(~ifelse(is.na(.), "missing", .))
# data_list[[2]] <- data_list[[2]] %>%
#   mutate_all(~ifelse(is.na(.), "missing", .))
# data_list[[3]] <- data_list[[3]] %>%
#   mutate_all(~ifelse(is.na(.), "missing", .))

#check the class of the data frame 
class(data_list[[1]])
class(data_list[[2]])
# class(data_list[[3]])

# Filter out rows from data_list[[2]] that have study_id's already in data_list[[1]]
data_list[[6]] <- anti_join(data_list[[2]], data_list[[1]], by = "study_id")

# Now perform the full_join with the filtered data_list[[2]]
summary_table <- full_join(data_list[[1]], data_list[[6]], by = "study_id")

class(summary_table)
colnames(summary_table)
nrow(summary_table)

# Remove leading and trailing spaces from study_id
summary_table <- summary_table %>%
  mutate(study_id = trimws(study_id))

view(summary_table)
#------------------risk of bias numbers-------------------------------
# Count occurrences in "risk_of_bias.x"
count_risk_x <- summary_table %>%
  group_by(risk_of_bias.x) %>%
  summarise(count = n())

# Count occurrences in "risk_of_bias.y"
count_risk_y <- summary_table %>%
  group_by(risk_of_bias.y) %>%
  summarise(count = n())

# View the results
count_risk_x
count_risk_y

# Find rows with "" in risk_of_bias.x
missing_risk_x <- summary_table %>%
  filter(risk_of_bias.x == "")

# View the study IDs with missing risk_of_bias.x
missing_study_ids <- missing_risk_x$study_id

# Print the study IDs
missing_study_ids
#54OH was missing the ROB which was high when I checked the excel file and corrected it 

#----------------sample size range and categories-------------------------------------
#use summary_table which has the unique study_id for 159 included studes
#however I need to exclude the study id's which end in P because that is my indicator that only percent prevalence provided and no N

# Categorize and count for num_older_wom.x
categorized_counts_x <- summary_table %>%
  filter(!str_detect(study_id, "P")) %>%
  mutate(
    category_x = case_when(
      num_older_wom.x < 100 ~ "<100",
      num_older_wom.x >= 100 & num_older_wom.x <= 300 ~ "100-300",
      num_older_wom.x >= 301 & num_older_wom.x <= 500 ~ "301-500",
      num_older_wom.x >= 501 & num_older_wom.x <= 999 ~ "501-999",
      num_older_wom.x >= 1000 ~ ">1000",
      TRUE ~ "Other"
    )
  ) %>%
  count(category_x)

# Categorize and count for num_older_wom.y
categorized_counts_y <- summary_table %>%
  filter(!str_detect(study_id, "P")) %>%
  mutate(
    category_y = case_when(
      num_older_wom.y < 100 ~ "<100",
      num_older_wom.y >= 100 & num_older_wom.y <= 300 ~ "100-300",
      num_older_wom.y >= 301 & num_older_wom.y <= 500 ~ "301-500",
      num_older_wom.y >= 501 & num_older_wom.y <= 999 ~ "501-999",
      num_older_wom.y >= 1000 ~ ">1000",
      TRUE ~ "Other"
    )
  ) %>%
  count(category_y)

# View the categorized counts
categorized_counts_x
categorized_counts_y

#--------------------------# TYPE specific quantitative analysis--------------------------------------------------------------------
# Count unique study_ids containing 'G' but not containing 'P' for summary_table 
count_study_id_G_no_P_summary<- summary_table %>%
  filter(str_detect(study_id, "G") & !str_detect(study_id, "P")) %>%  # Filter study_id containing 'G' and not 'P'
  distinct(study_id) %>%                                             # Get unique study_ids
  nrow()                                                             # Count the number of unique study_ids

# Print the count which is 45
count_study_id_G_no_P_summary

#-------------list of study ids which are included in the TYPE specific quantitative analysis
study_ids_with_G_no_P_any <- summary_table %>%
  filter(str_detect(study_id, "G") & !str_detect(study_id, "P")) %>%  # Filter study_id containing 'G' and not 'P'
  select(study_id)                                                   # Select the study_id column

# Print the study_ids
study_ids_with_G_no_P


#--------------------------# anytype quantitative analysis--------------------------------------------------------------------
# Count unique study_ids containing 'G' but not containing 'P' for summary_table 
count_study_id_O_no_P_summary<- summary_table %>%
  filter(str_detect(study_id, "O") & !str_detect(study_id, "P")) %>%  # Filter study_id containing 'G' and not 'P'
  distinct(study_id) %>%                                             # Get unique study_ids
  nrow()                                                             # Count the number of unique study_ids

# Print the count which is 83
count_study_id_O_no_P_summary

#-------------list of study ids which are included in the anytype quantitative analysis
study_ids_with_O_no_P_any <- summary_table %>%
  filter(str_detect(study_id, "O") & !str_detect(study_id, "P")) %>%  # Filter study_id containing 'G' and not 'P'
  select(study_id)                                                   # Select the study_id column

# Print the study_ids
study_ids_with_O_no_P

#--------------------------# HR quantitative analysis--------------------------------------------------------------------
# Count unique study_ids containing 'G' but not containing 'P' for summary_table 
count_study_id_H_no_P_summary<- summary_table %>%
  filter(str_detect(study_id, "H") & !str_detect(study_id, "P")) %>%  # Filter study_id containing 'G' and not 'P'
  distinct(study_id) %>%                                             # Get unique study_ids
  nrow()                                                             # Count the number of unique study_ids

# Print the count which is 
count_study_id_H_no_P_summary

#-------------list of study ids which are included in the HR quantitative analysis
study_ids_with_O_no_P_any <- summary_table %>%
  filter(str_detect(study_id, "H") & !str_detect(study_id, "P")) %>%  # Filter study_id containing 'G' and not 'P'
  select(study_id)                                                   # Select the study_id column

# Print the study_ids
study_ids_with_H_no_P


#Qualitative analysis studies have a study id which ends in P-------------------------------

# # Remove leading and trailing spaces and convert back to numeric
# data_list[[1]] <- data_list[[1]] %>%
#   mutate(num_older_wom = as.numeric(trimws(as.character(num_older_wom))))
# 
# # Check the structure to confirm the changes
# str(data_list[[1]])
# 
# # Remove leading and trailing spaces and convert back to numeric
# data_list[[2]] <- data_list[[2]] %>%
#   mutate(num_older_wom = as.numeric(trimws(as.character(num_older_wom))))
# 
# # Check the structure to confirm the changes
# str(data_list[[2]])

# Count study_ids containing 'P' in summary_table
count_study_id_P <- summary_table %>%
  filter(str_detect(study_id, "P")) %>%  # Filter study_id containing 'P'
  nrow()                                # Count the number of study_ids

# Select study_ids containing 'P'
study_ids_with_P <- summary_table %>%
  filter(str_detect(study_id, "P")) %>%  # Filter study_id containing 'P'
  select(study_id)                      # Select the study_id column

# Print the list of study_ids which P which are for the  qualitative analysis ie prevalence range
study_ids_with_P

#17 studies in total 
# 4 studies in any type
# 15 studies in HR type
# 8 studies in type specifics
# combine

#__________________Descriptive Characteristics of TABLE 1____________________________####


cols_to_check <- c("study_design.x", "world_region.x" , "world_subregion.x", "focus.x" , "pap_method.x",
                   "recruitment_setting.x", "hpv_test.x", "study_design.y", "world_region.y",
                   "world_subregion.y", "focus.y" , "pap_method.y", "recruitment_setting.y", "hpv_test.y")

count_unique_values <- function(data, column_name) {data %>%
    count(!!sym(column_name)) %>%
    arrange(desc(n))}

results <- lapply(cols_to_check, function(col) count_unique_values(summary_table, col))
names(results) <- cols_to_check

results

# data clean up for the summary table numbers

summary_table <- summary_table %>%
  mutate(
    study_design.x = case_when(
      study_design.x == "case control " ~ "case-control",
      study_design.x == "prospective " ~ "prospective",
      TRUE ~ study_design.x  # This line means "in all other cases, leave the value as is"
    )
  )

summary_table <- summary_table %>%
  mutate(
    study_design.y = case_when(
      study_design.y == "case control " ~ "case-control",
      study_design.y == "prospective " ~ "prospective",
      TRUE ~ study_design.y  # This line means "in all other cases, leave the value as is"
    )
  )

summary_table <- summary_table %>%
  mutate(
    world_region.x = case_when(
      world_region.x == "Asia " ~ "Asia",
      world_region.x == "Americas " ~ "Americas",
      world_region.x == "Africa " ~ "Africa",
      TRUE ~ world_region.x  # This line means "in all other cases, leave the value as is"
    )
  )

summary_table <- summary_table %>%
  mutate(
    world_region.y = case_when(
      world_region.y == "Africe" ~ "Africa",
      world_region.y == "Africa " ~ "Africa",
      world_region.y == "America" ~ "Americas",
      world_region.y == "Americas " ~ "Americas",
      world_region.y == "Asia " ~ "Asia",
      world_region.y == "Europe " ~ "Europe",
      TRUE ~ world_region.y  # This line means "in all other cases, leave the value as is"
    )
  )

summary_table <- summary_table %>%
  mutate(
    focus.x = case_when(
      focus.x == "both " ~ "Both",
      focus.x == "Both " ~ "Both",
      TRUE ~ focus.x  # This line means "in all other cases, leave the value as is"
    )
  )

summary_table <- summary_table %>%
  mutate(
    focus.y = case_when(
      focus.y == "Both " ~ "Both",
      TRUE ~ focus.y  # This line means "in all other cases, leave the value as is"
    )
  )

summary_table <- summary_table %>%
  mutate(
    world_subregion.y = case_when(
      world_subregion.y == "Northern America" ~ "North America",
      world_subregion.y == "Austraiia and New Zealand" ~ "Australia and New Zealand",
      world_subregion.y == "Middle Adrica" ~ "Middle Africa ",
      world_subregion.y == "South-eastern Asia" ~ "Southeastern asia",
      world_subregion.y == "Southeastern asia" ~ "Southeastern Asia",
      world_subregion.y == "Middle Africa " ~ "Middle Africa",
      TRUE ~ world_subregion.y  # This line means "in all other cases, leave the value as is"
    )
  )

summary_table <- summary_table %>%
  mutate(
    world_subregion.x = case_when(
      world_subregion.x == "Southeastern asia" ~ "Southeastern Asia",
      world_subregion.x == "Middle Africa " ~ "Middle Africa",
      TRUE ~ world_subregion.x  # This line means "in all other cases, leave the value as is"
    )
  )

summary_table <- summary_table %>%
  mutate(
    pap_method.x = case_when(
      pap_method.x == "conventional " ~ "conventional",
      TRUE ~ pap_method.x  # This line means "in all other cases, leave the value as is"
    )
  )

summary_table <- summary_table %>%
  mutate(
    pap_method.y = case_when(
      pap_method.y == "conventional " ~ "conventional",
      TRUE ~ pap_method.y  # This line means "in all other cases, leave the value as is"
    )
  )


summary_table <- summary_table %>%
  mutate(
    recruitment_setting.x = case_when(
      recruitment_setting.x == "screening " ~ "screening",
      recruitment_setting.x == "clinical setting " ~ "clinical setting",
      TRUE ~ recruitment_setting.x  # This line means "in all other cases, leave the value as is"
    )
  )

summary_table <- summary_table %>%
  mutate(
    recruitment_setting.y = case_when(
      recruitment_setting.y == "screening " ~ "screening",
      recruitment_setting.y == "Screening" ~ "screening",
      TRUE ~ recruitment_setting.y  # This line means "in all other cases, leave the value as is"
    )
  )

summary_table <- summary_table %>%
  mutate(
    hpv_test.x = case_when(
      hpv_test.x == "PGMY09/11 " ~ "PGMY09/11",
      hpv_test.x == "MY09/MY11 " ~ "MY09/11",
      hpv_test.x == "MY09/11 and GP5+/6+" ~ "MY09/11 and GP5+/GP6+",
      hpv_test.x == "GP5+/6+" ~ "GP5+/GP6+",
      hpv_test.x == "GP5+/6+" ~ "GP5+/GP6+",
      TRUE ~ hpv_test.x  # This line means "in all other cases, leave the value as is"
    )
  )

summary_table <- summary_table %>%
  mutate(
    hpv_test.y = case_when(
      hpv_test.y == "PGMY09/11 " ~ "PGMY09/11",
      hpv_test.y == "Other PCR" ~ "other PCR",
      hpv_test.y == "cobas" ~ "Cobas",
      hpv_test.y == "PCR L1 or E1\n" ~ "PCR L1 or E1",
      hpv_test.y == "other PCR\n" ~ "other PCR",
      hpv_test.y == "multiple\n" ~ "multiple",
      hpv_test.y == "Mulitple\n" ~ "multiple",
      hpv_test.y == "MY09/MY11 " ~ "MY09/11",
      hpv_test.y == "MY09/11 and GP5+/6+" ~ "MY09/11 and GP5+/GP6+",
      hpv_test.y == "GP5+/6+" ~ "GP5+/GP6+",
      hpv_test.y == "GP5+/6+" ~ "GP5+/GP6+",
      TRUE ~ hpv_test.y  # This line means "in all other cases, leave the value as is"
    )
  )

count_unique_values <- function(data, column_name) {data %>%
    count(!!sym(column_name)) %>%
    arrange(desc(n))}

results <- lapply(cols_to_check, function(col) count_unique_values(summary_table, col))
names(results) <- cols_to_check

results

# First, standardize column names and remove NAs
standardized_results <- lapply(results, function(df) {
  colnames(df) <- c("Value", "Count") # Set a standard column name
  filter(df, !is.na(Value)) # Remove NA values
})

# Now, combine all data frames into one
combined_results <- bind_rows(standardized_results)

merged_results <- combined_results %>%
  group_by(Value) %>%             # Group by the unique values
  summarise(TotalCount = sum(Count))  # Sum the counts for each group

# View the merged results
merged_results
view(merged_results)

library(dplyr)

merged_results <- merged_results %>%
  mutate(Value = case_when(
    Value == "America" ~ "Americas",
    Value == "Latin America" ~ "Latin America and the Caribbean",
    Value == "South America" ~ "Latin America and the Caribbean",
    Value == "Caribbean" ~ "Latin America and the Caribbean",
    Value == "North America" ~ "Northern America",
    Value == "Central America" ~ "Latin America and the Caribbean",
    Value == "Southeastern asia" ~ "Southeastern Asia",
    Value == "Middle Africa " ~ "Middle Africa",
    Value == "cross-sectional " ~ "cross-sectional",
    
    # Add more conditions as needed
    TRUE ~ Value  # Keeps the value unchanged if none of the above conditions are met
  ))
view(merged_results)

# Add a new column for grouping
grouped_results <- merged_results %>%
  mutate(Group = case_when(
    Value %in% c("cross-sectional", "retrospective", "prospective", "case-control") ~ "Study Design",
    Value %in% c("clinical setting", "screening", "general population", "controls") ~ "Recruitment Setting",
    Value %in% c("LBC", "conventional", "LBC and conventional", "missing", "VIA") ~ "Pap Method",
    Value %in% c("Both", "Older") ~ "Study Population Focus",
    Value %in% c("other PCR", "GP5+/GP6+", "PGMY09/11", "MY09/11", "HC2", "LA","PCR L1 or E1", "PGMY09/11 and GP5+/GP6+",
                 "MY09/11 and GP5+/GP6+", "PCR E6/E7", "Cobas", "multiple", "Multiple") ~ "HPV Testing Method",
    Value %in% c("Asia", "Americas", "Europe", "Africa", "Oceania") ~ "World Region",
    Value %in% c("Eastern Asia", "Western Asia", "Southeastern Asia", "Southern Asia", "Northern Europe",
                 "Southern Europe", "Western Europe", "Northern America", "Northern Africa", "Western Africa", 
                 "Eastern Africa", "Southern Africa", "Australia and New Zealand","Melanesia", "Central Asia", 
                 "Middle Africa", "Sub-Saharan Africa", "Latin America and the Caribbean") ~ "World Subregion",
    # Add more conditions as needed
    TRUE ~ "Other"  # Default group
  ))

view(grouped_results)




#-----study end date ranges---------------------------------------------------------------------

library(readr)
HR_years <- read_csv("HR_years.csv")

overall_years <- read_csv("overall_years.csv")

str(overall_years)
str(HR_years)

intersect(overall_years$study_id, HR_years$study_id)

library(dplyr)
unique_HR_years <- anti_join(HR_years, overall_years, by = "study_id")
# Perform a full join to combine overall_years with the unique rows from HR_years
combined_years <- full_join(overall_years, unique_HR_years, by = "study_id")

nrow(combined_years)

# View column names of the combined_years data frame
column_names <- names(combined_years)

# Print the column names
column_names

# Combine the ending_year.x and ending_year.y columns into a single vector
combined_ending_year <- c(combined_years$ending_year.x, combined_years$ending_year.y)

# Calculate the minimum and maximum, excluding NA values
min_ending_year <- min(combined_ending_year, na.rm = TRUE)
max_ending_year <- max(combined_ending_year, na.rm = TRUE)

# Print the results
min_ending_year
max_ending_year


# Function to assign year to a range
year_range_label <- function(year) {
  case_when(
    year >= 1989 & year <= 1999 ~ "1989-1999",
    year >= 2000 & year <= 2005 ~ "2000-2005",
    year >= 2006 & year <= 2010 ~ "2006-2010",
    year >= 2011 & year <= 2015 ~ "2011-2015",
    year >= 2016 & year <= 2019 ~ "2016-2019",
    TRUE ~ "Other"  # For years outside the specified range, if any
  )
}

# Apply the function to both year columns
combined_years <- combined_years %>%
  mutate(
    year_range_x = year_range_label(ending_year.x),
    year_range_y = year_range_label(ending_year.y)
  )

# Count the number of entries in each range for both columns
count_by_range_x <- combined_years %>% count(year_range_x)
count_by_range_y <- combined_years %>% count(year_range_y)

# View the results
count_by_range_x
count_by_range_y



