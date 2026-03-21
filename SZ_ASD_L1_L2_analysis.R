# Arthur's thesis analysis code

# To do in R ----

# == Data cleaning ==
# Sort out AoA data
# - including cols that were deleted from Qualtrics:
#   dlq_e_a4_1, dlq_e_a4_2, dlq_e_a4_3, dlq_e_a4_4, same thing with "f" not "e"
# Sort out QC


# Color scheme ----

color_sz  <- "skyblue3"
color_asd <- "mediumpurple1"

color_l1  <- "darkolivegreen3"
color_l2  <- "lightsalmon2"

color_lt   <- "indianred2"
color_ca_c <- "seagreen3"
color_sm   <- "steelblue2"
color_p    <- "plum3"


# Set working directory and load libraries ----

setwd("C:/Users/Arthur/OneDrive/Documents/R/PhD_Thesis/SZ_ASD_L1_L2_analysis/Data")

library(tidyverse)
library(janitor)
library(Hmisc)
library(GGally)     # For correlograms
library(MatchIt)    # For propensity score matching (PSM)
library(plotrix)    # For std.error function
library(gridExtra)  # For a gallery of plots
library(patchwork)  # Newer library for a gallery of plots
library(lme4)       # For LMEs
library(car)        # To test for collinearity
library(emmeans)    # To visualize LMEs
library(ggeffects)  # Also to visualize LMEs
library(ggbreak)    # For a broken axis on a histogram

# library(boot)       # For bootstrap confidence intervals
# library(partR2)     # For partial R-squared
# library(plspm)
# library(quantreg)   # For quantile regression

# This library/function can be used for assert statements
# https://www.rdocumentation.org/packages/testit/versions/0.13/topics/assert

# traceback() to get error call stack


# ====== Create functions ====== ----

# Create functions for creating totals ----

# Creating function to get the sum of a block (returning NA if there are no values)
# my_data is a df, my_row is a number, my_block is the name of the block we're summing (a string)
# Note: This gives the number of correct multiple-choice responses (i.e., the number marked as "1"), not the actual sum
sum_block <- function(my_data, my_row, my_block)
{
  first_col <- block_col_nums[block_col_nums$block == my_block, "first_col"]
  last_col  <- block_col_nums[block_col_nums$block == my_block, "last_col"] 
  col_count = last_col - first_col + 1
  
  result = NA
  if(sum(is.na(my_data[my_row, first_col:last_col])) < col_count)
  {
    result = sum(my_data[my_row, first_col:last_col] == 1, na.rm=TRUE)
  }
  
  out <- result
}

# Creating function to create L1 or L2 column for a task
# The is_for_l1 parameter should be set to 0 for L2 and 1 for L1
convert_totals_to_l1_or_l2 <- function(my_data, is_for_l1, en_column, fr_column)
{
  out <- ifelse(my_data[["is_l1_en"]] == is_for_l1, my_data[[en_column]], my_data[[fr_column]])
}

# Creating function to add a single column for a question (merging EN and FR responses)
# This can only be used after the "is_l1_en" column is created
add_merged_col <- function(my_df, en_col_name, fr_col_name, new_col_name)
{
  my_df$new_col <- NA
  colnames(my_df)[which(colnames(my_df)=="new_col")] <- new_col_name
  for (my_row in 1:nrow(my_df))
  {
    if (my_df[my_row, "is_l1_en"] == 1)
      my_df[my_row, new_col_name] <- my_df[my_row, en_col_name]
    else
      my_df[my_row, new_col_name] <- my_df[my_row, fr_col_name]
  }
  out <- my_df
}


# Wrapper functions to create plots in consistent style ----

# Correlogram
create_correlogram <- function(my_df, my_cols, my_title, my_labels)
{
  my_reduced_df <- my_df[,my_cols]
  ggpairs(my_reduced_df, title=my_title,
          upper = list(continuous = "cor", combo = "box_no_facet", discrete = "facetbar", na = "na"),
          lower = list(continuous = wrap("smooth", color='darkred', alpha=0.05), combo = "facethist", discrete = "facetbar", na = "na"),
          diag = list(continuous = "blankDiag", discrete = "barDiag", na = "naDiag"),
          columnLabels = my_labels, axisLabels = "none") +
          theme_bw() + theme(plot.title=element_text(hjust=0.5))
}

# Histogram
create_histogram <- function(my_df, my_var, my_title,
                             my_xlab, my_ylab, my_color,
                             x_min = NULL, x_max = NULL)
{
ggplot(data=my_df, aes(x={{my_var}}, rm.na = TRUE)) +
       geom_histogram(fill=my_color, color='gray60', binwidth=1, na.rm=TRUE) +
       xlab(my_xlab) + ylab(my_ylab) +
       ggtitle(my_title) +
       theme_bw() + theme(plot.title=element_text(hjust=0.5)) +
       coord_cartesian(xlim = c(x_min, x_max))
}


# ====== Clean data and calculate totals ====== ----

# Loading and cleaning data ----

# Loading Qualtrics data
my_data <- read.csv(file = "QualtricsData_2025_12_02.csv", encoding="UTF-8")
# my_data <- read.csv(file = "QualtricsData_2025_06_25.csv", encoding="UTF-8")

# Cleaning data frame with the janitor package
my_data <- clean_names(my_data)

# # Remove first two rows, which are just info from Qualtrics
my_data <- my_data[-c(1, 2), ]

# Switch duration in seconds to numeric value
my_data$duration_in_seconds <- as.numeric(my_data$duration_in_seconds)

# Standardizing source names
my_data <- within(my_data, source[source == "cogsci_sona"] <- "sona_cogsci")

# Correcting the source for two participants
# I entered the sources manually due to technical problems, but got the sources wrong for these ones and only noticed after they started
my_data <- within(my_data, source[random_id == 33534377 & source == "sona_cogsci"] <- "sona_psych")
my_data <- within(my_data, source[random_id == 98996051 & source == "sona_psych"]  <- "sona_cogsci")


# First round of checks to remove bad data (bots, test runs, etc.) ----
# Move this below the import of Pavlovia data to unify all QC?

# Filter by consent and age
my_data <- my_data[my_data$give_consent == "I consent to participate. / Je donne mon consentement à participer.", ]
my_data <- subset(my_data, my_data$eligibility_age == "Yes / Oui" | my_data$eligibility_age == "")

# Filter by language eligibility
                           # 2023 data - SONA
my_data <- subset(my_data, my_data$eligibility_language == "Yes / Oui" |
                           # 2023 data - Prolific
                           (my_data$eligibility_language == "" & (my_data$source == "prolific_l1e" | my_data$source == "prolific_l1f")) |
                           # Newer data - bilinguals
                           (my_data$eligibility_english == "Yes / Oui" & my_data$eligibility_french == "Yes / Oui") |
                           # Newer data - monolinguals
                           (my_data$eligibility_english == "Yes / Oui" & my_data$eligibility_french == "No / Non" & eligibility_other_lan == "No / Non")
                 )

my_data <- subset(my_data, my_data$first_language == "English / Anglais" | my_data$first_language == "French / Français")
my_data <- subset(my_data, my_data$l1 == "English / Anglais" | my_data$l1 == "French / Français")
# Qualtrics already blocks participants who do not meet eligibility criteria
# Accepting blank eligibility responses because Prolific required us to do that externally so that subset weren't asked
# Should modify the above to only accept blank eligibility responses from Prolific participants (just in case, though should make no difference)


# Filter to confirm the entry comes from a legitimate source / recruitment avenue
my_data <- my_data[my_data$status == "IP Address", ]
my_data <- subset(my_data, my_data$source == "sona_cogsci" | my_data$source == "sona_psych" | my_data$source == "prolific_l1e" | my_data$source == "prolific_l1f")

# Filter by completion
my_data <- my_data[my_data$finished == "True", ]
# my_data <- my_data[my_data$duration_in_seconds >= 1200, ]   # This is now below

# Remove rows that were practice trials
# I believe they are all excluded under other checks above, but just in case
qualtrics_test_runs <- c("R_ex2UqHGk3N4rr4l", "R_0p64Pt4Pn1xe5fb", "R_2YbRN1GFsE1FhdG",
                         "R_1k17QY4T3sbgvwZ", "R_ywfTKwA3Qcc1lXH", "R_0p64Pt4Pn1xe5fb",
                         "R_65tMJaBcWfyAVQf", "R_3CrHVKG4Db5xIgB", "R_7aP80Int6z55TNk",
                         "R_5kjMPUh4nf7uvIz")
my_data <- my_data[!(my_data$response_id %in% qualtrics_test_runs), ]

# Keep only the first row for each unique random ID
my_data <- my_data %>% distinct(random_id, .keep_all = TRUE)
# Keep only the first row for each 6-digit SONA ID (note that false starts will already have been removed above)
my_data <- my_data %>% 
  group_by(id, source) %>%
  filter(nchar(id) != 6 | 
         !(source %in% c("sona_cogsci", "sona_psych")) | 
         row_number() == 1) %>%
  ungroup()

# # Adding new data to test the code above - should remove only one of the last two (works at last test)
# my_data <- add_row(my_data, id = "654321", source = "sona_cogsci")
# my_data <- add_row(my_data, id = "654321", source = "sona_psych")
# my_data <- add_row(my_data, id = "987654", source = "outside")
# my_data <- add_row(my_data, id = "987654", source = "outside")
# my_data <- add_row(my_data, id = "23467890", source = "sona_cogsci")
# my_data <- add_row(my_data, id = "23467890", source = "sona_cogsci")
# my_data <- add_row(my_data, id = "098764", source = "sona_cogsci")
# my_data <- add_row(my_data, id = "098764", source = "sona_cogsci")


# Extra line to exclude the 2023 participants (if this was wanted)
# my_data_prev <- my_data_prev[substr(my_data_prev$start_date, 1, 4) != "2023", ]



# Load and clean data from old Qualtrics survey version (up to end of 2023) ----
my_data_prev <- read.csv(file = "QualtricsData_2023_12_11.csv", encoding="UTF-8")
my_data_prev <- clean_names(my_data_prev)
my_data_prev <- my_data_prev[-c(1, 2), ]
my_data_prev$duration_in_seconds <- as.numeric(my_data_prev$duration_in_seconds)

# Filter by consent and eligibility
my_data_prev <- my_data_prev[my_data_prev$give_consent == "I consent to participate. / Je donne mon consentement à participer.", ]
my_data_prev <- subset(my_data_prev, my_data_prev$eligibility_age == "Yes / Oui" | 
                                     my_data_prev$eligibility_age == "")
my_data_prev <- subset(my_data_prev, my_data_prev$eligibility_language == "Yes / Oui" |
                                     my_data_prev$eligibility_language == "")
my_data_prev <- subset(my_data_prev, my_data_prev$first_language == "English / Anglais" | 
                                     my_data_prev$first_language == "French / Français")
my_data_prev <- subset(my_data_prev, my_data_prev$l1 == "English / Anglais" | 
                                     my_data_prev$l1 == "French / Français")


# Filter to confirm the entry comes from a legitimate source / recruitment avenue
my_data_prev <- my_data_prev[my_data_prev$status == "IP Address", ]
my_data_prev <- subset(my_data_prev, my_data_prev$source == "sona_cogsci" | my_data_prev$source == "sona_psych" | my_data_prev$source == "prolific_l1e" | my_data_prev$source == "prolific_l1f")

# Filter by completion and completion time
my_data_prev <- my_data_prev[my_data_prev$finished == "True", ]
my_data_prev <- my_data_prev[my_data_prev$duration_in_seconds >= 1200, ]

# Remove rows that were practice trials
# I believe they are all excluded under other checks above, but just in case
# qualtrics_test_runs is specified above
my_data_prev <- my_data_prev[!(my_data_prev$response_id %in% qualtrics_test_runs), ]

# Reduce to only the cols needed for merging to main df
my_data_prev <- my_data_prev[,c("random_id","dlq_e_1","dlq_e_11","dlq_f_1","dlq_f_11")]


# # Export selected cols
# my_data_prev_exp <- my_data_prev[,c(8,363:375)]
# setwd("C:/Users/Arthur/OneDrive/Documents/R/PhD_Thesis/SZ_ASD_L1_L2_analysis/Output")
# write.csv(my_data_prev_exp, "AH_Thesis_2024_ForQC.csv")
# setwd("C:/Users/Arthur/OneDrive/Documents/R/PhD_Thesis/SZ_ASD_L1_L2_analysis/Data")


# Create cols with language status variables ----

# Create col with L1 as 1 (English) or 0 (French)
my_data[,"is_l1_en"] <- ifelse(my_data[,"first_language"] == "English / Anglais", 1, 0)

# Create col with monolingualism as 1 (monolingual) or 0 (bilingual)
my_data[,"is_monoling"] <- ifelse(my_data[,"eligibility_french"] == "No / Non", 1, 0)


# Move attention check cols to end ----

my_data <- my_data %>% relocate(sm_e1_attention_check, .after = last_col())
my_data <- my_data %>% relocate(sm_e2_attention_check, .after = last_col())
my_data <- my_data %>% relocate(sm_f1_attention_check, .after = last_col())
my_data <- my_data %>% relocate(sm_f2_attention_check, .after = last_col())
my_data <- my_data %>% relocate( p_e1_attention_check, .after = last_col())
my_data <- my_data %>% relocate( p_e2_attention_check, .after = last_col())
my_data <- my_data %>% relocate( p_f1_attention_check, .after = last_col())
my_data <- my_data %>% relocate( p_f2_attention_check, .after = last_col())


# Binarize participant responses as correct or incorrect ----

# Loading answer key -- we will check participant responses against this below
my_key <- read.csv(file = "QualtricsAnswerKeyNoParticipants_2023_07_12.csv", encoding="UTF-8")
# Cleaning answer key with the janitor package so that the column names match the main data frame
my_key <- clean_names(my_key)

# Simplifying 4-point ASQ scale to 2-point scale (agree/disagree)
col_range <- which(colnames(my_data)=="asq_e_1"):which(colnames(my_data)=="asq_e_50")
for (my_col in col_range)
{
  my_data[,my_col] <- ifelse(my_data[,my_col] == "",
                             "",
                             ifelse(my_data[,my_col] == "Definitely agree" |
                                    my_data[,my_col] == "Slightly agree",
                                    "Agree",
                                    "Disagree")
                             )
}
rm(col_range)

# Same as above but for French
# Calls to the 3rd row of my_key are the French spellings to get the special characters
col_range <- which(colnames(my_data)=="asq_f_1"):which(colnames(my_data)=="asq_f_50")
for (my_col in col_range)
{
  my_data[,my_col] <- ifelse(my_data[,my_col] == "",
                             "",
                             ifelse(my_data[,my_col] == my_key[3,1] |
                                    my_data[,my_col] == my_key[3,2],
                                    my_key[3,5],
                                    my_key[3,6])
                             )
}
rm(col_range)

# Convert values to correct/incorrect (1 or 0) based on the answer key
for (my_col in 1:ncol(my_data))
{
  col_name <- colnames(my_data)[my_col]
  if(any(colnames(my_key) == col_name))   # Check if the column is in the key
  {
  # Note to self: may eventually want an expanded key
    if (my_key[2,col_name] != "")
    {
      my_data[,col_name] <- ifelse(my_data[,col_name] == "",
                                   NA,
                                   ifelse(my_data[,col_name] == my_key[2,col_name],
                                          1,
                                          0)
                                  )
    }
  }
}


# Load coded DLQ scores and merge into main df ----

# Import .csv files as dfs
dlq_3_data     <- read.csv(file = "DLQ_3_CodedByJenn.csv",     encoding="UTF-8")
dlq_4_5_6_data <- read.csv(file = "DLQ_4_5_6_CodedByJenn.csv", encoding="UTF-8")
dlq_7_data     <- read.csv(file = "DLQ_7_CodedByJenn.csv",     encoding="UTF-8")
dlq_8_data     <- read.csv(file = "DLQ_8_CodedByJenn.csv",     encoding="UTF-8")

# Clean names for the DLQ dfs
dlq_3_data     <- clean_names(dlq_3_data)
dlq_4_5_6_data <- clean_names(dlq_4_5_6_data)
dlq_7_data     <- clean_names(dlq_7_data)
dlq_8_data     <- clean_names(dlq_8_data)

# Remove remove non-data rows and cols
dlq_3_data     <- dlq_3_data[-c(1:2),c(1,7:11)]
dlq_4_5_6_data <- dlq_4_5_6_data[-c(1:3),c(1,9:14)]
dlq_7_data     <- dlq_7_data[-c(1:2),c(1,9:14)]
dlq_8_data     <- dlq_8_data[-c(1:2),c(1,11:18)]

# Adjust column names
colnames(dlq_3_data) <- c("random_id",
                          "en_flu_rank", "fr_flu_rank",
                          "en_acq_rank", "fr_acq_rank",
                          "num_langs")
colnames(dlq_4_5_6_data) <- c("response_id",
                              "en_exposure", "en_reading", "en_speaking",
                              "fr_exposure", "fr_reading", "fr_speaking")
colnames(dlq_7_data) <- c("random_id",
                          "en_country_yrs",     "fr_country_yrs",
                          "en_family_yrs",      "fr_family_yrs",
                          "en_school_work_yrs", "fr_school_work_yrs")
colnames(dlq_8_data) <- c("response_id",
                          "prof_e_speak",   "prof_f_speak",
                          "prof_e_underst", "prof_f_underst",
                          "prof_e_write",   "prof_f_write",
                          "prof_e_read",    "prof_f_read")

# Convert cols to numeric
for (my_col in 2:6)
{
  dlq_3_data[,my_col] <- as.numeric(dlq_3_data[,my_col])
}
for (my_col in 2:7)
{
  dlq_4_5_6_data[,my_col] <- as.numeric(dlq_4_5_6_data[,my_col])
}
for (my_col in 2:7)
{
  dlq_7_data[,my_col] <- as.numeric(dlq_7_data[,my_col])
}
for (my_col in 2:9)
{
  dlq_8_data[,my_col] <- as.numeric(dlq_8_data[,my_col])
}

# Merge into main df
my_data <- left_join(my_data, dlq_3_data,     by="random_id")
my_data <- left_join(my_data, dlq_4_5_6_data, by="response_id")
my_data <- left_join(my_data, dlq_7_data,     by="random_id")
my_data <- left_join(my_data, dlq_8_data,     by="response_id")

# Cleaning environment now that these are merged into main df
rm(dlq_3_data, dlq_4_5_6_data, dlq_7_data, dlq_8_data)


# Merge in values for LEAP-Q questions from more recent participants ----

merge_in_en_fr_vals <- function(my_df, dest_col, en_col, fr_col)
{
  for (cur_partic in 1:nrow(my_df))
  {
    if (is.na(my_df[cur_partic, dest_col]) == TRUE)
    {
      if (my_df[cur_partic, en_col] != "")
        my_df[cur_partic, dest_col] <- as.numeric(my_df[cur_partic, en_col])
      if (my_df[cur_partic, fr_col] != "")
        my_df[cur_partic, dest_col] <- as.numeric(my_df[cur_partic, fr_col])
      #my_df[cur_partic, dest_col] <- as.numeric(my_df[cur_partic, dest_col])
    }
  }
  return(my_df)
}

my_data <- merge_in_en_fr_vals(my_data, "en_flu_rank",    "dlq_e_3a_r_1", "dlq_f_3a_r_1")
my_data <- merge_in_en_fr_vals(my_data, "fr_flu_rank",    "dlq_e_3a_r_2", "dlq_f_3a_r_2")
my_data <- merge_in_en_fr_vals(my_data, "en_acq_rank",    "dlq_e_3b_r_1", "dlq_f_3b_r_1")
my_data <- merge_in_en_fr_vals(my_data, "fr_acq_rank",    "dlq_e_3b_r_2", "dlq_f_3b_r_2")
my_data <- merge_in_en_fr_vals(my_data, "en_exposure",    "dlq_e_4r_1",   "dlq_f_4r_1")
my_data <- merge_in_en_fr_vals(my_data, "fr_exposure",    "dlq_e_4r_2",   "dlq_f_4r_2")
my_data <- merge_in_en_fr_vals(my_data, "en_reading",     "dlq_e_5r_1",   "dlq_f_5r_1")
my_data <- merge_in_en_fr_vals(my_data, "fr_reading",     "dlq_e_5r_2",   "dlq_f_5r_2")
my_data <- merge_in_en_fr_vals(my_data, "en_speaking",    "dlq_e_6r_1",   "dlq_f_6r_1")
my_data <- merge_in_en_fr_vals(my_data, "fr_speaking",    "dlq_e_6r_2",   "dlq_f_6r_2")
my_data <- merge_in_en_fr_vals(my_data, "en_country_yrs", "dlq_e_7r_1",   "dlq_f_7r_1")
my_data <- merge_in_en_fr_vals(my_data, "fr_country_yrs", "dlq_e_7r_4",   "dlq_f_7r_4")
my_data <- merge_in_en_fr_vals(my_data, "en_family_yrs",  "dlq_e_7r_2",   "dlq_f_7r_2")
my_data <- merge_in_en_fr_vals(my_data, "fr_family_yrs",  "dlq_e_7r_5",   "dlq_f_7r_5")
my_data <- merge_in_en_fr_vals(my_data, "en_school_work_yrs", "dlq_e_7r_3", "dlq_f_7r_3")
my_data <- merge_in_en_fr_vals(my_data, "fr_school_work_yrs", "dlq_e_7r_6", "dlq_f_7r_6")
my_data <- merge_in_en_fr_vals(my_data, "prof_e_speak",   "dlq_e_8r_1",   "dlq_f_8r_5")
my_data <- merge_in_en_fr_vals(my_data, "prof_f_speak",   "dlq_e_8r_5",   "dlq_f_8r_1")
my_data <- merge_in_en_fr_vals(my_data, "prof_e_underst", "dlq_e_8r_2",   "dlq_f_8r_6")
my_data <- merge_in_en_fr_vals(my_data, "prof_f_underst", "dlq_e_8r_6",   "dlq_f_8r_2")
my_data <- merge_in_en_fr_vals(my_data, "prof_e_write",   "dlq_e_8r_3",   "dlq_f_8r_7")
my_data <- merge_in_en_fr_vals(my_data, "prof_f_write",   "dlq_e_8r_7",   "dlq_f_8r_3")
my_data <- merge_in_en_fr_vals(my_data, "prof_e_read",    "dlq_e_8r_4",   "dlq_f_8r_8")
my_data <- merge_in_en_fr_vals(my_data, "prof_f_read",    "dlq_e_8r_8",   "dlq_f_8r_4")


# WROTE THIS OUT TO KEEP THINGS STRAIGHT IN MY MIND
# English - order of fluency (1-8)  -> dlq_e_3a_r_1 OR dlq_f_3a_r_1
# French - order of fluency (1-8)   -> dlq_e_3a_r_2 OR dlq_f_3a_r_2
# English - order of acq (1-8)      -> dlq_e_3b_r_1 OR dlq_f_3b_r_1
# French - order of acq (1-8)       -> dlq_e_3b_r_2 OR dlq_f_3b_r_2
# English - cur exposure (%)        -> dlq_e_4r_1   OR dlq_f_4r_1
# French - cur exposure (%)         -> dlq_e_4r_2   OR dlq_f_4r_2
# English - choose to read text (%) -> dlq_e_5r_1   OR dlq_f_5r_1
# French - choose to read text (%)  -> dlq_e_5r_2   OR dlq_f_5r_2
# English - choose to speak (%)     -> dlq_e_6r_1   OR dlq_f_6r_1
# French - choose to speak (%)      -> dlq_e_6r_2   OR dlq_f_6r_2
# English - country (yrs)           -> dlq_e_7r_1   OR dlq_f_7r_1
# French - country (yrs)            -> dlq_e_7r_4   OR dlq_f_7r_4
# English - family (yrs)            -> dlq_e_7r_2   OR dlq_f_7r_2
# French - family (yrs)             -> dlq_e_7r_5   OR dlq_f_7r_5
# English - school or work (yrs)    -> dlq_e_7r_3   OR dlq_f_7r_3
# French - school or work (yrs)     -> dlq_e_7r_6   OR dlq_f_7r_6
# English - prof. spoken (0-10)     -> dlq_e_8r_1   OR dlq_f_8r_5
# English - prof. understand (0-10) -> dlq_e_8r_2   OR dlq_f_8r_6
# English - prof. written (0-10)    -> dlq_e_8r_3   OR dlq_f_8r_7
# English - prof. read (0-10)       -> dlq_e_8r_4   OR dlq_f_8r_8
# French - prof. spoken (0-10)      -> dlq_e_8r_5   OR dlq_f_8r_1
# French - prof. understand (0-10)  -> dlq_e_8r_6   OR dlq_f_8r_2
# French - prof. written (0-10)     -> dlq_e_8r_7   OR dlq_f_8r_3
# French - prof. reading (0-10)     -> dlq_e_8r_8   OR dlq_f_8r_4


# Merge in data from my_data_prev for age and diagnoses ----

for (cur_partic in 1:nrow(my_data_prev))
{
  # Determine row in up-to-date df with same random_id
  if (my_data_prev[cur_partic,"random_id"] != "")
    cur_id <- my_data_prev[cur_partic,"random_id"]
  # If there is one
  if (cur_id %in% my_data$random_id)
    # Assign it a value to call below
    corresp_row <- which(my_data$random_id == cur_id)

  # If no age (EN or FR) in the up-to-date df
  if (my_data[corresp_row,"dlq_e_1r_1"] == "" &
      my_data[corresp_row,"dlq_f_1r_1"] == "")
  {
    # If age (EN) in prev df
    if (my_data_prev[cur_partic,"dlq_e_1"] != "")
      # Overwrite age (EN) in up-to-date df with age (EN) from prev df
      my_data[corresp_row,"dlq_e_1r_1"] <- my_data_prev[cur_partic,"dlq_e_1"]
    # If age (FR) in prev df
    if (my_data_prev[cur_partic,"dlq_f_1"] != "")
      # Overwrite age (FR) in up-to-date df with age (FR) from prev df
      my_data[corresp_row,"dlq_f_1r_1"] <- my_data_prev[cur_partic,"dlq_f_1"]
  }

  # Repeat the above (for age) for the diagnoses col
  if (my_data[corresp_row,"dlq_e_11r"] == "" &
      my_data[corresp_row,"dlq_f_11r"] == "")
  {
    if (my_data_prev[cur_partic,"dlq_e_11"] != "")
      my_data[corresp_row,"dlq_e_11r"] <- my_data_prev[cur_partic,"dlq_e_11"]
    if (my_data_prev[cur_partic,"dlq_f_11"] != "")
      my_data[corresp_row,"dlq_f_11r"] <- my_data_prev[cur_partic,"dlq_f_11"]
  }
}


# Further processing of DLQ data (age, gender, diagnoses) ----

# Age
# my_data <- add_merged_col(my_data, "dlq_e_1", "dlq_f_1", "age")  # Old version of cols

# Somebody wrote "18 years old" -- could add code to turn that into 18

my_data <- add_merged_col(my_data, "dlq_e_1r_1", "dlq_f_1r_1", "age")
my_data$age <- as.numeric(my_data$age)

# Gender
my_data <- add_merged_col(my_data, "dlq_e_2", "dlq_f_2", "gender")
my_data[(my_data$gender == "Female" | my_data$gender == "Femme"), "gender"] <- "0"
my_data[(my_data$gender == "Male" | my_data$gender == "Homme"), "gender"] <- "1"
my_data[(my_data$gender == "Other (please specify)" | my_data$gender == "Autre (veuillez préciser)"), "gender"] <- "2"
my_data$gender <- as.numeric(my_data$gender)

# Diagnoses
# my_data <- add_merged_col(my_data, "dlq_e_11", "dlq_f_11", "diagnoses")   # Old version of cols
my_data <- add_merged_col(my_data, "dlq_e_11r", "dlq_f_11r", "diagnoses")

# There are also some additional options for participants, need to check the questions

# Creating columns by specific diagnosis
my_data$diag_asd   <- NA
my_data$diag_adhd  <- NA
my_data$diag_sli   <- NA
my_data$diag_dysl  <- NA
my_data$diag_epil  <- NA
my_data$diag_sz    <- NA
my_data$diag_mood  <- NA
my_data$diag_anx   <- NA

# Set column value to TRUE or FALSE for each diagnosis
for (my_row in 1:nrow(my_data))
{
  # My memory says just put | within the quotation marks to check for either of two strings, but can ChatGPT (already asked it)
  my_data[my_row, "diag_asd"]  <- grepl("Autism spectrum disorder|Trouble du spectre", my_data[my_row, "diagnoses"])
  my_data[my_row, "diag_adhd"] <- grepl("ADHD|TDAH", my_data[my_row, "diagnoses"])
  my_data[my_row, "diag_sli"]  <- grepl("Specific Language Impairment|cifique du langage", my_data[my_row, "diagnoses"])
  my_data[my_row, "diag_dysl"] <- grepl("Dyslexia|Dyslexie", my_data[my_row, "diagnoses"])
  my_data[my_row, "diag_epil"] <- grepl("Epilepsy|pilepsie", my_data[my_row, "diagnoses"])
  my_data[my_row, "diag_sz"]   <- grepl("Schizophrenia|trouble psychotique", my_data[my_row, "diagnoses"])
  my_data[my_row, "diag_mood"] <- grepl("Mood disorder|pression ou trouble bipolaire", my_data[my_row, "diagnoses"])
  my_data[my_row, "diag_anx"]  <- grepl("Anxiety disorder|Trouble anxieux", my_data[my_row, "diagnoses"])
}

# Remove participants whose reported age is too young
my_data <- my_data[(my_data$age >= 18 | is.na(my_data$age)),]


# Process files from Pavlovia and merge into main df ----

# Load all .csv files from the Pavlovia data folder
# folder_path_pavlovia_data <- ("C:/Users/Arthur/OneDrive/Documents/R/PhD_Thesis/SZ_ASD_L1_L2_analysis/Data/PavloviaData_2023_12_11")
folder_path_pavlovia_data <- ("C:/Users/Arthur/OneDrive/Documents/R/PhD_Thesis/SZ_ASD_L1_L2_analysis/Data/PavloviaData_2025_12_02")
pavlovia_files <- list.files(path = folder_path_pavlovia_data,
                             pattern = "\\.csv$",
                             full.names = TRUE)

# Create a function to check if a .csv file is empty
is_file_empty <- function(file_name)
{
  tryCatch(expr={my_csv <- read.csv(file_name, encoding="UTF-8", header=FALSE, nrows=20)
                 (nrow(my_csv) < 4 | ncol(my_csv) < 4)},
           error = function(e) {return(TRUE)}
          )
}

# Create vectors of non-empty and empty csv files
non_empty_pavlovia_files <- pavlovia_files[!sapply(pavlovia_files, is_file_empty)]
empty_pavlovia_files <- pavlovia_files[sapply(pavlovia_files, is_file_empty)]

# Create data frame to combine the info from many .csv files
pavlovia_df_cols <- c("random_id",
                      "dsct_total", "dsct_incorr",
                      "sst_1", "sst_2", "sst_3", "sst_total",
                      "fhat_1", "fhat_2", "fhat_3", "fhat_4", "fhat_5", "fhat_6", "fhat_7", "fhat_8", "fhat_9", "fhat_total")
pavlovia_data <- data.frame(matrix(ncol=length(pavlovia_df_cols), nrow=0))
colnames(pavlovia_data) <- pavlovia_df_cols
rm(pavlovia_df_cols)

# Correct responses for the FHAT task
fhat_answer_key <- c(3,1,3,2,1,1,2,3,2)
# Start and end cols to guarantee the correct cols are selected
start_col_fhat <- which(colnames(pavlovia_data)=="fhat_1")
end_col_fhat   <- which(colnames(pavlovia_data)=="fhat_9")

for (file_num in non_empty_pavlovia_files)
{
  # Load in participant's .csv file
  my_csv <- read.csv(file = file_num, encoding="UTF-8")
  my_csv <- clean_names(my_csv)

  # Remove files with no participant ID (mostly test runs)
  if (!("random_id" %in% colnames(my_csv)))
    next

  # Get participant ID
  pavlovia_data[file_num,"random_id"] <- my_csv[1,"random_id"]

  # Get DSCT info
  if ("resp_corr" %in% colnames(my_csv))
  {
    # Extract a vector of the DSCT responses (corr/incorr)
    non_na_groups_dsct <- rle(is.na(my_csv$resp_corr))
    # Note: There is -1 at the end because of a problem in the PsychoPy code, because of which they always get one more response after the timer runs out
    dsct_resps <- my_csv[(non_na_groups_dsct[["lengths"]][1] + 1):(non_na_groups_dsct[["lengths"]][1] + non_na_groups_dsct[["lengths"]][2] - 1),"resp_corr"]
    # Add DSCT scores to df
    pavlovia_data[file_num,"dsct_total"]   <- sum(dsct_resps == 1, na.rm=TRUE)
    pavlovia_data[file_num,"dsct_incorr"] <- sum(dsct_resps == 0, na.rm=TRUE)
  }

  # Get SST info
  if ("output_sst" %in% colnames(my_csv))
  {
    # Find groups of blank and non-blank rows in the column showing what they typed in the SST
    non_na_groups_sst <- rle(my_csv$output_sst != "")
    # This gives the maximum number of squares they correctly remembered by trial
    sst_correct <- c(non_na_groups_sst[["lengths"]][2], # first trial
                     non_na_groups_sst[["lengths"]][4], # second trial
                     non_na_groups_sst[["lengths"]][6]) # third trial
    # This gives the number of stages they correctly passed
    # I believe this is how this will be scored but this should be examined further
    pavlovia_data[file_num,"sst_1"] <- sst_correct[1] - 1
    pavlovia_data[file_num,"sst_2"] <- sst_correct[2] - 1
    pavlovia_data[file_num,"sst_3"] <- sst_correct[3] - 1
    pavlovia_data[file_num,"sst_total"] <- sum(pavlovia_data[file_num,c("sst_1","sst_2","sst_3")])
  }

  # Get FHAT info
  if ("fhat_resp_keys" %in% colnames(my_csv))
  {
    # Vector of the participant's responses for each FHAT trial
    non_na_groups_fhat <- rle(is.na(my_csv$fhat_resp_keys))
    fhat_resps <- my_csv[non_na_groups_fhat[["lengths"]][1] + 1:9,"fhat_resp_keys"]
    # Compare answer key to participant to determine if they were correct
    fhat_data <- data.frame(fhat_answer_key, fhat_resps)
    fhat_data$fhat_corr <- NA
    for (i in 1:9)
    {
      fhat_data[i,"fhat_corr"] <- ifelse(fhat_data[i,"fhat_answer_key"] == fhat_data[i,"fhat_resps"],
                                         1,
                                         0)
    }
    # Add FHAT scores by question to df
    pavlovia_data[file_num,start_col_fhat:end_col_fhat] <- fhat_data$fhat_corr
    # Add total FHAT score to df
    pavlovia_data[file_num,"fhat_total"] <- sum(pavlovia_data[file_num,start_col_fhat:end_col_fhat])
  }
}

# Keep only the first row for each unique ID value
pavlovia_data <- pavlovia_data %>% distinct(random_id, .keep_all = TRUE)

# # Before filtering bad entries from Pavlovia data, copy col to QC df
# my_qc <- left_join(my_qc,
#                    pavlovia_data[,c("random_id",
#                                     "dsct_incorr")],
#                    by="random_id")

# # Filter bad entries from Pavlovia data -- this needs expansion
# # Remove participants using the strategy where they just run their finger along random keys
# pavlovia_data <- pavlovia_data[pavlovia_data$dsct_incorr < 50,]

# Delete rows with no random ID
pavlovia_data <- pavlovia_data %>% filter(!is.na(random_id))

# Merge Pavlovia data into Qualtrics df
my_data <- left_join(my_data, pavlovia_data, by="random_id")

# Keeping environment clean
rm(fhat_data, non_na_groups_dsct, non_na_groups_sst, non_na_groups_fhat,
   dsct_resps, start_col_fhat, end_col_fhat, fhat_resps, sst_correct)

# Calculate total scores by participant for each task or questionnaire ----

# Start and end cols to guarantee the correct cols are selected
block_col_nums <- data.frame(block = character(),
                             first_col = numeric(),
                             last_col = numeric())

# Function to add a new block to the data frame
add_block_cols <- function(block_num, block_name, first_col_name, last_col_name)
{
  # Must use <<- instead of <- or it only modifies a local copy within the function
  block_col_nums[block_num,"block"]     <<- block_name
  block_col_nums[block_num,"first_col"] <<- which(colnames(my_data)==first_col_name)
  block_col_nums[block_num,"last_col"]  <<- which(colnames(my_data)==last_col_name)
}

# Add the column numbers for the blocks
add_block_cols(1, "lt_e",   "lt_e_1",   "lt_e_30")
add_block_cols(2, "lt_f",   "lt_f_1",   "lt_f_30")
add_block_cols(3, "ca_c_e1","ca_c_e1_1","ca_c_e1_16")
add_block_cols(4, "ca_c_f2","ca_c_f2_1","ca_c_f2_16")
add_block_cols(5, "ca_c_e2","ca_c_e2_1","ca_c_e2_16")
add_block_cols(6, "ca_c_f1","ca_c_f1_1","ca_c_f1_16")
add_block_cols(7, "sm_e1",  "sm_e1_1",  "sm_e1_16")
add_block_cols(8, "sm_f2",  "sm_f2_1",  "sm_f2_16")
add_block_cols(9, "sm_e2",  "sm_e2_1",  "sm_e2_16")
add_block_cols(10,"sm_f1",  "sm_f1_1",  "sm_f1_16")
add_block_cols(11,"p_e1",   "p_e1_1",   "p_e1_18")
add_block_cols(12,"p_f2",   "p_f2_1",   "p_f2_18")
add_block_cols(13,"p_e2",   "p_e2_1",   "p_e2_18")
add_block_cols(14,"p_f1",   "p_f1_1",   "p_f1_18")
add_block_cols(15,"mss_b_e","mss_b_e_1","mss_b_e_38")
add_block_cols(16,"asq_e",  "asq_e_1",  "asq_e_50")
add_block_cols(17,"mss_b_f","mss_b_f_1","mss_b_f_38")
add_block_cols(18,"asq_f",  "asq_f_1",  "asq_f_50")

# Total for each block in the experiment by participant
for (my_row in 1:nrow(my_data))
{
  # Calculating total score by block
  # sum_block is my function and really counts occurrences
  my_data[my_row, "lt_e_total"]    <- sum_block(my_data, my_row, "lt_e")
  my_data[my_row, "lt_f_total"]    <- sum_block(my_data, my_row, "lt_f")
  my_data[my_row, "ca_c_e1_total"] <- sum_block(my_data, my_row, "ca_c_e1")
  my_data[my_row, "ca_c_f2_total"] <- sum_block(my_data, my_row, "ca_c_f2")
  my_data[my_row, "ca_c_e2_total"] <- sum_block(my_data, my_row, "ca_c_e2")
  my_data[my_row, "ca_c_f1_total"] <- sum_block(my_data, my_row, "ca_c_f1")
  my_data[my_row, "sm_e1_total"]   <- sum_block(my_data, my_row, "sm_e1")
  my_data[my_row, "sm_f2_total"]   <- sum_block(my_data, my_row, "sm_f2")
  my_data[my_row, "sm_e2_total"]   <- sum_block(my_data, my_row, "sm_e2")
  my_data[my_row, "sm_f1_total"]   <- sum_block(my_data, my_row, "sm_f1")
  my_data[my_row, "p_e1_total"]    <- sum_block(my_data, my_row, "p_e1")
  my_data[my_row, "p_f2_total"]    <- sum_block(my_data, my_row, "p_f2")
  my_data[my_row, "p_e2_total"]    <- sum_block(my_data, my_row, "p_e2")
  my_data[my_row, "p_f1_total"]    <- sum_block(my_data, my_row, "p_f1")
  my_data[my_row, "mss_b_e_total"] <- sum_block(my_data, my_row, "mss_b_e")
  my_data[my_row, "asq_e_total"]   <- sum_block(my_data, my_row, "asq_e")
  my_data[my_row, "mss_b_f_total"] <- sum_block(my_data, my_row, "mss_b_f")
  my_data[my_row, "asq_f_total"]   <- sum_block(my_data, my_row, "asq_f")

  # Calculating total scores for self-rated proficiency
  my_data[my_row, "prof_e_total"] <- sum(my_data[my_row, c("prof_e_speak","prof_e_underst","prof_e_write","prof_e_read")])
  my_data[my_row, "prof_f_total"] <- sum(my_data[my_row, c("prof_f_speak","prof_f_underst","prof_f_write","prof_f_read")])
}

for (my_row in 1:nrow(my_data))
{
  # Calculating MSS-B subscale scores (English)
  my_data[my_row,"mss_b_e_neg"] <- sum(as.numeric(as.matrix(c(
                  my_data[my_row,"mss_b_e_1"],  my_data[my_row,"mss_b_e_4"],  my_data[my_row,"mss_b_e_7"],
                  my_data[my_row,"mss_b_e_10"], my_data[my_row,"mss_b_e_13"], my_data[my_row,"mss_b_e_16"],
                  my_data[my_row,"mss_b_e_19"], my_data[my_row,"mss_b_e_22"], my_data[my_row,"mss_b_e_25"],
                  my_data[my_row,"mss_b_e_28"], my_data[my_row,"mss_b_e_31"], my_data[my_row,"mss_b_e_34"],
                  my_data[my_row,"mss_b_e_37"]), na.rm=TRUE)))
  my_data[my_row,"mss_b_e_pos"] <- sum(as.numeric(as.matrix(c(
                  my_data[my_row,"mss_b_e_2"],  my_data[my_row,"mss_b_e_5"],  my_data[my_row,"mss_b_e_8"],
                  my_data[my_row,"mss_b_e_11"], my_data[my_row,"mss_b_e_14"], my_data[my_row,"mss_b_e_17"],
                  my_data[my_row,"mss_b_e_20"], my_data[my_row,"mss_b_e_23"], my_data[my_row,"mss_b_e_26"],
                  my_data[my_row,"mss_b_e_29"], my_data[my_row,"mss_b_e_32"], my_data[my_row,"mss_b_e_35"],
                  my_data[my_row,"mss_b_e_38"]), na.rm=TRUE)))
  my_data[my_row,"mss_b_e_dis"] <- sum(as.numeric(as.matrix(c(
                  my_data[my_row,"mss_b_e_3"],  my_data[my_row,"mss_b_e_6"],  my_data[my_row,"mss_b_e_9"],
                  my_data[my_row,"mss_b_e_12"], my_data[my_row,"mss_b_e_15"], my_data[my_row,"mss_b_e_18"],
                  my_data[my_row,"mss_b_e_21"], my_data[my_row,"mss_b_e_24"], my_data[my_row,"mss_b_e_27"],
                  my_data[my_row,"mss_b_e_30"], my_data[my_row,"mss_b_e_33"], my_data[my_row,"mss_b_e_36"]),
                  na.rm=TRUE)))

  # Calculating MSS-B subscale scores (French)
  my_data[my_row,"mss_b_f_neg"] <- sum(as.numeric(as.matrix(c(
                 my_data[my_row,"mss_b_f_1"],  my_data[my_row,"mss_b_f_4"],  my_data[my_row,"mss_b_f_7"],
                 my_data[my_row,"mss_b_f_10"], my_data[my_row,"mss_b_f_13"], my_data[my_row,"mss_b_f_16"],
                 my_data[my_row,"mss_b_f_19"], my_data[my_row,"mss_b_f_22"], my_data[my_row,"mss_b_f_25"],
                 my_data[my_row,"mss_b_f_28"], my_data[my_row,"mss_b_f_31"], my_data[my_row,"mss_b_f_34"],
                 my_data[my_row,"mss_b_f_37"]), na.rm=TRUE)))
  my_data[my_row,"mss_b_f_pos"] <- sum(as.numeric(as.matrix(c(
                 my_data[my_row,"mss_b_f_2"],  my_data[my_row,"mss_b_f_5"],  my_data[my_row,"mss_b_f_8"],
                 my_data[my_row,"mss_b_f_11"], my_data[my_row,"mss_b_f_14"], my_data[my_row,"mss_b_f_17"],
                 my_data[my_row,"mss_b_f_20"], my_data[my_row,"mss_b_f_23"], my_data[my_row,"mss_b_f_26"],
                 my_data[my_row,"mss_b_f_29"], my_data[my_row,"mss_b_f_32"], my_data[my_row,"mss_b_f_35"],
                 my_data[my_row,"mss_b_f_38"]), na.rm=TRUE)))
  my_data[my_row,"mss_b_f_dis"] <- sum(as.numeric(as.matrix(c(
                 my_data[my_row,"mss_b_f_3"],  my_data[my_row,"mss_b_f_6"],  my_data[my_row,"mss_b_f_9"],
                 my_data[my_row,"mss_b_f_12"], my_data[my_row,"mss_b_f_15"], my_data[my_row,"mss_b_f_18"],
                 my_data[my_row,"mss_b_f_21"], my_data[my_row,"mss_b_f_24"], my_data[my_row,"mss_b_f_27"],
                 my_data[my_row,"mss_b_f_30"], my_data[my_row,"mss_b_f_33"], my_data[my_row,"mss_b_f_36"]),
                 na.rm=TRUE)))

  # Calculating AQ subscale scores (English)
  my_data[my_row,"asq_e_comm"]       <- sum(as.numeric(as.matrix(c(
                  my_data[my_row,"asq_e_7"],  my_data[my_row,"asq_e_17"], my_data[my_row,"asq_e_18"],
                  my_data[my_row,"asq_e_26"], my_data[my_row,"asq_e_27"], my_data[my_row,"asq_e_31"],
                  my_data[my_row,"asq_e_33"], my_data[my_row,"asq_e_35"], my_data[my_row,"asq_e_38"],
                  my_data[my_row,"asq_e_39"]), na.rm=TRUE)))
  my_data[my_row,"asq_e_social"]     <- sum(as.numeric(as.matrix(c(
                  my_data[my_row,"asq_e_1"],  my_data[my_row,"asq_e_11"], my_data[my_row,"asq_e_13"],
                  my_data[my_row,"asq_e_15"], my_data[my_row,"asq_e_22"], my_data[my_row,"asq_e_36"],
                  my_data[my_row,"asq_e_44"], my_data[my_row,"asq_e_45"], my_data[my_row,"asq_e_47"],
                  my_data[my_row,"asq_e_48"]), na.rm=TRUE)))
  my_data[my_row,"asq_e_imagine"]    <- sum(as.numeric(as.matrix(c(
                  my_data[my_row,"asq_e_3"],  my_data[my_row,"asq_e_8"],  my_data[my_row,"asq_e_14"],
                  my_data[my_row,"asq_e_20"], my_data[my_row,"asq_e_21"], my_data[my_row,"asq_e_24"],
                  my_data[my_row,"asq_e_40"], my_data[my_row,"asq_e_41"], my_data[my_row,"asq_e_42"],
                  my_data[my_row,"asq_e_50"]), na.rm=TRUE)))
  my_data[my_row,"asq_e_detail"]     <- sum(as.numeric(as.matrix(c(
                  my_data[my_row,"asq_e_5"],  my_data[my_row,"asq_e_6"],  my_data[my_row,"asq_e_9"],
                  my_data[my_row,"asq_e_12"], my_data[my_row,"asq_e_19"], my_data[my_row,"asq_e_23"],
                  my_data[my_row,"asq_e_28"], my_data[my_row,"asq_e_29"], my_data[my_row,"asq_e_30"],
                  my_data[my_row,"asq_e_49"]), na.rm=TRUE)))
  my_data[my_row,"asq_e_att_switch"] <- sum(as.numeric(as.matrix(c(
                  my_data[my_row,"asq_e_2"],  my_data[my_row,"asq_e_4"],  my_data[my_row,"asq_e_10"],
                  my_data[my_row,"asq_e_16"], my_data[my_row,"asq_e_25"], my_data[my_row,"asq_e_32"],
                  my_data[my_row,"asq_e_34"], my_data[my_row,"asq_e_37"], my_data[my_row,"asq_e_43"],
                  my_data[my_row,"asq_e_46"]), na.rm=TRUE)))

  # Calculating AQ subscale scores (French)
  my_data[my_row,"asq_f_comm"]       <- sum(as.numeric(as.matrix(c(
                  my_data[my_row,"asq_f_7"],  my_data[my_row,"asq_f_17"], my_data[my_row,"asq_f_18"],
                  my_data[my_row,"asq_f_26"], my_data[my_row,"asq_f_27"], my_data[my_row,"asq_f_31"],
                  my_data[my_row,"asq_f_33"], my_data[my_row,"asq_f_35"], my_data[my_row,"asq_f_38"],
                  my_data[my_row,"asq_f_39"]), na.rm=TRUE)))
  my_data[my_row,"asq_f_social"]     <- sum(as.numeric(as.matrix(c(
                  my_data[my_row,"asq_f_1"],  my_data[my_row,"asq_f_11"], my_data[my_row,"asq_f_13"],
                  my_data[my_row,"asq_f_15"], my_data[my_row,"asq_f_22"], my_data[my_row,"asq_f_36"],
                  my_data[my_row,"asq_f_44"], my_data[my_row,"asq_f_45"], my_data[my_row,"asq_f_47"],
                  my_data[my_row,"asq_f_48"]), na.rm=TRUE)))
  my_data[my_row,"asq_f_imagine"]    <- sum(as.numeric(as.matrix(c(
                  my_data[my_row,"asq_f_3"],  my_data[my_row,"asq_f_8"],  my_data[my_row,"asq_f_14"],
                  my_data[my_row,"asq_f_20"], my_data[my_row,"asq_f_21"], my_data[my_row,"asq_f_24"],
                  my_data[my_row,"asq_f_40"], my_data[my_row,"asq_f_41"], my_data[my_row,"asq_f_42"],
                  my_data[my_row,"asq_f_50"]), na.rm=TRUE)))
  my_data[my_row,"asq_f_detail"]     <- sum(as.numeric(as.matrix(c(
                  my_data[my_row,"asq_f_5"],  my_data[my_row,"asq_f_6"],  my_data[my_row,"asq_f_9"],
                  my_data[my_row,"asq_f_12"], my_data[my_row,"asq_f_19"], my_data[my_row,"asq_f_23"],
                  my_data[my_row,"asq_f_28"], my_data[my_row,"asq_f_29"], my_data[my_row,"asq_f_30"],
                  my_data[my_row,"asq_f_49"]), na.rm=TRUE)))
  my_data[my_row,"asq_f_att_switch"] <- sum(as.numeric(as.matrix(c(
                  my_data[my_row,"asq_f_2"],  my_data[my_row,"asq_f_4"],  my_data[my_row,"asq_f_10"],
                  my_data[my_row,"asq_f_16"], my_data[my_row,"asq_f_25"], my_data[my_row,"asq_f_32"],
                  my_data[my_row,"asq_f_34"], my_data[my_row,"asq_f_37"], my_data[my_row,"asq_f_43"],
                  my_data[my_row,"asq_f_46"]), na.rm=TRUE)))
}

for (my_row in 1:nrow(my_data))
{
  # Merged cols for blocks of same task and language by participant
  my_data[my_row,"ca_c_e_total"] <- sum(c(
      as.numeric(as.matrix(my_data[my_row,"ca_c_e1_total"])),
      as.numeric(as.matrix(my_data[my_row,"ca_c_e2_total"]))), na.rm=TRUE)
  my_data[my_row,"ca_c_f_total"] <- sum(c(
      as.numeric(as.matrix(my_data[my_row,"ca_c_f1_total"])),
      as.numeric(as.matrix(my_data[my_row,"ca_c_f2_total"]))), na.rm=TRUE)
  my_data[my_row,"sm_e_total"] <- sum(c(
      as.numeric(as.matrix(my_data[my_row,"sm_e1_total"])),
      as.numeric(as.matrix(my_data[my_row,"sm_e2_total"]))), na.rm=TRUE)
  my_data[my_row,"sm_f_total"] <- sum(c(
      as.numeric(as.matrix(my_data[my_row,"sm_f1_total"])),
      as.numeric(as.matrix(my_data[my_row,"sm_f2_total"]))), na.rm=TRUE)
  my_data[my_row,"p_e_total"] <- sum(c(
      as.numeric(as.matrix(my_data[my_row,"p_e1_total"])),
      as.numeric(as.matrix(my_data[my_row,"p_e2_total"]))), na.rm=TRUE)
  my_data[my_row,"p_f_total"] <- sum(c(
      as.numeric(as.matrix(my_data[my_row,"p_f1_total"])),
      as.numeric(as.matrix(my_data[my_row,"p_f2_total"]))), na.rm=TRUE)

  # Merged cols for English and French symptom scores by participant
  my_data[my_row,"mss_b_total"] <- sum(c(
      as.numeric(as.matrix(my_data[my_row,"mss_b_e_total"])),
      as.numeric(as.matrix(my_data[my_row,"mss_b_f_total"]))), na.rm=TRUE)
  my_data[my_row,"mss_b_neg"] <- sum(c(
      as.numeric(as.matrix(my_data[my_row,"mss_b_e_neg"])),
      as.numeric(as.matrix(my_data[my_row,"mss_b_f_neg"]))), na.rm=TRUE)
  my_data[my_row,"mss_b_pos"] <- sum(c(
      as.numeric(as.matrix(my_data[my_row,"mss_b_e_pos"])),
      as.numeric(as.matrix(my_data[my_row,"mss_b_f_pos"]))), na.rm=TRUE)
  my_data[my_row,"mss_b_dis"] <- sum(c(
      as.numeric(as.matrix(my_data[my_row,"mss_b_e_dis"])),
      as.numeric(as.matrix(my_data[my_row,"mss_b_f_dis"]))), na.rm=TRUE)
  my_data[my_row,"asq_total"] <- sum(c(
      as.numeric(as.matrix(my_data[my_row,"asq_e_total"])),
      as.numeric(as.matrix(my_data[my_row,"asq_f_total"]))), na.rm=TRUE)
  my_data[my_row,"asq_comm"] <- sum(c(
      as.numeric(as.matrix(my_data[my_row,"asq_e_comm"])),
      as.numeric(as.matrix(my_data[my_row,"asq_f_comm"]))), na.rm=TRUE)
  my_data[my_row,"asq_social"] <- sum(c(
      as.numeric(as.matrix(my_data[my_row,"asq_e_social"])),
      as.numeric(as.matrix(my_data[my_row,"asq_f_social"]))), na.rm=TRUE)
  my_data[my_row,"asq_imagine"] <- sum(c(
      as.numeric(as.matrix(my_data[my_row,"asq_e_imagine"])),
      as.numeric(as.matrix(my_data[my_row,"asq_f_imagine"]))), na.rm=TRUE)
  my_data[my_row,"asq_detail"] <- sum(c(
      as.numeric(as.matrix(my_data[my_row,"asq_e_detail"])),
      as.numeric(as.matrix(my_data[my_row,"asq_f_detail"]))), na.rm=TRUE)
  my_data[my_row,"asq_att_switch"] <- sum(c(
      as.numeric(as.matrix(my_data[my_row,"asq_e_att_switch"])),
      as.numeric(as.matrix(my_data[my_row,"asq_f_att_switch"]))), na.rm=TRUE)
}


# Alternate LexTALE scoring ----

# Obtain vectors of "word" and "nonword" cols in LexTALE (English)
lt_e_word_cols    <- character(0)
lt_e_nonword_cols <- character(0)
first_col <- which(colnames(my_key)=="lt_e_1")
last_col  <- which(colnames(my_key)=="lt_e_30")
for (my_col in first_col:last_col)
{
  if (my_key[2,my_col] == "Yes")
    lt_e_word_cols    <- c(lt_e_word_cols, names(my_key)[my_col])
  if (my_key[2,my_col] == "No")
    lt_e_nonword_cols <- c(lt_e_nonword_cols, names(my_key)[my_col])
}
rm(first_col, last_col)

# Same for LexTALE (French)
lt_f_word_cols    <- character(0)
lt_f_nonword_cols <- character(0)
first_col <- which(colnames(my_key)=="lt_f_1")
last_col  <- which(colnames(my_key)=="lt_f_30")
for (my_col in first_col:last_col)
{
  if (my_key[2,my_col] == "Oui")
    lt_f_word_cols    <- c(lt_f_word_cols, names(my_key)[my_col])
  if (my_key[2,my_col] == "Non")
    lt_f_nonword_cols <- c(lt_f_nonword_cols, names(my_key)[my_col])
}
rm(first_col, last_col)

# Calculate subtotals for LexTALE (words / non-words), then weighted score as a percentage
for (my_row in 1:nrow(my_data))
{
  # Subtotals
  my_data[my_row,"lt_e_words"]    <- sum(my_data[my_row, lt_e_word_cols] == 1,
                                         na.rm=TRUE)
  my_data[my_row,"lt_e_nonwords"] <- sum(my_data[my_row, lt_e_nonword_cols] == 1,
                                         na.rm=TRUE)
  my_data[my_row,"lt_f_words"]    <- sum(my_data[my_row, lt_f_word_cols] == 1,
                                         na.rm=TRUE)
  my_data[my_row,"lt_f_nonwords"] <- sum(my_data[my_row, lt_f_nonword_cols] == 1,
                                         na.rm=TRUE)
  
  # Weighted scores
  my_data[my_row,"lt_e_weighted"] <- (my_data[my_row,"lt_e_words"]    * (50/20)) +
                                     (my_data[my_row,"lt_e_nonwords"] * (50/10))
  my_data[my_row,"lt_f_weighted"] <- (my_data[my_row,"lt_f_words"]    * (50/20)) +
                                     (my_data[my_row,"lt_f_nonwords"] * (50/10))  
}


# Determine is_l1pr_en for each participant ----

# Assigns self-declared L1 by default
my_data$is_l1pr_en <- my_data$is_l1_en
# Then overwrites it with language with higher self-declared proficiency (unless they're equal)
for (my_row in 1:nrow(my_data))
{
  if (!is.na(my_data[my_row,"prof_e_total"]) && !is.na(my_data[my_row,"prof_f_total"]))
  {
    if (my_data[my_row,"prof_e_total"] > my_data[my_row,"prof_f_total"])
      my_data[my_row,"is_l1pr_en"] <- 1
    if (my_data[my_row,"prof_f_total"] > my_data[my_row,"prof_e_total"])
      my_data[my_row,"is_l1pr_en"] <- 0
  }
}


# ====== Additional QC ====== ----

# Ideas for screening ----

# High number of wrong answers for DSCT (currently below)
# Check for many blank responses
# Check for accuracy around chance on tasks (1/2, 1/4, 1/2, 1/4)
# Flag those who indicated unlikely age (too high or low)
# Flag those who indicated implausible education level
# Other implausible DLQ scores to flag?
# First language in DLQ does not match first language indicated at start


# Adding cols to help with QC ----

# # Add a col to check whether under 20 mins
# my_data$under_20_mins <- FALSE
# my_data[my_data$duration_in_seconds < 1200, "under_20_mins"] <- TRUE
# my_data <- my_data %>% relocate(under_20_mins, .before = sm_e1_attention_check)

# Add truncated col for time spent
my_data$dur_seconds_trunc <- my_data$duration_in_seconds
my_data[my_data$dur_seconds_trunc > 7200, "dur_seconds_trunc"] <- 7200

# Prepare to run loop and rename attention check cols
attn_check_cols <- c("sm_e1_attention_check", "sm_e2_attention_check",
                     "sm_f1_attention_check", "sm_f2_attention_check",
                     "p_e1_attention_check",  "p_e2_attention_check",
                     "p_f1_attention_check",  "p_f2_attention_check")
failed_attn_names <- c("failed_attn_sm_e1", "failed_attn_sm_e2",
                       "failed_attn_sm_f1", "failed_attn_sm_f2",
                       "failed_attn_p_e1",  "failed_attn_p_e2",
                       "failed_attn_p_f1",  "failed_attn_p_f2")

# Set attn check cols TRUE or FALSE based on whether they failed the check
# And rename cols to "failed_attn..."
for (i in 1:length(attn_check_cols))
{
  my_col <- attn_check_cols[i]
  my_new_col <- failed_attn_names[i]
  my_data[[my_col]] <- ifelse(my_data[[my_col]] != "", TRUE, FALSE)
  my_data <- my_data %>% rename(!!my_new_col := all_of(my_col))
}

# Calculate total number of failed attn checks
my_data <- my_data %>% mutate(failed_attn_checks = rowSums(across(all_of(failed_attn_names)) == TRUE,
                                                       na.rm = TRUE))

# Keeping environment clean
rm(attn_check_cols, failed_attn_names)


# Second round of filtering of bad entries (mainly aimed at participants not taking it seriously) ----

# Filter by completion time
my_data <- my_data[(my_data$duration_in_seconds >= 1200 | 
                      is.na(my_data$duration_in_seconds)), ]

# Filtering those who used a "cheat" strategy on DSCT
my_data <- my_data[(my_data$dsct_incorr < 50 | 
                      is.na(my_data$dsct_incorr)), ]

# Removing participants who failed at least 2 attention checks
my_data <- my_data[(my_data$failed_attn_checks <= 1 | 
                     is.na(my_data$failed_attn_checks)), ]

# # Harsher criterion - unsure whether to use this
# # Could impute a "failed_attn_checks" score instead
# my_data <- my_data[my_data$sm_e_total > 4, ]
# my_data <- my_data[my_data$sm_f_total > 4, ]
# my_data <- my_data[my_data$p_e_total  > 4.5, ]
# my_data <- my_data[my_data$p_f_total  > 4.5, ]



# Excluding monolings, SZ, ASD ----

# Filter by bilingual status
my_data <- my_data[(my_data$is_monoling == 0 | 
                      is.na(my_data$is_monoling)), ]
# Filter by diagnoses
my_data <- my_data[(my_data$diag_sz == FALSE | 
                      is.na(my_data$diag_sz)), ]
my_data <- my_data[(my_data$diag_asd == FALSE | 
                      is.na(my_data$diag_asd)), ]


# ====== More processing of data ====== ----

# Create synthetic columns for analyzing L1/L2 rather than En/Fr ----

# Create columns for L1 and L2 totals for each linguistic task
my_data$lt_l1_total   <- convert_totals_to_l1_or_l2(
    my_data, 1,   "lt_e_total",   "lt_f_total")
my_data$lt_l2_total   <- convert_totals_to_l1_or_l2(
    my_data, 0,   "lt_e_total",   "lt_f_total")
my_data$ca_c_l1_total <- convert_totals_to_l1_or_l2(
    my_data, 1, "ca_c_e_total", "ca_c_f_total")
my_data$ca_c_l2_total <- convert_totals_to_l1_or_l2(
    my_data, 0, "ca_c_e_total", "ca_c_f_total")
my_data$sm_l1_total   <- convert_totals_to_l1_or_l2(
    my_data, 1,   "sm_e_total",   "sm_f_total")
my_data$sm_l2_total   <- convert_totals_to_l1_or_l2(
    my_data, 0,   "sm_e_total",   "sm_f_total")
my_data$p_l1_total    <- convert_totals_to_l1_or_l2(
    my_data, 1,    "p_e_total",    "p_f_total")
my_data$p_l2_total    <- convert_totals_to_l1_or_l2(
    my_data, 0,    "p_e_total",    "p_f_total")

# Same for weighted LexTALE scores
my_data$lt_l1_weighted <- convert_totals_to_l1_or_l2(
    my_data, 1, "lt_e_weighted", "lt_f_weighted")
my_data$lt_l2_weighted <- convert_totals_to_l1_or_l2(
    my_data, 0, "lt_e_weighted", "lt_f_weighted")

# Convert scores to L1/L2 for proficiency scores
my_data$prof_l1_speak   <- convert_totals_to_l1_or_l2(
    my_data, 1, "prof_e_speak",   "prof_f_speak")
my_data$prof_l2_speak   <- convert_totals_to_l1_or_l2(
    my_data, 0, "prof_e_speak",   "prof_f_speak")
my_data$prof_l1_underst <- convert_totals_to_l1_or_l2(
    my_data, 1, "prof_e_underst", "prof_f_underst")
my_data$prof_l2_underst <- convert_totals_to_l1_or_l2(
    my_data, 0, "prof_e_underst", "prof_f_underst")
my_data$prof_l1_write   <- convert_totals_to_l1_or_l2(
    my_data, 1, "prof_e_write",   "prof_f_write")
my_data$prof_l2_write   <- convert_totals_to_l1_or_l2(
    my_data, 0, "prof_e_write",   "prof_f_write")
my_data$prof_l1_read    <- convert_totals_to_l1_or_l2(
    my_data, 1, "prof_e_read",    "prof_f_read")
my_data$prof_l2_read    <- convert_totals_to_l1_or_l2(
    my_data, 0, "prof_e_read",    "prof_f_read")
my_data$prof_l1_total   <- convert_totals_to_l1_or_l2(
    my_data, 1, "prof_e_total",    "prof_f_total")
my_data$prof_l2_total   <- convert_totals_to_l1_or_l2(
    my_data, 0, "prof_e_total",    "prof_f_total")

# Convert scores to L1/L2 for other LEAP-Q data
my_data$exposure_l1   <- convert_totals_to_l1_or_l2(
    my_data, 1, "en_exposure",   "fr_exposure")
my_data$exposure_l2   <- convert_totals_to_l1_or_l2(
    my_data, 0, "en_exposure",   "fr_exposure")
my_data$reading_choice_l1   <- convert_totals_to_l1_or_l2(
    my_data, 1, "en_reading",   "fr_reading")
my_data$reading_choice_l2   <- convert_totals_to_l1_or_l2(
    my_data, 0, "en_reading",   "fr_reading")
my_data$speaking_choice_l1   <- convert_totals_to_l1_or_l2(
    my_data, 1, "en_speaking",   "fr_speaking")
my_data$speaking_choice_l2   <- convert_totals_to_l1_or_l2(
    my_data, 0, "en_speaking",   "fr_speaking")
my_data$country_yrs_l1   <- convert_totals_to_l1_or_l2(
    my_data, 1, "en_country_yrs",   "fr_country_yrs")
my_data$country_yrs_l2   <- convert_totals_to_l1_or_l2(
    my_data, 0, "en_country_yrs",   "fr_country_yrs")
my_data$family_yrs_l1   <- convert_totals_to_l1_or_l2(
    my_data, 1, "en_family_yrs",   "fr_family_yrs")
my_data$family_yrs_l2   <- convert_totals_to_l1_or_l2(
    my_data, 0, "en_family_yrs",   "fr_family_yrs")
my_data$school_work_yrs_l1   <- convert_totals_to_l1_or_l2(
    my_data, 1, "en_school_work_yrs",   "fr_school_work_yrs")
my_data$school_work_yrs_l2   <- convert_totals_to_l1_or_l2(
    my_data, 0, "en_school_work_yrs",   "fr_school_work_yrs")

# Convert cols to numeric
my_data <- my_data %>%
  mutate(
    lt_l1_total    = as.numeric(as.character(lt_l1_total)),
    lt_l2_total    = as.numeric(as.character(lt_l2_total)),
    ca_c_l1_total  = as.numeric(as.character(ca_c_l1_total)),
    ca_c_l2_total  = as.numeric(as.character(ca_c_l2_total)),
    sm_l1_total    = as.numeric(as.character(sm_l1_total)),
    sm_l2_total    = as.numeric(as.character(sm_l2_total)),
    p_l1_total     = as.numeric(as.character(p_l1_total)),
    p_l2_total     = as.numeric(as.character(p_l2_total)),
    lt_l1_weighted = as.numeric(as.character(lt_l1_weighted)),
    lt_l2_weighted = as.numeric(as.character(lt_l2_weighted))
  )

# Create columns for L1/L2 differences for each linguistic task
my_data$lt_l1_l2_difference   <- my_data$lt_l1_total   - my_data$lt_l2_total
my_data$ca_c_l1_l2_difference <- my_data$ca_c_l1_total - my_data$ca_c_l2_total
my_data$sm_l1_l2_difference   <- my_data$sm_l1_total   - my_data$sm_l2_total
my_data$p_l1_l2_difference    <- my_data$p_l1_total    - my_data$p_l2_total
my_data$prof_l1_l2_diff       <- my_data$prof_l1_total - my_data$prof_l2_total

# Create columns for L1/L2 ratios for linguistic tasks
my_data$lt_l1_l2_ratio   <- my_data$lt_l2_total   / my_data$lt_l1_total
my_data$ca_c_l1_l2_ratio <- my_data$ca_c_l2_total / my_data$ca_c_l1_total
my_data$sm_l1_l2_ratio   <- my_data$sm_l2_total   / my_data$sm_l1_total
my_data$p_l1_l2_ratio    <- my_data$p_l2_total    / my_data$p_l1_total


# # Create column for total L1/L2 difference across linguistic tasks
# # Need to double check why I made two versions and which is correct
# my_data$total_l1_l2_difference <- sum(c(my_data$lt_l1_l2_difference,
#                                         my_data$ca_c_l1_l2_difference,
#                                         my_data$sm_l1_l2_difference,
#                                         my_data$p_l1_l2_difference),
#                                         rm.na=FALSE)
# my_data$total_l1_l2_difference <- rowSums(my_data[,c("lt_l1_l2_difference",
#                                                      "ca_c_l1_l2_difference",
#                                                      "sm_l1_l2_difference",
#                                                      "p_l1_l2_difference")])


# # Remove outliers ----
# 
# Not sure this is best way for my data
# May be better to determine a participant not trying based on various indicators
#
# # Create list of cols to apply removal to
# lt_e_col <- which(colnames(my_data)=="lt_e_total")
# lt_f_col <- which(colnames(my_data)=="lt_f_total")
# first_outl_col <- which(colnames(my_data)=="ca_c_e_total")
# last_outl_col  <- which(colnames(my_data)=="p_f_total")
# cols_for_outl_removal <- c(lt_e_col, lt_f_col, first_outl_col:last_outl_col)
# 
# # Iterate over the list of cols and remove outliers using Tukey's criterion (lower bound only)
# for (my_col in cols_for_outl_removal)
# {
#   my_first_quart <- quantile(my_data[[my_col]], probs = .25, na.rm = FALSE)
#   my_IQR <- IQR(my_data[[my_col]])
#   my_data[my_data[[my_col]] < (my_first_quart - (1.5*my_IQR)), my_col] <- NA
# }
# 
# # Clean environment
# rm(lt_e_col, lt_f_col, first_outl_col, last_outl_col, cols_for_outl_removal)


# Calculate synthetic scores (e.g., z-scores) from block totals ----

for (my_row in 1:nrow(my_data))
{
  # Calculate task scores as percentages from totals
  # Excepted LexTALE, which is converted to % above as part of the weighting
  my_data[my_row,"ca_c_l1_pct"] <- my_data[my_row,"ca_c_l1_total"] * (100/16)
  my_data[my_row,"ca_c_l2_pct"] <- my_data[my_row,"ca_c_l2_total"] * (100/16)
  my_data[my_row,  "sm_l1_pct"] <- my_data[my_row,  "sm_l1_total"] * (100/16)
  my_data[my_row,  "sm_l2_pct"] <- my_data[my_row,  "sm_l2_total"] * (100/16)
  my_data[my_row,   "p_l1_pct"] <- my_data[my_row,   "p_l1_total"] * (100/18)
  my_data[my_row,   "p_l2_pct"] <- my_data[my_row,   "p_l2_total"] * (100/18)
  
  # # Calculate z-scores for linguistic tasks
  # # df$z_score <- (df$score - mean(df$score, na.rm = TRUE)) / sd(df$score, na.rm = TRUE)
  # my_data[my_row, "lt_e_z"]   <- (my_data[my_row, "lt_e_total"] - mean(my_data$lt_e_total, na.rm = TRUE)) / sd(my_data$lt_e_total, na.rm = TRUE)
  # my_data[my_row, "lt_f_z"]   <- (my_data[my_row, "lt_f_total"] - mean(my_data$lt_f_total, na.rm = TRUE)) / sd(my_data$lt_f_total, na.rm = TRUE)
  # my_data[my_row, "ca_c_e_z"] <- (my_data[my_row, "ca_c_e_total"] - mean(my_data$ca_c_e_total, na.rm = TRUE)) / sd(my_data$ca_c_e_total, na.rm = TRUE)
  # my_data[my_row, "ca_c_f_z"] <- (my_data[my_row, "ca_c_f_total"] - mean(my_data$ca_c_f_total, na.rm = TRUE)) / sd(my_data$ca_c_f_total, na.rm = TRUE)
  # my_data[my_row, "sm_e_z"]   <- (my_data[my_row, "sm_e_total"] - mean(my_data$sm_e_total, na.rm = TRUE)) / sd(my_data$sm_e_total, na.rm = TRUE)
  # my_data[my_row, "sm_f_z"]   <- (my_data[my_row, "sm_f_total"] - mean(my_data$sm_f_total, na.rm = TRUE)) / sd(my_data$sm_f_total, na.rm = TRUE)
  # my_data[my_row, "p_e_z"]    <- (my_data[my_row, "p_e_total"] - mean(my_data$p_e_total, na.rm = TRUE)) / sd(my_data$p_e_total, na.rm = TRUE)
  # my_data[my_row, "p_f_z"]    <- (my_data[my_row, "p_f_total"] - mean(my_data$p_f_total, na.rm = TRUE)) / sd(my_data$p_f_total, na.rm = TRUE)
}

# # Create columns for L1 and L2 z-scores for each linguistic task
# my_data$lt_l1_z   <- convert_totals_to_l1_or_l2(my_data, 1,   "lt_e_z",   "lt_f_z")
# my_data$lt_l2_z   <- convert_totals_to_l1_or_l2(my_data, 0,   "lt_e_z",   "lt_f_z")
# my_data$ca_c_l1_z <- convert_totals_to_l1_or_l2(my_data, 1, "ca_c_e_z", "ca_c_f_z")
# my_data$ca_c_l2_z <- convert_totals_to_l1_or_l2(my_data, 0, "ca_c_e_z", "ca_c_f_z")
# my_data$sm_l1_z   <- convert_totals_to_l1_or_l2(my_data, 1,   "sm_e_z",   "sm_f_z")
# my_data$sm_l2_z   <- convert_totals_to_l1_or_l2(my_data, 0,   "sm_e_z",   "sm_f_z")
# my_data$p_l1_z    <- convert_totals_to_l1_or_l2(my_data, 1,    "p_e_z",    "p_f_z")
# my_data$p_l2_z    <- convert_totals_to_l1_or_l2(my_data, 0,    "p_e_z",    "p_f_z")
# 
# # Create columns for diffs. b/w L1/L2 z-scores for linguistic tasks
# my_data$lt_l1_l2_z_diff   <- my_data$lt_l1_z   - my_data$lt_l2_z
# my_data$ca_c_l1_l2_z_diff <- my_data$ca_c_l1_z - my_data$ca_c_l2_z
# my_data$sm_l1_l2_z_diff   <- my_data$sm_l1_z   - my_data$sm_l2_z
# my_data$p_l1_l2_z_diff    <- my_data$p_l1_z    - my_data$p_l2_z
# my_data$all_l1_l2_z_diff  <- (my_data$lt_l1_l2_z_diff + 
#                              my_data$ca_c_l1_l2_z_diff + 
#                              my_data$sm_l1_l2_z_diff +
#                              my_data$p_l1_l2_z_diff)


# ====== Summarizing cleaned data ====== ----
#
# # Export a .csv with all or selected columns ----
#
# # Create df for export
# full_df_to_export <- my_data
# summary_df_to_export <- my_data[,c(492,493,503,508,518,519,
#                                    552:581)]
# 
# # Export selected columns
# selected_cols_to_export <- my_data[complete.cases(my_data[,c("mss_b_pos",
#                                                              "mss_b_neg",
#                                                              "mss_b_dis",
#                                                              "prof_l1_l2_diff",
#                                                              "lt_l1_l2_difference")]),
#                                    c("mss_b_pos", "mss_b_neg", "mss_b_dis",
#                                      "prof_l1_l2_diff", "lt_l1_l2_difference")]
# 
# # Do the export
# setwd("C:/Users/Arthur/OneDrive/Documents/R/PhD_Thesis/SZ_ASD_L1_L2_analysis/Output")
# # write.csv(full_df_to_export, "AH-Thesis-Data-Processed.csv")
# # write.csv(summary_df_to_export, "AH-Thesis-Data-Processed-Summary.csv")
# write.csv(selected_cols_to_export, "AH-Thesis-Power-Sim.csv")
# setwd("C:/Users/Arthur/OneDrive/Documents/R/PhD_Thesis/SZ_ASD_L1_L2_analysis/Data")
#
# # Another version
# en_cols_to_export <- my_data[my_data$l1 == "English / Anglais",
#                              c(8, 363, 365:377)]
# fr_cols_to_export <- my_data[my_data$l1 == "French / Français",
#                              c(8, 471, 473:485)]
# 
# # Do the export
# setwd("C:/Users/Arthur/OneDrive/Documents/R/PhD_Thesis/SZ_ASD_L1_L2_analysis/Output")
# write.csv(en_cols_to_export, "AH_Thesis_screening_cols_EN.csv")
# write.csv(fr_cols_to_export, "AH_Thesis_screening_cols_FR.csv")
# setwd("C:/Users/Arthur/OneDrive/Documents/R/PhD_Thesis/SZ_ASD_L1_L2_analysis/Data")


# Legend for interpreting column titles ----
#
# mss-b = schizotypy
# asq = autistic load
# lt = LexTale
# ca_c = Camel and Cactus
# sm = Syntactic Modification
# p = Pragmatics
# dsct = non-verbal processing speed
# sst = non-verbal working memory
# fhat = non-verbal theory of mind
# dlq = Demographic and Linguistic Questionnaire (mostly from LEAP-Q)



# Create list containing results ----

# Remove participants without covariates
my_data <- 
  my_data[complete.cases(my_data[, c("dsct_total",
                                     "sst_total",
                                     "fhat_total")]), ]

# Create the empty list
my_results <- list()

# Add sample size by characteristics
my_results[["sample_counts"]] <- as.data.frame(matrix(NA, nrow = 11, ncol = 2))
colnames(my_results[["sample_counts"]]) <- c("group", "N")
my_results[["sample_counts"]]$group <- c("Total",
                                         "Female",  "Male",  "OtherGender", "GenderNA",
                                         "L1_En",   "L1_Fr",
                                         "Monolingual", "Bilingual",
                                         "SONA", "Prolific")
my_results[["sample_counts"]][1,"N"] <- nrow(my_data)
my_results[["sample_counts"]][2,"N"] <- sum(my_data$gender == 0,
                                            na.rm = TRUE)
my_results[["sample_counts"]][3,"N"] <- sum(my_data$gender == 1,
                                            na.rm = TRUE)
my_results[["sample_counts"]][4,"N"] <- sum(my_data$gender == 2,
                                            na.rm = TRUE)
my_results[["sample_counts"]][5,"N"] <- sum(is.na(my_data$gender))
my_results[["sample_counts"]][6,"N"] <- sum(my_data$is_l1_en == 1,
                                            na.rm = TRUE)
my_results[["sample_counts"]][7,"N"] <- sum(my_data$is_l1_en == 0,
                                            na.rm = TRUE)
my_results[["sample_counts"]][8,"N"] <- sum(my_data$is_monoling == 1,
                                            na.rm = TRUE)
my_results[["sample_counts"]][9,"N"] <- sum(my_data$is_monoling == 0,
                                            na.rm = TRUE)
my_results[["sample_counts"]][10,"N"] <- 
    sum(my_data$source == "sona_cogsci" | 
          my_data$source == "sona_psych",
        na.rm = TRUE)
my_results[["sample_counts"]][11,"N"] <- 
    sum(my_data$source == "prolific_l1e" | 
          my_data$source == "prolific_l1f",
        na.rm = TRUE)

# Add diagnosis counts to my_results
my_results[["diag_counts"]] <- as.data.frame(matrix(NA, nrow = 8, ncol = 2))
colnames(my_results[["diag_counts"]]) <- c("diagnosis", "N")
my_results[["diag_counts"]]$diagnosis <- c("ASD",          "ADHD",
                                           "SLI",          "Dyslexia",
                                           "Epilepsy",     "Schizophrenia",
                                           "MoodDisorder", "AnxietyDisorder")
my_results[["diag_counts"]][1,"N"] <- sum(my_data$diag_asd  == TRUE,
                                          na.rm = TRUE)
my_results[["diag_counts"]][2,"N"] <- sum(my_data$diag_adhd == TRUE,
                                          na.rm = TRUE)
my_results[["diag_counts"]][3,"N"] <- sum(my_data$diag_sli  == TRUE,
                                          na.rm = TRUE)
my_results[["diag_counts"]][4,"N"] <- sum(my_data$diag_dysl == TRUE,
                                          na.rm = TRUE)
my_results[["diag_counts"]][5,"N"] <- sum(my_data$diag_epil == TRUE,
                                          na.rm = TRUE)
my_results[["diag_counts"]][6,"N"] <- sum(my_data$diag_sz   == TRUE,
                                          na.rm = TRUE)
my_results[["diag_counts"]][7,"N"] <- sum(my_data$diag_mood == TRUE,
                                          na.rm = TRUE)
my_results[["diag_counts"]][8,"N"] <- sum(my_data$diag_anx  == TRUE,
                                          na.rm = TRUE)

# Add means and SDs of sample characteristics
my_results[["means_sds"]] <- as.data.frame(matrix(NA, nrow = 21, ncol = 3))
colnames(my_results[["means_sds"]]) <- c("Variable", "M", "SD")

# Function to fill means/SDs table
add_entry_means_sds <- function(row_num, my_name, my_col)
{
  # Must use <<- instead of <- or it only modifies a local copy within the function
  my_results[["means_sds"]][row_num,"Variable"] <<- my_name
  my_results[["means_sds"]][row_num,"M"]        <<- mean(my_data[[my_col]],
                                                      na.rm = TRUE)
  my_results[["means_sds"]][row_num,"SD"]       <<- sd(my_data[[my_col]],
                                                      na.rm = TRUE)
}

# Populate the means/SDs df
add_entry_means_sds(1, "Age", "age")
add_entry_means_sds(2, "L1_Prof_Total", "prof_l1_total")
add_entry_means_sds(3, "L2_Prof_Total", "prof_l2_total")
add_entry_means_sds(4, "L1_Prof_Speak", "prof_l1_speak")
add_entry_means_sds(5, "L2_Prof_Speak", "prof_l2_speak")
add_entry_means_sds(6, "L1_Prof_Underst", "prof_l1_underst")
add_entry_means_sds(7, "L2_Prof_Underst", "prof_l2_underst")
add_entry_means_sds(8, "L1_Prof_Write", "prof_l1_write")
add_entry_means_sds(9, "L2_Prof_Write", "prof_l2_write")
add_entry_means_sds(10, "L1_Prof_Read", "prof_l1_read")
add_entry_means_sds(11, "L2_Prof_Read", "prof_l2_read")
add_entry_means_sds(12, "MSS-B_Total", "mss_b_total")
add_entry_means_sds(13, "MSS-B_Pos", "mss_b_pos")
add_entry_means_sds(14, "MSS-B_Neg", "mss_b_neg")
add_entry_means_sds(15, "MSS-B_Dis", "mss_b_dis")
add_entry_means_sds(16, "AQ_Total", "asq_total")
add_entry_means_sds(17, "AQ_Comm", "asq_comm")
add_entry_means_sds(18, "AQ_Social", "asq_social")
add_entry_means_sds(19, "AQ_Imagine", "asq_imagine")
add_entry_means_sds(20, "AQ_Detail", "asq_detail")
add_entry_means_sds(21, "AQ_Att_Switch", "asq_att_switch")
add_entry_means_sds(22, "L1_Cur_Exposure", "exposure_l1")
add_entry_means_sds(23, "L2_Cur_Exposure", "exposure_l2")
add_entry_means_sds(24, "L1_Reading_Choice", "reading_choice_l1")
add_entry_means_sds(25, "L2_Reading_Choice", "reading_choice_l2")
add_entry_means_sds(26, "L1_Speaking_Choice", "speaking_choice_l1")
add_entry_means_sds(27, "L2_Speaking_Choice", "speaking_choice_l2")
add_entry_means_sds(28, "L1_Country_Yrs", "country_yrs_l1")
add_entry_means_sds(29, "L2_Country_Yrs", "country_yrs_l2")
add_entry_means_sds(30, "L1_Family_Yrs", "family_yrs_l1")
add_entry_means_sds(31, "L2_Family_Yrs", "family_yrs_l2")
add_entry_means_sds(32, "L1_School_Work_Yrs", "school_work_yrs_l1")
add_entry_means_sds(33, "L2_School_Work_Yrs", "school_work_yrs_l2")

# Get means and SDs for task performance
my_results[["task_means_sds"]] <- as.data.frame(matrix(NA, nrow = 11, ncol = 3))
colnames(my_results[["task_means_sds"]]) <- c("Task", "M", "SD")

# Extra col for SST average score
my_data$sst_average <- my_data$sst_total / 3

# Function to fill task means/SDs table
add_entry_task_means_sds <- function(row_num, my_name, my_col)
{
  # Must use <<- instead of <- or it only modifies a local copy within the function
  my_results[["task_means_sds"]][row_num,"Task"] <<- my_name
  my_results[["task_means_sds"]][row_num,"M"]    <<- mean(my_data[[my_col]],
                                                          na.rm = TRUE)
  my_results[["task_means_sds"]][row_num,"SD"]   <<- sd(my_data[[my_col]],
                                                        na.rm = TRUE)
}

# Populate the task means/SDs df
add_entry_task_means_sds(1, "LexTALE_L1", "lt_l1_weighted")
add_entry_task_means_sds(2, "LexTALE_L2", "lt_l2_weighted")
add_entry_task_means_sds(3, "C&C_L1", "ca_c_l1_pct")
add_entry_task_means_sds(4, "C&C_L2", "ca_c_l2_pct")
add_entry_task_means_sds(5, "SM_L1", "sm_l1_pct")
add_entry_task_means_sds(6, "SM_L2", "sm_l2_pct")
add_entry_task_means_sds(7, "Pragmatics_L1", "p_l1_pct")
add_entry_task_means_sds(8, "Pragmatics_L2", "p_l2_pct")
add_entry_task_means_sds(9, "DSCT", "dsct_total")
add_entry_task_means_sds(10, "SST", "sst_average")
add_entry_task_means_sds(11, "FHAT", "fhat_total")


# ====== Descriptives ====== ----

# Histograms ----

# 1. Histogram for L1/L2 current exposure
exposure_long <- my_data %>%
  select(exposure_l1, exposure_l2) %>%
  pivot_longer(
    cols = everything(),
    names_to = "language",
    values_to = "exposure"
  )
hist_exposure <- ggplot(
    data=exposure_long, aes(x=exposure, fill=language, na.rm = TRUE)) +
    geom_histogram(position="identity", color="gray60",
                   binwidth=10, boundary=0,
                   alpha=0.5, na.rm=TRUE) +
    scale_fill_manual(
      values = c(exposure_l1 = color_l1,
                 exposure_l2 = color_l2),
      labels = c(exposure_l1 = "L1",
                 exposure_l2 = "L2"),
      name = "Language") +
    xlab("Current Exposure (%)") + ylab("Count") +
    # coord_cartesian(xlim = c(0, 40)) +
    scale_x_continuous(breaks = seq(0,100,20)) +
    theme_bw() +
    theme(plot.title=element_text(size=11, hjust=0.5),
          axis.title.x=element_text(size=9),
          axis.title.y=element_text(size=9),
          legend.position = "none")
rm(exposure_long)

# 2. Histogram for L1/L2 reading choice
reading_choice_long <- my_data %>%
  select(reading_choice_l1, reading_choice_l2) %>%
  pivot_longer(
    cols = everything(),
    names_to = "language",
    values_to = "reading_choice"
  )
# This line is a bit ad hoc, should be addressed above
reading_choice_long <- 
  reading_choice_long[reading_choice_long$reading_choice <= 100,]
# Return to plotting
hist_reading_choice <- ggplot(
    data=reading_choice_long, 
    aes(x=reading_choice, fill=language, na.rm = TRUE)) +
    geom_histogram(position="identity", color="gray60",
                   binwidth=10, boundary=0,
                   alpha=0.5, na.rm=TRUE) +
    scale_fill_manual(
      values = c(reading_choice_l1 = color_l1,
                 reading_choice_l2 = color_l2),
      labels = c(reading_choice_l1 = "L1",
                 reading_choice_l2 = "L2"),
      name = "Language") +
    xlab("Choice to Read In (%)") + ylab(NULL) +
    # coord_cartesian(xlim = c(0, 40)) +
    scale_x_continuous(breaks = seq(0,100,20)) +
    theme_bw() +
    theme(plot.title=element_text(size=11, hjust=0.5),
          axis.title.x=element_text(size=9),
          axis.title.y=element_text(size=9),
          legend.position = "none")
rm(reading_choice_long)

# 3. Histogram for L1/L2 speaking choice
speaking_choice_long <- my_data %>%
  select(speaking_choice_l1, speaking_choice_l2) %>%
  pivot_longer(
    cols = everything(),
    names_to = "language",
    values_to = "speaking_choice"
  )
hist_speaking_choice <- ggplot(
    data=speaking_choice_long, 
    aes(x=speaking_choice, fill=language, na.rm = TRUE)) +
    geom_histogram(position="identity", color="gray60",
                   binwidth=10, boundary=0,
                   alpha=0.5, na.rm=TRUE) +
    scale_fill_manual(
      values = c(speaking_choice_l1 = color_l1,
                 speaking_choice_l2 = color_l2),
      labels = c(speaking_choice_l1 = "L1",
                 speaking_choice_l2 = "L2"),
      name = "Language") +
    xlab("Choice to Speak In (%)") + ylab(NULL) +
    # ggtitle("Subjective Proficiency") +
    # coord_cartesian(xlim = c(0, 40)) +
    scale_x_continuous(breaks = seq(0,100,20)) +
    theme_bw() +
    theme(plot.title=element_text(size=11, hjust=0.5),
          axis.title.x=element_text(size=9),
          axis.title.y=element_text(size=9),
          legend.position = "none")
rm(speaking_choice_long)

# 4. Histogram for years lived in country
country_yrs_long <- my_data %>%
  select(country_yrs_l1, country_yrs_l2) %>%
  pivot_longer(
    cols = everything(),
    names_to = "language",
    values_to = "country_yrs"
  )
# This line is a bit ad hoc, should be addressed above
country_yrs_long <- 
  country_yrs_long[country_yrs_long$country_yrs < 60,]
# Return to plotting
hist_country_yrs <- ggplot(
    data=country_yrs_long, 
    aes(x=country_yrs, fill=language, na.rm = TRUE)) +
    geom_histogram(position="identity", color="gray60",
                   binwidth=2, boundary=0,
                   alpha=0.5, na.rm=TRUE) +
    scale_fill_manual(
      values = c(country_yrs_l1 = color_l1,
                 country_yrs_l2 = color_l2),
      labels = c(country_yrs_l1 = "L1",
                 country_yrs_l2 = "L2"),
      name = "Language") +
    xlab("Time in Country (Yrs.)") + ylab("Count") +
    # ggtitle("Subjective Proficiency") +
    # coord_cartesian(xlim = c(0, 40)) +
    # scale_x_continuous(breaks = seq(0,100,20)) +
    theme_bw() +
    theme(plot.title=element_text(size=11, hjust=0.5),
          axis.title.x=element_text(size=9),
          axis.title.y=element_text(size=9),
          legend.position = "none")
rm(country_yrs_long)

# 5. Histogram for years lived in family
family_yrs_long <- my_data %>%
  select(family_yrs_l1, family_yrs_l2) %>%
  pivot_longer(
    cols = everything(),
    names_to = "language",
    values_to = "family_yrs"
  )
# This line is a bit ad hoc, should be addressed above
family_yrs_long <- 
  family_yrs_long[family_yrs_long$family_yrs < 60,]
# Return to plotting
hist_family_yrs <- ggplot(
    data=family_yrs_long, 
    aes(x=family_yrs, fill=language, na.rm = TRUE)) +
    geom_histogram(position="identity", color="gray60",
                   binwidth=2, boundary=0,
                   alpha=0.5, na.rm=TRUE) +
    scale_fill_manual(
      values = c(family_yrs_l1 = color_l1,
                 family_yrs_l2 = color_l2),
      labels = c(family_yrs_l1 = "L1",
                 family_yrs_l2 = "L2"),
      name = "Language") +
    xlab("Time in Family (Yrs.)") + ylab(NULL) +
    # ggtitle("Subjective Proficiency") +
    # coord_cartesian(xlim = c(0, 40)) +
    # scale_x_continuous(breaks = seq(0,100,20)) +
    theme_bw() +
    theme(plot.title=element_text(size=11, hjust=0.5),
          axis.title.x=element_text(size=9),
          axis.title.y=element_text(size=9),
          legend.position = "none")
rm(family_yrs_long)

# 6. Histogram for years lived in school/work environment
school_work_yrs_long <- my_data %>%
  select(school_work_yrs_l1, school_work_yrs_l2) %>%
  pivot_longer(
    cols = everything(),
    names_to = "language",
    values_to = "school_work_yrs"
  )
# This line is a bit ad hoc, should be addressed above
school_work_yrs_long <- 
  school_work_yrs_long[school_work_yrs_long$school_work_yrs < 60,]
# Return to plotting
hist_school_work_yrs <- ggplot(
    data=school_work_yrs_long, 
    aes(x=school_work_yrs, fill=language, na.rm = TRUE)) +
    geom_histogram(position="identity", color="gray60",
                   binwidth=2, boundary=0,
                   alpha=0.5, na.rm=TRUE) +
    scale_fill_manual(
      values = c(school_work_yrs_l1 = color_l1,
                 school_work_yrs_l2 = color_l2),
      labels = c(school_work_yrs_l1 = "L1",
                 school_work_yrs_l2 = "L2"),
      name = "Language") +
    xlab("Time in School/Job (Yrs.)") + ylab(NULL) +
    # coord_cartesian(xlim = c(0, 40)) +
    # scale_x_continuous(breaks = seq(0,100,20)) +
    theme_bw() +
    theme(plot.title=element_text(size=11, hjust=0.5),
          axis.title.x=element_text(size=9),
          axis.title.y=element_text(size=9),
          legend.position = "none")
rm(school_work_yrs_long)

# 7. Histogram for L2 Proficiency
profic_long <- my_data %>%
  select(prof_l1_total, prof_l2_total) %>%
  pivot_longer(
    cols = everything(),
    names_to = "language",
    values_to = "proficiency"
  )
hist_profic <- ggplot(
    data=profic_long, aes(x=proficiency, fill=language, na.rm = TRUE)) +
    geom_histogram(position="identity", color="gray60",
                   binwidth=4, boundary=0,
                   alpha=0.5, na.rm=TRUE) +
    scale_fill_manual(
      values = c(prof_l1_total = color_l1,
                 prof_l2_total = color_l2),
      labels = c(prof_l1_total = "L1",
                 prof_l2_total = "L2"),
      name = "Language") +
    xlab("Subjective Proficiency") + ylab("Count") +
    coord_cartesian(xlim = c(0, 40)) +
    # scale_y_break(c(60, 560)) +
    # scale_y_continuous(
    #   limits = c(0, 800),
    #   breaks = c(0, 20, 40, 60, 560, 580, 600)) +
    theme_bw() +
    theme(plot.title=element_text(size=11, hjust=0.5),
          axis.title.x=element_text(size=9),
          axis.title.y=element_text(size=9),
          axis.text.y.right  = element_blank(),
          axis.ticks.y.right = element_blank(),
          axis.line.y.right  = element_blank())#,
          #legend.position = "none")
rm(profic_long)

# 8. Histogram for MSS-B
hist_mss_b <- ggplot(
    data=my_data, aes(x=mss_b_total, na.rm = TRUE)) +
    geom_histogram(fill=color_sz, color="gray60",
                   binwidth=2, boundary=0,
                   na.rm=TRUE) +
    xlab("MSS-B Total Score") + ylab(NULL) +
    ggtitle("Schizotypy") +
    theme_bw() +
    theme(plot.title=element_text(size=11, hjust=0.5),
          axis.title.x=element_text(size=9),
          axis.title.y=element_text(size=9)) +
    coord_cartesian(xlim = c(0, 38))

# 9. Histogram for AQ
hist_asq <- ggplot(
    data=my_data, aes(x=asq_total, na.rm = TRUE)) +
    geom_histogram(fill=color_asd, color="gray60",
                   binwidth=2, boundary=0,
                   na.rm=TRUE) +
    xlab("AQ Total Score") + ylab(NULL) +
    ggtitle("Autistic Traits") +
    theme_bw() +
    theme(plot.title=element_text(size=11, hjust=0.5),
          axis.title.x=element_text(size=9),
          axis.title.y=element_text(size=9)) +
    coord_cartesian(xlim = c(0, 50))

# Design for patchwork
design <- "ABC\nDEF\n.GG"
# Create gallery of seven histograms
# It displays the first one in the bottom row, I can't figure out why but oh well it works with this code
hist_profic + hist_exposure + hist_reading_choice + hist_speaking_choice + 
  hist_country_yrs + hist_family_yrs + hist_school_work_yrs +
  plot_layout(design=design, guides = "collect") &
    theme(legend.position = "right")

# Create gallery of seven histograms
(hist_exposure | hist_reading_choice | hist_speaking_choice) /
(hist_country_yrs | hist_family_yrs | hist_school_work_yrs) /
(hist_profic | plot_spacer() | plot_spacer())  +
  plot_layout(widths = unit(c(1,1,1), "null"),
              guides = "collect") &
  theme(legend.position = "right")

# Two histograms together
hist_mss_b + hist_asq


# Violin plot of task scores in L1 and L2 ----

# Create smaller df for this analysis
my_data_violin <- my_data[,c("random_id",
                             "lt_l1_total", "ca_c_l1_total", "sm_l1_total", "p_l1_total",
                             "lt_l2_total", "ca_c_l2_total", "sm_l2_total", "p_l2_total",
                             "is_l1pr_en")]
# Change scores to percentages
my_data_violin$lt_l1_total   <- my_data_violin$lt_l1_total   * 100/30
my_data_violin$ca_c_l1_total <- my_data_violin$ca_c_l1_total * 100/16
my_data_violin$sm_l1_total   <- my_data_violin$sm_l1_total   * 100/16
my_data_violin$p_l1_total    <- my_data_violin$p_l1_total    * 100/18
my_data_violin$lt_l2_total   <- my_data_violin$lt_l2_total   * 100/30
my_data_violin$ca_c_l2_total <- my_data_violin$ca_c_l2_total * 100/16
my_data_violin$sm_l2_total   <- my_data_violin$sm_l2_total   * 100/16
my_data_violin$p_l2_total    <- my_data_violin$p_l2_total    * 100/18

# Switch to long format
violin_data <- tidyr::pivot_longer(my_data_violin,
                                   cols = c("lt_l1_total", "ca_c_l1_total", "sm_l1_total", "p_l1_total",
                                            "lt_l2_total", "ca_c_l2_total", "sm_l2_total", "p_l2_total"),
                                   names_to = "task_and_lang",
                                   values_to = "score")
violin_data$Task     <- rep(c("LexTALE", "C&C", "S.M.", "Pragmatics"),
                            nrow(violin_data) / 4)
violin_data$Language <- rep(c(rep("L1", 4), rep("L2", 4)),
                            nrow(violin_data) / 8)

# Make task an ordered factor
violin_data$Task <- factor(violin_data$Task,
                           unique(violin_data$Task),
                           ordered=TRUE)

# Violin plot of the data by task and language
ggplot(violin_data, aes(x = Task, y = score, fill = as.factor(Language))) +
       geom_violin(position = position_dodge(width = 0.8), adjust = 0.5, bw = 6) +
       facet_wrap( ~ Task, scales = "free_x", ncol = 4) +
       scale_fill_manual(values = c(color_l1, color_l2), name = "Language") +
       labs(x= NULL, y = "Mean Total Score (%)", title = element_blank()) +
       theme_minimal() +
       #theme(axis.text.x = element_text(angle = 45, hjust = 1),
       theme(#axis.text.x = element_text(hjust = 1),
             strip.text = element_blank())


# ====== Research Question #1 ====== ----

# # LexTALE ----
# 
# # Initialize the data frame
# df_rq1_lt <- my_data[,c("random_id","lt_l1_total","lt_l2_total",
#                         "dsct_total","sst_total","fhat_total")]
# df_rq1_lt <- df_rq1_lt %>% rename(
#     L1 = lt_l1_total,
#     L2 = lt_l2_total)
# df_rq1_lt <- na.omit(df_rq1_lt)
# 
# # Pivot to long format
# df_rq1_lt <- tidyr::pivot_longer(df_rq1_lt,
#                                  cols = c("L1", "L2"),
#                                  names_to = "task_language",
#                                  values_to = "score")
# 
# # Create the LME model
# mod_rq1_lt <- lmer(score ~ task_language + dsct_total + 
#                    sst_total + fhat_total + (1|random_id),
#                    data=df_rq1_lt, REML = FALSE)
# summary(mod_rq1_lt)
# mod_rq1_lt_null <- lmer(score ~ dsct_total + sst_total + 
#                         fhat_total + (1|random_id),
#                         data=df_rq1_lt, REML = FALSE)
# summary(mod_rq1_lt_null)
# anova(mod_rq1_lt, mod_rq1_lt_null)
# 
# # Assumption 1: Linearity
# 
# checks_df_rq1_lt <- df_rq1_lt %>% 
#   mutate(fitted = fitted(mod_rq1_lt),
#          resid = resid(mod_rq1_lt))
# 
# ggplot(checks_df_rq1_lt,
#        aes(x = fitted, y = resid, color=task_language)) +
#   geom_point(alpha = 0.1) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   theme_bw() +
#   guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
#   labs(
#     x = "Fitted values",
#     y = "Residuals",
#     title = "Residuals vs Fitted Values"
#   )
# 
# # Assumption 2: Absence of colinearity
# 
# # Using an equivalent fixed-effects regression
# lm_rq1_lt <- lm(score ~ task_language + dsct_total + 
#                    sst_total + fhat_total,
#                    data=df_rq1_lt)
# # Rule of thumb: < 3 is not worrisome
# vif(lm_rq1_lt)
# 
# # Assumption 3: Homoskedasticity
# # Same as the ggplot from Assumption 1 with "abs" added
# ggplot(checks_df_rq1_lt,
#        aes(x = fitted, y = abs(resid), color=task_language)) +
#   geom_point(alpha = 0.1) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   theme_bw() +
#   guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
#   labs(
#     x = "Fitted values",
#     y = "Residuals",
#     title = "Residuals vs Fitted Values"
#   )
# 
# # Assumption 4: Normality of residuals
# # Not generally important given large sample
# hist(resid(mod_rq1_lt))
# qqnorm(resid(mod_rq1_lt))
# 
# # Assumption 5: Absence of influential data points
# 
# # Values > 1 are likely influential, values < 0.5 are generally safe
# 
# cooks_d_rq1_lt <- cooks.distance(lm_rq1_lt)
# checks_df_rq1_lt$cooks_d <- cooks_d_rq1_lt
# checks_df_rq1_lt <- checks_df_rq1_lt %>% 
#   mutate(obs_id = row_number())
# 
# ggplot(checks_df_rq1_lt, aes(x = obs_id, y = cooks_d)) +
#   geom_segment(aes(xend = obs_id, yend = 0)) +
#   geom_hline(
#     yintercept = 0.5,
#     linetype = "dashed",
#     color = "red") +
#   theme_bw() +
#   labs(x = "Observation index", y = "Cook's distance")
# 
# # Extract Estimated Marginal Means for figure
# emm_rq1_lt <- emmeans(mod_rq1_lt, ~ task_language)
# emm_df_rq1_lt <- as.data.frame(emm_rq1_lt)
# 
# # Plot estimated marginal means
# ggplot(emm_df_rq1_lt, aes(x = task_language, y = emmean)) +
#   geom_point(size = 3) +
#   geom_errorbar(
#     aes(ymin = lower.CL, ymax = upper.CL),
#     width = 0.1) +
#   theme_bw() +
#   labs(x = "Task language", y = "Estimated mean score")
# 
# 
# # Camel and Cactus ----
# 
# # Initialize the data frame
# df_rq1_ca_c <- my_data[,c("random_id","ca_c_l1_total","ca_c_l2_total",
#                         "dsct_total","sst_total","fhat_total")]
# df_rq1_ca_c <- df_rq1_ca_c %>% rename(
#     L1 = ca_c_l1_total,
#     L2 = ca_c_l2_total)
# df_rq1_ca_c <- na.omit(df_rq1_ca_c)
# 
# # Pivot to long format
# df_rq1_ca_c <- tidyr::pivot_longer(df_rq1_ca_c,
#                                  cols = c("L1", "L2"),
#                                  names_to = "task_language",
#                                  values_to = "score")
# 
# # Create the LME model
# mod_rq1_ca_c <- lmer(score ~ task_language + dsct_total + 
#                    sst_total + fhat_total + (1|random_id),
#                    data=df_rq1_ca_c, REML = FALSE)
# summary(mod_rq1_ca_c)
# mod_rq1_ca_c_null <- lmer(score ~ dsct_total + sst_total + 
#                         fhat_total + (1|random_id),
#                         data=df_rq1_ca_c, REML = FALSE)
# summary(mod_rq1_ca_c_null)
# anova(mod_rq1_ca_c, mod_rq1_ca_c_null)
# 
# # Assumption 1: Linearity
# 
# checks_df_rq1_ca_c <- df_rq1_ca_c %>% 
#   mutate(fitted = fitted(mod_rq1_ca_c),
#          resid = resid(mod_rq1_ca_c))
# 
# ggplot(checks_df_rq1_ca_c,
#        aes(x = fitted, y = resid, color=task_language)) +
#   geom_point(alpha = 0.1) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   theme_bw() +
#   guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
#   labs(
#     x = "Fitted values",
#     y = "Residuals",
#     title = "Residuals vs Fitted Values"
#   )
# 
# # Assumption 2: Absence of colinearity
# 
# # Using an equivalent fixed-effects regression
# lm_rq1_ca_c <- lm(score ~ task_language + dsct_total + 
#                    sst_total + fhat_total,
#                    data=df_rq1_ca_c)
# # Rule of thumb: < 3 is not worrisome
# vif(lm_rq1_ca_c)
# 
# # Assumption 3: Homoskedasticity
# # Same as the ggplot from Assumption 1 with "abs" added
# ggplot(checks_df_rq1_ca_c,
#        aes(x = fitted, y = abs(resid), color=task_language)) +
#   geom_point(alpha = 0.1) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   theme_bw() +
#   guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
#   labs(
#     x = "Fitted values",
#     y = "Residuals",
#     title = "Residuals vs Fitted Values"
#   )
# 
# # Assumption 4: Normality of residuals
# # Not generally important given large sample
# hist(resid(mod_rq1_ca_c))
# qqnorm(resid(mod_rq1_ca_c))
# 
# # Assumption 5: Absence of influential data points
# 
# # Values > 1 are likely influential, values < 0.5 are generally safe
# 
# cooks_d_rq1_ca_c <- cooks.distance(lm_rq1_ca_c)
# checks_df_rq1_ca_c$cooks_d <- cooks_d_rq1_ca_c
# checks_df_rq1_ca_c <- checks_df_rq1_ca_c %>% 
#   mutate(obs_id = row_number())
# 
# ggplot(checks_df_rq1_ca_c, aes(x = obs_id, y = cooks_d)) +
#   geom_segment(aes(xend = obs_id, yend = 0)) +
#   geom_hline(
#     yintercept = 0.5,
#     linetype = "dashed",
#     color = "red") +
#   theme_bw() +
#   labs(x = "Observation index", y = "Cook's distance")
# 
# 
# # Syntactic Modification ----
# 
# # Initialize the data frame
# df_rq1_sm <- my_data[,c("random_id","sm_l1_total","sm_l2_total",
#                         "dsct_total","sst_total","fhat_total")]
# df_rq1_sm <- df_rq1_sm %>% rename(
#     L1 = sm_l1_total,
#     L2 = sm_l2_total)
# df_rq1_sm <- na.omit(df_rq1_sm)
# 
# # Pivot to long format
# df_rq1_sm <- tidyr::pivot_longer(df_rq1_sm,
#                                  cols = c("L1", "L2"),
#                                  names_to = "task_language",
#                                  values_to = "score")
# 
# # Create the LME model
# mod_rq1_sm <- lmer(score ~ task_language + dsct_total + 
#                    sst_total + fhat_total + (1|random_id),
#                    data=df_rq1_sm, REML = FALSE)
# summary(mod_rq1_sm)
# mod_rq1_sm_null <- lmer(score ~ dsct_total + sst_total + 
#                         fhat_total + (1|random_id),
#                         data=df_rq1_sm, REML = FALSE)
# summary(mod_rq1_sm_null)
# anova(mod_rq1_sm, mod_rq1_sm_null)
# 
# # Assumption 1: Linearity
# 
# checks_df_rq1_sm <- df_rq1_sm %>% 
#   mutate(fitted = fitted(mod_rq1_sm),
#          resid = resid(mod_rq1_sm))
# 
# ggplot(checks_df_rq1_sm,
#        aes(x = fitted, y = resid, color=task_language)) +
#   geom_point(alpha = 0.1) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   theme_bw() +
#   guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
#   labs(
#     x = "Fitted values",
#     y = "Residuals",
#     title = "Residuals vs Fitted Values"
#   )
# 
# # Assumption 2: Absence of colinearity
# 
# # Using an equivalent fixed-effects regression
# lm_rq1_sm <- lm(score ~ task_language + dsct_total + 
#                    sst_total + fhat_total,
#                    data=df_rq1_sm)
# # Rule of thumb: < 3 is not worrisome
# vif(lm_rq1_sm)
# 
# # Assumption 3: Homoskedasticity
# # Same as the ggplot from Assumption 1 with "abs" added
# ggplot(checks_df_rq1_sm,
#        aes(x = fitted, y = abs(resid), color=task_language)) +
#   geom_point(alpha = 0.1) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   theme_bw() +
#   guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
#   labs(
#     x = "Fitted values",
#     y = "Residuals",
#     title = "Residuals vs Fitted Values"
#   )
# 
# # Assumption 4: Normality of residuals
# # Not generally important given large sample
# hist(resid(mod_rq1_sm))
# qqnorm(resid(mod_rq1_sm))
# 
# # Assumption 5: Absence of influential data points
# 
# # Values > 1 are likely influential, values < 0.5 are generally safe
# 
# cooks_d_rq1_sm <- cooks.distance(lm_rq1_sm)
# checks_df_rq1_sm$cooks_d <- cooks_d_rq1_sm
# checks_df_rq1_sm <- checks_df_rq1_sm %>% 
#   mutate(obs_id = row_number())
# 
# ggplot(checks_df_rq1_sm, aes(x = obs_id, y = cooks_d)) +
#   geom_segment(aes(xend = obs_id, yend = 0)) +
#   geom_hline(
#     yintercept = 0.5,
#     linetype = "dashed",
#     color = "red") +
#   theme_bw() +
#   labs(x = "Observation index", y = "Cook's distance")
# 
# 
# # Pragmatics ----
# 
# # Initialize the data frame
# df_rq1_p <- my_data[,c("random_id","p_l1_total","p_l2_total",
#                         "dsct_total","sst_total","fhat_total")]
# df_rq1_p <- df_rq1_p %>% rename(
#     L1 = p_l1_total,
#     L2 = p_l2_total)
# df_rq1_p <- na.omit(df_rq1_p)
# 
# # Pivot to long format
# df_rq1_p <- tidyr::pivot_longer(df_rq1_p,
#                                  cols = c("L1", "L2"),
#                                  names_to = "task_language",
#                                  values_to = "score")
# 
# # Create the LME model
# mod_rq1_p <- lmer(score ~ task_language + dsct_total + 
#                    sst_total + fhat_total + (1|random_id),
#                    data=df_rq1_p, REML = FALSE)
# summary(mod_rq1_p)
# mod_rq1_p_null <- lmer(score ~ dsct_total + sst_total + 
#                         fhat_total + (1|random_id),
#                         data=df_rq1_p, REML = FALSE)
# summary(mod_rq1_p_null)
# anova(mod_rq1_p, mod_rq1_p_null)
# 
# # Assumption 1: Linearity
# 
# checks_df_rq1_p <- df_rq1_p %>% 
#   mutate(fitted = fitted(mod_rq1_p),
#          resid = resid(mod_rq1_p))
# 
# ggplot(checks_df_rq1_p,
#        aes(x = fitted, y = resid, color=task_language)) +
#   geom_point(alpha = 0.1) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   theme_bw() +
#   guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
#   labs(
#     x = "Fitted values",
#     y = "Residuals",
#     title = "Residuals vs Fitted Values"
#   )
# 
# # Assumption 2: Absence of collinearity
# 
# # Using an equivalent fixed-effects regression
# lm_rq1_p <- lm(score ~ task_language + dsct_total + 
#                    sst_total + fhat_total,
#                    data=df_rq1_p)
# # Rule of thumb: < 3 is not worrisome
# vif(lm_rq1_p)
# 
# # Assumption 3: Homoskedasticity
# # Same as the ggplot from Assumption 1 with "abs" added
# ggplot(checks_df_rq1_p,
#        aes(x = fitted, y = abs(resid), color=task_language)) +
#   geom_point(alpha = 0.1) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   theme_bw() +
#   guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
#   labs(
#     x = "Fitted values",
#     y = "Residuals",
#     title = "Residuals vs Fitted Values"
#   )
# 
# # Assumption 4: Normality of residuals
# # Not generally important given large sample
# hist(resid(mod_rq1_p))
# qqnorm(resid(mod_rq1_p))
# 
# # Assumption 5: Absence of influential data points
# 
# # Values > 1 are likely influential, values < 0.5 are generally safe
# 
# cooks_d_rq1_p <- cooks.distance(lm_rq1_p)
# checks_df_rq1_p$cooks_d <- cooks_d_rq1_p
# checks_df_rq1_p <- checks_df_rq1_p %>% 
#   mutate(obs_id = row_number())
# 
# ggplot(checks_df_rq1_p, aes(x = obs_id, y = cooks_d)) +
#   geom_segment(aes(xend = obs_id, yend = 0)) +
#   geom_hline(
#     yintercept = 0.5,
#     linetype = "dashed",
#     color = "red") +
#   theme_bw() +
#   labs(x = "Observation index", y = "Cook's distance")


# # Combined model of four tasks ----
# 
# # Initialize the data frame
# df_rq1_all <- my_data[,c("random_id",
#                         "lt_l1_total",  "lt_l2_total",
#                         "ca_c_l1_total","ca_c_l2_total",
#                         "sm_l1_total",  "sm_l2_total",
#                         "p_l1_total",   "p_l2_total",
#                         "dsct_total","sst_total","fhat_total")]
# df_rq1_all <- na.omit(df_rq1_all)
# 
# # Change scores to percentages
# df_rq1_all$lt_l1_total   <- df_rq1_all$lt_l1_total   * 100/30
# df_rq1_all$ca_c_l1_total <- df_rq1_all$ca_c_l1_total * 100/16
# df_rq1_all$sm_l1_total   <- df_rq1_all$sm_l1_total   * 100/16
# df_rq1_all$p_l1_total    <- df_rq1_all$p_l1_total    * 100/18
# df_rq1_all$lt_l2_total   <- df_rq1_all$lt_l2_total   * 100/30
# df_rq1_all$ca_c_l2_total <- df_rq1_all$ca_c_l2_total * 100/16
# df_rq1_all$sm_l2_total   <- df_rq1_all$sm_l2_total   * 100/16
# df_rq1_all$p_l2_total    <- df_rq1_all$p_l2_total    * 100/18
# 
# # Pivot to long format
# df_rq1_all <- tidyr::pivot_longer(df_rq1_all,
#                                  cols = c(
#                                    "lt_l1_total",  "lt_l2_total",
#                                    "ca_c_l1_total","ca_c_l2_total",
#                                    "sm_l1_total",  "sm_l2_total",
#                                    "p_l1_total",   "p_l2_total",),
#                                  names_to = "task_and_language",
#                                  values_to = "score")
# # Add task and language cols
# df_rq1_all$task     <- rep(c(rep("LexTALE", 2), rep("C&C", 2),
#                              rep("S.M.", 2), rep("Pragmatics", 2)),
#                            nrow(df_rq1_all) / 8)
# df_rq1_all$lang <- rep(c("L1","L2"),
#                        nrow(df_rq1_all) / 2)
# 
# # Create the LME model
# mod_rq1_all <- lmer(score ~ task * lang + dsct_total + 
#                    sst_total + fhat_total + (1|random_id),
#                    data=df_rq1_all, REML = FALSE)
# summary(mod_rq1_all)
# mod_rq1_all_null <- lmer(score ~ task + dsct_total + sst_total + 
#                         fhat_total + (1|random_id),
#                         data=df_rq1_all, REML = FALSE)
# summary(mod_rq1_all_null)
# anova(mod_rq1_all, mod_rq1_all_null)
# 
# 
# # Assumption 1: Linearity
# 
# checks_df_rq1_all <- df_rq1_all %>% 
#   mutate(fitted = fitted(mod_rq1_all),
#          resid = resid(mod_rq1_all))
# 
# ggplot(checks_df_rq1_all,
#        aes(x = fitted, y = resid, color=task)) +
#   geom_point(alpha = 0.1) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   theme_bw() +
#   guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
#   labs(
#     x = "Fitted values",
#     y = "Residuals",
#     title = "Residuals vs Fitted Values"
#   )
# 
# # Assumption 2: Absence of colinearity
# 
# # Using an equivalent fixed-effects regression
# lm_rq1_all <- lm(score ~ task * lang + dsct_total + 
#                    sst_total + fhat_total,
#                    data=df_rq1_all)
# # Rule of thumb: < 3 is not worrisome
# vif(lm_rq1_all, type="predictor")
# 
# # Assumption 3: Homoskedasticity
# # Same as the ggplot from Assumption 1 with "abs" added
# ggplot(checks_df_rq1_all,
#        aes(x = fitted, y = abs(resid), color=task)) +
#   geom_point(alpha = 0.1) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   theme_bw() +
#   guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
#   labs(
#     x = "Fitted values",
#     y = "Residuals",
#     title = "Residuals vs Fitted Values"
#   )
# 
# # Assumption 4: Normality of residuals
# # Not generally important given large sample
# hist(resid(mod_rq1_all))
# qqnorm(resid(mod_rq1_all))
# 
# # Assumption 5: Absence of influential data points
# 
# # Values > 1 are likely influential, values < 0.5 are generally safe
# 
# cooks_d_rq1_all <- cooks.distance(lm_rq1_all)
# checks_df_rq1_all$cooks_d <- cooks_d_rq1_all
# checks_df_rq1_all <- checks_df_rq1_all %>% 
#   mutate(obs_id = row_number())
# 
# ggplot(checks_df_rq1_all, aes(x = obs_id, y = cooks_d)) +
#   geom_segment(aes(xend = obs_id, yend = 0)) +
#   geom_hline(
#     yintercept = 0.5,
#     linetype = "dashed",
#     color = "red") +
#   theme_bw() +
#   labs(x = "Observation index", y = "Cook's distance")
# 
# # Visualize the findings
# 
# # Extract the Estimated Marginal Means as a df
# emm_rq1_all <- emmeans(mod_rq1_all, ~ lang | task,
#                        pbkrtest.limit = 8000,
#                        lmerTest.limit = 8000)
# emm_df_rq1_all <- as.data.frame(emm_rq1_all)
# # Reorder the df to have the tasks in the correct order
# emm_df_rq1_all$task <- factor(
#   emm_df_rq1_all$task,
#   levels = c("LexTALE", "C&C",
#              "S.M.", "Pragmatics"))
# # Create a plot of the EMMs
# ggplot(emm_df_rq1_all, aes(x = task, y = emmean, color = lang)) +
#   geom_point(size = 2) +
#   geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
#                 width = 0.15) +
#   labs(x="Task", y="Estimated Marginal Means (%)",
#        color="Language", title=element_blank()) +
#   theme_bw() +
#   theme(panel.grid.major.y = element_line(color = "grey80"),
#         panel.grid.major.x = element_blank()) +
#   scale_color_manual(name   = "Language",
#                      values = c("L1" = color_l1,
#                                 "L2" = color_l2)) +
#   coord_cartesian(ylim = c(49, 91))
# 
# # Get estimates and p-values for L1/L2 contrasts by task
# contrasts_rq1_all <- contrast(
#   emm_rq1_all, method="pairwise", by="task",
#   adjust="holm")
# summary(contrasts_rq1_all, infer=TRUE)
# 
# # Now get Cohen's d for each contrast
# cohen_d_rq1_all <- eff_size(emm_rq1_all,
#                             sigma = sigma(mod_rq1_all),
#                             edf = df.residual(mod_rq1_all),
#                             method = "pairwise")
# summary(cohen_d_rq1_all)
# 
# 
# # ====== Research Question #2 ====== ----
# 
# # Initialize the data frame ----
# 
# df_rq2_all <- my_data[,c("random_id",
#                          "lt_l1_total", "ca_c_l1_total",
#                          "sm_l1_total", "p_l1_total",
#                          "mss_b_total", "asq_total",
#                          "dsct_total",  "sst_total", "fhat_total")]
# df_rq2_all <- na.omit(df_rq2_all)
# 
# # Change scores to percentages
# df_rq2_all$lt_l1_total   <- df_rq2_all$lt_l1_total   * 100/30
# df_rq2_all$ca_c_l1_total <- df_rq2_all$ca_c_l1_total * 100/16
# df_rq2_all$sm_l1_total   <- df_rq2_all$sm_l1_total   * 100/16
# df_rq2_all$p_l1_total    <- df_rq2_all$p_l1_total    * 100/18
# 
# # Standardize MSS-B and ASQ scores
# df_rq2_all$mss_b_total <- as.numeric(scale(df_rq2_all$mss_b_total,
#                                            center=TRUE,
#                                            scale=TRUE))
# df_rq2_all$asq_total <- as.numeric(scale(df_rq2_all$asq_total,
#                                          center=TRUE,
#                                          scale=TRUE))
# 
# # Pivot to long format
# df_rq2_all <- tidyr::pivot_longer(df_rq2_all,
#                                   cols = c(
#                                     "lt_l1_total", "ca_c_l1_total",
#                                     "sm_l1_total", "p_l1_total"),
#                                   names_to = "task",
#                                   values_to = "score")
# 
# # LME model for MSS-B ----
# 
# # Note: There are additional models / model comparisons here that I did in meeting with Olessia
# # Fitting the model
# mod_rq2_mss_b <- lmer(score ~ task * mss_b_total + 
#                         dsct_total + sst_total + fhat_total + 
#                         (1|random_id),
#                    data=df_rq2_all, REML = FALSE)
# summary(mod_rq2_mss_b)
# # Simplified model
# mod_rq2_mss_b_sim <- lmer(score ~ task + mss_b_total + 
#                         dsct_total + sst_total + fhat_total + 
#                         (1|random_id),
#                    data=df_rq2_all, REML = FALSE)
# summary(mod_rq2_mss_b_sim)
# # Comparing to a null model
# mod_rq2_mss_b_null <- lmer(score ~ task + dsct_total + sst_total + 
#                         fhat_total + (1|random_id),
#                         data=df_rq2_all, REML = FALSE)
# summary(mod_rq2_mss_b_null)
# # Comparing to a null model
# mod_rq2_mss_b_null2 <- lmer(score ~ mss_b_total + dsct_total + sst_total + 
#                         fhat_total + (1|random_id),
#                         data=df_rq2_all, REML = FALSE)
# summary(mod_rq2_mss_b_null2)
# # Comparing to a null model
# mod_rq2_mss_b_null3 <- lmer(score ~ dsct_total + sst_total + 
#                         fhat_total + (1|random_id),
#                         data=df_rq2_all, REML = FALSE)
# summary(mod_rq2_mss_b_null3)
# anova(mod_rq2_mss_b, mod_rq2_mss_b_null)
# anova(mod_rq2_mss_b, mod_rq2_mss_b_sim)
# anova(mod_rq2_mss_b_sim, mod_rq2_mss_b_null)
# anova(mod_rq2_mss_b_sim, mod_rq2_mss_b_null2)
# anova(mod_rq2_mss_b_null2, mod_rq2_mss_b_null3)
# 
# 
# # mod_rq2_mss_b_alt <- lmer(score ~ task * mss_b_total + 
# #                         (1|random_id),
# #                    data=df_rq2_all, REML = FALSE)
# # summary(mod_rq2_mss_b_alt)
# 
# # Assumption 1: Linearity
# 
# checks_df_rq2_mss_b <- df_rq2_all %>% 
#   mutate(fitted = fitted(mod_rq2_mss_b),
#          resid = resid(mod_rq2_mss_b))
# 
# ggplot(checks_df_rq2_mss_b,
#        aes(x = fitted, y = resid, color=task)) +
#   geom_point(alpha = 0.05) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   theme_bw() +
#   guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
#   labs(
#     x = "Fitted values",
#     y = "Residuals",
#     title = "Residuals vs Fitted Values"
#   )
# 
# # Assumption 2: Absence of collinearity
# 
# # Using an equivalent fixed-effects regression
# lm_rq2_mss_b <- lm(score ~ task * mss_b_total + dsct_total + 
#                    sst_total + fhat_total,
#                    data=df_rq2_all)
# # Rule of thumb: < 3 is not worrisome
# vif(lm_rq2_mss_b, type="predictor")
# 
# # Assumption 3: Homoskedasticity
# # Same as the ggplot from Assumption 1 with "abs" added
# ggplot(checks_df_rq2_mss_b,
#        aes(x = fitted, y = abs(resid), color=task)) +
#   geom_point(alpha = 0.05) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   theme_bw() +
#   guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
#   labs(
#     x = "Fitted values",
#     y = "Residuals",
#     title = "Residuals vs Fitted Values"
#   )
# 
# # Assumption 4: Normality of residuals
# # Not generally important given large sample
# hist(resid(mod_rq2_mss_b))
# qqnorm(resid(mod_rq2_mss_b))
# 
# # Assumption 5: Absence of influential data points
# 
# # Values > 1 are likely influential, values < 0.5 are generally safe
# 
# cooks_d_rq2_mss_b <- cooks.distance(lm_rq2_mss_b)
# checks_df_rq2_mss_b$cooks_d <- cooks_d_rq2_mss_b
# checks_df_rq2_mss_b <- checks_df_rq2_mss_b %>% 
#   mutate(obs_id = row_number())
# 
# ggplot(checks_df_rq2_mss_b, aes(x = obs_id, y = cooks_d)) +
#   geom_segment(aes(xend = obs_id, yend = 0)) +
#   geom_hline(
#     yintercept = 0.5,
#     linetype = "dashed",
#     color = "red") +
#   theme_bw() +
#   labs(x = "Observation index", y = "Cook's distance")
# 
# # Obtain p-values for the model-adjusted slopes
# # And mss_b_total.trend is the beta coefficients, meaning:
# # 1 SD in MSS-B corresponds to beta SDs on that task
# slopes_rq2_mss_b <- emtrends(
#   mod_rq2_mss_b, ~ task, var = "mss_b_total",
#   pbkrtest.limit = 4000, lmerTest.limit = 4000)
# summary(slopes_rq2_mss_b, infer=c(TRUE, TRUE), adjust="holm")
# 
# # # Calculate partial R-squared per effect
# # partR2(mod_rq2_mss_b,
# #        partvars = c("mss_b_total", "task:mss_b_total"),
# #        data = df_rq2_all)
# 
# # To test which tasks were sig. diff. from each other in the effect of MSS-B on score
# pairs(slopes_rq2_mss_b)
# 
# 
# # LME model for ASQ ----
# 
# mod_rq2_asq <- lmer(score ~ task * asq_total + 
#                         dsct_total + sst_total + fhat_total + 
#                         (1|random_id),
#                    data=df_rq2_all, REML = FALSE)
# summary(mod_rq2_asq)
# mod_rq2_asq_null <- lmer(score ~ task + dsct_total + sst_total + 
#                         fhat_total + (1|random_id),
#                         data=df_rq2_all, REML = FALSE)
# summary(mod_rq2_asq_null)
# anova(mod_rq2_asq, mod_rq2_asq_null)
# 
# # Assumption 1: Linearity
# 
# checks_df_rq2_asq <- df_rq2_all %>% 
#   mutate(fitted = fitted(mod_rq2_asq),
#          resid = resid(mod_rq2_asq))
# 
# ggplot(checks_df_rq2_asq,
#        aes(x = fitted, y = resid, color=task)) +
#   geom_point(alpha = 0.05) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   theme_bw() +
#   guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
#   labs(
#     x = "Fitted values",
#     y = "Residuals",
#     title = "Residuals vs Fitted Values"
#   )
# 
# # Assumption 2: Absence of collinearity
# 
# # Using an equivalent fixed-effects regression
# lm_rq2_asq <- lm(score ~ task * asq_total + dsct_total + 
#                    sst_total + fhat_total,
#                    data=df_rq2_all)
# # Rule of thumb: < 3 is not worrisome
# vif(lm_rq2_asq, type="predictor")
# 
# # Assumption 3: Homoskedasticity
# # Same as the ggplot from Assumption 1 with "abs" added
# ggplot(checks_df_rq2_asq,
#        aes(x = fitted, y = abs(resid), color=task)) +
#   geom_point(alpha = 0.05) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   theme_bw() +
#   guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
#   labs(
#     x = "Fitted values",
#     y = "Residuals",
#     title = "Residuals vs Fitted Values"
#   )
# 
# # Assumption 4: Normality of residuals
# # Not generally important given large sample
# hist(resid(mod_rq2_asq))
# qqnorm(resid(mod_rq2_asq))
# 
# # Assumption 5: Absence of influential data points
# 
# # Values > 1 are likely influential, values < 0.5 are generally safe
# 
# cooks_d_rq2_asq <- cooks.distance(lm_rq2_asq)
# checks_df_rq2_asq$cooks_d <- cooks_d_rq2_asq
# checks_df_rq2_asq <- checks_df_rq2_asq %>% 
#   mutate(obs_id = row_number())
# 
# ggplot(checks_df_rq2_asq, aes(x = obs_id, y = cooks_d)) +
#   geom_segment(aes(xend = obs_id, yend = 0)) +
#   geom_hline(
#     yintercept = 0.5,
#     linetype = "dashed",
#     color = "red") +
#   theme_bw() +
#   labs(x = "Observation index", y = "Cook's distance")
# 
# # Obtain p-values for the model-adjusted slopes
# # And asq_total.trend is the beta coefficients, meaning:
# # 1 SD in ASQ corresponds to beta SDs on that task
# slopes_rq2_asq <- emtrends(
#   mod_rq2_asq, ~ task, var = "asq_total",
#   pbkrtest.limit = 4000, lmerTest.limit = 4000)
# summary(slopes_rq2_asq, infer=c(TRUE, TRUE), adjust="holm")
# 
# 
# # Visualizing model-implied slopes for figure ----
# 
# # Obtain the predicted values for MSS-B model
# pred_rq2_mss_b <- ggpredict(mod_rq2_mss_b,
#                             terms = c("mss_b_total", "task"))
# pred_df_rq2_mss_b <- as.data.frame(pred_rq2_mss_b)
# # Reorder the df to have the tasks in the correct order
# pred_df_rq2_mss_b$group <- factor(
#   pred_df_rq2_mss_b$group,
#   levels = c("lt_l1_total", "ca_c_l1_total",
#              "sm_l1_total", "p_l1_total"))
# 
# # Obtain the predicted values for ASQ model
# pred_rq2_asq <- ggpredict(mod_rq2_asq,
#                             terms = c("asq_total", "task"))
# pred_df_rq2_asq <- as.data.frame(pred_rq2_asq)
# # Reorder the df to have the tasks in the correct order
# pred_df_rq2_asq$group <- factor(
#   pred_df_rq2_asq$group,
#   levels = c("lt_l1_total", "ca_c_l1_total",
#              "sm_l1_total", "p_l1_total"))
# 
# # Create plot of model-implied slopes for MSS-B model
# rq2_mss_b_plot <- ggplot(
#   pred_df_rq2_mss_b,
#   aes(x = x, y = predicted, color = group)) +
#   geom_line(linewidth = 1) +
#   geom_ribbon(aes(ymin = conf.low,
#                   ymax = conf.high,
#                   fill = group),
#               alpha = 0.2, color = NA) +
#   theme_bw() +
#   theme(panel.grid.major.y = element_line(color = "grey80"),
#         panel.grid.major.x = element_blank(),
#         legend.position = "none") +
#   labs(x = "MSS-B (Z-Scores)",
#        y = "Predicted Score (%)",
#        color = "Task", fill  = "Task") +
#   coord_cartesian(ylim = c(49, 91)) +
#   scale_color_manual(
#   values = c(
#     "ca_c_l1_total" = color_ca_c,
#     "lt_l1_total"   = color_lt,
#     "p_l1_total"    = color_p,
#     "sm_l1_total"   = color_sm),
#   labels = c(
#     "ca_c_l1_total" = "C&C",
#     "lt_l1_total"   = "LexTALE",
#     "p_l1_total"    = "Pragmatics",
#     "sm_l1_total"   = "S.M.")) +
#   scale_fill_manual(
#   values = c(
#     "ca_c_l1_total" = color_ca_c,
#     "lt_l1_total"   = color_lt,
#     "p_l1_total"    = color_p,
#     "sm_l1_total"   = color_sm),
#   labels = c(
#     "ca_c_l1_total" = "C&C",
#     "lt_l1_total"   = "LexTALE",
#     "p_l1_total"    = "Pragmatics",
#     "sm_l1_total"   = "S.M."))
#   
#   
# # Create plot of model-implied slopes for ASQ model
# rq2_asq_plot <- ggplot(
#   pred_df_rq2_asq,
#   aes(x = x, y = predicted, color = group)) +
#   geom_line(linewidth = 1) +
#   geom_ribbon(aes(ymin = conf.low,
#                   ymax = conf.high,
#                   fill = group),
#               alpha = 0.2, color = NA) +
#   theme_bw() +
#   theme(panel.grid.major.y = element_line(color = "grey80"),
#         panel.grid.major.x = element_blank()) +
#   labs(x = "AQ (Z-Scores)",
#        y = NULL,
#        color = "Task", fill  = "Task")  +
#   coord_cartesian(ylim = c(49, 91)) +
#   scale_color_manual(
#   values = c(
#     "ca_c_l1_total" = color_ca_c,
#     "lt_l1_total"   = color_lt,
#     "p_l1_total"    = color_p,
#     "sm_l1_total"   = color_sm),
#   labels = c(
#     "ca_c_l1_total" = "C&C",
#     "lt_l1_total"   = "LexTALE",
#     "p_l1_total"    = "Pragmatics",
#     "sm_l1_total"   = "S.M.")) +
#   scale_fill_manual(
#   values = c(
#     "ca_c_l1_total" = color_ca_c,
#     "lt_l1_total"   = color_lt,
#     "p_l1_total"    = color_p,
#     "sm_l1_total"   = color_sm),
#   labels = c(
#     "ca_c_l1_total" = "C&C",
#     "lt_l1_total"   = "LexTALE",
#     "p_l1_total"    = "Pragmatics",
#     "sm_l1_total"   = "S.M."))
# 
# # Display the two plots together as a figure
# rq2_mss_b_plot + rq2_asq_plot
# 
# 
# # ====== Research Question #3 ====== ----
# 
# # Initialize the data frame ----
# 
# df_rq3_all <- my_data[,c("random_id",
#                          "lt_l1_l2_difference",
#                          "ca_c_l1_l2_difference",
#                          "sm_l1_l2_difference",
#                          "p_l1_l2_difference",
#                          "mss_b_total", "asq_total",
#                          "prof_l1_l2_diff",
#                          "dsct_total",  "sst_total", "fhat_total")]
# # Remove any participants with missing data
# df_rq3_all <- na.omit(df_rq3_all)
# 
# # Change scores to percentages
# df_rq3_all$lt_l1_l2_difference   <- df_rq3_all$lt_l1_l2_difference   * 100/30
# df_rq3_all$ca_c_l1_l2_difference <- df_rq3_all$ca_c_l1_l2_difference * 100/16
# df_rq3_all$sm_l1_l2_difference   <- df_rq3_all$sm_l1_l2_difference   * 100/16
# df_rq3_all$p_l1_l2_difference    <- df_rq3_all$p_l1_l2_difference    * 100/18
# 
# # Standardize MSS-B and ASQ scores
# df_rq3_all$mss_b_total <- as.numeric(scale(df_rq3_all$mss_b_total,
#                                            center=TRUE,
#                                            scale=TRUE))
# df_rq3_all$asq_total   <- as.numeric(scale(df_rq3_all$asq_total,
#                                            center=TRUE,
#                                            scale=TRUE))
# 
# # Pivot to long format
# df_rq3_all <- tidyr::pivot_longer(df_rq3_all,
#                                   cols = c(
#                                     "lt_l1_l2_difference",
#                                     "ca_c_l1_l2_difference",
#                                     "sm_l1_l2_difference",
#                                     "p_l1_l2_difference"),
#                                   names_to = "task",
#                                   values_to = "score")
# 
# 
# # LME model for MSS-B ----
# 
# mod_rq3_mss_b <- lmer(score ~ task * mss_b_total + 
#                         prof_l1_l2_diff +
#                         dsct_total + sst_total + fhat_total + 
#                         (1|random_id),
#                       data=df_rq3_all, REML = FALSE)
# summary(mod_rq3_mss_b)
# mod_rq3_mss_b_null <- lmer(score ~ task + prof_l1_l2_diff + 
#                              dsct_total + sst_total + 
#                              fhat_total + (1|random_id),
#                            data=df_rq3_all, REML = FALSE)
# summary(mod_rq3_mss_b_null)
# anova(mod_rq3_mss_b, mod_rq3_mss_b_null)
# 
# # Assumption 1: Linearity
# 
# checks_df_rq3_mss_b <- df_rq3_all %>% 
#   mutate(fitted = fitted(mod_rq3_mss_b),
#          resid = resid(mod_rq3_mss_b))
# 
# ggplot(checks_df_rq3_mss_b,
#        aes(x = fitted, y = resid, color=task)) +
#   geom_point(alpha = 0.05) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   theme_bw() +
#   guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
#   labs(
#     x = "Fitted values",
#     y = "Residuals",
#     title = "Residuals vs Fitted Values"
#   )
# 
# # Assumption 2: Absence of collinearity
# 
# # Using an equivalent fixed-effects regression
# lm_rq3_mss_b <- lm(score ~ task * mss_b_total + 
#                    prof_l1_l2_diff + dsct_total + 
#                    sst_total + fhat_total,
#                    data=df_rq3_all)
# # Rule of thumb: < 3 is not worrisome
# vif(lm_rq3_mss_b, type="predictor")
# 
# # Assumption 3: Homoskedasticity
# # Same as the ggplot from Assumption 1 with "abs" added
# ggplot(checks_df_rq3_mss_b,
#        aes(x = fitted, y = abs(resid), color=task)) +
#   geom_point(alpha = 0.05) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   theme_bw() +
#   guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
#   labs(
#     x = "Fitted values",
#     y = "Residuals",
#     title = "Residuals vs Fitted Values"
#   )
# 
# # Assumption 4: Normality of residuals
# # Not generally important given large sample
# hist(resid(mod_rq3_mss_b))
# qqnorm(resid(mod_rq3_mss_b))
# 
# # Assumption 5: Absence of influential data points
# 
# # Values > 1 are likely influential, values < 0.5 are generally safe
# 
# cooks_d_rq3_mss_b <- cooks.distance(lm_rq3_mss_b)
# checks_df_rq3_mss_b$cooks_d <- cooks_d_rq3_mss_b
# checks_df_rq3_mss_b <- checks_df_rq3_mss_b %>% 
#   mutate(obs_id = row_number())
# 
# ggplot(checks_df_rq3_mss_b, aes(x = obs_id, y = cooks_d)) +
#   geom_segment(aes(xend = obs_id, yend = 0)) +
#   geom_hline(
#     yintercept = 0.5,
#     linetype = "dashed",
#     color = "red") +
#   theme_bw() +
#   labs(x = "Observation index", y = "Cook's distance")
# 
# # Obtain p-values for the model-adjusted slopes
# # And mss_b_total.trend is the beta coefficients, meaning:
# # 1 SD in MSS-B corresponds to beta SDs on that task
# slopes_rq3_mss_b <- emtrends(
#   mod_rq3_mss_b, ~ task, var = "mss_b_total",
#   pbkrtest.limit = 4000, lmerTest.limit = 4000)
# summary(slopes_rq3_mss_b, infer=c(TRUE, TRUE), adjust="holm")
# 
# 
# # LME model for ASQ ----
# 
# mod_rq3_asq <- lmer(score ~ task * asq_total + 
#                         prof_l1_l2_diff +
#                         dsct_total + sst_total + fhat_total + 
#                         (1|random_id),
#                       data=df_rq3_all, REML = FALSE)
# summary(mod_rq3_asq)
# mod_rq3_asq_null <- lmer(score ~ task + prof_l1_l2_diff + 
#                              dsct_total + sst_total + 
#                              fhat_total + (1|random_id),
#                            data=df_rq3_all, REML = FALSE)
# summary(mod_rq3_asq_null)
# anova(mod_rq3_asq, mod_rq3_asq_null)
# 
# # Assumption 1: Linearity
# 
# checks_df_rq3_asq <- df_rq3_all %>% 
#   mutate(fitted = fitted(mod_rq3_asq),
#          resid = resid(mod_rq3_asq))
# 
# ggplot(checks_df_rq3_asq,
#        aes(x = fitted, y = resid, color=task)) +
#   geom_point(alpha = 0.05) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   theme_bw() +
#   guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
#   labs(
#     x = "Fitted values",
#     y = "Residuals",
#     title = "Residuals vs Fitted Values"
#   )
# 
# # Assumption 2: Absence of collinearity
# 
# # Using an equivalent fixed-effects regression
# lm_rq3_asq <- lm(score ~ task * asq_total + 
#                    prof_l1_l2_diff + dsct_total + 
#                    sst_total + fhat_total,
#                    data=df_rq3_all)
# # Rule of thumb: < 3 is not worrisome
# vif(lm_rq3_asq, type="predictor")
# 
# # Assumption 3: Homoskedasticity
# # Same as the ggplot from Assumption 1 with "abs" added
# ggplot(checks_df_rq3_asq,
#        aes(x = fitted, y = abs(resid), color=task)) +
#   geom_point(alpha = 0.05) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   theme_bw() +
#   guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
#   labs(
#     x = "Fitted values",
#     y = "Residuals",
#     title = "Residuals vs Fitted Values"
#   )
# 
# # Assumption 4: Normality of residuals
# # Not generally important given large sample
# hist(resid(mod_rq3_asq))
# qqnorm(resid(mod_rq3_asq))
# 
# # Assumption 5: Absence of influential data points
# 
# # Values > 1 are likely influential, values < 0.5 are generally safe
# 
# cooks_d_rq3_asq <- cooks.distance(lm_rq3_asq)
# checks_df_rq3_asq$cooks_d <- cooks_d_rq3_asq
# checks_df_rq3_asq <- checks_df_rq3_asq %>% 
#   mutate(obs_id = row_number())
# 
# ggplot(checks_df_rq3_asq, aes(x = obs_id, y = cooks_d)) +
#   geom_segment(aes(xend = obs_id, yend = 0)) +
#   geom_hline(
#     yintercept = 0.5,
#     linetype = "dashed",
#     color = "red") +
#   theme_bw() +
#   labs(x = "Observation index", y = "Cook's distance")
# 
# # Obtain p-values for the model-adjusted slopes
# # And asq_total.trend is the beta coefficients, meaning:
# # 1 SD in ASQ corresponds to beta SDs on that task
# slopes_rq3_asq <- emtrends(
#   mod_rq3_asq, ~ task, var = "asq_total",
#   pbkrtest.limit = 4000, lmerTest.limit = 4000)
# summary(slopes_rq3_asq, infer=c(TRUE, TRUE), adjust="holm")
# 
# 
# # Visualizing model-implied slopes for figure ----
# 
# # Obtain the predicted values for MSS-B model
# pred_rq3_mss_b <- ggpredict(mod_rq3_mss_b,
#                             terms = c("mss_b_total", "task"))
# pred_df_rq3_mss_b <- as.data.frame(pred_rq3_mss_b)
# # Reorder the df to have the tasks in the correct order
# pred_df_rq3_mss_b$group <- factor(
#   pred_df_rq3_mss_b$group,
#   levels = c("lt_l1_l2_difference", "ca_c_l1_l2_difference",
#              "sm_l1_l2_difference", "p_l1_l2_difference"))
# 
# # Obtain the predicted values for ASQ model
# pred_rq3_asq <- ggpredict(mod_rq3_asq,
#                             terms = c("asq_total", "task"))
# pred_df_rq3_asq <- as.data.frame(pred_rq3_asq)
# # Reorder the df to have the tasks in the correct order
# pred_df_rq3_asq$group <- factor(
#   pred_df_rq3_asq$group,
#   levels = c("lt_l1_l2_difference", "ca_c_l1_l2_difference",
#              "sm_l1_l2_difference", "p_l1_l2_difference"))
# 
# # Create plot of model-implied slopes for MSS-B model
# rq3_mss_b_plot <- ggplot(
#   pred_df_rq3_mss_b,
#   aes(x = x, y = predicted, color = group)) +
#   geom_line(linewidth = 1) +
#   geom_ribbon(aes(ymin = conf.low,
#                   ymax = conf.high,
#                   fill = group),
#               alpha = 0.2, color = NA) +
#   theme_bw() +
#   theme(panel.grid.major.y = element_line(color = "grey80"),
#         panel.grid.major.x = element_blank(),
#         legend.position = "none") +
#   labs(x = "MSS-B (Z-Scores)",
#        y = "Predicted L1/L2 Diff. (%)",
#        color = "Task", fill  = "Task") +
#   coord_cartesian(ylim = c(-1, 31)) +
#   scale_color_manual(
#   values = c(
#     "ca_c_l1_l2_difference" = color_ca_c,
#     "lt_l1_l2_difference"   = color_lt,
#     "p_l1_l2_difference"    = color_p,
#     "sm_l1_l2_difference"   = color_sm),
#   labels = c(
#     "ca_c_l1_l2_difference" = "C&C",
#     "lt_l1_l2_difference"   = "LexTALE",
#     "p_l1_l2_difference"    = "Pragmatics",
#     "sm_l1_l2_difference"   = "S.M.")) +
#   scale_fill_manual(
#   values = c(
#     "ca_c_l1_l2_difference" = color_ca_c,
#     "lt_l1_l2_difference"   = color_lt,
#     "p_l1_l2_difference"    = color_p,
#     "sm_l1_l2_difference"   = color_sm),
#   labels = c(
#     "ca_c_l1_l2_difference" = "C&C",
#     "lt_l1_l2_difference"   = "LexTALE",
#     "p_l1_l2_difference"    = "Pragmatics",
#     "sm_l1_l2_difference"   = "S.M."))
# 
# # Create plot of model-implied slopes for ASQ model
# rq3_asq_plot <- ggplot(
#   pred_df_rq3_asq,
#   aes(x = x, y = predicted, color = group)) +
#   geom_line(linewidth = 1) +
#   geom_ribbon(aes(ymin = conf.low,
#                   ymax = conf.high,
#                   fill = group),
#               alpha = 0.2, color = NA) +
#   theme_bw() +
#   theme(panel.grid.major.y = element_line(color = "grey80"),
#         panel.grid.major.x = element_blank()) +
#   labs(x = "AQ (Z-Scores)",
#        y = NULL,
#        color = "Task", fill  = "Task")  +
#   coord_cartesian(ylim = c(-1, 31)) +
#   scale_color_manual(
#   values = c(
#     "ca_c_l1_l2_difference" = color_ca_c,
#     "lt_l1_l2_difference"   = color_lt,
#     "p_l1_l2_difference"    = color_p,
#     "sm_l1_l2_difference"   = color_sm),
#   labels = c(
#     "ca_c_l1_l2_difference" = "C&C",
#     "lt_l1_l2_difference"   = "LexTALE",
#     "p_l1_l2_difference"    = "Pragmatics",
#     "sm_l1_l2_difference"   = "S.M.")) +
#   scale_fill_manual(
#   values = c(
#     "ca_c_l1_l2_difference" = color_ca_c,
#     "lt_l1_l2_difference"   = color_lt,
#     "p_l1_l2_difference"    = color_p,
#     "sm_l1_l2_difference"   = color_sm),
#   labels = c(
#     "ca_c_l1_l2_difference" = "C&C",
#     "lt_l1_l2_difference"   = "LexTALE",
#     "p_l1_l2_difference"    = "Pragmatics",
#     "sm_l1_l2_difference"   = "S.M."))
# 
# # Display the two plots together as a figure
# rq3_mss_b_plot + rq3_asq_plot


# ====== Redoing analyses with combined LME models ====== ----

# Initialize the data frame ----

df_lmes <- my_data[,c("random_id",
                      "lt_l1_weighted", "lt_l2_weighted",
                      "ca_c_l1_total", "ca_c_l2_total",
                      "sm_l1_total", "sm_l2_total",
                      "p_l1_total", "p_l2_total",
                      "mss_b_total", "asq_total", "prof_l1_l2_diff",
                      "dsct_total",  "sst_total", "fhat_total")]
df_lmes <- na.omit(df_lmes)

# Change scores to percentages (except LexTALE - already adjusted during weighting)
df_lmes$ca_c_l1_total <- df_lmes$ca_c_l1_total * 100/16
df_lmes$sm_l1_total   <- df_lmes$sm_l1_total   * 100/16
df_lmes$p_l1_total    <- df_lmes$p_l1_total    * 100/18
df_lmes$ca_c_l2_total <- df_lmes$ca_c_l2_total * 100/16
df_lmes$sm_l2_total   <- df_lmes$sm_l2_total   * 100/16
df_lmes$p_l2_total    <- df_lmes$p_l2_total    * 100/18

# Standardize MSS-B and ASQ scores
df_lmes$mss_b_total <- as.numeric(scale(df_lmes$mss_b_total,
                                           center=TRUE,
                                           scale=TRUE))
df_lmes$asq_total <- as.numeric(scale(df_lmes$asq_total,
                                         center=TRUE,
                                         scale=TRUE))

# Pivot to long format
df_lmes <- tidyr::pivot_longer(df_lmes,
                               cols = c(
                                  "lt_l1_weighted", "lt_l2_weighted",
                                  "ca_c_l1_total", "ca_c_l2_total",
                                  "sm_l1_total", "sm_l2_total",
                                  "p_l1_total", "p_l2_total"),
                               names_to = "task_and_lang",
                               values_to = "score")

# Add task and language cols
df_lmes$task     <- rep(c(rep("LexTALE", 2), rep("C&C", 2),
                          rep("S.M.", 2), rep("Pragmatics", 2)),
                        nrow(df_lmes) / 8)
df_lmes$lang <- rep(c("L1","L2"),
                    nrow(df_lmes) / 2)


# LME model for MSS-B ----

# Create the model itself
mod_mss_b <- lmer(score ~ task*lang*mss_b_total + 
                        prof_l1_l2_diff +
                        dsct_total + sst_total + fhat_total + 
                        (1|random_id),
                      data=df_lmes, REML = FALSE)
summary(mod_mss_b)

# Test sig. of 3-way interaction
mod_mss_b_1 <- lmer(score ~ task*lang + task*mss_b_total + lang*mss_b_total +
                        prof_l1_l2_diff +
                        dsct_total + sst_total + fhat_total + 
                        (1|random_id),
                      data=df_lmes, REML = FALSE)
summary(mod_mss_b_1)
anova(mod_mss_b, mod_mss_b_1)

# Three-way interaction was significant so test simple interactions
# This is for RQ on L1/L2 differences
slopes_mss_b_1 <- emtrends(mod_mss_b, ~ lang | task,
                         var = "mss_b_total",
                         pbkrtest.limit = 8000,
                         lmerTest.limit = 8000)
# First level contrasts (L1 - L2)
inter_mss_b_1 <- contrast(slopes_mss_b_1, interaction = "pairwise")
summary(inter_mss_b_1, adjust="none", by=NULL)
# summary(inter_mss_b_1, adjust="holm", by=NULL)
# Planned comparisons of simple interactions
pairs(inter_mss_b_1, adjust="holm", by = NULL)

# Now testing a different set of simple interactions
# This is for RQ on symptoms and task performance
slopes_mss_b_2 <- emtrends(mod_mss_b, ~ task | lang,
                         var = "mss_b_total",
                         pbkrtest.limit = 8000,
                         lmerTest.limit = 8000)
test(slopes_mss_b_2, adjust="none")
# test(slopes_mss_b_2, adjust="holm")
# First level contrasts (differences between tasks)
inter_mss_b_2 <- contrast(slopes_mss_b_2, interaction = "pairwise")
summary(inter_mss_b_2, adjust="holm", by=NULL)


# Assumption checks for MSS-B model ----

# Assumption 1: Linearity

checks_df_mss_b <- df_lmes %>% 
  mutate(fitted = fitted(mod_mss_b),
         resid = resid(mod_mss_b))

ggplot(checks_df_mss_b,
       aes(x = fitted, y = resid, color=task)) +
  geom_point(alpha = 0.05) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  labs(
    x = "Fitted values",
    y = "Residuals",
    title = "Residuals vs Fitted Values"
  )

# Assumption 2: Absence of collinearity

# Using an equivalent fixed-effects regression
lm_mss_b <- lm(score ~ task*mss_b_total*lang + 
                   prof_l1_l2_diff + dsct_total + 
                   sst_total + fhat_total,
                   data=df_lmes)
# Rule of thumb: < 3 is not worrisome
vif(lm_mss_b, type="predictor")

# Assumption 3: Homoskedasticity
# Same as the ggplot from Assumption 1 with "abs" added
ggplot(checks_df_mss_b,
       aes(x = fitted, y = abs(resid), color=task)) +
  geom_point(alpha = 0.05) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  labs(
    x = "Fitted values",
    y = "Residuals",
    title = "Residuals vs Fitted Values"
  )

# Assumption 4: Normality of residuals
# Not generally important given large sample
hist(resid(mod_mss_b))
qqnorm(resid(mod_mss_b))

# Assumption 5: Absence of influential data points

# Values > 1 are likely influential, values < 0.5 are generally safe

cooks_d_mss_b <- cooks.distance(lm_mss_b)
checks_df_mss_b$cooks_d <- cooks_d_mss_b
checks_df_mss_b <- checks_df_mss_b %>% 
  mutate(obs_id = row_number())

ggplot(checks_df_mss_b, aes(x = obs_id, y = cooks_d)) +
  geom_segment(aes(xend = obs_id, yend = 0)) +
  geom_hline(
    yintercept = 0.5,
    linetype = "dashed",
    color = "red") +
  theme_bw() +
  labs(x = "Observation index", y = "Cook's distance")


# LME model for ASQ ----

# Create the model itself
mod_asq <- lmer(score ~ task*lang*asq_total + 
                        prof_l1_l2_diff +
                        dsct_total + sst_total + fhat_total + 
                        (1|random_id),
                      data=df_lmes, REML = FALSE)
summary(mod_asq)

# Test sig. of 3-way interaction
mod_asq_1 <- lmer(score ~ task*lang + task*asq_total + lang*asq_total +
                        prof_l1_l2_diff +
                        dsct_total + sst_total + fhat_total + 
                        (1|random_id),
                      data=df_lmes, REML = FALSE)
summary(mod_asq_1)
anova(mod_asq, mod_asq_1)

# Not significant so we move on and test the two-way interactions
# Test sig. of lang:asq_total interaction
mod_asq_2 <- lmer(score ~ task*lang + task*asq_total +
                        prof_l1_l2_diff +
                        dsct_total + sst_total + fhat_total +
                        (1|random_id),
                      data=df_lmes, REML = FALSE)
anova(mod_asq_1, mod_asq_2)
# Test sig. of task:asq_total interaction
mod_asq_3 <- lmer(score ~ task*lang + lang*asq_total +
                        prof_l1_l2_diff +
                        dsct_total + sst_total + fhat_total +
                        (1|random_id),
                      data=df_lmes, REML = FALSE)
anova(mod_asq_1, mod_asq_3)
# Test sig. of task:lang interaction
mod_asq_4 <- lmer(score ~ task*asq_total + lang*asq_total +
                        prof_l1_l2_diff +
                        dsct_total + sst_total + fhat_total +
                        (1|random_id),
                      data=df_lmes, REML = FALSE)
anova(mod_asq_1, mod_asq_4)

# Test simple interactions
# This is for RQ on symptoms and task performance
slopes_asq <- emtrends(mod_asq, ~ task,
                         var = "asq_total",
                         pbkrtest.limit = 8000,
                         lmerTest.limit = 8000)
test(slopes_asq, adjust="none")
# test(slopes_asq, adjust="holm")
# First level contrasts (differences between tasks)
inter_asq <- contrast(slopes_asq, interaction = "pairwise")
summary(inter_asq, adjust="holm", by=NULL)


# Assumption checks for ASQ model ----

# Assumption 1: Linearity

checks_df_asq <- df_lmes %>% 
  mutate(fitted = fitted(mod_asq),
         resid = resid(mod_asq))

ggplot(checks_df_asq,
       aes(x = fitted, y = resid, color=task)) +
  geom_point(alpha = 0.05) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  labs(
    x = "Fitted values",
    y = "Residuals",
    title = "Residuals vs Fitted Values"
  )

# Assumption 2: Absence of collinearity

# Using an equivalent fixed-effects regression
lm_asq <- lm(score ~ task*asq_total*lang + 
                   prof_l1_l2_diff + dsct_total + 
                   sst_total + fhat_total,
                   data=df_lmes)
# Rule of thumb: < 3 is not worrisome
vif(lm_asq, type="predictor")

# Assumption 3: Homoskedasticity
# Same as the ggplot from Assumption 1 with "abs" added
ggplot(checks_df_asq,
       aes(x = fitted, y = abs(resid), color=task)) +
  geom_point(alpha = 0.05) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  labs(
    x = "Fitted values",
    y = "Residuals",
    title = "Residuals vs Fitted Values"
  )

# Assumption 4: Normality of residuals
# Not generally important given large sample
hist(resid(mod_asq))
qqnorm(resid(mod_asq))

# Assumption 5: Absence of influential data points

# Values > 1 are likely influential, values < 0.5 are generally safe

cooks_d_asq <- cooks.distance(lm_asq)
checks_df_asq$cooks_d <- cooks_d_asq
checks_df_asq <- checks_df_asq %>% 
  mutate(obs_id = row_number())

ggplot(checks_df_asq, aes(x = obs_id, y = cooks_d)) +
  geom_segment(aes(xend = obs_id, yend = 0)) +
  geom_hline(
    yintercept = 0.5,
    linetype = "dashed",
    color = "red") +
  theme_bw() +
  labs(x = "Observation index", y = "Cook's distance")


# Plots for performance by task by language ----

# Prepare df for MSS-B plot
mss_b_seq <- with(
  df_lmes,
  seq(min(mss_b_total, na.rm = TRUE),
      max(mss_b_total, na.rm = TRUE),
      length.out = 50)
)
mss_b_emm_grid <- emmeans(
  mod_mss_b,
  ~ lang * mss_b_total * task,
  at = list(mss_b_total = mss_b_seq),
  pbkrtest.limit = 8000, lmerTest.limit = 8000
)
mss_b_emm_df <- as.data.frame(mss_b_emm_grid)

# Order the tasks correctly
mss_b_emm_df$task <- factor(
  mss_b_emm_df$task,
  levels = c("LexTALE", "C&C", "S.M.", "Pragmatics")
)

# Create MSS-B plot
mss_b_langs_plot <- ggplot(
  mss_b_emm_df,
  aes(
    x = mss_b_total,
    y = emmean,
    linetype = lang,
    color = task,
    group = interaction(task, lang)
  )
) +
  geom_line(linewidth = 1) +
  labs(
    x = "Schizotypy (MSS-B Z-Score)",
    y = "Predicted Score (%)",
    linetype = "Language",
    color = "Task"
  ) +
  theme_bw(base_size=12) +
  # theme(
  #   legend.position = "right",
  #   legend.box = "vertical"
  # ) +
  theme(legend.position = "none") +     # Because the legend is on the other plot
  coord_cartesian(ylim = c(35, 90)) +
  scale_color_manual(
  values = c(
    "C&C" = color_ca_c,
    "LexTALE"   = color_lt,
    "Pragmatics"    = color_p,
    "S.M."   = color_sm)
  )

# Prepare df for ASQ plot
asq_seq <- with(
  df_lmes,
  seq(min(asq_total, na.rm = TRUE),
      max(asq_total, na.rm = TRUE),
      length.out = 50)
)
asq_emm_grid <- emmeans(
  mod_asq,
  ~ lang * asq_total * task,
  at = list(asq_total = asq_seq),
  pbkrtest.limit = 8000, lmerTest.limit = 8000
)
asq_emm_df <- as.data.frame(asq_emm_grid)

# Order the tasks correctly
asq_emm_df$task <- factor(
  asq_emm_df$task,
  levels = c("LexTALE", "C&C", "S.M.", "Pragmatics")
)

# Create ASQ plot
asq_langs_plot <- ggplot(
  asq_emm_df,
  aes(
    x = asq_total,
    y = emmean,
    linetype = lang,
    color = task,
    group = interaction(task, lang)
  )
) +
  geom_line(linewidth = 1) +
  labs(
    x = "Autistic Traits (AQ Z-Score)",
    y = NULL,     # On the other plot
    linetype = "Language",
    color = "Task"
  ) +
  theme_bw(base_size=12) +
  theme(
    legend.position = "right",
    legend.box = "vertical" 
  ) +
  coord_cartesian(ylim = c(35, 90)) +
  scale_color_manual(
  values = c(
    "C&C"        = color_ca_c,
    "LexTALE"    = color_lt,
    "Pragmatics" = color_p,
     "S.M."      = color_sm)
  )

# Joint plot
mss_b_langs_plot + asq_langs_plot


# Plot of L1/L2 differences ----

# Compute the data for L1/L2 differences (with MSS-B)
mss_b_diff <- contrast(
  mss_b_emm_grid,
  method = "pairwise",   # L1 − L2
  by = c("task", "mss_b_total"),
  adjust = "none"        # no multiplicity adjustment for plotting
)
mss_b_diff_df <- as.data.frame(
  summary(mss_b_diff, infer=TRUE))

# Order the tasks correctly
mss_b_diff_df$task <- factor(
  mss_b_diff_df$task,
  levels = c("LexTALE", "C&C", "S.M.", "Pragmatics")
)

# Create the MSS-B L1/L2 differences plot
mss_b_l1_l2_plot <- ggplot(mss_b_diff_df, aes(
  x = mss_b_total, y = estimate,
  color = task, fill = task)) +
  geom_hline(yintercept = 0, linewidth = 0.6) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = lower.CL, ymax = upper.CL, fill = task),
    alpha = 0.12,
    color = NA, show.legend = FALSE) +
  labs(
    x = "Schizotypy (MSS-B Z-Score)",
    y = "Predicted Difference (L1 − L2)",
    color = "Task") +
  theme_bw(base_size = 12) +
  # theme(legend.position = "right") +
  theme(legend.position = "none") +     # Because the legend is on the other plot
  coord_cartesian(ylim = c(-1, 31)) +
  scale_color_manual(
  values = c(
    "C&C"        = color_ca_c,
    "LexTALE"    = color_lt,
    "Pragmatics" = color_p,
     "S.M."      = color_sm)) + 
  scale_fill_manual(
    values = c(
      "C&C" = color_ca_c,
      "LexTALE" = color_lt,
      "Pragmatics" = color_p,
      "S.M." = color_sm)) +
    guides(fill = "none")

# Compute the data for L1/L2 differences (with ASQ)
asq_diff <- contrast(
  asq_emm_grid,
  method = "pairwise",   # L1 − L2
  by = c("task", "asq_total"),
  adjust = "none"        # no multiplicity adjustment for plotting
)
asq_diff_df <- as.data.frame(
  summary(asq_diff, infer=TRUE))

# Order the tasks correctly
asq_diff_df$task <- factor(
  asq_diff_df$task,
  levels = c("LexTALE", "C&C", "S.M.", "Pragmatics")
)

# Create the ASQ L1/L2 differences plot
asq_l1_l2_plot <- ggplot(asq_diff_df, aes(
  x = asq_total, y = estimate,
  color = task, fill = task)) +
  geom_hline(yintercept = 0, linewidth = 0.6) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = lower.CL, ymax = upper.CL, fill = task),
    alpha = 0.12,
    color = NA, show.legend = FALSE) +
  labs(
    x = "Autistic Traits (AQ Z-Score)",
    # y = "Predicted Difference (L1 − L2)",
    y = NULL,
    color = "Task") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right") +
  coord_cartesian(ylim = c(-1, 31)) +
  scale_color_manual(
  values = c(
    "C&C"        = color_ca_c,
    "LexTALE"    = color_lt,
    "Pragmatics" = color_p,
     "S.M."      = color_sm)) + 
  scale_fill_manual(
    values = c(
      "C&C" = color_ca_c,
      "LexTALE" = color_lt,
      "Pragmatics" = color_p,
      "S.M." = color_sm)) +
    guides(fill = "none")

# Create a joint plot
mss_b_l1_l2_plot + asq_l1_l2_plot

stop()


# ====== For backpocket slides ====== ----

# # More histograms ----
#
# # Histogram using my function
# create_histogram(my_df=my_data, my_var=age,
#                  my_xlab="Age", my_ylab="Count",
#                  my_color='mediumpurple1'
#                  x_min=18, x_max=30,
#                  my_title="Participant Ages")
#
# 
# # Switch to long format (L1 prof and L2 prof separate)
# profic_data <- my_data[,c("prof_l1_total","prof_l2_total")]
# colnames(profic_data) <- c("L1","L2")
# profic_data_long <- tidyr::pivot_longer(profic_data,
#                                         cols = c("L1","L2"),
#                                         names_to = "Language",
#                                         values_to = "score")
# # Histogram of L1 and L2 proficiency scores
# # Ended up editing this one in paint to shorten the last L1 bar (and mark jump on y-axis)
# ggplot(data=profic_data_long,
#        aes(x=score, fill=Language, colour=Language, rm.na = TRUE)) +
#        geom_histogram(binwidth=2, alpha=0.3, position="identity") +
#        scale_y_continuous(breaks=seq(0,360,by=28), "") +
#        xlab("Self-Rated Proficiency") + #ylab("") + ylim(0,100) +
#        #scale_fill_discrete(name = "Language", labels = c("L1", "L2")) +
#        #ggtitle("Participant Second-Language Proficiency") +
#        theme_bw() + theme(plot.title=element_text(hjust=0.5),
#                           axis.ticks.y=element_blank(),
#                           axis.text.y=element_blank(),
#                           #panel.grid.major.y = element_blank(),
#                           #panel.grid.major.x = element_blank(),
#                           panel.grid.minor = element_blank()
#                           )
# 
# # Histograms by L1 linguistic task scores
# ggplot(data=my_data, aes(x=lt_l1_total, rm.na = TRUE)) +
#        geom_histogram(fill='blue', color='lightgrey') +
#        xlab("LexTale score") + ylab("Number of participants") +
#        ggtitle("LexTale scores (L1)") +
#        theme_bw() + theme(plot.title=element_text(hjust=0.5))
# ggplot(data=my_data, aes(x=ca_c_l1_total, rm.na = TRUE)) +
#        geom_histogram(fill='blue', color='lightgrey') +
#        xlab("Camel and Cactus task score") + ylab("Number of participants") +
#        ggtitle("Camel and Cactus scores (L1)") +
#        theme_bw() + theme(plot.title=element_text(hjust=0.5))
# ggplot(data=my_data, aes(x=sm_l1_total, rm.na = TRUE)) +
#        geom_histogram(fill='blue', color='lightgrey') +
#        xlab("Syntactic modification task score") + ylab("Number of participants") +
#        ggtitle("Syntactic modification scores (L1)") +
#        theme_bw() + theme(plot.title=element_text(hjust=0.5))
# ggplot(data=my_data, aes(x=p_l1_total, rm.na = TRUE)) +
#        geom_histogram(fill='blue', color='lightgrey') +
#        xlab("Pragmatics task scores") + ylab("Number of participants") +
#        ggtitle("Pragmatics scores (L1)") +
#        theme_bw() + theme(plot.title=element_text(hjust=0.5))
# 
# # Histograms for scores of MSS-B and subscales
# ggplot(data=my_data, aes(x=mss_b_total, rm.na = TRUE)) +
#        geom_histogram(fill='blue', color='lightgrey', binwidth=1) +
#        xlab("MSS-B score") + ylab("Number of participants") +
#        ggtitle("Distribution of Schizotypy Scores") +
#        #facet_wrap("gender") + # labeller=as_labeller(my_facet_labels)) +
#        theme_bw() + theme(plot.title=element_text(hjust=0.5))
# # Faceted by gender
# gender_labels_facets <- c("0" = "Female", "1" = "Male")
# ggplot(data=my_data[my_data$gender %in% c(0,1),],
#        aes(x=mss_b_total, rm.na = TRUE)) +
#        geom_histogram(aes(y=after_stat(density)), fill='blue',
#                       color='lightgrey', binwidth=1) +
#        xlab("MSS-B Score") + ylab("Proportion") +
#        ggtitle("Distribution of Schizotypy Scores by Gender") +
#        facet_wrap("gender", labeller=as_labeller(gender_labels_facets)) +
#        theme_bw() + theme(plot.title=element_text(hjust=0.5))
#
# ggplot(data=my_data, aes(x=mss_b_pos, rm.na = TRUE)) +
#        geom_histogram(fill='blue', color='lightgrey', binwidth=1) +
#        xlab("MSS-B (positive subscale) score") + ylab("Number of participants") +
#        ggtitle("Distribution of MSS-B (positive) scores") +
#        facet_wrap("gender") +
#        theme_bw() + theme(plot.title=element_text(hjust=0.5))
# ggplot(data=my_data, aes(x=mss_b_neg, rm.na = TRUE)) +
#        geom_histogram(fill='blue', color='lightgrey', binwidth=1) +
#        xlab("MSS-B (negative subscale) score") + ylab("Number of participants") +
#        ggtitle("Distribution of MSS-B (negative) scores") +
#        facet_wrap("gender") +
#        theme_bw() + theme(plot.title=element_text(hjust=0.5))
# ggplot(data=my_data, aes(x=mss_b_dis, rm.na = TRUE)) +
#        geom_histogram(fill='blue', color='lightgrey', binwidth=1) +
#        xlab("MSS-B (disorganized subscale) score") + ylab("Number of participants") +
#        ggtitle("Distribution of MSS-B (disorganized) scores") +
#        facet_wrap("gender") +
#        theme_bw() + theme(plot.title=element_text(hjust=0.5))
# 
# # Histograms for scores of ASQ and subscales (overall and facet by gender)
# 
# gender_labels_facets <- c("0" = "Female", "1" = "Male")
# 
# ggplot(data=my_data, aes(x=asq_total, rm.na = TRUE)) +
#        geom_histogram(fill='blue', color='lightgrey', binwidth=1) +
#        xlab("ASQ score") + ylab("Number of participants") +
#        ggtitle("Distribution of autistic traits scores") +
#        theme_bw() + theme(plot.title=element_text(hjust=0.5))
# ggplot(data=my_data[my_data$gender %in% c(0,1),],
#        aes(x=asq_total, rm.na = TRUE)) +
#        geom_histogram(aes(y=after_stat(density)), fill='blue',
#                       color='lightgrey', binwidth=1) +
#        xlab("ASQ score") + ylab("Density") +
#        ggtitle("Distribution of autistic traits scores") +
#        facet_wrap("gender", labeller=as_labeller(gender_labels_facets)) +
#        theme_bw() + theme(plot.title=element_text(hjust=0.5))
# 
# ggplot(data=my_data, aes(x=asq_comm, rm.na = TRUE)) +
#        geom_histogram(fill='blue', color='lightgrey', binwidth=1) +
#        xlab("ASQ (communication subscale) score") + ylab("Number of participants") +
#        ggtitle("Distribution of ASQ (communication) scores") +
#        theme_bw() + theme(plot.title=element_text(hjust=0.5))
# ggplot(data=my_data[my_data$gender %in% c(0,1),],
#        aes(x=asq_comm, rm.na = TRUE)) +
#        geom_histogram(aes(y=after_stat(density)), fill='blue',
#                       color='lightgrey', binwidth=1) +
#        xlab("ASQ (communication subscale) score") + ylab("Density") +
#        ggtitle("Distribution of ASQ (communication) scores") +
#        facet_wrap("gender", labeller=as_labeller(gender_labels_facets)) +
#        theme_bw() + theme(plot.title=element_text(hjust=0.5))
# 
# ggplot(data=my_data, aes(x=asq_social, rm.na = TRUE)) +
#        geom_histogram(fill='blue', color='lightgrey', binwidth=1) +
#        xlab("ASQ (social subscale) score") + ylab("Number of participants") +
#        ggtitle("Distribution of ASQ (social) scores") +
#        theme_bw() + theme(plot.title=element_text(hjust=0.5))
# ggplot(data=my_data[my_data$gender %in% c(0,1),],
#        aes(x=asq_social, rm.na = TRUE)) +
#        geom_histogram(aes(y=after_stat(density)), fill='blue',
#                       color='lightgrey', binwidth=1) +
#        xlab("ASQ (social subscale) score") + ylab("Density") +
#        ggtitle("Distribution of ASQ (social) scores") +
#        facet_wrap("gender", labeller=as_labeller(gender_labels_facets)) +
#        theme_bw() + theme(plot.title=element_text(hjust=0.5))
# 
# ggplot(data=my_data, aes(x=asq_imagine, rm.na = TRUE)) +
#        geom_histogram(fill='blue', color='lightgrey', binwidth=1) +
#        xlab("ASQ (imagination subscale) score") + ylab("Number of participants") +
#        ggtitle("Distribution of ASQ (imagination) scores") +
#        theme_bw() + theme(plot.title=element_text(hjust=0.5))
# ggplot(data=my_data[my_data$gender %in% c(0,1),],
#        aes(x=asq_imagine, rm.na = TRUE)) +
#        geom_histogram(aes(y=after_stat(density)), fill='blue',
#                       color='lightgrey', binwidth=1) +
#        xlab("ASQ (imagination subscale) score") + ylab("Density") +
#        ggtitle("Distribution of ASQ (imagination) scores") +
#        facet_wrap("gender", labeller=as_labeller(gender_labels_facets)) +
#        theme_bw() + theme(plot.title=element_text(hjust=0.5))
# 
# ggplot(data=my_data, aes(x=asq_detail, rm.na = TRUE)) +
#        geom_histogram(fill='blue', color='lightgrey', binwidth=1) +
#        xlab("ASQ (detail subscale) score") + ylab("Number of participants") +
#        ggtitle("Distribution of ASQ (detail) scores") +
#        theme_bw() + theme(plot.title=element_text(hjust=0.5))
# ggplot(data=my_data[my_data$gender %in% c(0,1),],
#        aes(x=asq_detail, rm.na = TRUE)) +
#        geom_histogram(aes(y=after_stat(density)), fill='blue',
#                       color='lightgrey', binwidth=1) +
#        xlab("ASQ (detail subscale) score") + ylab("Density") +
#        ggtitle("Distribution of ASQ (detail) scores") +
#        facet_wrap("gender", labeller=as_labeller(gender_labels_facets)) +
#        theme_bw() + theme(plot.title=element_text(hjust=0.5))
# 
# ggplot(data=my_data, aes(x=asq_att_switch, rm.na = TRUE)) +
#        geom_histogram(fill='blue', color='lightgrey', binwidth=1) +
#        xlab("ASQ (attention-switching subscale) score") + ylab("Number of participants") +
#        ggtitle("Distribution of ASQ (attention-switching) scores") +
#        theme_bw() + theme(plot.title=element_text(hjust=0.5))
# ggplot(data=my_data[my_data$gender %in% c(0,1),],
#        aes(x=asq_att_switch, rm.na = TRUE)) +
#        geom_histogram(aes(y=after_stat(density)), fill='blue',
#                       color='lightgrey', binwidth=1) +
#        xlab("ASQ (attention-switching subscale) score") + ylab("Density") +
#        ggtitle("Distribution of ASQ (attention-switching) scores") +
#        facet_wrap("gender", labeller=as_labeller(gender_labels_facets)) +
#        theme_bw() + theme(plot.title=element_text(hjust=0.5))


# ====== Correlational analyses ====== ----

# # Test and plot relationship between MSS-B and ASQ scores ----
# 
# # # Test correlation between ASQ and MSS-B scores
# # symptoms_cor <- cor.test(my_data$asq_total, my_data$mss_b_total)
# # symptoms_cor
# 
# # Create scatterplot of MSS-B vs AQ
# ggplot(data=my_data, aes(x=mss_b_total, y=asq_total)) +
#    geom_vline(xintercept = 0) + geom_vline(xintercept = 38) +
#    geom_hline(yintercept = 0) + geom_hline(yintercept = 50) +
#    geom_point(color='darkred', alpha=0.4) + geom_smooth(method='lm', formula=y~x, se=TRUE, color='black') +
#    xlab("MSS-B Score") + ylab("AQ Score") + ggtitle("Relationship Between Schizotypal and Autistic Symptoms") +
#    scale_x_continuous(limits = c(0, 38)) + scale_y_continuous(limits = c(0, 50)) +
#    theme_minimal() + theme(plot.title=element_text(hjust=0.5))
# 
# 
#
#
# grid.arrange(lt_l1_plot, ca_c_l1_plot, sm_l1_plot, p_l1_plot,
#              lt_l2_plot, ca_c_l2_plot, sm_l2_plot, p_l2_plot,
#              nrow = 2)


# # Correlograms for cognitive tasks ----
# 
# # Linguistic tasks in L1
# create_correlogram(my_df     = my_data,
#                    my_cols   = c("lt_l1_total", "ca_c_l1_total", "sm_l1_total", "p_l1_total"),
#                    my_title  = "Correlogram of Linguistic Tasks (L1 Only)",
#                    my_labels = c("LT", "C&C", "SM", "Prag.")
#                   )
# 
# # Linguistic tasks in L2
# create_correlogram(my_df     = my_data,
#                    my_cols   = c("lt_l2_total", "ca_c_l2_total", "sm_l2_total", "p_l2_total"),
#                    my_title  = "Correlogram of Linguistic Tasks (L2 Only)",
#                    my_labels = c("LT", "C&C", "SM", "Prag.")
#                   )
# 
# # Non-linguistic tasks
# create_correlogram(my_df     = my_data,
#                    my_cols   = c("dsct_total", "sst_total", "fhat_total"),
#                    my_title  = "Correlogram of Non-Linguistic Tasks",
#                    my_labels = c("DSC (PS)", "SS (WM)", "FHA (ToM)")
#                   )
# 
# # Linguistic tasks in L1 and non-linguistic tasks
# create_correlogram(my_df     = my_data,
#                    my_cols   = c("lt_l1_total", "ca_c_l1_total", "sm_l1_total", "p_l1_total",
#                                  "dsct_total", "sst_total", "fhat_total"),
#                    my_title  = "Correlogram of Cognitive Tasks (L1 Only)",
#                    my_labels = c("LT", "C&C", "SM", "Prag.", "DSC", "SS", "FHA")
#                   )
# 
# # Linguistic tasks in L2 and non-linguistic tasks
# create_correlogram(my_df     = my_data,
#                    my_cols   = c("lt_l2_total", "ca_c_l2_total", "sm_l2_total", "p_l2_total",
#                                  "dsct_total", "sst_total", "fhat_total"),
#                    my_title  = "Correlogram of Non-Linguistic and L2 Linguistic Tasks",
#                    my_labels = c("LT", "C&C", "SM", "Prag.", "DSC", "SS", "FHA")
#                   )
# 
# # L1/L2 differences on linguistic tasks, and non-linguistic tasks
# create_correlogram(my_df     = my_data,
#                    my_cols   = c("lt_l1_l2_difference", "ca_c_l1_l2_difference", "sm_l1_l2_difference", "p_l1_l2_difference",
#                                  "dsct_total", "sst_total", "fhat_total"),
#                    my_title  = "Correlogram of L1/L2 Diffs. and Non-Linguistic Tasks",
#                    my_labels = c("LT", "C&C", "SM", "Prag.", "DSC", "SS", "FHA")
#                   )
# 
# 
# # Correlograms for symptoms and cognitive tasks ----
# 
# # Correlations among symptoms
# create_correlogram(my_df     = my_data,
#                    my_cols   = c("mss_b_total", "mss_b_neg", "mss_b_pos", "mss_b_dis", "asq_total"),
#                    my_title  = "Correlogram of Symptom Scales and Sub-Scales",
#                    my_labels = c("MSS-B", "MSS-B (Neg.)", "MSS-B (Pos.)", "MSS-B (Dis.)", "ASQ")
#                   )
# 
# # Linguistic tasks in L1 and symptoms
# create_correlogram(my_df     = my_data,
#                    my_cols   = c("lt_l1_total", "ca_c_l1_total", "sm_l1_total", "p_l1_total",
#                                  "mss_b_total", "asq_total"),
#                    my_title  = "Correlogram of L1 Scores and Symptoms",
#                    my_labels = c("LT", "C&C", "SM", "Prag.",
#                                  "MSS-B", "ASQ")
#                   )
# 
# # Linguistic tasks in L1 and symptoms
# create_correlogram(my_df     = my_data,
#                    my_cols   = c("lt_l1_total", "ca_c_l1_total", "sm_l1_total", "p_l1_total",
#                                  "mss_b_total", "asq_total"),
#                    my_title  = "Correlogram of L1 Scores and Symptoms",
#                    my_labels = c("LT", "C&C", "SM", "Prag.",
#                                  "MSS-B", "ASQ")
#                   )
# 
# # L1 ling task scores and MSS-B subscales
# create_correlogram(my_df     = my_data,
#                    my_cols   = c("lt_l1_total", "ca_c_l1_total", "sm_l1_total", "p_l1_total",
#                                  "mss_b_total", "mss_b_pos", "mss_b_neg", "mss_b_dis"),
#                    my_title  = "L1 Tasks vs. MSS-B Subscales",
#                    my_labels = c("LT", "C&C", "SM", "Prag.",
#                                  "SZ-Total", "SZ-Pos", "SZ-Neg", "SZ-Dis")
#                   )
# 
# # L1 ling task scores and ASQ subscales
# create_correlogram(my_df     = my_data,
#                    my_cols   = c("lt_l1_total", "ca_c_l1_total", "sm_l1_total", "p_l1_total",
#                                  "asq_total", "asq_comm", "asq_social", "asq_imagine", "asq_att_switch", "asq_detail"),
#                    my_title  = "L1 Tasks vs. AQ Subscales",
#                    my_labels = c("LT", "C&C", "SM", "Prag.",
#                                  "AQ-Total", "AQ-Comm", "AQ-Social", "AQ-Imagine", "AQ-Att-Switch", "AQ-Detail")
#                   )
# 
# # Linguistic tasks in L2 and symptoms
# create_correlogram(my_df     = my_data,
#                    my_cols   = c("lt_l2_total", "ca_c_l2_total", "sm_l2_total", "p_l2_total",
#                                  "mss_b_total", "asq_total"),
#                    my_title  = "Correlogram of L2 Scores and Symptoms",
#                    my_labels = c("LT", "C&C", "SM", "Prag.",
#                                  "MSS-B", "ASQ")
#                   )
# 
# # L2 ling task scores and MSS-B subscales
# create_correlogram(my_df     = my_data,
#                    my_cols   = c("lt_l2_total", "ca_c_l2_total", "sm_l2_total", "p_l2_total",
#                                  "mss_b_total", "mss_b_pos", "mss_b_neg", "mss_b_dis"),
#                    my_title  = "L2 Tasks vs. MSS-B Subscales",
#                    my_labels = c("LT", "C&C", "SM", "Prag.",
#                                  "SZ-Total", "SZ-Pos", "SZ-Neg", "SZ-Dis")
#                   )
# 
# # L2 ling task scores and ASQ subscales
# create_correlogram(my_df     = my_data,
#                    my_cols   = c("lt_l2_total", "ca_c_l2_total", "sm_l2_total", "p_l2_total",
#                                  "asq_total", "asq_comm", "asq_social", "asq_imagine", "asq_att_switch", "asq_detail"),
#                    my_title  = "L2 Tasks vs. AQ Subscales",
#                    my_labels = c("LT", "C&C", "SM", "Prag.",
#                                  "AQ-Total", "AQ-Comm", "AQ-Social", "AQ-Imagine", "AQ-Att-Switch", "AQ-Detail")
#                   )
# 
# # L1/L2 differences and symptoms
# create_correlogram(my_df     = my_data,
#                    my_cols   = c("lt_l1_l2_difference", "ca_c_l1_l2_difference", "sm_l1_l2_difference", "p_l1_l2_difference",
#                                  "mss_b_total", "asq_total"),
#                    my_title  = "Correlogram of L1/L2 Differences and Symptoms",
#                    my_labels = c("LT Diff.", "C&C Diff.", "SM Diff.", "Prag. Diff.",
#                                  "MSS-B", "ASQ")
#                   )
# 
# # L1/L2 diffs and MSS-B subscales
# create_correlogram(my_df     = my_data,
#                    my_cols   = c("lt_l1_l2_difference", "ca_c_l1_l2_difference", "sm_l1_l2_difference", "p_l1_l2_difference",
#                                  "mss_b_total", "mss_b_pos", "mss_b_neg", "mss_b_dis"),
#                    my_title  = "L1/L2 Diffs. vs. MSS-B Subscales",
#                    my_labels = c("LT", "C&C", "SM", "Prag.",
#                                  "SZ-Total", "SZ-Pos", "SZ-Neg", "SZ-Dis")
#                   )
# 
# # L1/L2 ratios and MSS-B subscales
# create_correlogram(my_df     = my_data,
#                    my_cols   = c("lt_l1_l2_ratio", "ca_c_l1_l2_ratio", "sm_l1_l2_ratio", "p_l1_l2_ratio",
#                                  "mss_b_total", "mss_b_pos", "mss_b_neg", "mss_b_dis"),
#                    my_title  = "L1/L2 Ratios vs. MSS-B Subscales",
#                    my_labels = c("LT", "C&C", "SM", "Prag.",
#                                  "SZ-Total", "SZ-Pos", "SZ-Neg", "SZ-Dis")
#                   )
# 
# # L1/L2 z-score diffs. and MSS-B subscales
# create_correlogram(my_df     = my_data,
#                    my_cols   = c("lt_l1_l2_z_diff", "ca_c_l1_l2_z_diff", "sm_l1_l2_z_diff", "p_l1_l2_z_diff", "all_l1_l2_z_diff",
#                                  "mss_b_total", "mss_b_pos", "mss_b_neg", "mss_b_dis"),
#                    my_title  = "L1/L2 Z-Score Diffs. vs. MSS-B Subscales",
#                    my_labels = c("LT", "C&C", "SM", "Prag.", "All",
#                                  "SZ-Total", "SZ-Pos", "SZ-Neg", "SZ-Dis")
#                   )
# 
# # L1/L2 diffs and ASQ subscales
# create_correlogram(my_df     = my_data,
#                    my_cols   = c("lt_l1_l2_difference", "ca_c_l1_l2_difference", "sm_l1_l2_difference", "p_l1_l2_difference",
#                                  "asq_total", "asq_comm", "asq_social", "asq_imagine", "asq_att_switch", "asq_detail"),
#                    my_title  = "L1/L2 Diffs. vs. AQ Subscales",
#                    my_labels = c("LT", "C&C", "SM", "Prag.",
#                                  "AQ-Total", "AQ-Comm", "AQ-Social", "AQ-Imagine", "AQ-Att-Switch", "AQ-Detail")
#                   )
# 
# # Symptoms and non-linguistic tasks
# create_correlogram(my_df     = my_data,
#                    my_cols   = c("mss_b_total", "asq_total",
#                                  "dsct_total", "sst_total", "fhat_total"),
#                    my_title  = "Correlogram of Symptoms and Non-Linguistic Task Scores",
#                    my_labels = c("MSS-B", "ASQ",
#                                  "DSC", "SS", "FHA")
#                   )
# 
# 
# # Correlograms exploring the impact of L2 proficiency ----
# 
# # L2 proficiency and symptoms
# create_correlogram(my_df     = my_data,
#                    my_cols   = c("prof_l2_total", "prof_l2_speak", "prof_l2_underst", "prof_l2_write", "prof_l2_read",
#                                  "mss_b_total", "asq_total"),
#                    my_title  = "Correlogram of L2 Prof. and Symptoms",
#                    my_labels = c("L2 Prof.", "L2 Speak", "L2 Under.", "L2 Write", "L2 Read",
#                                  "MSS-B", "ASQ")
#                   )
# 
# # L1 task scores and L2 proficiency
# create_correlogram(my_df     = my_data,
#                    my_cols   = c("lt_l1_total", "ca_c_l1_total", "sm_l1_total", "p_l1_total",
#                                  "prof_l2_total", "prof_l2_speak", "prof_l2_underst", "prof_l2_write", "prof_l2_read"),
#                    my_title  = "Correlogram of L1 Tasks and L2 Prof.",
#                    my_labels = c("LT", "C&C", "SM", "Prag.",
#                                  "L2 Prof.", "L2 Speak", "L2 Under.", "L2 Write", "L2 Read")
#                   )
# 
# # L2 task scores and L2 proficiency
# create_correlogram(my_df     = my_data,
#                    my_cols   = c("lt_l2_total", "ca_c_l2_total", "sm_l2_total", "p_l2_total",
#                                  "prof_l2_total", "prof_l2_speak", "prof_l2_underst", "prof_l2_write", "prof_l2_read"),
#                    my_title  = "Correlogram of L2 Tasks and L2 Prof.",
#                    my_labels = c("LT", "C&C", "SM", "Prag.",
#                                  "L2 Prof.", "L2 Speak", "L2 Under.", "L2 Write", "L2 Read")
#                   )
# 
# # L1/L2 differences in task scores and L2 proficiency
# create_correlogram(my_df     = my_data,
#                    my_cols   = c("lt_l1_l2_difference", "ca_c_l1_l2_difference", "sm_l1_l2_difference", "p_l1_l2_difference",
#                                  "prof_l2_total", "prof_l2_speak", "prof_l2_underst", "prof_l2_write", "prof_l2_read"),
#                    my_title  = "Correlogram of Task L1/L2 Diffs. and L2 Prof.",
#                    my_labels = c("LT", "C&C", "SM", "Prag.",
#                                  "L2 Prof.", "L2 Speak", "L2 Under.", "L2 Write", "L2 Read")
#                   )
# 
# # Fr vars and symptoms (L1_En only)
# create_correlogram(my_df     = my_data[my_data$l1 == "English / Anglais",],
#                    my_cols   = c("fr_exposure", "fr_reading", "fr_speaking",
#                                  "fr_country_yrs", "fr_family_yrs", "fr_school_work_yrs",
#                                  "mss_b_total", "asq_total"),
#                    my_title  = "Correlogram of Fr. Vars. and Symptoms",
#                    my_labels = c("Expo.", "Reading", "Speaking",
#                                  "Country", "Family", "Sch/Work",
#                                  "MSS-B", "ASQ")
#                   )
# 
# # # L2 proficiency and non-linguistic tasks
# # create_correlogram(my_df     = my_data,
# #                    my_cols   = c("prof_l2_total", "prof_l2_speak", "prof_l2_underst", "prof_l2_write", "prof_l2_read",
# #                                  "dsct_total", "sst_total", "fhat_total"),
# #                    my_title  = "Correlogram of L2 Prof. and Non-Linguistic Tasks",
# #                    my_labels = c("L2 Prof.", "L2 Speak", "L2 Under.", "L2 Write", "L2 Read",
# #                                  "DSC", "SS", "FHA")
# #                   )
# 
# # Correlogram for symptom subscales ----
# 
# # L2 proficiency and symptoms
# create_correlogram(my_df     = my_data,
#                    my_cols   = c("mss_b_pos", "mss_b_neg", "mss_b_dis",
#                                  "asq_comm", "asq_social", "asq_imagine", "asq_att_switch", "asq_detail"),
#                    my_title  = "Correlogram of Symptom Subscales",
#                    my_labels = c("SZ-Pos", "SZ-Neg", "SZ-Dis",
#                                  "AQ-Comm", "AQ-Social", "AQ-Imagine", "AQ-Att-Switch", "AQ-Detail")
#                   )
#
#


# ====== Analyses with PSM ====== ----

# # Setup for PSM ----
# # Mainly relied on this: https://www.r-bloggers.com/2016/06/how-to-use-r-for-matching-samples-propensity-score/
# 
# # Adding columns for arbitrary high/low split in MSS-B and ASQ
# my_data$has_high_mss_b <- FALSE
# my_data[my_data$mss_b_total >= 13,]$has_high_mss_b <- TRUE
# my_data$has_high_asq <- FALSE
# my_data[my_data$asq_total >= 23,]$has_high_asq <- TRUE
# 
# # Creating data frame for matching (only includes participants who provided their age and gender (and only included male and female for purposes of matching)
# df_for_matching <- my_data[!is.na(my_data$age) & !is.na(my_data$gender) & !is.na(my_data$dsct_total) & !is.na(my_data$sst_total),]
# 
# # Analyses with PSM - high vs. low MSS-B ----
# # match_it_mss_b <- matchit(has_high_mss_b ~ age + gender + dsct_total + sst_total, data = df_for_matching, method="nearest", ratio=1)
# match_it_mss_b <- matchit(has_high_mss_b ~ age + gender, data = df_for_matching, method="nearest", ratio=1)
# summary(match_it_mss_b)
# df_matched_mss_b <- match.data(match_it_mss_b)[1:ncol(my_data)]
# 
# # # Create data frame of means by high vs. low MSS-B (L1 and L2)
# # mss_b_group    <- c(rep("high", 4), rep("low", 4), rep("high", 4), rep("low", 4))
# # Language       <- c(rep("L1", 8), rep("L2", 8))
# # task           <- rep(c("LT", "C&C", "SM", "P"), 4)
# # my_mss_b_means <- c(mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == TRUE,]$lt_l1_total,    na.rm=TRUE),
# #                     mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == TRUE,]$ca_c_l1_total,  na.rm=TRUE),
# #                     mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == TRUE,]$sm_l1_total,    na.rm=TRUE),
# #                     mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == TRUE,]$p_l1_total,     na.rm=TRUE),
# #                     mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == FALSE,]$lt_l1_total,   na.rm=TRUE),
# #                     mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == FALSE,]$ca_c_l1_total, na.rm=TRUE),
# #                     mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == FALSE,]$sm_l1_total,   na.rm=TRUE),
# #                     mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == FALSE,]$p_l1_total,    na.rm=TRUE),
# #                     mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == TRUE,]$lt_l2_total,    na.rm=TRUE),
# #                     mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == TRUE,]$ca_c_l2_total,  na.rm=TRUE),
# #                     mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == TRUE,]$sm_l2_total,    na.rm=TRUE),
# #                     mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == TRUE,]$p_l2_total,     na.rm=TRUE),
# #                     mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == FALSE,]$lt_l2_total,   na.rm=TRUE),
# #                     mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == FALSE,]$ca_c_l2_total, na.rm=TRUE),
# #                     mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == FALSE,]$sm_l2_total,   na.rm=TRUE),
# #                     mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == FALSE,]$p_l2_total,    na.rm=TRUE))
# # df_mss_b_bars <- data.frame(mss_b_group, Language, task, my_mss_b_means)
# # 
# # my_facet_labels <- c(`high`="High", `low`="Low")
# # 
# # # Make them ordered factors
# # df_mss_b_bars$task <- factor(df_mss_b_bars$task,
# #                              unique(df_mss_b_bars$task),
# #                              ordered=TRUE)
# # 
# # # Plot of high vs. low MSS-B groups
# # ggplot(df_mss_b_bars, aes(fill=Language, y=my_mss_b_means, x=task)) +
# #     geom_bar(position="dodge", stat="identity") +
# #     xlab("Linguistic Task") + ylab("Mean Task Score") +
# #     ggtitle("Task Scores by MSS-B Group") +
# #     facet_wrap("mss_b_group", labeller=as_labeller(my_facet_labels)) +
# #     theme_bw() + theme(plot.title=element_text(hjust=0.5)) +
# #     scale_fill_discrete(name = "Language", labels = c("L1", "L2")) +
# #     scale_fill_manual(values=c("darkred", "salmon"))
# 
# # Convert task scores to percentages
# df_matched_mss_b$lt_l1_l2_diff_perc   <- df_matched_mss_b$lt_l1_l2_difference   * 100/30
# df_matched_mss_b$ca_c_l1_l2_diff_perc <- df_matched_mss_b$ca_c_l1_l2_difference * 100/16
# df_matched_mss_b$sm_l1_l2_diff_perc   <- df_matched_mss_b$sm_l1_l2_difference   * 100/16
# df_matched_mss_b$p_l1_l2_diff_perc    <- df_matched_mss_b$p_l1_l2_difference    * 100/18
# 
# # Create data frame of means by high vs. low MSS-B (L1 and L2)
# mss_b_group    <- rep(c("Low", "High"), 4)
# task           <- c(rep("LT", 2), rep("C&C", 2), rep("SM", 2), rep("P", 2))
# my_mss_b_means <- c(mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == FALSE,]$lt_l1_l2_diff_perc,   na.rm=TRUE),
#                     mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == TRUE,]$lt_l1_l2_diff_perc,    na.rm=TRUE),
#                     mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == FALSE,]$ca_c_l1_l2_diff_perc, na.rm=TRUE),
#                     mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == TRUE,]$ca_c_l1_l2_diff_perc,  na.rm=TRUE),
#                     mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == FALSE,]$sm_l1_l2_diff_perc,   na.rm=TRUE),
#                     mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == TRUE,]$sm_l1_l2_diff_perc,    na.rm=TRUE),
#                     mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == FALSE,]$p_l1_l2_diff_perc,    na.rm=TRUE),
#                     mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == TRUE,]$p_l1_l2_diff_perc,     na.rm=TRUE))
# # Standard errors for error bars
# my_mss_b_errors <- c(std.error(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == FALSE,]$lt_l1_l2_diff_perc,   na.rm=TRUE),
#                      std.error(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == TRUE,]$lt_l1_l2_diff_perc,    na.rm=TRUE),
#                      std.error(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == FALSE,]$ca_c_l1_l2_diff_perc, na.rm=TRUE),
#                      std.error(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == TRUE,]$ca_c_l1_l2_diff_perc,  na.rm=TRUE),
#                      std.error(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == FALSE,]$sm_l1_l2_diff_perc,   na.rm=TRUE),
#                      std.error(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == TRUE,]$sm_l1_l2_diff_perc,    na.rm=TRUE),
#                      std.error(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == FALSE,]$p_l1_l2_diff_perc,    na.rm=TRUE),
#                      std.error(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == TRUE,]$p_l1_l2_diff_perc,     na.rm=TRUE))
# 
# df_mss_b_bars <- data.frame(mss_b_group, task, my_mss_b_means, my_mss_b_errors)
# 
# # Make them ordered factors
# df_mss_b_bars$task <- factor(df_mss_b_bars$task,
#                              unique(df_mss_b_bars$task),
#                              ordered=TRUE)
# df_mss_b_bars$mss_b_group <- factor(df_mss_b_bars$mss_b_group,
#                                     unique(df_mss_b_bars$mss_b_group),
#                                     ordered=TRUE)
# 
# # Plot of high vs. low MSS-B groups
# ggplot(df_mss_b_bars, aes(fill=mss_b_group, y=my_mss_b_means, x=task)) +
#     geom_bar(position="dodge", stat="identity") +
#     geom_errorbar(aes(ymin = my_mss_b_means - my_mss_b_errors,
#                       ymax = my_mss_b_means + my_mss_b_errors),
#                   width = 0.2, position = position_dodge(0.9)) +
#     xlab("Linguistic Task") + ylab("Mean Diff. (% of Total)") + #ylim(0,30) +
#     ggtitle("L1/L2 Differences by Level of Schizotypy") +
#     theme_bw() + theme(plot.title=element_text(hjust=0.5),
#                        panel.grid.major.x = element_blank(),
#                        panel.grid.minor = element_blank()) +
#     scale_fill_manual(values=c("darkred", "salmon"), name="Schizotypy") +
#     scale_y_continuous(breaks = seq(0, 30, by = 5), limits = c(0,30))
# 
# # Test for difference in level of L1/L2 difference between high- vs. low-MSS-B samples
# mss_b_t_test_lt   <- t.test(
#     x=df_matched_mss_b[df_matched_mss_b$has_high_mss_b == FALSE,]$lt_l1_l2_difference,
#     y=df_matched_mss_b[df_matched_mss_b$has_high_mss_b == TRUE,]$lt_l1_l2_difference,
#     alternative="two.sided", paired=FALSE, var.equal=FALSE)
# mss_b_t_test_ca_c <- t.test(
#     x=df_matched_mss_b[df_matched_mss_b$has_high_mss_b == FALSE,]$ca_c_l1_l2_difference,
#     y=df_matched_mss_b[df_matched_mss_b$has_high_mss_b == TRUE,]$ca_c_l1_l2_difference,
#     alternative="two.sided", paired=FALSE, var.equal=FALSE)
# mss_b_t_test_sm   <- t.test(
#     x=df_matched_mss_b[df_matched_mss_b$has_high_mss_b == FALSE,]$sm_l1_l2_difference,
#     y=df_matched_mss_b[df_matched_mss_b$has_high_mss_b == TRUE,]$sm_l1_l2_difference,
#     alternative="two.sided", paired=FALSE, var.equal=FALSE)
# mss_b_t_test_p    <- t.test(
#     x=df_matched_mss_b[df_matched_mss_b$has_high_mss_b == FALSE,]$p_l1_l2_difference,
#     y=df_matched_mss_b[df_matched_mss_b$has_high_mss_b == TRUE,]$p_l1_l2_difference,
#     alternative="two.sided", paired=FALSE, var.equal=FALSE)
# 
# mss_b_t_tests_p_values <- c(mss_b_t_test_lt["p.value"],
#                             mss_b_t_test_ca_c["p.value"],
#                             mss_b_t_test_sm["p.value"],
#                             mss_b_t_test_p["p.value"]
#                            )
# mss_b_t_tests_corr_p <- p.adjust(mss_b_t_tests_p_values,
#                                  method = "bonferroni",
#                                  n = length(mss_b_t_tests_p_values))
# # cohens_d <- abs(diff(mss_b_t_test_lt$estimate) / sqrt(mss_b_t_test_lt$parameter))

# # Analyses with PSM - high vs. low ASQ ----
# match_it_asq <- matchit(has_high_asq ~ age + gender + dsct_total + sst_total, data = df_for_matching, method="nearest", ratio=1)
# summary(match_it_asq)
# df_matched_asq <- match.data(match_it_asq)[1:ncol(my_data)]
# 
# # Create data frame of means by high vs. low ASQ
# asq_group    <- c(rep("high", 4), rep("low", 4), rep("high", 4), rep("low", 4))
# Language       <- c(rep("L1", 8), rep("L2", 8))
# task           <- rep(c("LT", "C&C", "SM", "P"), 4)
# my_asq_means <- c(mean(df_matched_asq[df_matched_asq$has_high_asq == TRUE,]$lt_l1_total,    na.rm=TRUE),
#                   mean(df_matched_asq[df_matched_asq$has_high_asq == TRUE,]$ca_c_l1_total,  na.rm=TRUE),
#                   mean(df_matched_asq[df_matched_asq$has_high_asq == TRUE,]$sm_l1_total,    na.rm=TRUE),
#                   mean(df_matched_asq[df_matched_asq$has_high_asq == TRUE,]$p_l1_total,     na.rm=TRUE),
#                   mean(df_matched_asq[df_matched_asq$has_high_asq == FALSE,]$lt_l1_total,   na.rm=TRUE),
#                   mean(df_matched_asq[df_matched_asq$has_high_asq == FALSE,]$ca_c_l1_total, na.rm=TRUE),
#                   mean(df_matched_asq[df_matched_asq$has_high_asq == FALSE,]$sm_l1_total,   na.rm=TRUE),
#                   mean(df_matched_asq[df_matched_asq$has_high_asq == FALSE,]$p_l1_total,    na.rm=TRUE),
#                   mean(df_matched_asq[df_matched_asq$has_high_asq == TRUE,]$lt_l2_total,    na.rm=TRUE),
#                   mean(df_matched_asq[df_matched_asq$has_high_asq == TRUE,]$ca_c_l2_total,  na.rm=TRUE),
#                   mean(df_matched_asq[df_matched_asq$has_high_asq == TRUE,]$sm_l2_total,    na.rm=TRUE),
#                   mean(df_matched_asq[df_matched_asq$has_high_asq == TRUE,]$p_l2_total,     na.rm=TRUE),
#                   mean(df_matched_asq[df_matched_asq$has_high_asq == FALSE,]$lt_l2_total,   na.rm=TRUE),
#                   mean(df_matched_asq[df_matched_asq$has_high_asq == FALSE,]$ca_c_l2_total, na.rm=TRUE),
#                   mean(df_matched_asq[df_matched_asq$has_high_asq == FALSE,]$sm_l2_total,   na.rm=TRUE),
#                   mean(df_matched_asq[df_matched_asq$has_high_asq == FALSE,]$p_l2_total,    na.rm=TRUE))
# df_asq_bars <- data.frame(asq_group, Language, task, my_asq_means)
# 
# # Make them ordered factors
# df_asq_bars$task <- factor(df_asq_bars$task,
#                            unique(df_asq_bars$task),
#                            ordered=TRUE)
# 
# # Need to add the my_facet_labels code from an older version -- got removed by accident at some point
#
# # Plot of high vs. low ASQ groups
# ggplot(df_asq_bars, aes(fill=Language, y=my_asq_means, x=task)) +
#        geom_bar(position="dodge", stat="identity") +
#        xlab("Linguistic Task") + ylab("Mean Task Score") +
#        ggtitle("Task Scores by ASQ Group") +
#        facet_wrap("asq_group", labeller=as_labeller(my_facet_labels)) +
#        theme_bw() + theme(plot.title=element_text(hjust=0.5)) +
#        scale_fill_discrete(name = "Language", labels = c("L1", "L2")) +
#        scale_fill_manual(values=c("darkred", "salmon"))

# # Test for difference in level of L1/L2 difference between high- vs. low-ASQ samples
# t.test(x=df_matched_asq[df_matched_asq$has_high_asq == FALSE,]$lt_l1_l2_difference,
#        y=df_matched_asq[df_matched_asq$has_high_asq == TRUE,]$lt_l1_l2_difference,
#        alternative="two.sided", paired=FALSE, var.equal=FALSE)
# t.test(x=df_matched_asq[df_matched_asq$has_high_asq == FALSE,]$ca_c_l1_l2_difference,
#        y=df_matched_asq[df_matched_asq$has_high_asq == TRUE,]$ca_c_l1_l2_difference,
#        alternative="two.sided", paired=FALSE, var.equal=FALSE)
# t.test(x=df_matched_asq[df_matched_asq$has_high_asq == FALSE,]$sm_l1_l2_difference,
#        y=df_matched_asq[df_matched_asq$has_high_asq == TRUE,]$sm_l1_l2_difference,
#        alternative="two.sided", paired=FALSE, var.equal=FALSE)
# t.test(x=df_matched_asq[df_matched_asq$has_high_asq == FALSE,]$p_l1_l2_difference,
#        y=df_matched_asq[df_matched_asq$has_high_asq == TRUE,]$p_l1_l2_difference,
#        alternative="two.sided", paired=FALSE, var.equal=FALSE)


# ====== Other analyses ====== ----

# # Regression ----
# 
# lm_lt_covar1   <- lm(data = my_data,
#                     formula = lt_l1_l2_difference ~ mss_b_total)
# summary(lm_lt_covar1)
# 
# lm_lt_covar2   <- lm(data = my_data,
#                     formula = lt_l1_l2_difference ~ mss_b_total + prof_l1_l2_diff)
# summary(lm_lt_covar2)
# 
# lm_lt_covar3   <- lm(data = my_data,
#                     formula = lt_l1_l2_difference ~ mss_b_pos + mss_b_neg + mss_b_dis + prof_l1_l2_diff)
# summary(lm_lt_covar3)
# 
# lm_lt_covar4   <- lm(data = my_data,
#                     formula = lt_l1_l2_difference ~ fr_speaking + mss_b_pos + mss_b_neg + mss_b_dis + prof_l1_l2_diff)
# summary(lm_lt_covar4)
# 
#
# lm_lt_sz_each <- lm(data = my_data,
#                     formula = lt_l1_l2_difference ~ mss_b_e_1 + mss_b_e_2 + mss_b_e_3 + mss_b_e_4 +
#                                                     mss_b_e_5 + mss_b_e_6 + mss_b_e_7 + mss_b_e_8 +
#                                                     mss_b_e_9 + mss_b_e_10 + mss_b_e_11 + mss_b_e_12 +
#                                                     mss_b_e_13 + mss_b_e_14 + mss_b_e_15 + mss_b_e_16 +
#                                                     mss_b_e_17 + mss_b_e_18 + mss_b_e_19 + mss_b_e_20 +
#                                                     mss_b_e_21 + mss_b_e_22 + mss_b_e_23 + mss_b_e_24 +
#                                                     mss_b_e_25 + mss_b_e_26 + mss_b_e_27 + mss_b_e_28 +
#                                                     mss_b_e_29 + mss_b_e_30 + mss_b_e_31 + mss_b_e_32 +
#                                                     mss_b_e_33 + mss_b_e_34 + mss_b_e_35 + mss_b_e_36 +
#                                                     mss_b_e_37 + mss_b_e_38)
# summary(lm_lt_sz_each)
# 
# lm_lt_sz_each <- lm(data = my_data,
#                     formula = lt_l1_l2_difference ~ mss_b_e_1  + mss_b_e_6  + mss_b_e_15 + mss_b_e_21 + 
#                                                     mss_b_e_23 + mss_b_e_24 + mss_b_e_25 + mss_b_e_28 +
#                                                     mss_b_e_36)
# summary(lm_lt_sz_each)
# 
# lm_ca_c_covar   <- lm(data = my_data,
#                     formula = ca_c_l1_l2_difference ~ mss_b_pos + mss_b_neg + mss_b_dis + prof_l1_l2_diff)
# summary(lm_ca_c_covar)
# 
# lm_sm_covar   <- lm(data = my_data,
#                     formula = sm_l1_l2_difference ~ mss_b_pos + mss_b_neg + mss_b_dis + prof_l1_l2_diff)
# summary(lm_sm_covar)
# 
# lm_p_covar   <- lm(data = my_data,
#                     formula = p_l1_l2_difference ~ mss_b_pos + mss_b_neg + mss_b_dis + prof_l1_l2_diff)
# summary(lm_p_covar)


# # Trying correlations and regressions that include L2 proficiency ----
# 
# # Correlations between L2 proficiency and performance on the linguistic tasks in L2
# cor.test(my_data$l2_prof_total, my_data$lt_l2_total)
# cor.test(my_data$l2_prof_total, my_data$ca_c_l2_total)
# cor.test(my_data$l2_prof_total, my_data$sm_l2_total)
# cor.test(my_data$l2_prof_total, my_data$p_l2_total)
# # Correlations between L2 proficiency and difference in L1 and L2 scores on the linguistic tasks
# cor.test(my_data$l2_prof_total, my_data$lt_l1_l2_difference)
# cor.test(my_data$l2_prof_total, my_data$ca_c_l1_l2_difference)
# cor.test(my_data$l2_prof_total, my_data$sm_l1_l2_difference)
# cor.test(my_data$l2_prof_total, my_data$p_l1_l2_difference)
# 
# # Attempting to predict L1/L2 differences in scores using L2 proficiency and symptoms
# lm_lt   <- lm(data = my_data, formula = lt_l1_l2_difference ~ mss_b_total + asq_total + l2_prof_total)
# lm_ca_c <- lm(data = my_data, formula = ca_c_l1_l2_difference ~ mss_b_total + asq_total + l2_prof_total)
# lm_sm   <- lm(data = my_data, formula = sm_l1_l2_difference ~ mss_b_total + asq_total + l2_prof_total)
# lm_p    <- lm(data = my_data, formula = p_l1_l2_difference ~ mss_b_total + asq_total + l2_prof_total)
# summary(lm_lt)
# summary(lm_ca_c)
# summary(lm_sm)
# summary(lm_p)
# 
# # Attempting to predict L1/L2 differences in scores from symptoms in a low-L2-proficiency subset of participants
# lm2_lt   <- lm(data = my_data, subset = (l2_prof_total <= median(my_data$l2_prof_total, na.rm = TRUE)),
#                                          formula = lt_l1_l2_difference ~ mss_b_total + asq_total)
# lm2_ca_c <- lm(data = my_data, subset = (l2_prof_total <= median(my_data$l2_prof_total, na.rm = TRUE)),
#                                          formula = ca_c_l1_l2_difference ~ mss_b_total + asq_total)
# lm2_sm   <- lm(data = my_data, subset = (l2_prof_total <= median(my_data$l2_prof_total, na.rm = TRUE)),
#                                          formula = sm_l1_l2_difference ~ mss_b_total + asq_total)
# lm2_p    <- lm(data = my_data, subset = (l2_prof_total <= median(my_data$l2_prof_total, na.rm = TRUE)),
#                                          formula = p_l1_l2_difference ~ mss_b_total + asq_total)
# summary(lm2_lt)
# summary(lm2_ca_c)
# summary(lm2_sm)
# summary(lm2_p)






