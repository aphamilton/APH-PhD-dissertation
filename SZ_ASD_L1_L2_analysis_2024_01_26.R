# Setting working directory and loading libraries ----

setwd("C:/Users/Arthur/OneDrive - Carleton University/Documents/R/PhD_Thesis/SZ_ASD_L1_L2_analysis/Data")

library(tidyverse)
library(janitor)
library(Hmisc)
library(GGally)    # For correlograms
library(MatchIt)   # For propensity score matching (PSM)

# Creating functions ----

# Creating function to get the sum of a block (returning NA if there are no values)
# Note: This gives the number of correct multiple-choice responses (i.e., the number marked as "1"), not the actual sum
sum_block <- function(my_data, my_row, first_col, last_col)
{
  col_count = last_col - first_col + 1
  
  result = NA
  if(sum(is.na(my_data[my_row, first_col:last_col])) < col_count)
    result = sum(my_data[my_row, first_col:last_col] == 1, na.rm=TRUE)

  out <- result
}

# Creating function to create L1 or L2 column for a task
convert_totals_to_l1_or_l2 <- function(my_data, is_for_l1, en_column, fr_column)
{
  out <- ifelse(my_data[,"is_l1_en"] == is_for_l1, my_data[,en_column], my_data[,fr_column])
}

# Creating function to add a single column for a question (merging EN and FR responses)
# This can only be used after the "is_l1_en" column is created
add_merged_col <- function(my_data, en_col_name, fr_col_name)
{
  new_col_vector <- vector("character", nrow(my_data))
  for (my_row in 1:nrow(my_data))
  {
    if (my_data[my_row, "is_l1_en"] == 1)
      new_col_vector[my_row] <- my_data[my_row, en_col_name]
    else
      new_col_vector[my_row] <- my_data[my_row, fr_col_name]
  out <- new_col_vector
  }
}


# Loading and cleaning data ----

# Loading Qualtrics data
my_data <- read.csv(file = "QualtricsData_2023_12_11.csv", encoding="UTF-8")
# my_data <- read.csv(file = "QualtricsData_2023_12_09.csv", encoding="UTF-8")

# Cleaning data frame with the janitor package
my_data <- clean_names(my_data)

# Remove first two rows, which are just info from Qualtrics
my_data <- my_data[-c(1, 2), ]

# Switch duration in seconds to numeric value
my_data$duration_in_seconds <- as.numeric(my_data$duration_in_seconds)

# Removing rows that are not good data
my_data <- my_data[my_data$status == "IP Address", ]
my_data <- my_data[my_data$finished == "True", ]
my_data <- my_data[my_data$give_consent == "I consent to participate. / Je donne mon consentement à participer.", ]
my_data <- subset(my_data, my_data$eligibility_age == "Yes / Oui" | my_data$eligibility_age == "")
my_data <- subset(my_data, my_data$eligibility_language == "Yes / Oui" | my_data$eligibility_language == "")
my_data <- subset(my_data, my_data$first_language == "English / Anglais" | my_data$first_language == "French / Français")
my_data <- subset(my_data, my_data$l1 == "English / Anglais" | my_data$l1 == "French / Français")
my_data <- my_data[my_data$duration_in_seconds >= 1200, ]
# Should modify the above to only accept blank eligibility responses from Prolific participants

# Removing rows that were practice trials
qualtrics_test_runs <- c("R_ex2UqHGk3N4rr4l", "R_0p64Pt4Pn1xe5fb", "R_2YbRN1GFsE1FhdG", "R_1k17QY4T3sbgvwZ", "R_ywfTKwA3Qcc1lXH", "R_0p64Pt4Pn1xe5fb")
my_data <- my_data[!(my_data$response_id %in% qualtrics_test_runs), ]

# Correcting the source for two participants
# I entered the sources manually due to technical problems, but got the sources wrong and only noticed after they started
my_data <- within(my_data, source[random_id == 33534377 & source == "sona_cogsci"] <- "sona_psych")
my_data <- within(my_data, source[random_id == 98996051 & source == "sona_psych"] <- "sona_cogsci")

# Standardizing source names
my_data <- within(my_data, source[source == "cogsci_sona"] <- "sona_cogsci")

# Excluding participants whose source does not correspond to a recruitment avenue
my_data <- subset(my_data, my_data$source == "sona_cogsci" | my_data$source == "sona_psych" | my_data$source == "prolific_l1e" | my_data$source == "prolific_l1f")

# Loading answer key -- we will check participant responses against this below
my_key <- read.csv(file = "QualtricsAnswerKeyNoParticipants_2023_07_12.csv", encoding="UTF-8")

# Cleaning answer key with the janitor package so that the column names match the main data frame
my_key <- clean_names(my_key)


# Binarize participant responses as correct or incorrect ----

# Simplifying 4-point ASQ scale to 2-point scale (agree/disagree)
for (my_col in 313:362)
{
  my_data[,my_col] <- ifelse(my_data[,my_col] == "",
                             "",
                             ifelse(my_data[,my_col] == "Definitely agree" |
                                    my_data[,my_col] == "Slightly agree",
                                    "Agree",
                                    "Disagree")
                             )
}

# Same as above but for French
# Calls to the 3rd row of my_key are the French spellings to get the special characters
for (my_col in 421:470)
{
  my_data[,my_col] <- ifelse(my_data[,my_col] == "",
                             "",
                             ifelse(my_data[,my_col] == my_key[3,1] |
                                    my_data[,my_col] == my_key[3,2],
                                    my_key[3,5],
                                    my_key[3,6])
                             )
}

# Convert values to correct/incorrect (1 or 0)
for (my_col in 1:ncol(my_data))
{
  if (my_key[2,my_col] != "")
  {
    my_data[,my_col] <- ifelse(my_data[,my_col] == "",
                               NA,
                               ifelse(my_data[,my_col] == my_key[2,my_col],
                                      1,
                                      0)
                               )
  }
}

# Create column with L1 as 1 (English) or 0 (French)
my_data[,"is_l1_en"] <- ifelse(my_data[,"first_language"] == "English / Anglais", 1, 0)

# Calculate total scores by participant for each task or questionnaire ----

# Total for each block in the experiment by participant
for (my_row in 1:nrow(my_data))
{
  # sum_block really counts occurrences
  my_data[my_row,   "lt_e_total"]  <- sum_block(my_data, my_row,  15,  44)
  my_data[my_row,   "lt_f_total"]  <- sum_block(my_data, my_row,  45,  74)
  my_data[my_row, "ca_c_e1_total"] <- sum_block(my_data, my_row,  75,  90)
  my_data[my_row, "ca_c_f2_total"] <- sum_block(my_data, my_row,  91, 106)
  my_data[my_row, "ca_c_e2_total"] <- sum_block(my_data, my_row, 107, 122)
  my_data[my_row, "ca_c_f1_total"] <- sum_block(my_data, my_row, 123, 138)
  my_data[my_row,   "sm_e1_total"] <- sum_block(my_data, my_row, 139, 154)
  my_data[my_row,   "sm_f2_total"] <- sum_block(my_data, my_row, 155, 170)
  my_data[my_row,   "sm_e2_total"] <- sum_block(my_data, my_row, 171, 186)
  my_data[my_row,   "sm_f1_total"] <- sum_block(my_data, my_row, 187, 202)
  my_data[my_row,    "p_e1_total"] <- sum_block(my_data, my_row, 203, 220)
  my_data[my_row,    "p_f2_total"] <- sum_block(my_data, my_row, 221, 238)
  my_data[my_row,    "p_e2_total"] <- sum_block(my_data, my_row, 239, 256)
  my_data[my_row,    "p_f1_total"] <- sum_block(my_data, my_row, 257, 274)
  my_data[my_row, "mss_b_e_total"] <- sum_block(my_data, my_row, 275, 312)
  my_data[my_row,   "asq_e_total"] <- sum_block(my_data, my_row, 313, 362)
  my_data[my_row, "mss_b_f_total"] <- sum_block(my_data, my_row, 383, 420)
  my_data[my_row,   "asq_f_total"] <- sum_block(my_data, my_row, 421, 470)
  
  # Calculating scores for symptom sub-scales

  my_data[my_row,"mss_b_e_neg"] <- sum(c(my_data[my_row,"mss_b_e_1"],  my_data[my_row,"mss_b_e_4"],  my_data[my_row,"mss_b_e_7"],
                                         my_data[my_row,"mss_b_e_10"], my_data[my_row,"mss_b_e_13"], my_data[my_row,"mss_b_e_16"],
                                         my_data[my_row,"mss_b_e_19"], my_data[my_row,"mss_b_e_22"], my_data[my_row,"mss_b_e_25"],
                                         my_data[my_row,"mss_b_e_28"], my_data[my_row,"mss_b_e_31"], my_data[my_row,"mss_b_e_34"],
                                         my_data[my_row,"mss_b_e_37"]), na.rm=TRUE)
  my_data[my_row,"mss_b_e_pos"] <- sum(c(my_data[my_row,"mss_b_e_2"],  my_data[my_row,"mss_b_e_5"],  my_data[my_row,"mss_b_e_8"],
                                         my_data[my_row,"mss_b_e_11"], my_data[my_row,"mss_b_e_14"], my_data[my_row,"mss_b_e_17"],
                                         my_data[my_row,"mss_b_e_20"], my_data[my_row,"mss_b_e_23"], my_data[my_row,"mss_b_e_26"],
                                         my_data[my_row,"mss_b_e_29"], my_data[my_row,"mss_b_e_32"], my_data[my_row,"mss_b_e_35"],
                                         my_data[my_row,"mss_b_e_38"]), na.rm=TRUE)
  my_data[my_row,"mss_b_e_dis"] <- sum(c(my_data[my_row,"mss_b_e_3"],  my_data[my_row,"mss_b_e_6"],  my_data[my_row,"mss_b_e_9"],
                                         my_data[my_row,"mss_b_e_12"], my_data[my_row,"mss_b_e_15"], my_data[my_row,"mss_b_e_18"],
                                         my_data[my_row,"mss_b_e_21"], my_data[my_row,"mss_b_e_24"], my_data[my_row,"mss_b_e_27"],
                                         my_data[my_row,"mss_b_e_30"], my_data[my_row,"mss_b_e_33"], my_data[my_row,"mss_b_e_36"]),
                                         na.rm=TRUE)
  
  my_data[my_row,"mss_b_f_neg"] <- sum(c(my_data[my_row,"mss_b_f_1"],  my_data[my_row,"mss_b_f_4"],  my_data[my_row,"mss_b_f_7"],
                                         my_data[my_row,"mss_b_f_10"], my_data[my_row,"mss_b_f_13"], my_data[my_row,"mss_b_f_16"],
                                         my_data[my_row,"mss_b_f_19"], my_data[my_row,"mss_b_f_22"], my_data[my_row,"mss_b_f_25"],
                                         my_data[my_row,"mss_b_f_28"], my_data[my_row,"mss_b_f_31"], my_data[my_row,"mss_b_f_34"],
                                         my_data[my_row,"mss_b_f_37"]), na.rm=TRUE)
  my_data[my_row,"mss_b_f_pos"] <- sum(c(my_data[my_row,"mss_b_f_2"],  my_data[my_row,"mss_b_f_5"],  my_data[my_row,"mss_b_f_8"],
                                         my_data[my_row,"mss_b_f_11"], my_data[my_row,"mss_b_f_14"], my_data[my_row,"mss_b_f_17"],
                                         my_data[my_row,"mss_b_f_20"], my_data[my_row,"mss_b_f_23"], my_data[my_row,"mss_b_f_26"],
                                         my_data[my_row,"mss_b_f_29"], my_data[my_row,"mss_b_f_32"], my_data[my_row,"mss_b_f_35"],
                                         my_data[my_row,"mss_b_f_38"]), na.rm=TRUE)
  my_data[my_row,"mss_b_f_dis"] <- sum(c(my_data[my_row,"mss_b_f_3"],  my_data[my_row,"mss_b_f_6"],  my_data[my_row,"mss_b_f_9"],
                                         my_data[my_row,"mss_b_f_12"], my_data[my_row,"mss_b_f_15"], my_data[my_row,"mss_b_f_18"],
                                         my_data[my_row,"mss_b_f_21"], my_data[my_row,"mss_b_f_24"], my_data[my_row,"mss_b_f_27"],
                                         my_data[my_row,"mss_b_f_30"], my_data[my_row,"mss_b_f_33"], my_data[my_row,"mss_b_f_36"]),
                                         na.rm=TRUE)
  
  # Merging blocks for same task and language by participant
  my_data[my_row,"ca_c_e_total"] <- sum(c(my_data[my_row,"ca_c_e1_total"], my_data[my_row,"ca_c_e2_total"]), na.rm=TRUE)
  my_data[my_row,"ca_c_f_total"] <- sum(c(my_data[my_row,"ca_c_f1_total"], my_data[my_row,"ca_c_f2_total"]), na.rm=TRUE)
  my_data[my_row,  "sm_e_total"] <- sum(c(my_data[my_row,  "sm_e1_total"], my_data[my_row,  "sm_e2_total"]), na.rm=TRUE)
  my_data[my_row,  "sm_f_total"] <- sum(c(my_data[my_row,  "sm_f1_total"], my_data[my_row,  "sm_f2_total"]), na.rm=TRUE)
  my_data[my_row,   "p_e_total"] <- sum(c(my_data[my_row,   "p_e1_total"], my_data[my_row,   "p_e2_total"]), na.rm=TRUE)
  my_data[my_row,   "p_f_total"] <- sum(c(my_data[my_row,   "p_f1_total"], my_data[my_row,   "p_f2_total"]), na.rm=TRUE)
  
  # Merging symptom scores by participant
  my_data[my_row,"mss_b_total"] <- sum(c(my_data[my_row,"mss_b_e_total"], my_data[my_row,"mss_b_f_total"]), na.rm=TRUE)
  my_data[my_row,  "mss_b_neg"] <- sum(c(my_data[my_row,  "mss_b_e_neg"], my_data[my_row,  "mss_b_f_neg"]), na.rm=TRUE)
  my_data[my_row,  "mss_b_pos"] <- sum(c(my_data[my_row,  "mss_b_e_pos"], my_data[my_row,  "mss_b_f_pos"]), na.rm=TRUE)
  my_data[my_row,  "mss_b_dis"] <- sum(c(my_data[my_row,  "mss_b_e_dis"], my_data[my_row,  "mss_b_f_dis"]), na.rm=TRUE)
  my_data[my_row,  "asq_total"] <- sum(c(my_data[my_row,  "asq_e_total"], my_data[my_row,  "asq_f_total"]), na.rm=TRUE)
}

# This code below here can definitely be simplified

# my_data$dlq_1 <- add_merged_col(my_data, "dlq_e_1", "dlq_f_1")
# my_data$dlq_2 <- add_merged_col(my_data, "dlq_e_2", "dlq_f_2")

# Creating column for age
my_data$dlq_1 <- NA
for (my_row in 1:nrow(my_data))
{
  if (my_data[my_row, "is_l1_en"] == 1)
    my_data[my_row, "dlq_1"] <- my_data[my_row, "dlq_e_1"]
  else
    my_data[my_row, "dlq_1"] <- my_data[my_row, "dlq_f_1"]
}
my_data$dlq_1 <- as.numeric(my_data$dlq_1)

# Creating column for gender
my_data$dlq_2 <- NA
for (my_row in 1:nrow(my_data))
{
  if (my_data[my_row, "is_l1_en"] == 1)
    my_data[my_row, "dlq_2"] <- my_data[my_row, "dlq_e_2"]
  else
    my_data[my_row, "dlq_2"] <- my_data[my_row, "dlq_f_2"]
}
my_data[(my_data$dlq_2 == "Female" | my_data$dlq_2 == "Femme"), "dlq_2"] <- 0
my_data[(my_data$dlq_2 == "Male" | my_data$dlq_2 == "Homme"), "dlq_2"] <- 1
my_data$dlq_2 <- as.numeric(my_data$dlq_2)


# # Creating merged column for a DLQ result
# my_data$dlq_1 <- NA
# for (my_row in 1:nrow(my_data))
# {
# if (my_data[my_row, "is_l1_en"] == 1)
#   my_data[my_row, "dlq_1"] <- my_data[my_row, "dlq_e_1"]
# else
#   my_data[my_row, "dlq_1"] <- my_data[my_row, "dlq_f_1"]
# }


# Create synthetic columns for analyzing L1/L2 rather than En/Fr ----

# Create columns for L1 and L2 totals for each linguistic task
my_data$lt_l1_total   <- convert_totals_to_l1_or_l2(my_data, 1,   "lt_e_total",   "lt_f_total")
my_data$lt_l2_total   <- convert_totals_to_l1_or_l2(my_data, 0,   "lt_e_total",   "lt_f_total")
my_data$ca_c_l1_total <- convert_totals_to_l1_or_l2(my_data, 1, "ca_c_e_total", "ca_c_f_total")
my_data$ca_c_l2_total <- convert_totals_to_l1_or_l2(my_data, 0, "ca_c_e_total", "ca_c_f_total")
my_data$sm_l1_total   <- convert_totals_to_l1_or_l2(my_data, 1,   "sm_e_total",   "sm_f_total")
my_data$sm_l2_total   <- convert_totals_to_l1_or_l2(my_data, 0,   "sm_e_total",   "sm_f_total")
my_data$p_l1_total    <- convert_totals_to_l1_or_l2(my_data, 1,    "p_e_total",    "p_f_total")
my_data$p_l2_total    <- convert_totals_to_l1_or_l2(my_data, 0,    "p_e_total",    "p_f_total")

# Create columns for L1/L2 differences for each linguistic task
my_data$lt_l1_l2_difference   <- my_data$lt_l1_total   - my_data$lt_l2_total
my_data$ca_c_l1_l2_difference <- my_data$ca_c_l1_total - my_data$ca_c_l2_total
my_data$sm_l1_l2_difference   <- my_data$sm_l1_total   - my_data$sm_l2_total
my_data$p_l1_l2_difference    <- my_data$p_l1_total    - my_data$p_l2_total

# Load hand-coded results and merge into main data frame ----

# Loading Pavlovia data
my_data_pavlovia <- read.csv(file = "PavloviaCombinedData_2023_07_14.csv", encoding="UTF-8")
# Merge Pavlovia data into main data frame
my_data <- merge(my_data, my_data_pavlovia, by="random_id", all.x = TRUE, all.y = TRUE)

# Loading coded DLQ data
my_data_dlq <- read.csv(file = "DLQCodedData_2023_10_28.csv", encoding="UTF-8")
# Cleaning coded DLQ data with janitor package
my_data_dlq <- clean_names(my_data_dlq)
# Remove first two rows, which are just info from Qualtrics
my_data_dlq <- my_data_dlq[-c(1, 2), ]
# Keep only relevant columns
my_data_dlq <- my_data_dlq[ ,c("random_id","l2_speaking","l2_understanding","l2_writing","l2_reading")]
# Merge coded DLQ data data into main data frame
my_data <- merge(my_data, my_data_dlq, by="random_id", all.x = TRUE, all.y = TRUE)
#Calculating total across L2 proficiency scores
for (my_row in 1:nrow(my_data))
{
  if (!any(is.na(my_data[my_row, c("l2_speaking","l2_understanding","l2_writing","l2_reading")])))
    my_data[my_row, "l2_prof_total"] <- sum(my_data[my_row, c("l2_speaking","l2_understanding","l2_writing","l2_reading")], na.rm=TRUE)
  else
    my_data[my_row, "l2_prof_total"] <- NA
}


# Histograms ----

# Histogram by time to complete in Qualtrics
# mean(my_data[my_data$duration_in_seconds > 1200 & my_data$duration_in_seconds < 10600,]$duration_in_seconds)
# median(my_data[my_data$duration_in_seconds > 1200 & my_data$duration_in_seconds < 10600,]$duration_in_seconds)
# median(my_data$duration_in_seconds)
ggplot(data=my_data[my_data$duration_in_seconds > 1200 & my_data$duration_in_seconds < 10600,], aes(x=duration_in_seconds, rm.na = TRUE)) +
       geom_histogram(fill='blue', color='lightgrey') +
       xlab("Time (seconds)") + ylab("Number of participants") +
       ggtitle("Distribution of Qualtrics Completion Times") +
       theme_bw() + theme(plot.title=element_text(hjust=0.5))

# Histograms by L1 linguistic task scores
ggplot(data=my_data, aes(x=lt_l1_total, rm.na = TRUE)) +
       geom_histogram(fill='blue', color='lightgrey') +
       xlab("LexTale score") + ylab("Number of participants") +
       ggtitle("LexTale scores (L1)") +
       theme_bw() + theme(plot.title=element_text(hjust=0.5))
ggplot(data=my_data, aes(x=ca_c_l1_total, rm.na = TRUE)) +
       geom_histogram(fill='blue', color='lightgrey') +
       xlab("Camel and Cactus task score") + ylab("Number of participants") +
       ggtitle("Camel and Cactus scores (L1)") +
       theme_bw() + theme(plot.title=element_text(hjust=0.5))
ggplot(data=my_data, aes(x=sm_l1_total, rm.na = TRUE)) +
       geom_histogram(fill='blue', color='lightgrey') +
       xlab("Syntactic modification task score") + ylab("Number of participants") +
       ggtitle("Syntactic modification scores (L1)") +
       theme_bw() + theme(plot.title=element_text(hjust=0.5))
ggplot(data=my_data, aes(x=p_l1_total, rm.na = TRUE)) +
       geom_histogram(fill='blue', color='lightgrey') +
       xlab("Pragmatics task scores") + ylab("Number of participants") +
       ggtitle("Pragmatics scores (L1)") +
       theme_bw() + theme(plot.title=element_text(hjust=0.5))

# Histograms by MSS-B and ASQ scores
ggplot(data=my_data, aes(x=mss_b_total, rm.na = TRUE)) +
       geom_histogram(fill='blue', color='lightgrey') +
       xlab("MSS-B score") + ylab("Number of participants") +
       ggtitle("Distribution of schizotypy scores") +
       theme_bw() + theme(plot.title=element_text(hjust=0.5))
ggplot(data=my_data, aes(x=asq_total, rm.na = TRUE)) +
       geom_histogram(fill='blue', color='lightgrey') +
       xlab("ASQ score") + ylab("Number of participants") +
       ggtitle("Distribution of autistic traits scores") +
       theme_bw() + theme(plot.title=element_text(hjust=0.5))


# Matching ----
# Mainly relied on this: https://www.r-bloggers.com/2016/06/how-to-use-r-for-matching-samples-propensity-score/

# Adding columns for arbitrary high/low split in MSS-B and ASQ
my_data$has_high_mss_b <- FALSE
my_data[my_data$mss_b_total >= 13,]$has_high_mss_b <- TRUE
my_data$has_high_asq <- FALSE
my_data[my_data$asq_total >= 30,]$has_high_asq <- TRUE

# Creating data frame for matching (only includes participants who provided their age and gender (and only included male and female for purposes of matching)
df_for_matching <- my_data[!is.na(my_data$dlq_1) & !is.na(my_data$dlq_2),]

# High vs. low MSS-B
match_it_mss_b <- matchit(has_high_mss_b ~ dlq_1 + dlq_2, data = df_for_matching, method="nearest", ratio=1)
summary(match_it_mss_b)
df_matched_mss_b <- match.data(match_it_mss_b)[1:ncol(my_data)]

# Create data frame of means by high vs. low MSS-B
mss_b_group <- c(rep("high", 4), rep("low", 4), rep("high", 4), rep("low", 4))
language    <- c(rep("l1", 8), rep("l2", 8))
task        <- rep(c("lt", "ca_c", "sm", "p"), 4)
my_means    <- c(mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == TRUE,]$lt_l1_total,    na.rm=TRUE),
                 mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == TRUE,]$ca_c_l1_total,  na.rm=TRUE),
                 mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == TRUE,]$sm_l1_total,    na.rm=TRUE),
                 mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == TRUE,]$p_l1_total,     na.rm=TRUE),
                 mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == FALSE,]$lt_l1_total,   na.rm=TRUE),
                 mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == FALSE,]$ca_c_l1_total, na.rm=TRUE),
                 mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == FALSE,]$sm_l1_total,   na.rm=TRUE),
                 mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == FALSE,]$p_l1_total,    na.rm=TRUE),
                 mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == TRUE,]$lt_l2_total,    na.rm=TRUE),
                 mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == TRUE,]$ca_c_l2_total,  na.rm=TRUE),
                 mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == TRUE,]$sm_l2_total,    na.rm=TRUE),
                 mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == TRUE,]$p_l2_total,     na.rm=TRUE),
                 mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == FALSE,]$lt_l2_total,   na.rm=TRUE),
                 mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == FALSE,]$ca_c_l2_total, na.rm=TRUE),
                 mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == FALSE,]$sm_l2_total,   na.rm=TRUE),
                 mean(df_matched_mss_b[df_matched_mss_b$has_high_mss_b == FALSE,]$p_l2_total,    na.rm=TRUE))
df_mss_b_bars <- data.frame(mss_b_group, language, task, my_means)

ggplot(df_mss_b_bars, aes(fill=language, y=my_means, x=task)) + 
    geom_bar(position="dodge", stat="identity") +
    ggtitle("Task scores by MSS-B group") +
    facet_wrap("mss_b_group") +
    theme_bw() + theme(plot.title=element_text(hjust=0.5))

# Test for difference in level of L1/L2 difference between high- vs. low-MSS-B samples
t.test(x=df_matched_mss_b[df_matched_mss_b$has_high_mss_b == FALSE,]$lt_l1_l2_difference,
       y=df_matched_mss_b[df_matched_mss_b$has_high_mss_b == TRUE,]$lt_l1_l2_difference,
       alternative="two.sided", paired=FALSE, var.equal=FALSE)
t.test(x=df_matched_mss_b[df_matched_mss_b$has_high_mss_b == FALSE,]$ca_c_l1_l2_difference,
       y=df_matched_mss_b[df_matched_mss_b$has_high_mss_b == TRUE,]$ca_c_l1_l2_difference,
       alternative="two.sided", paired=FALSE, var.equal=FALSE)
t.test(x=df_matched_mss_b[df_matched_mss_b$has_high_mss_b == FALSE,]$sm_l1_l2_difference,
       y=df_matched_mss_b[df_matched_mss_b$has_high_mss_b == TRUE,]$sm_l1_l2_difference,
       alternative="two.sided", paired=FALSE, var.equal=FALSE)
t.test(x=df_matched_mss_b[df_matched_mss_b$has_high_mss_b == FALSE,]$p_l1_l2_difference,
       y=df_matched_mss_b[df_matched_mss_b$has_high_mss_b == TRUE,]$p_l1_l2_difference,
       alternative="two.sided", paired=FALSE, var.equal=FALSE)

# High vs. low ASQ
match_it_asq <- matchit(has_high_asq ~ dlq_1 + dlq_2, data = df_for_matching, method="nearest", ratio=1)
summary(match_it_asq)
df_matched_asq <- match.data(match_it_asq)[1:ncol(my_data)]

# Should add a barplot for high vs. low ASQ analogous to the MSS-B barplot above

# Test for difference in level of L1/L2 difference between high- vs. low-ASQ samples
t.test(x=df_matched_asq[df_matched_asq$has_high_asq == FALSE,]$lt_l1_l2_difference,
       y=df_matched_asq[df_matched_asq$has_high_asq == TRUE,]$lt_l1_l2_difference,
       alternative="two.sided", paired=FALSE, var.equal=FALSE)
t.test(x=df_matched_asq[df_matched_asq$has_high_asq == FALSE,]$ca_c_l1_l2_difference,
       y=df_matched_asq[df_matched_asq$has_high_asq == TRUE,]$ca_c_l1_l2_difference,
       alternative="two.sided", paired=FALSE, var.equal=FALSE)
t.test(x=df_matched_asq[df_matched_asq$has_high_asq == FALSE,]$sm_l1_l2_difference,
       y=df_matched_asq[df_matched_asq$has_high_asq == TRUE,]$sm_l1_l2_difference,
       alternative="two.sided", paired=FALSE, var.equal=FALSE)
t.test(x=df_matched_asq[df_matched_asq$has_high_asq == FALSE,]$p_l1_l2_difference,
       y=df_matched_asq[df_matched_asq$has_high_asq == TRUE,]$p_l1_l2_difference,
       alternative="two.sided", paired=FALSE, var.equal=FALSE)

# Create reduced data frames for analysis and compute correlation matrices ----

# Create data frame with only columns for analysis
my_reduced_data <- subset(my_data, select = c(mss_b_total:fhat_total_score))
# Create correlation matrix of that data frame
cor_matrix <- rcorr(as.matrix(my_reduced_data),type="pearson")
cor_matrix

# Create data frame with symptom scale and sub-scale totals
all_symptoms_data <- subset(my_data, select = c(mss_b_total, mss_b_neg, mss_b_pos, mss_b_dis, asq_total))

# Create data frame with only columns for analysis
all_l1_tasks_data <- subset(my_data, select = c(lt_l1_total, ca_c_l1_total, sm_l1_total, p_l1_total,
                                                dsct_total_score, sst_total_score, fhat_total_score))
# Create correlation matrix of that data frame
cor_matrix_l1_tasks <- rcorr(as.matrix(all_l1_tasks_data),type="pearson")
cor_matrix_l1_tasks

# Create data frame with only columns for analysis
my_reduced_data_l1en <- subset(my_data, select = c(l1,mss_b_total:fhat_total_score))
my_reduced_data_l1en <- subset(my_reduced_data_l1en, my_reduced_data_l1en$l1 == "English / Anglais")
my_reduced_data_l1en <- my_reduced_data_l1en[ , !names(my_reduced_data_l1en) %in% c("l1")]
# Create correlation matrix of that data frame
cor_matrix_l1en <- rcorr(as.matrix(my_reduced_data_l1en),type="pearson")
cor_matrix_l1en

# Create data frame with only columns for analysis
my_reduced_data_l1fr <- subset(my_data, select = c(l1,mss_b_total:fhat_total_score))
my_reduced_data_l1fr <- subset(my_reduced_data_l1fr, my_reduced_data_l1fr$l1 != "English / Anglais")
my_reduced_data_l1fr <- my_reduced_data_l1fr[ , !names(my_reduced_data_l1fr) %in% c("l1")]
# Create correlation matrix of that data frame
cor_matrix_l1fr <- rcorr(as.matrix(my_reduced_data_l1fr),type="pearson")
cor_matrix_l1fr

# Create data frame with L1 scores and symptom totals
l1_scores_symptoms_data <- subset(my_data, select = c(lt_l1_total, ca_c_l1_total, sm_l1_total, p_l1_total,
                                                        mss_b_total, mss_b_neg, mss_b_pos, mss_b_dis, asq_total))

# Create data frame with L2 scores and symptom totals
l2_scores_symptoms_data <- subset(my_data, select = c(lt_l2_total, ca_c_l2_total, sm_l2_total, p_l2_total,
                                                        mss_b_total, mss_b_neg, mss_b_pos, mss_b_dis, asq_total))

# Create data frame with L2 scores and self-rated L2 proficiency
l2_scores_proficiency_data <- subset(my_data, select = c(lt_l2_total, ca_c_l2_total, sm_l2_total, p_l2_total,
                                                        l2_speaking, l2_understanding, l2_writing, l2_reading, l2_prof_total))

# Create data frame with L1/L2 differences and symptom totals
l1_l2_diffs_symptoms_data <- subset(my_data, select = c(lt_l1_l2_difference, ca_c_l1_l2_difference,
                                                        sm_l1_l2_difference, p_l1_l2_difference,
                                                        mss_b_total, mss_b_neg, mss_b_pos, mss_b_dis, asq_total))

# Test correlation between ASQ and MSS-B scores
symptoms_cor <- cor.test(my_data$asq_total, my_data$mss_b_total)
symptoms_cor

# Create visualizations ----

# Create scatterplot of MSS-B vs AQ
ggplot(data=my_data, aes(x=mss_b_total, y=asq_total)) +
   geom_vline(xintercept = 0) + geom_vline(xintercept = 38) +
   geom_hline(yintercept = 0) + geom_hline(yintercept = 50) +
   geom_point(color='darkred', alpha=0.25) + geom_smooth(method='lm', formula=y~x, se=TRUE, color='black') +
   xlab("MSS-B Score") + ylab("AQ Score") + ggtitle("Relationship Between Schizotypal and Autistic Symptoms") +
   scale_x_continuous(limits = c(0, 38)) + scale_y_continuous(limits = c(0, 50)) +
   theme_minimal() + theme(plot.title=element_text(hjust=0.5))

#Create correlogram of symptom scale and sub-scale totals
ggpairs(all_symptoms_data, title="Correlogram of Symptom Scales and Sub-Scales",
        upper = list(continuous = "cor", combo = "box_no_facet", discrete = "facetbar", na = "na"),
        lower = list(continuous = wrap("smooth", color='darkred', alpha=0.15), combo = "facethist", discrete = "facetbar", na = "na"),
        diag = list(continuous = "blankDiag", discrete = "barDiag", na = "naDiag"),
        columnLabels = c("MSS-B", "MSS-B (Neg.)", "MSS-B (Pos.)", "MSS-B (Dis.)", "ASQ"), axisLabels = "none") +
        theme_bw() + theme(plot.title=element_text(hjust=0.5))

#Create correlogram of L1 tasks (incl. visual tasks)
ggpairs(all_l1_tasks_data, title="Correlogram of Cognitive Tasks (L1 Only)",
        upper = list(continuous = "cor", combo = "box_no_facet", discrete = "facetbar", na = "na"),
        lower = list(continuous = wrap("smooth", color='darkred', alpha=0.15), combo = "facethist", discrete = "facetbar", na = "na"),
        diag = list(continuous = "blankDiag", discrete = "barDiag", na = "naDiag"),
        columnLabels = c("LT", "C&C", "SM", "Prag.", "DSC", "SS", "FHA"), axisLabels = "none") +
        theme_bw() + theme(plot.title=element_text(hjust=0.5))

#Create correlogram of L1 scores and symptom scales
ggpairs(l1_scores_symptoms_data, title="Correlogram of L1 Scores and Symptoms",
        upper = list(continuous = "cor", combo = "box_no_facet", discrete = "facetbar", na = "na"),
        lower = list(continuous = wrap("smooth", color='darkred', alpha=0.15), combo = "facethist", discrete = "facetbar", na = "na"),
        diag = list(continuous = "blankDiag", discrete = "barDiag", na = "naDiag"),
        columnLabels = c("LT", "C&C", "SM", "Prag.", "MSS-B", "MSS-B (Neg.)", "MSS-B (Pos.)", "MSS-B (Dis.)", "ASQ"), axisLabels = "none") +
        theme_bw() + theme(plot.title=element_text(hjust=0.5))

#Create correlogram of L2 scores and symptom scales
ggpairs(l2_scores_symptoms_data, title="Correlogram of L2 Scores and Symptoms",
        upper = list(continuous = "cor", combo = "box_no_facet", discrete = "facetbar", na = "na"),
        lower = list(continuous = wrap("smooth", color='darkred', alpha=0.15), combo = "facethist", discrete = "facetbar", na = "na"),
        diag = list(continuous = "blankDiag", discrete = "barDiag", na = "naDiag"),
        columnLabels = c("LT", "C&C", "SM", "Prag.", "MSS-B", "MSS-B (Neg.)", "MSS-B (Pos.)", "MSS-B (Dis.)", "ASQ"), axisLabels = "none") +
        theme_bw() + theme(plot.title=element_text(hjust=0.5))

#Create correlogram of L2 scores and self-rated proficiency
ggpairs(l2_scores_proficiency_data, title="Correlogram of L2 Task and Proficiency Scores",
        upper = list(continuous = "cor", combo = "box_no_facet", discrete = "facetbar", na = "na"),
        lower = list(continuous = wrap("smooth", color='darkred', alpha=0.15), combo = "facethist", discrete = "facetbar", na = "na"),
        diag = list(continuous = "blankDiag", discrete = "barDiag", na = "naDiag"),
        columnLabels = c("LT", "C&C", "SM", "Prag.", "Speaking", "Listening", "Writing", "Reading", "Prof. Total"), axisLabels = "none") +
        theme_bw() + theme(plot.title=element_text(hjust=0.5))

#Create correlogram of L1/L2 differences and symptom scales
ggpairs(l1_l2_diffs_symptoms_data, title="Correlogram of L1/L2 Differences and Symptoms",
        upper = list(continuous = "cor", combo = "box_no_facet", discrete = "facetbar", na = "na"),
        lower = list(continuous = wrap("smooth", color='darkred', alpha=0.15), combo = "facethist", discrete = "facetbar", na = "na"),
        diag = list(continuous = "blankDiag", discrete = "barDiag", na = "naDiag"),
        columnLabels = c("LT Diff.", "C&C Diff.", "SM Diff.", "Prag. Diff.", "MSS-B", "MSS-B (Neg.)", "MSS-B (Pos.)", "MSS-B (Dis.)", "ASQ"), axisLabels = "none") +
        theme_bw() + theme(plot.title=element_text(hjust=0.5))

# Testing impact of gender
#
# t.test(x=my_data[(my_data$dlq_e_2 == "Male" | my_data$dlq_f_2 == "Homme"),]$lt_l1_total, 
#        y=my_data[(my_data$dlq_e_2 == "Female" | my_data$dlq_f_2 == "Femme"),]$lt_l1_total, 
#        alternative="two.sided", paired=FALSE)
# 
# t.test(x=my_data[(my_data$dlq_e_2 == "Male" | my_data$dlq_f_2 == "Homme"),]$ca_c_l1_total, 
#        y=my_data[(my_data$dlq_e_2 == "Female" | my_data$dlq_f_2 == "Femme"),]$ca_c_l1_total, 
#        alternative="two.sided", paired=FALSE)
# 
# t.test(x=my_data[(my_data$dlq_e_2 == "Male" | my_data$dlq_f_2 == "Homme"),]$sm_l1_total, 
#        y=my_data[(my_data$dlq_e_2 == "Female" | my_data$dlq_f_2 == "Femme"),]$sm_l1_total, 
#        alternative="two.sided", paired=FALSE)
# 
# t.test(x=my_data[(my_data$dlq_e_2 == "Male" | my_data$dlq_f_2 == "Homme"),]$p_l1_total, 
#        y=my_data[(my_data$dlq_e_2 == "Female" | my_data$dlq_f_2 == "Femme"),]$p_l1_total, 
#        alternative="two.sided", paired=FALSE)
# 
# df_male   <- my_data[(my_data$dlq_e_2 == "Male" | my_data$dlq_f_2 == "Homme"),]
# df_female <- my_data[(my_data$dlq_e_2 == "Female" | my_data$dlq_f_2 == "Femme"),]
# 
# cor.test(df_male$lt_l1_total, df_male$mss_b_total)
# cor.test(df_female$lt_l1_total, df_female$mss_b_total)
# cor.test(df_male$ca_c_l1_total, df_male$mss_b_total)
# cor.test(df_female$ca_c_l1_total, df_female$mss_b_total)
# cor.test(df_male$sm_l1_total, df_male$mss_b_total)
# cor.test(df_female$sm_l1_total, df_female$mss_b_total)
# cor.test(df_male$p_l1_total, df_male$mss_b_total)
# cor.test(df_female$p_l1_total, df_female$mss_b_total)


# Compare L1 vs L2 performance on each linguistic task ----

t.test(x=my_data$lt_l1_total, y=my_data$lt_l2_total, alternative="two.sided", paired=TRUE)
t.test(x=my_data$ca_c_l1_total, y=my_data$ca_c_l2_total, alternative="two.sided", paired=TRUE)
t.test(x=my_data$sm_l1_total, y=my_data$sm_l2_total, alternative="two.sided", paired=TRUE)
t.test(x=my_data$p_l1_total, y=my_data$p_l2_total, alternative="two.sided", paired=TRUE)

# Doing analyses that take into account L2 proficiency ----

# Correlations between L2 proficiency and performance on the linguistic tasks in L2
cor.test(my_data$l2_prof_total, my_data$lt_l2_total)
cor.test(my_data$l2_prof_total, my_data$ca_c_l2_total)
cor.test(my_data$l2_prof_total, my_data$sm_l2_total)
cor.test(my_data$l2_prof_total, my_data$p_l2_total)
# Correlations between L2 proficiency and difference in L1 and L2 scores on the linguistic tasks
cor.test(my_data$l2_prof_total, my_data$lt_l1_l2_difference)
cor.test(my_data$l2_prof_total, my_data$ca_c_l1_l2_difference)
cor.test(my_data$l2_prof_total, my_data$sm_l1_l2_difference)
cor.test(my_data$l2_prof_total, my_data$p_l1_l2_difference)

# Attempting to predict L1/L2 differences in scores using L2 proficiency and symptoms
lm_lt   <- lm(data = my_data, formula = lt_l1_l2_difference ~ mss_b_total + asq_total + l2_prof_total)
lm_ca_c <- lm(data = my_data, formula = ca_c_l1_l2_difference ~ mss_b_total + asq_total + l2_prof_total)
lm_sm   <- lm(data = my_data, formula = sm_l1_l2_difference ~ mss_b_total + asq_total + l2_prof_total)
lm_p    <- lm(data = my_data, formula = p_l1_l2_difference ~ mss_b_total + asq_total + l2_prof_total)
summary(lm_lt)
summary(lm_ca_c)
summary(lm_sm)
summary(lm_p)

# Attempting to predict L1/L2 differences in scores from symptoms in a low-L2-proficiency subset of participants
lm2_lt   <- lm(data = my_data, subset = (l2_prof_total <= median(my_data$l2_prof_total, na.rm = TRUE)),
                                         formula = lt_l1_l2_difference ~ mss_b_total + asq_total)
lm2_ca_c <- lm(data = my_data, subset = (l2_prof_total <= median(my_data$l2_prof_total, na.rm = TRUE)),
                                         formula = ca_c_l1_l2_difference ~ mss_b_total + asq_total)
lm2_sm   <- lm(data = my_data, subset = (l2_prof_total <= median(my_data$l2_prof_total, na.rm = TRUE)),
                                         formula = sm_l1_l2_difference ~ mss_b_total + asq_total)
lm2_p    <- lm(data = my_data, subset = (l2_prof_total <= median(my_data$l2_prof_total, na.rm = TRUE)),
                                         formula = p_l1_l2_difference ~ mss_b_total + asq_total)
summary(lm2_lt)
summary(lm2_ca_c)
summary(lm2_sm)
summary(lm2_p)

# # Attempting (unsuccessfully so far) to create a gallery of eight scatterplots ----
#
# symptom_task_cor_df <- data.frame(matrix(NA, nrow = 2, ncol = 4))
# for (my_symptom in 1:length(symptoms))
# {
#   for (my_task in 1:length(tasks_l1_l2_difference))
#   {
#     symptom_index = which(colnames(my_data) == symptoms[my_symptom])
#     task_index    = which(colnames(my_data) == tasks_l1_l2_difference[my_task])
# 
#     # cat("my_symptom: ", symptoms[my_symptom], "my_task: ", tasks_l1_l2_difference[my_task], "\n")
#     # cat("symptom_index: ", symptom_index, "task_index: ", task_index, "\n")
#     result <- cor(my_data[symptom_index], my_data[task_index])
#     print(result)
#     symptom_task_cor_df[my_symptom, my_task] <- result
#   }
# }

# Create data frame with the main variables for analysis ----
my_data_main_variables <- subset(my_data, select = c(mss_b_total:l2_prof_total))
# mss-b = schizotypy
# asq = autistic load
# lt = LexTale
# ca_c = Camel and Cactus
# sm = Syntactic Modification
# p = Pragmatics
# dsct = non-verbal processing speed
# sst = non-verbal working memory
# fhat = non-verbal theory of mind
# last five columns are self-rated L2 proficiency; the last column is the total of the previous four


# Calculate column means and st deviations ----
# Keep at end so they aren't included in statistical tests as participants
for (my_col in 1:ncol(my_data))
{
  if(sum(is.character(my_data[,my_col])) == 0)
  {
    my_data["mean",my_col] <- mean(my_data[,my_col], na.rm=TRUE)
    my_data["sd",  my_col] <- sd  (my_data[,my_col], na.rm=TRUE)
  }
}



