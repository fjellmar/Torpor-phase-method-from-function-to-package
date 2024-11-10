#### Packages needed ####
library(dplyr)
library(testthat)
#library(tidyr)

# Libraries needed for plotting (after running the function)
library(ggplot2)
library(gridExtra)
library(ggplotify)


### R-CODES FOR DETERMINING TORPOR ENTRIES AND EXITS FROM SKIN- OR BODY-TEMPERATURE DATA ###

# This script contains descriptions and codes for applying the methods presented in the JTB article by
# Fjelldal et al. on datasets containing Tb- or Tskin-data of animals employing torpor. The function in
# this script simply adds data-columns to the dataset it is being used on, hopefully making the function
# output easy to interpret.

### PREPARING THE DATASET ###

# 1. Upload the dataset and name it "data_testing"
# 2. The dataset can have as many rows and columns as is preferred, but
#    it needs to contain the following two columns: one with the unique ID
#    of each animal (even if the dataset contains data from only one individual)
#    and one column containing continous skin- or body-temperature measures. The
#    ID-column should be named "ID" and the temperature column "Tb". Tb-data
#    for different individuals should be kept in the same column as opposed to separate
#    columns for each individual. Missing Tb-datapoints are allowed.
#
#   NOTE: Time and date columns are not needed (but can be kept in the dataset) as the codes assume
#   that within the ID-group each datapoint following the next have equal distance for all datapoints
#   (for the dataset used when creating this code the Tskin was recorded every 10 minutes).
#   If this is not true for the dataset, e.g. if an individual was measured on two separate
#   days without continuity, the ID-column should be unique for each individual AND period.
#   This is due to the lag-functions (e.g. line 107-108), which calculates the Tb-differences from
#   previous points and need the specification of staying within the group to not mix values
#   between individuals or periods.
#   Missing data periods can also be filled in with NA values in the Tb column so that each individual
#   keeps a consecutive data period from the first observation until the last (there were cases of
#   several days of missing data for some individuals in the dataset for which this code was built).


interval_length <- 10 # (The time between datapoints given in minutes, only necessary for graphing and analyses purposes after running the function)

# To run the function it is necessary to specify an initial torpor threshold value to be used in the first detection of torpor.
# For our data we used the torpor onset calculation presented by Willis (2007) and determined species-specific torpor onset values,
# but any method for deciding on a threshold value is possible to use.
# As some datasets might only need a single torpor onset value to be used across the whole dataset, while others have several
# values (species-dependent or varying across time and/or individuals) we have specified three alternatives for including such torpor
# onset values in this code.

# CHOOSE ONE OF THE FOLLOWING OPTIONS

### ----- OPTION 1 ----- ###

# If the dataset contains data from only one species, or only one single torpor threshold should be used throughout the dataset,
# specify the single torpor onset threshold with the following code.

data_testing$T_onset <- 30.1

### ----- OPTION 2 ----- ###

# If there is more than one species in the dataset, specify the species-specific torpor onset threshold (here shown as example with 3 species).
# To run this code, the dataset needs to contain an additional column named "Species". The species-info in the following code can then
# be updated to fit the data. If only two species, remove line 74 (or 75). If more than three species, replicate line 74 (or 75).

T_onset_threshold_Plecotus <- 30.1
T_onset_threshold_Enil <- 30.1
T_onset_threshold_Brandtii <- 29.9

data_testing$T_onset <- NA

data_testing <- data_testing %>%
  group_by(Species) %>%
  mutate(T_onset = ifelse(Species == "Plecotus", T_onset_threshold_Plecotus, T_onset)) %>%
  mutate(T_onset = ifelse(Species == "Enil", T_onset_threshold_Enil, T_onset)) %>%
  mutate(T_onset = ifelse(Species == "Brandtii", T_onset_threshold_Brandtii, T_onset))


### ----- OPTION 3 ----- ###

# For varying torpor onset thresholds (between days and/or individuals) the values should be included before uploading the dataset.
# The column containing the thresholds should be named "T_onset". See Supplementary Materials 2 for example.


# CONTINUE FROM HERE AFTER PERFORMING OPTION 1 OR 2 OR 3

data_testing$below_threshold <- data_testing$Tb < data_testing$T_onset

# Because field-data often contains missing data-points, the function presented in this script offers the possibility to
# specify a value for the number of consecutive NA-values that will be allowed into the determination of euthermic and torpor periods.
# This does not concern the phase-determination as phases with missing values will be assigned to the
# "Not complete" category and can thus be excluded from analysis.
# The "na_allow" (below) only applies for missing Tb-values either during stable torpor or during euthermic periods,
# and thus assigns the missing values as either "Torpor" or "Euthermic" if the consecutive number of missing values =< na_allow,
# AND ONLY if both the previous and following Tb-value recorded for the individual is below the T_onset (for torpor) or
# above (for euthermic periods). Code-details regarding this can be found on line 538-542 and 694-698.

na_allow <- 4

#### RUNNING THE FUNCTION ####

#torpor_phase_function <- function(sensitivity, within_bout_sensitivity) {
testdata <- readr::read_csv("testdata.csv")
testdata <- janitor::clean_names(testdata)
names(testdata)[names(testdata) == 'tb'] <- 'body_temp'
names(testdata)[names(testdata) == 't_onset'] <- 'torpor_onset'

#### 1. function --- checking that the dataset has the necessary columns needed and that temperature values are logic ####  

check_dataset <- function(.data) {
  if(!any(names(.data)=="body_temp")) {
    stop("column body_temp is missing")  # checks if the dataset has the column-name "body_temp" included
  }
  if(!any(names(.data)=="id")) {
    stop("column id is missing")  # checks if the dataset has the column-name "id" included
  }
  if(!any(names(.data)=="torpor_onset")) {
    stop("column torpor_onset is missing")  # checks if the dataset has the column-name "torpor_onset" included
  }
  if(max(.data$body_temp, na.rm=T) > 60) {
    stop("body temperature values above 60 degrees detected - check values") # checks if the temperature values are biologically sound
  }
  if(min(.data$body_temp, na.rm=T) < -20) {
    stop("body temperature values below -20 degrees detected - check values") # checks if the temperature values are biologically sound
  }
}

#### 2. function --- add three columns: body temperature with one or two lags, and column for whether body temperature is below threshold ####

add_temp_lag_cols <- function(.data) {
  test |>
    dplyr::group_by(id) |>
    dplyr::mutate(
      body_temp_diff = (body_temp - lag(body_temp, n = 1, default = NA)), # adds column with lag 1 temperatures
      body_temp_diff_2 = (body_temp - lag(body_temp, n = 2, default = NA)), # adds column with lag 2 temperatures
      below_threshold = ifelse(is.na(body_temp < torpor_onset), "No data", body_temp < torpor_onset) # makes new column with TRUE/FALSE based on whether body_temp is below or above the threshold, or no data for missing values
    )
}  

#### 3. function --- adding three columns: id of each section of unique below_threshold, length of each such section, and whether bat is torpid or not ####

myrleid <- function(x) { # small function for numbering each section
  x <- rle(x)$lengths
  rep(seq_along(x), times=x)
}

add_section_numbers <- function(.data) { # function for adding the three columns to the dataset
   .data |>
    dplyr::mutate(
      num = myrleid(below_threshold), # ading column with number-id of each section
      count_length = with(rle(below_threshold), rep(lengths, lengths)), # length of each section 
      torpor_col = ifelse(below_threshold == TRUE, "torpor", "not torpor") # adding column for torpor/non torpor
      )
}

#### 4. function --- adding the column phase_start for indicating the start of torpor entries and torpor exits ####

add_phase_start <- function(.data) {
.data |>
  dplyr::mutate(phase_start = case_when(
    below_threshold == F & lag(below_threshold == T) ~ "exiting",
    below_threshold == T & lag(below_threshold == F) ~ "entering",
    TRUE ~ "not phasestart"
  ))
}

#### 5. function --- ####
  my_vec <- ifelse(is.na(my_vec), "No", my_vec)
  torpor_dataset$phase_start <- my_vec

  Entering_rows <- which(torpor_dataset$phase_start == "Entering")
  Exiting_rows <- which(torpor_dataset$phase_start == "Exiting")

  my_vec_entering_1 <- numeric()
  my_vec_entering_1_notfull <- numeric()

  for (i in Entering_rows) {
    while (i < nrow(torpor_dataset)) {
      my_vec_entering_1[i] <- "Entering"
      # increment number by 1
      i <- i + 1

      if (!is.na(torpor_dataset[i, c("Tb_diff")]) & is.na(torpor_dataset[i + 1, c("Tb_diff")])) {
        if (torpor_dataset[i, c("Tb_diff")] > -sensitivity) {
          break
        } else {
          my_vec_entering_1[i] <- ifelse(torpor_dataset[i, c("Tb_diff")] <= -sensitivity, "Entering", NA)
          my_vec_entering_1_notfull[i] <- "Not complete"
        }
      } else {
        # break if following conditions are met
        if ((torpor_dataset[i, c("Tb_diff")] > -sensitivity & torpor_dataset[i + 1, c("Tb_diff")] > -sensitivity) |
          (torpor_dataset[i, c("Tb_diff")] > -sensitivity & torpor_dataset[i + 1, c("Tb_diff")] <= -sensitivity & torpor_dataset[i + 1, c("Tb_diff_2")] >= 0) |
          is.na(torpor_dataset[i, c("Tb_diff")])) {
          break
        }
      }
    }
  }

  my_vec_entering_2 <- numeric()
  my_vec_entering_2_notfull <- numeric()

  for (i in Entering_rows) {
    while (i <= nrow(torpor_dataset)) {
      my_vec_entering_2[i] <- "Entering"
      i <- i - 1

      if (!is.na(torpor_dataset[i, c("Tb_diff")]) & is.na(torpor_dataset[i - 1, c("Tb_diff")])) {
        if (torpor_dataset[i, c("Tb_diff")] > -sensitivity) {
          break
        } else {
          my_vec_entering_2[i] <- ifelse(torpor_dataset[i, c("Tb_diff")] <= -sensitivity, "Entering", NA)
          my_vec_entering_2_notfull[i] <- "Not complete"
        }
      } else {
        # break if following conditions are met
        if ((torpor_dataset[i, c("Tb_diff")] > -sensitivity & torpor_dataset[i - 1, c("Tb_diff")] > -sensitivity) |
          (torpor_dataset[i, c("Tb_diff")] > -sensitivity & torpor_dataset[i - 1, c("Tb_diff")] <= -sensitivity & torpor_dataset[i, c("Tb_diff_2")] >= 0) |
          is.na(torpor_dataset[i, c("Tb_diff")]) |
          (i == 1)) {
          break
        }
      }
    }
  }

  my_vec_exiting_1 <- numeric()
  my_vec_exiting_1_notfull <- numeric()

  for (i in Exiting_rows) {
    while (i < nrow(torpor_dataset)) {
      my_vec_exiting_1[i] <- "Exiting"
      i <- i + 1

      if (!is.na(torpor_dataset[i, c("Tb_diff")]) & is.na(torpor_dataset[i + 1, c("Tb_diff")])) {
        if (torpor_dataset[i, c("Tb_diff")] < sensitivity) {
          break
        } else {
          my_vec_exiting_1[i] <- ifelse(torpor_dataset[i, c("Tb_diff")] >= sensitivity, "Exiting", NA)
          my_vec_exiting_1_notfull[i] <- "Not complete"
        }
      } else {
        # break if following conditions are met
        if ((torpor_dataset[i, c("Tb_diff")] < sensitivity & torpor_dataset[i + 1, c("Tb_diff")] < sensitivity) |
          (torpor_dataset[i, c("Tb_diff")] < sensitivity & torpor_dataset[i + 1, c("Tb_diff")] >= sensitivity & torpor_dataset[i + 1, c("Tb_diff_2")] <= 0) |
          is.na(torpor_dataset[i, c("Tb_diff")])) {
          break
        }
      }
    }
  }

  my_vec_exiting_2 <- numeric()
  my_vec_exiting_2_notfull <- numeric()

  for (i in Exiting_rows) {
    while (i <= nrow(torpor_dataset)) {
      my_vec_exiting_2[i] <- "Exiting"
      i <- i - 1

      if (!is.na(torpor_dataset[i, c("Tb_diff")]) & is.na(torpor_dataset[i - 1, c("Tb_diff")])) {
        if (torpor_dataset[i, c("Tb_diff")] < sensitivity) {
          break
        } else {
          my_vec_exiting_2[i] <- ifelse(torpor_dataset[i, c("Tb_diff")] >= sensitivity, "Exiting", NA)
          my_vec_exiting_2_notfull[i] <- "Not complete"
        }
      } else {
        # break if following conditions are met
        if ((torpor_dataset[i, c("Tb_diff")] < sensitivity & torpor_dataset[i - 1, c("Tb_diff")] < sensitivity) |
          (torpor_dataset[i, c("Tb_diff")] < sensitivity & torpor_dataset[i - 1, c("Tb_diff")] >= sensitivity & torpor_dataset[i, c("Tb_diff_2")] <= 0) |
          is.na(torpor_dataset[i, c("Tb_diff")]) |
          (i == 1)) {
          break
        }
      }
    }
  }

  my_vec_entering_1 <- c(my_vec_entering_1, rep(NA, nrow(torpor_dataset) - length(my_vec_entering_1)))
  my_vec_entering_2 <- c(my_vec_entering_2, rep(NA, nrow(torpor_dataset) - length(my_vec_entering_2)))
  my_vec_exiting_1 <- c(my_vec_exiting_1, rep(NA, nrow(torpor_dataset) - length(my_vec_exiting_1)))
  my_vec_exiting_2 <- c(my_vec_exiting_2, rep(NA, nrow(torpor_dataset) - length(my_vec_exiting_2)))

  my_vec_entering_1_notfull <- c(my_vec_entering_1_notfull, rep(NA, nrow(torpor_dataset) - length(my_vec_entering_1_notfull)))
  my_vec_entering_2_notfull <- c(my_vec_entering_2_notfull, rep(NA, nrow(torpor_dataset) - length(my_vec_entering_2_notfull)))
  my_vec_exiting_1_notfull <- c(my_vec_exiting_1_notfull, rep(NA, nrow(torpor_dataset) - length(my_vec_exiting_1_notfull)))
  my_vec_exiting_2_notfull <- c(my_vec_exiting_2_notfull, rep(NA, nrow(torpor_dataset) - length(my_vec_exiting_2_notfull)))

  defaultW <- getOption("warn")
  options(warn = -1)

  phase_all1 <- ifelse(is.na(my_vec_entering_1) & !is.na(my_vec_entering_2), my_vec_entering_2, my_vec_entering_1)
  phase_all1 <- ifelse(is.na(phase_all1) & !is.na(my_vec_exiting_1), my_vec_exiting_1, phase_all1)
  phase_all1 <- ifelse(is.na(phase_all1) & !is.na(my_vec_exiting_2), my_vec_exiting_2, phase_all1)

  Not_full_phases <- ifelse(is.na(my_vec_entering_1_notfull) & !is.na(my_vec_entering_2_notfull), my_vec_entering_2_notfull, my_vec_entering_1_notfull)
  Not_full_phases <- ifelse(is.na(Not_full_phases) & !is.na(my_vec_exiting_1_notfull), my_vec_exiting_1_notfull, Not_full_phases)
  Not_full_phases <- ifelse(is.na(Not_full_phases) & !is.na(my_vec_exiting_2_notfull), my_vec_exiting_2_notfull, Not_full_phases)

  options(warn = defaultW)

  torpor_dataset$phase_all1 <- phase_all1
  torpor_dataset$Not_full_phases <- Not_full_phases

  # Including the within torpor bout entries and rewarmings #

  within_bout_entering_rows <- which(is.na(torpor_dataset$phase_all1) & torpor_dataset$Torpor == "Torpor" & torpor_dataset$Tb_diff < -within_bout_sensitivity)
  within_bout_exiting_rows <- which(is.na(torpor_dataset$phase_all1) & torpor_dataset$Torpor == "Torpor" & torpor_dataset$Tb_diff > within_bout_sensitivity)

  my_vec_entering_within_1 <- numeric()
  my_vec_entering_within_1_notfull <- numeric()

  for (i in within_bout_entering_rows) {
    while (i < nrow(torpor_dataset)) {
      my_vec_entering_within_1[i] <- "Entering_within"
      i <- i + 1

      if (!is.na(torpor_dataset[i, c("Tb_diff")]) & is.na(torpor_dataset[i + 1, c("Tb_diff")])) {
        if (torpor_dataset[i, c("Tb_diff")] > -sensitivity) {
          break
        } else {
          my_vec_entering_within_1[i] <- ifelse(torpor_dataset[i, c("Tb_diff")] <= -sensitivity, "Entering_within", NA)
          my_vec_entering_within_1_notfull[i] <- "Not complete"
        }
      } else {
        # break if following conditions are met
        if ((torpor_dataset[i, c("Tb_diff")] > -sensitivity & torpor_dataset[i + 1, c("Tb_diff")] > -sensitivity) |
          (torpor_dataset[i, c("Tb_diff")] > -sensitivity & torpor_dataset[i + 1, c("Tb_diff")] <= -sensitivity & torpor_dataset[i + 1, c("Tb_diff_2")] >= 0) |
          is.na(torpor_dataset[i, c("Tb_diff")])) {
          break
        }
      }
    }
  }

  my_vec_entering_within_2 <- numeric()
  my_vec_entering_within_2_notfull <- numeric()
  phase2 <- ifelse(is.na(torpor_dataset$phase_all1), 0, torpor_dataset$phase_all1)

  for (i in within_bout_entering_rows) {
    while (i <= nrow(torpor_dataset)) {
      my_vec_entering_within_2[i] <- "Entering_within"
      i <- i - 1

      if (!is.na(torpor_dataset[i, c("Tb_diff")]) & is.na(torpor_dataset[i - 1, c("Tb_diff")])) {
        if (torpor_dataset[i, c("Tb_diff")] > -sensitivity) {
          break
        } else {
          my_vec_entering_within_2[i] <- ifelse(torpor_dataset[i, c("Tb_diff")] <= -sensitivity, "Entering_within", NA)
          my_vec_entering_within_2_notfull[i] <- "Not complete"
        }
      } else {
        # break if following conditions are met
        if ((torpor_dataset[i, c("Tb_diff")] > -sensitivity & torpor_dataset[i - 1, c("Tb_diff")] > -sensitivity) |
          (torpor_dataset[i, c("Tb_diff")] > -sensitivity & torpor_dataset[i - 1, c("Tb_diff")] <= -sensitivity & torpor_dataset[i, c("Tb_diff_2")] >= 0) |
          (torpor_dataset[i, c("Tb_diff")] > -sensitivity & phase2[i - 1] == "Entering") |
          is.na(torpor_dataset[i, c("Tb_diff")]) |
          (i == 1)) {
          break
        }
      }
    }
  }

  my_vec_exiting_within_1 <- numeric()
  my_vec_exiting_within_1_notfull <- numeric()

  for (i in within_bout_exiting_rows) {
    while (i < nrow(torpor_dataset)) {
      my_vec_exiting_within_1[i] <- "Exiting_within"
      i <- i + 1

      if (!is.na(torpor_dataset[i, c("Tb_diff")]) & is.na(torpor_dataset[i + 1, c("Tb_diff")])) {
        if (torpor_dataset[i, c("Tb_diff")] < sensitivity) {
          break
        } else {
          my_vec_exiting_within_1[i] <- ifelse(torpor_dataset[i, c("Tb_diff")] >= sensitivity, "Exiting_within", NA)
          my_vec_exiting_within_1_notfull[i] <- "Not complete"
        }
      } else {
        # break if following conditions are met
        if ((torpor_dataset[i, c("Tb_diff")] < sensitivity & torpor_dataset[i + 1, c("Tb_diff")] < sensitivity) |
          (torpor_dataset[i, c("Tb_diff")] < sensitivity & torpor_dataset[i + 1, c("Tb_diff")] >= sensitivity & torpor_dataset[i + 1, c("Tb_diff_2")] <= 0) |
          is.na(torpor_dataset[i, c("Tb_diff")]) |
          (torpor_dataset[i, c("Tb_diff")] < sensitivity & phase2[i + 1] == "Exiting") |
          is.na(torpor_dataset[i + 1, c("Tb_diff")])) {
          break
        }
      }
    }
  }

  my_vec_exiting_within_2 <- numeric()
  my_vec_exiting_within_2_notfull <- numeric()

  for (i in within_bout_exiting_rows) {
    while (i <= nrow(torpor_dataset)) {
      my_vec_exiting_within_2[i] <- "Exiting_within"
      i <- i - 1

      if (!is.na(torpor_dataset[i, c("Tb_diff")]) & is.na(torpor_dataset[i - 1, c("Tb_diff")])) {
        if (torpor_dataset[i, c("Tb_diff")] < sensitivity) {
          break
        } else {
          my_vec_exiting_within_2[i] <- ifelse(torpor_dataset[i, c("Tb_diff")] >= sensitivity, "Exiting_within", NA)
          my_vec_exiting_within_2_notfull[i] <- "Not complete"
        }
      } else {
        # break if following conditions are met
        if ((torpor_dataset[i, c("Tb_diff")] < sensitivity & torpor_dataset[i - 1, c("Tb_diff")] < sensitivity) |
          (torpor_dataset[i, c("Tb_diff")] < sensitivity & torpor_dataset[i - 1, c("Tb_diff")] >= sensitivity & torpor_dataset[i, c("Tb_diff_2")] <= 0) |
          is.na(torpor_dataset[i, c("Tb_diff")]) |
          is.na(torpor_dataset[i - 1, c("Tb_diff")]) |
          (i == 1)) {
          break
        }
      }
    }
  }

  my_vec_entering_within_1 <- c(my_vec_entering_within_1, rep(NA, nrow(torpor_dataset) - length(my_vec_entering_within_1)))
  my_vec_entering_within_2 <- c(my_vec_entering_within_2, rep(NA, nrow(torpor_dataset) - length(my_vec_entering_within_2)))
  my_vec_exiting_within_1 <- c(my_vec_exiting_within_1, rep(NA, nrow(torpor_dataset) - length(my_vec_exiting_within_1)))
  my_vec_exiting_within_2 <- c(my_vec_exiting_within_2, rep(NA, nrow(torpor_dataset) - length(my_vec_exiting_within_2)))

  my_vec_entering_within_1_notfull <- c(my_vec_entering_within_1_notfull, rep(NA, nrow(torpor_dataset) - length(my_vec_entering_within_1_notfull)))
  my_vec_entering_within_2_notfull <- c(my_vec_entering_within_2_notfull, rep(NA, nrow(torpor_dataset) - length(my_vec_entering_within_2_notfull)))
  my_vec_exiting_within_1_notfull <- c(my_vec_exiting_within_1_notfull, rep(NA, nrow(torpor_dataset) - length(my_vec_exiting_within_1_notfull)))
  my_vec_exiting_within_2_notfull <- c(my_vec_exiting_within_2_notfull, rep(NA, nrow(torpor_dataset) - length(my_vec_exiting_within_2_notfull)))

  phase_all <- ifelse(is.na(my_vec_entering_1) & !is.na(my_vec_entering_2), my_vec_entering_2, my_vec_entering_1)
  phase_all <- ifelse(is.na(phase_all) & !is.na(my_vec_exiting_1), my_vec_exiting_1, phase_all)
  phase_all <- ifelse(is.na(phase_all) & !is.na(my_vec_exiting_2), my_vec_exiting_2, phase_all)
  phase_all <- ifelse(is.na(phase_all) & !is.na(my_vec_entering_within_1), my_vec_entering_within_1, phase_all)
  phase_all <- ifelse(is.na(phase_all) & !is.na(my_vec_entering_within_2), my_vec_entering_within_2, phase_all)
  phase_all <- ifelse(is.na(phase_all) & !is.na(my_vec_exiting_within_1), my_vec_exiting_within_1, phase_all)
  phase_all <- ifelse(is.na(phase_all) & !is.na(my_vec_exiting_within_2), my_vec_exiting_within_2, phase_all)

  Not_full_phases2 <- ifelse(is.na(my_vec_entering_1_notfull) & !is.na(my_vec_entering_2_notfull), my_vec_entering_2_notfull, my_vec_entering_1_notfull)
  Not_full_phases2 <- ifelse(is.na(Not_full_phases2) & !is.na(my_vec_exiting_1_notfull), my_vec_exiting_1_notfull, Not_full_phases2)
  Not_full_phases2 <- ifelse(is.na(Not_full_phases2) & !is.na(my_vec_exiting_2_notfull), my_vec_exiting_2_notfull, Not_full_phases2)
  Not_full_phases2 <- ifelse(is.na(Not_full_phases2) & !is.na(my_vec_entering_within_1_notfull), my_vec_entering_within_1_notfull, Not_full_phases2)
  Not_full_phases2 <- ifelse(is.na(Not_full_phases2) & !is.na(my_vec_entering_within_2_notfull), my_vec_entering_within_2_notfull, Not_full_phases2)
  Not_full_phases2 <- ifelse(is.na(Not_full_phases2) & !is.na(my_vec_exiting_within_1_notfull), my_vec_exiting_within_1_notfull, Not_full_phases2)
  Not_full_phases2 <- ifelse(is.na(Not_full_phases2) & !is.na(my_vec_exiting_within_2_notfull), my_vec_exiting_within_2_notfull, Not_full_phases2)

  torpor_dataset$phase_all <- phase_all
  torpor_dataset$Not_full_phases <- ifelse(is.na(torpor_dataset$Not_full_phases) & !is.na(Not_full_phases2), Not_full_phases2, torpor_dataset$Not_full_phases)

  # Label the different phases as "Complete", "Mixed", "Within bout" or "Not complete" #

  setDT(torpor_dataset)
  torpor_dataset <- as.data.frame(torpor_dataset[, phase_ID := rleid(phase_all)])

  torpor_dataset <- torpor_dataset %>%
    group_by(phase_ID) %>%
    mutate(phase_type = ifelse(phase_all == "Exiting" | phase_all == "Entering",
      "Complete", NA
    )) %>%
    mutate(phase_type = ifelse((phase_all == "Exiting" & any(Tb_diff < 0)) | (phase_all == "Entering" & any(Tb_diff > 0)),
      "Mixed", phase_type
    )) %>%
    mutate(phase_type = ifelse(phase_all == "Exiting_within" | phase_all == "Entering_within",
      "Within_bout", phase_type
    )) %>%
    mutate(phase_type = ifelse(any(!is.na(Not_full_phases)),
      "Not complete", phase_type
    ))


  my_vec2 <- numeric()

  for (i in 1:(nrow(torpor_dataset))) {
    if (i == 1) {
      my_vec2[i] <- NA
    } else if (i == nrow(torpor_dataset)) {
      my_vec2[i] <- NA
    } else {
      if ((!is.na(torpor_dataset[c(i), c("phase_all")]) & torpor_dataset[c(i), c("phase_type")] != "Not complete" & is.na(torpor_dataset[c(i + 1), c("Tb")])) |
        (!is.na(torpor_dataset[c(i), c("phase_all")]) & torpor_dataset[c(i), c("phase_type")] != "Not complete" & is.na(torpor_dataset[c(i - 1), c("Tb")]))) {
        my_vec2[i] <- "Not complete"
      } else {
        my_vec2[i] <- NA
      }
    }
  }

  torpor_dataset$my_vec2 <- my_vec2

  torpor_dataset <- torpor_dataset %>%
    group_by(phase_ID) %>%
    mutate(phase_type = ifelse(any(!is.na(my_vec2)),
      "Not complete", phase_type
    ))

  torpor_dataset$phase_ID <- ifelse(is.na(torpor_dataset$phase_all), NA, torpor_dataset$phase_ID)

  # Number each full torpor bout with corresponding torpor phases

  torpor_dataset <- torpor_dataset %>%
    group_by(ID) %>%
    mutate(last_noNA_Tb = Tb) %>%
    mutate(next_noNA_Tb = Tb) %>%
    fill(last_noNA_Tb, .direction = "down") %>%
    fill(next_noNA_Tb, .direction = "up")

  torpor_dataset$Torpor <- ifelse(torpor_dataset$below_threshold == "No data" &
    torpor_dataset$count_length <= na_allow & # Allowing up to x missing values within torpor bout to be identified as part of the bout
    torpor_dataset$last_noNA_Tb < torpor_dataset$T_onset &
    torpor_dataset$next_noNA_Tb < torpor_dataset$T_onset,
  "Torpor", torpor_dataset$Torpor
  )

  Entering_rows <- which(torpor_dataset$phase_start == "Entering")

  phase_all_new <- ifelse(is.na(torpor_dataset$phase_all), "none", torpor_dataset$phase_all)
  my_vec_full_num <- numeric()

  for (i in Entering_rows) {
    while (i < nrow(torpor_dataset)) {
      my_vec_full_num[i] <- "Torpor"
      i <- i + 1

      if (is.na(torpor_dataset[i, c("Torpor")]) |
        (torpor_dataset[i, c("Torpor")] == "Not torpor" & phase_all_new[i] == "none")) {
        break
      }
    }
  }


  my_vec_full_num2 <- numeric()

  for (i in Entering_rows) {
    while (i <= nrow(torpor_dataset)) {
      my_vec_full_num2[i] <- "Torpor"
      i <- i - 1

      if ((is.na(torpor_dataset[i, c("Torpor")])) |
        (torpor_dataset[i, c("Torpor")] == "Not torpor" & phase_all_new[i] == "none") |
        (phase_all_new[i] == "Exiting") |
        (i == 1)) {
        break
      }
    }
  }

  my_vec_full_num <- c(my_vec_full_num, rep(NA, nrow(torpor_dataset) - length(my_vec_full_num)))
  my_vec_full_num2 <- c(my_vec_full_num2, rep(NA, nrow(torpor_dataset) - length(my_vec_full_num2)))
  torpor2 <- ifelse(is.na(my_vec_full_num) & !is.na(my_vec_full_num2), my_vec_full_num2, my_vec_full_num)
  torpor_dataset$torpor2 <- torpor2

  setDT(torpor_dataset)
  torpor_dataset <- as.data.frame(torpor_dataset[, bout_ID := rleid(torpor2)])

  # Identifying back-to-back torpor bouts (repeating this step 4 times for up to 4 cases of back-to-back torpor bouts)

  torpor_dataset <- torpor_dataset %>%
    group_by(ID) %>%
    mutate(last_phase = lag(phase_all, n = 1, default = NA))

  double_bout_rows <- which(torpor_dataset$phase_all == "Entering" & torpor_dataset$last_phase == "Exiting")

  double_bout_vec <- numeric()
  double_bout_vec2 <- torpor_dataset$bout_ID

  for (i in double_bout_rows) {
    while (i < nrow(torpor_dataset)) {
      double_bout_vec[i] <- double_bout_vec2[i] + 0.1
      i <- i + 1

      # break if following conditions are met
      if ((double_bout_vec2[i] > double_bout_vec2[i - 1]) |
        is.na(torpor_dataset[i, c("bout_ID")]) |
        is.na(torpor_dataset[i + 1, c("bout_ID")])) {
        break
      }
    }
  }

  double_bout_vec <- c(double_bout_vec, rep(NA, nrow(torpor_dataset) - length(double_bout_vec)))
  torpor_dataset$bout_ID <- ifelse(!is.na(double_bout_vec), double_bout_vec, torpor_dataset$bout_ID)

  double_bout_rows2 <- which(torpor_dataset$phase_all == "Entering" & torpor_dataset$last_phase == "Exiting" &
    torpor_dataset$bout_ID == lag(torpor_dataset$bout_ID))

  double_bout_vec <- numeric()
  double_bout_vec2 <- torpor_dataset$bout_ID

  for (i in double_bout_rows2) {
    while (i < nrow(torpor_dataset)) {
      double_bout_vec[i] <- double_bout_vec2[i] + 0.1
      i <- i + 1

      # break if following conditions are met
      if ((double_bout_vec2[i] > double_bout_vec2[i - 1]) |
        is.na(torpor_dataset[i, c("bout_ID")]) |
        is.na(torpor_dataset[i + 1, c("bout_ID")])) {
        break
      }
    }
  }

  double_bout_vec <- c(double_bout_vec, rep(NA, nrow(torpor_dataset) - length(double_bout_vec)))
  torpor_dataset$bout_ID <- ifelse(!is.na(double_bout_vec), double_bout_vec, torpor_dataset$bout_ID)

  double_bout_rows3 <- which(torpor_dataset$phase_all == "Entering" & torpor_dataset$last_phase == "Exiting" &
    torpor_dataset$bout_ID == lag(torpor_dataset$bout_ID))

  double_bout_vec <- numeric()
  double_bout_vec2 <- torpor_dataset$bout_ID

  for (i in double_bout_rows3) {
    while (i < nrow(torpor_dataset)) {
      double_bout_vec[i] <- double_bout_vec2[i] + 0.1
      i <- i + 1

      # break if following conditions are met
      if ((double_bout_vec2[i] > double_bout_vec2[i - 1]) |
        is.na(torpor_dataset[i, c("bout_ID")]) |
        is.na(torpor_dataset[i + 1, c("bout_ID")])) {
        break
      }
    }
  }

  double_bout_vec <- c(double_bout_vec, rep(NA, nrow(torpor_dataset) - length(double_bout_vec)))
  torpor_dataset$bout_ID <- ifelse(!is.na(double_bout_vec), double_bout_vec, torpor_dataset$bout_ID)

  double_bout_rows4 <- which(torpor_dataset$phase_all == "Entering" & torpor_dataset$last_phase == "Exiting" &
    torpor_dataset$bout_ID == lag(torpor_dataset$bout_ID))

  double_bout_vec <- numeric()
  double_bout_vec2 <- torpor_dataset$bout_ID

  for (i in double_bout_rows4) {
    while (i < nrow(torpor_dataset)) {
      double_bout_vec[i] <- double_bout_vec2[i] + 0.1
      i <- i + 1

      # break if following conditions are met
      if ((double_bout_vec2[i] > double_bout_vec2[i - 1]) |
        is.na(torpor_dataset[i, c("bout_ID")]) |
        is.na(torpor_dataset[i + 1, c("bout_ID")])) {
        break
      }
    }
  }

  double_bout_vec <- c(double_bout_vec, rep(NA, nrow(torpor_dataset) - length(double_bout_vec)))
  torpor_dataset$bout_ID <- ifelse(!is.na(double_bout_vec), double_bout_vec, torpor_dataset$bout_ID)

  df3 <- torpor_dataset %>%
    filter(!is.na(bout_ID) & (phase_all == "Entering" | phase_all == "Exiting")) %>%
    group_by(bout_ID, ID) %>%
    count(phase_all, ID) %>%
    count(bout_ID, ID)

  count_phase <- df3$n[match(torpor_dataset$bout_ID, df3$bout_ID)]

  torpor_dataset$bout_ID <- ifelse(is.na(torpor_dataset$torpor2) | count_phase < 2, NA, torpor_dataset$bout_ID)

  torpor_dataset$state <- torpor_dataset$phase_all
  torpor_dataset$state <- ifelse(is.na(torpor_dataset$state) &
    ((!is.na(torpor_dataset$Tb) & torpor_dataset$Tb >= torpor_dataset$T_onset) |
      (is.na(torpor_dataset$Tb) & torpor_dataset$last_noNA_Tb >= torpor_dataset$T_onset &
        torpor_dataset$next_noNA_Tb >= torpor_dataset$T_onset & torpor_dataset$count_length <= na_allow)), # Allowing up to x missing values to be assigned as part of euthermic period
  "Euthermic", torpor_dataset$state
  )
  torpor_dataset$state <- ifelse(is.na(torpor_dataset$state) & torpor_dataset$Torpor == "Torpor", torpor_dataset$Torpor, torpor_dataset$state)
  torpor_dataset$state <- ifelse(is.na(torpor_dataset$state) & torpor_dataset$below_threshold == TRUE & torpor_dataset$Torpor == "Not torpor" &
    torpor_dataset$last_noNA_Tb >= torpor_dataset$T_onset & torpor_dataset$next_noNA_Tb >= torpor_dataset$T_onset,
  "Euthermic", torpor_dataset$state
  )


  # Clean up the dataset
  torpor_dataset$count_length <- NULL
  torpor_dataset$num <- NULL
  torpor_dataset$Tb_diff_2 <- NULL
  torpor_dataset$phase_start <- NULL
  torpor_dataset$Not_full_phases <- NULL
  torpor_dataset$phase_all1 <- NULL
  torpor_dataset$last_phase <- NULL
  torpor_dataset$last_noNA_Tb <- NULL
  torpor_dataset$next_noNA_Tb <- NULL
  torpor_dataset$my_vec2 <- NULL
  torpor_dataset$torpor2 <- NULL
  torpor_dataset$below_threshold <- NULL
  torpor_dataset$Torpor <- NULL

  return(torpor_dataset)
#}
