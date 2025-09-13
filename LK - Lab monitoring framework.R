#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$                                                                                    $#
#$                         Created by Lyubomir N. Kolev  II                           $#
#$                   SCDHEC, Division of Acute Disease Epidemiology                   $#
#$                                  January 30, 2024                                  $#
#$                                                                                    $#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

print(paste("Starting:", date()))

# installs and loads necessary packages automatically
if(!require(pacman)){install.packages("pacman")}
p_load(tictoc, odbc, magrittr, tidyverse, MMWRweek, RcppRoll)

output_path <-
  "C:/OneDrive/DIM Files"

print(paste("Step 1: Connecting to DWH and downloading data:", date()))
tic("Total runtime:")
tic("Downloading data")

DWH_connection <- 
  dbConnect(
    odbc(),
    dsn = "DPH_warehouse"
  )

labs_raw <-
  rbind(
    covid_raw,
    flu_raw,
    enteric_raw,
    hepA_raw,
    hepB_raw,
    hepC_raw,
    vbd_raw,
    vpd_raw
  )

toc()

print(paste("Step 2: Creating facility list, compute MMWR week and lag between dates:", date()))
tic("Facility list, MWMR, lag time")

labs_clean <-
  labs_raw %>% 
  mutate(
    MMWRweek(SENDING_DATE),
    'Week ending:' = MMWRweek2Date(MMWRyear = MMWRyear, MMWRweek = MMWRweek, MMWRday = 7),
    Facility_clean = coalesce(FACILITY, FACILITY__OTHER),
    Reporting_clean = coalesce(REPORTING_FACILITY, REPORTING_FACILITY_OTHER),
    'Lag between Specimen Date and Sending date' = SENDING_DATE - SPECIMEN_DATE
  ) %>% 
  relocate(MMWRyear:`Week ending:`, .after = EVENT_DATE)
toc()

print(paste("Step 3: Compute max date for sample, individual facility:", date()))
tic("Sample and facility sample dates")
latest_sample_date <-
  max(labs_clean$SENDING_DATE)

max_date_by_facility <-
  labs_clean %>% 
  group_by(Facility_clean, Reporting_clean, DISEASE) %>% 
  summarize(
    'Latest Facility Date' = max(SENDING_DATE)
  )

labs_clean %<>%
  left_join(max_date_by_facility, by = c("Facility_clean", "Reporting_clean", "DISEASE"))
toc()

print(paste("Step 4: Calculate # of days since each facility submitted a lab:", date()))
tic("Days since last submitted lab")
days_since_lastlab <-
  labs_clean %>% 
  group_by(Facility_clean, Reporting_clean, DISEASE) %>% 
  reframe(
    'Days since last submitted lab' = latest_sample_date - `Latest Facility Date`
  ) %>% 
  distinct()

labs_clean %<>%
  left_join(days_since_lastlab, by = c("Facility_clean", "Reporting_clean", "DISEASE"))
toc()

print(paste("Step 5: Compute average lag and SD lag by facility:", date()))
tic("Average lag, SD lag")
average_facility_lag <-
  labs_clean %>% 
  group_by(Facility_clean, Reporting_clean, DISEASE) %>% 
  summarize(
    'Average lag' = mean(`Lag between Specimen Date and Sending date`),
    'SD lag' = sd(`Lag between Specimen Date and Sending date`)
  )

labs_clean %<>%
  left_join(average_facility_lag, by = c("Facility_clean", "Reporting_clean", "DISEASE"))
toc()

print(paste("Step 6: Compute labs by MMWR week for each facility:", date()))
tic("Labs by week")
average_labs_by_MMWR <-
  labs_clean %>% 
  group_by(Facility_clean, Reporting_clean, DISEASE, MMWRweek) %>% 
  summarize(
    'Labs by week' = length(TEST)
  ) %>% 
  distinct()

labs_clean %<>%
  left_join(average_labs_by_MMWR, by = c("Facility_clean", "Reporting_clean", "DISEASE", "MMWRweek"))
toc()

print(paste("Step 7: Compute mean sample lag, rolling averages, volume and submission time flags:", date()))
tic("Rolling averages, volume, submission flags")
mean_sample_lag <-
  as.numeric(mean(labs_clean$`Lag between Specimen Date and Sending date`, na.rm = TRUE))

median_sample_lag <-
  as.numeric(median(labs_clean$`Lag between Specimen Date and Sending date`, na.rm = TRUE))

weekly_rolling_avg <-
  labs_clean %>%
  distinct(Facility_clean, Reporting_clean, DISEASE, MMWRyear, MMWRweek, `Labs by week`, `Days since last submitted lab`) %>%
  group_by(Facility_clean, Reporting_clean, DISEASE) %>%
  arrange(Facility_clean, Reporting_clean, desc(MMWRyear), desc(MMWRweek)) %>% 
  mutate(
    '4 week rolling average' = roll_mean(`Labs by week`, n = 4, fill = NA, align = "right"),
    '4 week rolling SD' = roll_sd(`Labs by week`, n = 4, fill = NA, align = "right"),
    'Lab volume flag - 1 SD' = case_when(
      `Labs by week` > (`4 week rolling average` + (1 * `4 week rolling SD`)) |
        `Labs by week` < (`4 week rolling average` - (1 * `4 week rolling SD`)) ~ "Flag", TRUE ~ "Normal"),
    'Lab volume flag - 2 SD' = case_when(
      `Labs by week` > (`4 week rolling average` + (2 * `4 week rolling SD`)) |
        `Labs by week` < (`4 week rolling average` - (2 * `4 week rolling SD`)) ~ "Flag", TRUE ~ "Normal"),
    'Lab submission time flag' = case_when(
      `Days since last submitted lab` > 7 ~ "Flag", TRUE ~ "Normal")
  )

labs_clean %<>%
  left_join(
    weekly_rolling_avg, by = c("Facility_clean", "Reporting_clean", "DISEASE", "MMWRyear", "MMWRweek", "Labs by week", "Days since last submitted lab")
  ) %>% 
  relocate(
    `4 week rolling average`, `4 week rolling SD`, .after = `Labs by week`
  )
toc()

print(paste("Step 8: Calculate maximum lag time by facility, create flag:", date()))
tic("Max lag and flag")
max_lag_by_facility <-
  labs_clean %>% 
  distinct(Facility_clean, Reporting_clean, DISEASE, MMWRyear, MMWRweek, `Lag between Specimen Date and Sending date`, .keep_all = TRUE) %>% 
  group_by(Facility_clean, Reporting_clean, DISEASE) %>%
  arrange(Facility_clean, Reporting_clean, desc(MMWRyear), desc(MMWRweek)) %>% 
  mutate(
    'Max lag' = max(`Lag between Specimen Date and Sending date`),
    'Abnormal lag' = case_when(
      `Max lag` > mean_sample_lag * 2 ~ "Flag", TRUE ~ "Normal")
  ) %>% 
  select(
    Facility_clean,
    Reporting_clean,
    DISEASE,
    MMWRyear,
    MMWRweek, 
    `Lag between Specimen Date and Sending date`, 
    `Max lag`, 
    `Abnormal lag`
  )

labs_clean %<>%
  left_join(
    max_lag_by_facility, by = c("Facility_clean", "Reporting_clean", "DISEASE", "MMWRyear", "MMWRweek", "Lag between Specimen Date and Sending date")
  ) %>% 
  relocate(
    `Max lag`, .after = `Average lag`
  ) %>% 
  select(-EVENT_DATE) %>% 
  mutate(
    Facility_clean = coalesce(Facility_clean, Reporting_clean, SENDING_FACILITY)
  )
toc()

print(paste("Step 9: Create facility aggregate and facilities:", date()))
tic("Aggregate and flagged facilities")
facilities_aggregated <-
  labs_clean %>%
  select(
    Facility_clean,
    Reporting_clean,
    DISEASE,
    TEST,
    SENDING_FACILITY,
    OBSERVATION_TYPE,
    MMWRyear,
    MMWRweek,
    `Week ending:`,
    `Latest Facility Date`,
    `Days since last submitted lab`,
    `Average lag`,
    `Max lag`,
    `SD lag`,
    `Labs by week`, 
    `4 week rolling average`, 
    `4 week rolling SD`,
    `Lab volume flag - 1 SD`,
    `Lab volume flag - 2 SD`,
    `Lab submission time flag`
  ) %>% 
  distinct(
    Facility_clean,
    Reporting_clean,
    DISEASE,
    MMWRyear,
    MMWRweek, 
    .keep_all = TRUE
  ) %>% 
  arrange(
    Facility_clean,
    Reporting_clean,
    desc(MMWRyear),
    desc(MMWRweek)
  )

facilities_flagged <-
  labs_clean %>%
  select(
    Facility_clean,
    Reporting_clean,
    DISEASE,
    TEST,
    SENDING_FACILITY,
    OBSERVATION_TYPE,
    MMWRyear,
    MMWRweek, 
    `Labs by week`, 
    `4 week rolling average`, 
    `4 week rolling SD`, 
    `Max lag`, 
    `Days since last submitted lab`, 
    `Lab volume flag - 1 SD`,
    `Lab volume flag - 2 SD`,
    `Abnormal lag`, 
    `Lab submission time flag`
  ) %>% 
  filter(
    `Lab volume flag - 1 SD` == 'Flag' |
      `Lab volume flag - 2 SD` == 'Flag' |
      `Lab submission time flag` == 'Flag' |
      `Abnormal lag` == 'Flag'
  ) %>% 
  distinct(
    Facility_clean, 
    Reporting_clean, 
    .keep_all =  TRUE) %>% 
  arrange(
    Facility_clean,
    Reporting_clean,
    desc(MMWRyear),
    desc(MMWRweek)
  )

labs_clean %<>%
  select(-c(`Latest Facility Date`:last_col()))

toc()

print(paste("Step 10: Saving:", date()))
tic("Saving")

write_csv(
  facilities_aggregated, paste0(output_path, "DIM facility aggregate, all diseases, past 180 days ", Sys.Date(),  ".csv")
)

write_csv(
  labs_clean, paste0(output_path, "DIM lab aggregate, all diseases, past 180 days ", Sys.Date(), ".csv")
)

toc()
toc()
print(paste("Finished:", date()))