#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$                                                                                    $#
#$                         Created by Lyubomir N. Kolev  II                           $#
#$                   SCDHEC, Division of Acute Disease Epidemiology                   $#
#$                                  November 22, 2023                                 $#
#$                                                                                    $#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

print(paste("Starting:", date()))
### install packages ####
# installs and loads necessary packages automatically
if(!require("pacman")){install.packages("pacman")}
p_load(tictoc, odbc, glue, MMWRweek, magrittr, tidyverse, readxl, writexl, dtplyr, plotly)

network_folder <-
  "network_path"

### set MMWR dates and years ####
current_year <- 
  year(Sys.Date())

MMWR <-
  MMWRweek(Sys.Date())

current_MMWR <- 
  MMWR$MMWRweek

previous_MMWR_end <- 
  MMWRweek2Date(
    MMWRyear = if_else(current_MMWR == 1, MMWR$MMWRyear - 1, MMWR$MMWRyear), # if_else for first week of new year, carries through file save
    MMWRweek = if_else(current_MMWR == 1, 52, MMWR$MMWRweek - 1),
    MMWRday = 7
  )

current_season_start <- 
  MMWRweek2Date(
    MMWRyear = if_else(current_MMWR >= 40, current_year, current_year - 1), 
    MMWRweek = 40, 
    MMWRday = 1)

years_back_to <- 
  if_else(between(current_MMWR, 40, 52), current_year - 5, current_year - 6)

years <- c(
  years_back_to, 
  years_back_to + 1, 
  years_back_to + 2, 
  years_back_to + 3, 
  years_back_to + 4)

year_5_start <- MMWRweek2Date(MMWRyear = years[[1]], MMWRweek = 40, MMWRday = 1)
year_4_start <- MMWRweek2Date(MMWRyear = years[[2]], MMWRweek = 40, MMWRday = 1)
year_3_start <- MMWRweek2Date(MMWRyear = years[[3]], MMWRweek = 40, MMWRday = 1)
year_2_start <- MMWRweek2Date(MMWRyear = years[[4]], MMWRweek = 40, MMWRday = 1)
year_1_start <- MMWRweek2Date(MMWRyear = years[[5]], MMWRweek = 40, MMWRday = 1)

historic_start <- 
  year_5_start

historic_end <- 
  current_season_start - 1

### connect to DWH and load historic data ####
tic("Total runtime:")
tic("Connecting and downloading")
print(paste("Step 1: Connecting to DWH and downloading data:", date()))

DWH_connection <-
  dbConnect(
    odbc(),
    dsn = "DPH_DWH"
  )

flu_historic <- # if beginning of a new month, run historic data first for accurate comparisonsk - otherwise this will error out
  read_excel(
    paste0(
      network_folder,
      pattern = paste0(
        "flulabshistoric",
        if_else(current_MMWR == 1, year(format(Sys.Date(), "%y-%m-%d")) - 1, year(format(Sys.Date(), "%y-%m-%d"))),
        if_else(current_MMWR == 1, "DEC", str_to_upper(month(Sys.Date(), label = TRUE))),
        ".xlsx"
        )
      )
  ) %>%
  as_tibble()

#SQL query ####
flu_dirty <-
  dbGetQuery(DWH_connection, glue("
  SELECT DISTINCT
	  UNIQUE_ID,
	  cases.PATIENT_ID,
	  cases.DISEASE,
	  cases.AGE_AT_ILLNESS,
	  COUNTY,
	  CURRENT_COUNTY,
	  cases.DISEASE_STATUS,
	  cases.CASE_EVENT_CLOSED,
	  NOTIFICATION_STATUS,
	  PRODUCT_CODE,
	  cases.EVENT_DATE,
	  cases.MMWR_WEEK,
	  cases.MMWR_YEAR,
	  labs.SPECIMEN_DATE,
	  labs.FACILITY,
	  labs.REPORTING_FACILITY,
	  labs.REPORTING_FACILITY_OTHER,
	  TEST,
	  RESULT,
	  RESULT_VALUE,
	  TEST_LOCAL_DESCRIPTION,
	  DEATH,
	  DEATH_DATE,
	  STATE,
	  labs.COMMENTS

  FROM
	  CASES_TABLE AS cases
	  INNER JOIN
	  LABS_TABLE AS labs
	  ON cases.PATIENT_ID = labs.PATIENT_ID

  WHERE
	  cases.EVENT_DATE BETWEEN '{current_season_start}' AND '{previous_MMWR_end}'
	  AND
	  STATE = 'SC'
  ")) %>%
  lazy_dt()
toc()

### group test types ####
tic("Test and result grouping")
print(paste("Step 2: Creating test and result groupings:", date()))

RNA_tests <-
  c(
    'FLUAV H1 2009 pandemic RNA XXX Ql PCR: ACnc: Pt: XXX: Ord: Prob.amp.tar',
    'FLUAV H1 HA gene Nph Ql NAA+probe',
    'FLUAV H1 RNA XXX Ql PCR: ACnc: Pt: XXX: Ord: Probe.amp.tar',
    'FLUAV H3 RNA Nph Ql NAA+probe',
    'FLUAV H3 RNA Upper resp Ql NAA+probe',
    'FLUAV RNA Lower resp Ql NAA+non-probe',
    'FLUAV RNA XXX QL NAA+probe',
    'FLUAV Subtyp XXX PCR: Prid: Pt: XXX: Nom: Probe.amp.tar',
    'FLUBV NS gene Nph Ql NAA+probe',
    'FLUBV RNA  XXX NAA+probe: NCnc: Pt: XXX: Qn: Probe.amp.tar',
    'FLUBV RNA Lower resp Ql NAA+non-probe',
    'FLUBV RNA XXX Ql PCR: ACnc: Pt: XXX: Ord: Probe.amp.tar',
    'Influenza virus A H1 2009 pandemic RNA  PrThr  Pt  Nph  Ord  Non-probe.amp.tar',
    'Influenza virus A H1 2009 pandemic RNA panel - Unspecified specimen by NAA with probe detection',
    'Influenza virus A H1 RNA  PrThr  Pt  Nph  Ord  Non-probe.amp.tar',
    'Influenza virus A H1 RNA  PrThr  Pt  Nph  Ord  Probe.amp.tar',
    'Influenza virus A H3 HA gene  PrThr  Pt  Nph  Ord  Probe.amp.tar',
    'Influenza virus A H3 RNA  PrThr  Pt  Nph  Ord  Non-probe.amp.tar',
    'Influenza virus A H3 RNA  PrThr  Pt  XXX  Ord  Probe.amp.tar',
    'Influenza virus A M gene  PrThr  Pt  Nph  Ord  Probe.amp.tar',
    'Influenza virus A RNA  PrThr  Pt  Nph  Ord  Non-probe.amp.tar',
    'Influenza virus A RNA  PrThr  Pt  Nph  Ord  Probe.amp.tar',
    'Influenza virus A RNA [Presence] in Respiratory specimen by NAA with probe detection',
    'Influenza virus A RNA [Presence] in Upper respiratory specimen by NAA with probe detection',
    'Influenza virus A RNA: ACnc: Pt: xxx: Ord: Probe.Amp.Tar',
    'Influenza virus A swine origin RNA  Prid  Pt  XXX  Nom  Probe.amp.tar',
    'Influenza virus B lineage RNA  Prid  Pt  XXX  Nom  Probe.amp.tar',
    'Influenza virus B RNA  PrThr  Pt  Nph  Ord  Non-probe.amp.tar',
    'Influenza virus B RNA [Presence] in Respiratory specimen by NAA with probe detection',
    'Influenza virus B RNA [Presence] in Upper respiratory specimen by NAA with probe detection',
    'Influenza virus B RNA PrThr Pt Nph Ord Probe.amp.tar',
    'FLUABV + SARS-CoV-2 Resp NAA+probe',
    'FLUABV+SARS-CoV-2+RSV Pnl Resp NAA+probe',
    'FLUAV+FLUBV RNA XXX PCR: Prid: Pt: XXX: Nom: Probe.amp.tar',
    'FLUAV+FLUBV+RSV RNA pnl Up resp NAA+pr',
    'Influenza virus A and B and SARS-CoV-2 (COVID-19) RNA panel NAA+probe (Resp)',
    'Influenza virus A+B RNA  PrThr  Pt  XXX  Ord  Probe.amp.tar'
  )

Culture_tests <-
  c(
    'Influenza virus identified: PrId: Pt: xxx: Nom: Shell vial culture',
    'FLUAV XXX Ql Cult: ACnc: Pt: XXX: Ord: Organism specific culture',
    'FLUV XXX Cult: ACnc: Pt: XXX: Nom: Organism specific culture',
    'Influenza virus B  PrThr  Pt  XXX  Ord  Organism specific culture'
  )

Ab_tests <-
  c(
    'FLUAV Ab Ser Ql',
    'FLUAV Ab Titr Ser CF: Titr: Pt: Ser: Qn: Comp fix',
    'FLUAV IgG Ser Ql',
    'FLUBV Ab Fld Ql',
    'FLUBV Ab Titr Ser CF: Titr: Pt: Ser: Qn: Comp fix'
  )

Ag_tests <-
  c(
    'FLUAV Ag Bronch Ql',
    'FLUAV Ag Nose Ql',
    'FLUAV Ag Nph Ql IA',
    'FLUAV Ag Spec Ql IA',
    'FLUAV Ag Upper resp Ql IA.rapid',
    'FLUAV and FLUBV Ag Nose EIA.rapid Prid: Pt: Node: Nom: EIA.rapid',
    'FLUAV+FLUBV Ag Nose Ql IA.rapid',
    'FLUAV+FLUBV Ag XXX Ql IF: ACnc: Pt: XXX: Ord: IF',
    'FLUBV Ag Spec Ql IA',
    'FLUBV Ag Throat Ql',
    'FLUBV Ag Throat Ql IF',
    'FLUBV Ag Upper resp Ql IA.rapid',
    'FLUBV Ag Nph Ql',
    'Influenza virus A Ag  PrThr  Pt  Thrt  Ord',
    'Influenza virus A Ag  PrThr  Pt  XXX  Ord  IF',
    'Influenza virus A Ag [Presence] in Nose by Immunoassay',
    'Influenza virus A Ag: ACnc: Pt: xxx: Ord:',
    'Influenza virus A+B Ag: ACnc: Pt: xxx: Ord:',
    'Influenza virus A+B+C Ag  PrThr  Pt  XXX  Ord',
    'Influenza virus B Ag  PrThr  Pt  Nose  Ord',
    'Influenza virus B Ag  PrThr  Pt  Nph  Ord  IA',
    'Influenza virus B Ag: ACnc: Pt: xxx: Ord:'
  )

Positive_results <-
  c(
    'Detected',
    'Influenza A',
    'Influenza A virus',
    'Influenza A virus (organism)',
    'Influenza A virus subtype H1 2009 pandemic',
    'Influenza A virus subtype H3',
    'Influenza B',
    'Influenza B virus',
    'Influenza B virus Victoria lineage',
    'Flu A',
    'Flu B',
    'Positive',
    'Reactive'
  )

toc()

### sum results by test type ####
tic("Summing positive results by test type")
print(paste("Step 3: Summing positive results by test type for each Case ID:", date()))
flu_clean <-
  flu_dirty %>%
  filter(
    TEST %in% c(RNA_tests, Ab_tests, Ag_tests, Culture_tests)
  ) %>%
  group_by(CASE_ID) %>%
  mutate(
    'Total tests' = length(TEST),
    'RV Positive' = grepl("detect|pos|lineage|H1|H3|Influenza A|Influenza A virus|Influenza B|Influenza B virus|Flu A|Flu B", RESULT_VALUE, ignore.case = TRUE) & !grepl("not|neg", RESULT_VALUE, ignore.case = TRUE), # made this change because some results were "Negative Influenza A" etc.
    'RNA Positives' = sum(case_when(
      TEST %in% RNA_tests & (RESULT %in% Positive_results | `RV Positive` == TRUE) ~ 1, TRUE ~ 0)),
    'Ag Positives' = sum(case_when(
      TEST %in% Ag_tests & (RESULT %in% Positive_results | `RV Positive` == TRUE) ~ 1, TRUE ~ 0)),
    'Ab Positives' = sum(case_when(
      TEST %in% Ab_tests & (RESULT %in% Positive_results | `RV Positive` == TRUE) ~ 1, TRUE ~ 0)),
    'Culture Positives' = sum(case_when(
      TEST %in% Culture_tests & (RESULT %in% Positive_results | `RV Positive` == TRUE) ~ 1, TRUE ~ 0))
  ) %>%
  relocate(
    `RV Positive`, .after = RESULT_VALUE
  )
toc()

### compute QC for case status ####
tic("Computing QC")
print(paste("Step 4: Computing QC variable:", date()))

flu_clean %<>%
  group_by(CASE_ID) %>% 
  mutate(
    'Final QC' = if_else(
      (DISEASE_STATUS == 'Confirmed' & DIED_FROM_DISEASE == 'Yes' & !is.na(DIED_FROM_DISEASE) |
         ((`RNA Positives` >= 1 | `Culture Positives` >= 1) & DISEASE_STATUS == 'Confirmed') |
         ((`RNA Positives` >= 0 | `Culture Positives` >= 0) & `Ag Positives` >= 1 & `Ab Positives` == 0 & (DISEASE_STATUS != 'Confirmed')) |
         ((`RNA Positives` >= 0 | `Culture Positives` >= 0) & `Ag Positives` == 0 & `Ab Positives` >= 1 & (DISEASE_STATUS != 'Confirmed')) |
         ((`RNA Positives` >= 0 | `Culture Positives` >= 0) & `Ag Positives` >= 1 & `Ab Positives` >= 1 & (DISEASE_STATUS != 'Confirmed')) |
         ((`RNA Positives` >= 0 | `Culture Positives` >= 0) & `Ag Positives` == 0 & `Ab Positives` == 0 & (DISEASE_STATUS != 'Confirmed')) |
         (STATE != 'SC' & DISEASE_STATUS == 'Not a case') |
         (`RNA Positives` >= 1 & DISEASE_STATUS == 'Suspect' & CASE_EVENT_CLOSED == 'No' & !(NOTIFICATION_STATUS %in% c('Submitted', 'Approved')))
      ), "Pass", "Fail"),
    'Final QC' = case_when(
      grepl('Yamagata', RESULT, ignore.case = TRUE) | grepl('Yamagata', RESULT_VALUE, ignore.case = TRUE) ~ 'Fail', TRUE ~ `Final QC`,
      (`RNA Positives` >= 1 & DISEASE_STATUS == 'Suspect' & CASE_EVENT_CLOSED == 'Yes' & is.na(NOTIFICATION_STATUS) & ymd(EVENT_DATE) < Sys.Date() - 32) ~ 'Fail', TRUE ~ 'Pass'), # adjust date to subtract from as needed
    'Time lag' = as_date(EVENT_DATE) - as_date(SPECIMEN_DATE)
  ) %>% 
  relocate(
    `Time lag`, .after = SPECIMEN_DATE) %>% 
  as_tibble()
toc()

### aggregate positive results ####
tic("Tallying positive results")
print(paste("Step 5: Summing all positive results:", date()))
positives_by_group <-
  flu_clean %>%
  filter(
    RESULT %in% Positive_results |
      `RV Positive` == TRUE
  ) %>% 
  mutate(
    'Test Group' = case_when(
      TEST %in% RNA_tests ~ "RNA",
      # TEST %in% Combo_tests ~ "Combo",
      TEST %in% Culture_tests ~ "Culture",
      TEST %in% Ab_tests ~ "Ab",
      TEST %in% Ag_tests ~ "Ag",
      TRUE ~ "Unknown")
  ) %>% 
  group_by(
    `Test Group`,
    TEST,
    RESULT,
    RESULT_VALUE,
    TEST_LOCAL_DESCRIPTION
  )

positives_by_group_summary <-
  positives_by_group %>% 
  summarize(
    'Count' = length(CASE_ID)
  )
toc()

### type and subtype ####
tic("Typing and subtyping")
print(paste("Step 6: Analyzing type and subtype:", date()))
flu_a_and_b <- # total for all flu a and b rna
  positives_by_group %>%
  filter(
    `Test Group` %in% c('RNA', 'Culture')
  ) %>% 
  mutate(
    'Flu Type' = case_when(
      (grepl("\\s*Influenza virus A\\s*|\\s*FLUAV\\s*", TEST, ignore.case = TRUE) & !grepl(" and |V\\+F|A\\+B|V\\+SARS|V \\+ SARS", TEST)) |
        grepl("Influenza A positive|Flu A positive|Flu A|Influenza A|Positive A|Positive flu A|Positive for flu A|Positive A virus", COMMENTS, ignore.case = TRUE) ~ "Flu A", # V + F covers FLUAV+FLUBV
      (grepl("\\s*Influenza virus B\\s*|\\s*FLUBV\\s*", TEST, ignore.case = TRUE) & !grepl(" and |V\\+F|A\\+B|V\\+SARS|V \\+ SARS", TEST)) |
        grepl("Influenza B positive|Flu B positive|Flu B|Influenza B|Positive B|Positive flu B|Positive for flu B|Positive B virus", COMMENTS, ignore.case = TRUE) ~ "Flu B",
      
      grepl("\\s*and\\s*|V\\+F|A\\+B|V\\s*\\+\\s*SARS|\\s*FLUV\\s*", TEST) & grepl("Influenza\\s*A", RESULT) |
        (grepl("Influenza\\s*A|Flu A|\\s*H1\\s*2009|\\s*H1-2009|\\s*H1\\s*(2009)\\s*|\\s*FLUAV\\s*", RESULT_VALUE, ignore.case = TRUE) |
           grepl("Influenza\\s*A|Flu A|H1\\s*2009|H1-2009|H1\\s*(2009)\\s*|\\s*FLUAV\\s*", TEST_LOCAL_DESCRIPTION, ignore.case = TRUE)) ~ "Flu A",
      
      grepl("\\s*and\\s*|V\\+F|A\\+B|V\\s*\\+\\s*SARS|\\s*FLUV\\s*", TEST) & grepl("Influenza\\s*B", RESULT) |
        (grepl("Influenza\\s*B|Flu\\s*B|\\s*FLUBV\\s*", RESULT_VALUE, ignore.case = TRUE) |
           grepl("Influenza\\s*B|Flu\\s*B|\\s*FLUBV\\s*", TEST_LOCAL_DESCRIPTION, ignore.case = TRUE)) ~ "Flu B",
      TRUE ~ "Unknown"),
    
    'Subtype' = case_when(
      grepl("H3", TEST) | grepl("H3", RESULT) | grepl("H3", RESULT_VALUE) | grepl("H3", TEST_LOCAL_DESCRIPTION) ~ "H3",
      grepl("H1 prior to 2010", RESULT) | grepl("H1 prior to 2010", RESULT_VALUE) | grepl("H1 prior to 2010", TEST_LOCAL_DESCRIPTION) ~ "H1 prior to 2010",
      grepl("H1 2009|H1-2009", TEST) | grepl("H1 2009|H1-2009", RESULT) | grepl("H1 2009|H1-2009", RESULT_VALUE) | grepl("H1 2009|H1-2009", TEST_LOCAL_DESCRIPTION) ~ "H1 2009",
      grepl("H1", TEST) ~ "H1",
      grepl("H1N1 pdm09", RESULT) | grepl("H1N1 pdm09", RESULT_VALUE) | grepl("H1N1 pdm09", TEST_LOCAL_DESCRIPTION) ~ "H1N1 pdm09",
      grepl("H2N2v", RESULT) | grepl("H2N2v", RESULT_VALUE) | grepl("H2N2v", TEST_LOCAL_DESCRIPTION) ~ "H2N2v",
      grepl("Victoria", RESULT) | grepl("Victoria", RESULT_VALUE) | grepl("Victoria", TEST_LOCAL_DESCRIPTION) ~ "Victoria",
      grepl("Yamagata", RESULT) | grepl("Yamagata", RESULT_VALUE) | grepl("Yamagata", TEST_LOCAL_DESCRIPTION) ~ "Yamagata",
      grepl("Unable", RESULT, ignore.case = TRUE) | 
        grepl("Unable", RESULT_VALUE, ignore.case = TRUE) | 
        grepl("Unable", TEST_LOCAL_DESCRIPTION, ignore.case = TRUE) & 
        !(grepl("invalid", RESULT, ignore.case = TRUE) |
            grepl("invalid", RESULT_VALUE, ignore.case = TRUE) |
            grepl("invalid", TEST_LOCAL_DESCRIPTION, ignore.case = TRUE)) ~ "Unable to Subtype",
      TRUE ~ NA),
    
    SCION_variable = case_when(
      `Flu Type` == "Flu A" & Subtype == "H1" ~ "Influenza A (H1N1)pdm09",
      `Flu Type` == "Flu A" & Subtype == "H1 2009" ~ "Influenza A (H1N1)pdm09",
      `Flu Type` == "Flu A" & Subtype == "H3" ~ "Influenza A (H3)",
      `Flu Type` == "Flu A" & is.na(Subtype) ~ "Influenza A (Subtyping Not Done)",
      `Flu Type` == "Flu B" & Subtype == "Victoria" ~ "Influenza B/Victoria lineage",
      `Flu Type` == "Flu B" & is.na(Subtype) ~ "Influenza B (Lineage Not Determined)",
      `Flu Type` == "Coinfection" ~ "Influenza virus co-infection",
      `Flu Type` == "Unknown" ~ "Unknown",
      Subtype == "Unable to Subtype" ~ "Influenza A (Unable to Subtype)")
  )

flu_a_and_b_join <- # created to properly join coinfections without expanding the rows
  flu_a_and_b %>%
  as_tibble() %>% 
  select(PARTY_UNID, CASE_ID, `Test Group`:last_col())

flu_clean %<>%
  left_join(flu_a_and_b_join, by = c("PARTY_UNID", "CASE_ID"))
toc()

### coinfections ####
tic("Coinfections")
print(paste("Step 7: Checking for coinfections:", date()))
coinfections_bydate <- 
  flu_a_and_b %>%
  group_by(CASE_ID) %>% 
  mutate(
    'Flu A count' = sum(case_when(
      `Flu Type` == "Flu A" ~ 1, TRUE ~ 0)),
    'Flu B count' = sum(case_when(
      `Flu Type` == "Flu B" ~ 1, TRUE ~ 0)),
    'Flu A date' = case_when(
      `Flu Type` == "Flu A" ~ SPECIMEN_DATE, TRUE ~ NA),
    'Flu B date' = case_when(
      `Flu Type` == "Flu A" ~ SPECIMEN_DATE, TRUE ~ NA)
  ) %>% 
  filter(
    `Flu A count` >= 1 & `Flu B count` >= 1,
    !is.na(`Flu A date`),
    !is.na(`Flu B date`)
  )
if (nrow(coinfections_bydate) > 0) { # added this 'if' because a mutate() with empty df causes error
  coinfections_bydate <-
    coinfections_bydate %>% 
    mutate(
      'Coinfection' = case_when(
        abs(difftime(`Flu A date`, `Flu B date`, units = "days")) <= 31 ~ "Yes",
        TRUE ~ NA
      )
    ) %>% 
    distinct(CASE_ID, .keep_all = TRUE) 
} else {
  print("No coinfections, skipping calculation.")
}

coinfections_comments <-
  flu_a_and_b %>% 
  group_by(CASE_ID) %>% 
  mutate(
    'Flu A count' = sum(case_when(
      `Flu Type` == "Flu A" ~ 1, TRUE ~ 0)),
    'Flu B count' = sum(case_when(
      `Flu Type` == "Flu B" ~ 1, TRUE ~ 0)),
    'Flu A date' = case_when(
      `Flu Type` == "Flu A" ~ SPECIMEN_DATE, TRUE ~ NA),
    'Flu B date' = case_when(
      `Flu Type` == "Flu A" ~ SPECIMEN_DATE, TRUE ~ NA),
    'Coinfection' = case_when(
      grepl("a and b|a&b|a & b|b&a|b & a|b and a|both|influenza a and|flu a and|influenza b and|flu b and", COMMENTS, ignore.case = TRUE) &
        !grepl("negative|rsv|cov", COMMENTS, ignore.case = TRUE) ~ "Yes",
      TRUE ~ NA
    )
  ) %>% 
  filter(Coinfection == "Yes") %>% 
  distinct(CASE_ID, .keep_all = TRUE)

coinfections <-
  rbind(coinfections_bydate, coinfections_comments) %>%
  select(PARTY_UNID, CASE_ID, Coinfection)

flu_clean %<>%
  left_join(coinfections, by = c("PARTY_UNID", "CASE_ID")) %>% 
  relocate(`Final QC`, .after = "Culture Positives") %>% 
  select(-Coinfection)

flu_failed <- # all failed labs, including RNA, Ag, Ab
  flu_clean %>% 
  filter(
    `Final QC` == 'Fail'
  )

failed_IDs <-
  flu_failed %>% 
  select(CASE_ID) %>% 
  distinct()

flu_a_and_b %<>%
  left_join(coinfections, by = c("PARTY_UNID", "CASE_ID")) %>% 
  mutate(
    `Flu Type` = case_when(
      Coinfection == "Yes" ~ "Coinfection", TRUE ~ `Flu Type`),
    SCEDSS_variable = case_when(
      Coinfection == "Yes" ~ "Influenza virus co-infection", TRUE ~ SCION_variable)
  ) %>% 
  select(-Coinfection)

flu_a_and_b_summary <- # for overall counts by type and subtype
  flu_a_and_b %>%
  group_by(
    `Flu Type`, `Subtype`
  ) %>%
  summarize(
    Count = length(CASE_ID)
  )  %>% 
  mutate(
    Count = case_when(
      `Flu Type` == "Coinfection" ~ nrow(coinfections), TRUE ~ Count
    )
  )
toc()

### compute historic comparisons ####
tic("Historic comparisons")
print(paste("Step 8: Comparing to historical data:", date()))
year_5 <-
  flu_historic %>% 
  filter(
    between(as_date(EVENT_DATE), year_5_start, (year_4_start - 1))
  )

year_5_mmwr_count <-
  year_5 %>% 
  select(CASE_ID, EVENT_DATE, MMWR_WEEK, MMWR_YEAR) %>% 
  arrange(EVENT_DATE) %>%
  group_by(MMWR_WEEK) %>% 
  mutate(
    'Weekly Count' = length(CASE_ID)
  ) %>% 
  ungroup() %>%
  select(-c(CASE_ID, EVENT_DATE)) %>%
  distinct(MMWR_YEAR, MMWR_WEEK, `Weekly Count`) %>% 
  mutate(
    'Cumulative Count' = cumsum(`Weekly Count`)
  ) %>% 
  filter(
    MMWR_WEEK == if_else(current_MMWR == 1, 52, MMWR$MMWRweek - 1)
  )

year_4 <-
  flu_historic %>% 
  filter(
    between(as_date(EVENT_DATE), year_4_start, (year_3_start - 1))
  )

year_4_mmwr_count <-
  year_4 %>% 
  select(CASE_ID, EVENT_DATE, MMWR_WEEK, MMWR_YEAR) %>% 
  arrange(EVENT_DATE) %>%
  group_by(MMWR_WEEK) %>% 
  mutate(
    'Weekly Count' = length(CASE_ID)
  ) %>% 
  ungroup() %>%
  select(-c(CASE_ID, EVENT_DATE)) %>%
  distinct(MMWR_YEAR, MMWR_WEEK, `Weekly Count`) %>% 
  mutate(
    'Cumulative Count' = cumsum(`Weekly Count`)
  ) %>% 
  filter(
    MMWR_WEEK == if_else(current_MMWR == 1, 52, MMWR$MMWRweek - 1)
  )

year_3 <-
  flu_historic %>% 
  filter(
    between(as_date(EVENT_DATE), year_3_start, (year_2_start - 1))
  )

year_3_mmwr_count <-
  year_3 %>% 
  select(CASE_ID, EVENT_DATE, MMWR_WEEK, MMWR_YEAR) %>% 
  arrange(EVENT_DATE) %>%
  group_by(MMWR_WEEK) %>% 
  mutate(
    'Weekly Count' = length(CASE_ID)
  ) %>% 
  ungroup() %>%
  select(-c(CASE_ID, EVENT_DATE)) %>%
  distinct(MMWR_YEAR, MMWR_WEEK, `Weekly Count`) %>% 
  mutate(
    'Cumulative Count' = cumsum(`Weekly Count`)
  ) %>% 
  filter(
    MMWR_WEEK == if_else(current_MMWR == 1, 52, MMWR$MMWRweek - 1)
  )

year_2 <-
  flu_historic %>% 
  filter(
    between(as_date(EVENT_DATE), year_2_start, (year_1_start - 1))
  )

year_2_mmwr_count <-
  year_2 %>% 
  select(CASE_ID, EVENT_DATE, MMWR_WEEK, MMWR_YEAR) %>% 
  arrange(EVENT_DATE) %>%
  group_by(MMWR_WEEK) %>% 
  mutate(
    'Weekly Count' = length(CASE_ID)
  ) %>% 
  ungroup() %>%
  select(-c(CASE_ID, EVENT_DATE)) %>%
  distinct(MMWR_YEAR, MMWR_WEEK, `Weekly Count`) %>% 
  mutate(
    'Cumulative Count' = cumsum(`Weekly Count`)
  ) %>% 
  filter(
    MMWR_WEEK == if_else(current_MMWR == 1, 52, MMWR$MMWRweek - 1)
  )

year_1 <-
  flu_historic %>% 
  filter(
    between(as_date(EVENT_DATE), year_1_start, (current_season_start - 1))
  )

year_1_mmwr_count <-
  year_1 %>% 
  select(CASE_ID, EVENT_DATE, MMWR_WEEK, MMWR_YEAR) %>% 
  arrange(EVENT_DATE) %>%
  group_by(MMWR_WEEK) %>% 
  mutate(
    'Weekly Count' = length(CASE_ID)
  ) %>% 
  ungroup() %>%
  select(-c(CASE_ID, EVENT_DATE)) %>%
  distinct(MMWR_YEAR, MMWR_WEEK, `Weekly Count`) %>% 
  mutate(
    'Cumulative Count' = cumsum(`Weekly Count`)
  ) %>% 
  filter(
    MMWR_WEEK == if_else(current_MMWR == 1, 52, MMWR$MMWRweek - 1)
  )

current_mmmwr_count <-
  flu_a_and_b %>%
  ungroup() %>% 
  arrange(EVENT_DATE) %>%
  group_by(MMWR_WEEK) %>% 
  mutate(
    'Weekly Count' = length(CASE_ID)
  ) %>% 
  ungroup() %>%
  select(-c(CASE_ID, EVENT_DATE)) %>%
  distinct(MMWR_YEAR, MMWR_WEEK, `Weekly Count`) %>% 
  mutate(
    'Cumulative Count' = cumsum(`Weekly Count`)
  ) %>% 
  filter(
    MMWR_WEEK == if_else(current_MMWR == 1, 52, MMWR$MMWRweek - 1)
  )

mmwr_comparison <- # will roll over with calendar year due to current MMWR week being in new calendar year
  rbind(
    year_5_mmwr_count, year_4_mmwr_count, year_3_mmwr_count, year_2_mmwr_count, year_1_mmwr_count, current_mmmwr_count
  )
toc()

### output files ####
tic("Saving")
print(paste("Step 9: Save files:", date())) # saves to csv then to xlsx due to parsing issue with a provider name which breaks xlsx

write_csv(
  flu_a_and_b,
  paste0(network_folder, "flulabsmmwr", if_else(current_MMWR == 1, year(format(Sys.Date(), "%y-%m-%d")) - 1, year(format(Sys.Date(), "%y-%m-%d"))), if_else(current_MMWR == 1, 52, current_MMWR - 1), ".csv")
)

csv_to_xlsx <-
  read_csv(
    paste0(
      network_folder, "flulabsmmwr", if_else(current_MMWR == 1, year(format(Sys.Date(), "%y-%m-%d")) - 1, year(format(Sys.Date(), "%y-%m-%d"))), if_else(current_MMWR == 1, 52, current_MMWR - 1), ".csv")
  )

write_xlsx(
  csv_to_xlsx,
  paste0(network_folder, "flulabsmmwr", if_else(current_MMWR == 1, year(format(Sys.Date(), "%y-%m-%d")) - 1, year(format(Sys.Date(), "%y-%m-%d"))), if_else(current_MMWR == 1, 52, current_MMWR - 1), ".xlsx")
)

file.remove(
  paste0(
    network_folder, "flulabsmmwr", if_else(current_MMWR == 1, year(format(Sys.Date(), "%y-%m-%d")) - 1, year(format(Sys.Date(), "%y-%m-%d"))), if_else(current_MMWR == 1, 52, current_MMWR - 1), ".csv")
)

trend_df <-
  flu_a_and_b %>%
  ungroup() %>% 
  mutate(
    EVENT_DATE = as.POSIXct(EVENT_DATE),
    PARTY_UNID = as.numeric(PARTY_UNID),
    CASE_ID = as.numeric(CASE_ID),
    SPECIMEN_DATE = as_datetime(SPECIMEN_DATE),
    `Time lag` = as.numeric(`Time lag`)
  ) %>% 
  rbind(year_1) %>% 
  group_by(
    Year = year(EVENT_DATE), Month = month(EVENT_DATE, label = TRUE, abbr = FALSE)
  ) %>% 
  summarize(
    'Count' = n()
  ) %>%
  ungroup() %>% 
  mutate(
    rownum = seq_along(1:n()),
    Month_Year = paste(Month, Year, sep = " ")
  )
  
trend_graph <-
  ggplot(
    data = trend_df,
    aes(
      x = rownum,
      y = Count
    )
  ) +
  geom_line() +
  geom_point() +
  theme_light() +
  scale_x_continuous(breaks = trend_df$rownum, labels = trend_df$Month_Year) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(ggplotly(trend_graph))

print(paste("Finished", date()))
toc()
toc()