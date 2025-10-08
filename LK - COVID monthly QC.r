#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$                                                                                    $#
#$                         Created by Lyubomir N. Kolev II                            $#
#$                   SCDHEC, Division of Acute Disease Epidemiology                   $#
#$                                  June 12, 2023                                     $#
#$                                                                                    $#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

print(paste("Starting at:", date()))
## install and load necessary packages automatically ####
if(!require(pacman)){install.packages("pacman")}
pacman::p_load(odbc, glue, magrittr, tidyverse, dplyr, writexl)

previous_month_start <- floor_date(Sys.Date() - months(1), unit = "month")
previous_month_end <- ceiling_date(Sys.Date() - months(1), unit = "month") - 1

user <- Sys.getenv("username")

## set file paths ####
save_path <-
  paste0('C:/Users/', user, '/OneDrive/', month(previous_month_start, label = TRUE, abbr = FALSE), " ", if_else(month(Sys.Date()) == 1, year(Sys.Date()) - 1, year(Sys.Date())), ' QC.xlsx')

## set up DWH connection####
print(paste("Connecting to DWH and downloading data:", date()))
DWH_connection <- 
  dbConnect(
    odbc(),
    dsn = "DPH_warehouse"
  )

## SQL query ####
QC_dirty <- 
  dbGetQuery(DWH_connection, glue("
	SELECT DISTINCT
	labs.PATIENT_ID,
	cases.DISEASE,
	cases.COUNTY,
	labs.DATE,
	labs.DISEASE_STATUS,
	cases.CASE_CLOSED,
	cases.NOTIFICATION_STATUS,
	labs.SPECIMEN_DATE,
	labs.SPECIMEN_ID,
	labs.SPECIMEN_RECEIVED_DATE,
	labs.SPECIMEN_SOURCE,
	labs.ACCESSION_NUMBER,
	labs.TEST,
	labs.RESULT,
	labs.RESULT_VALUE,
	labs.RESULT_DATE,
	
	/*Tally RNA Positive Results*/
	CASE WHEN 
		(
			(  
		   TEST IN ('2019 Novel Coronavirus E gene [Presence] in Unspecified specimen by NAA with probe detection', 
					'2019 Novel Coronavirus N gene [Presence] in Unspecified specimen by NAA with probe detection', 
					'2019 Novel Coronavirus N gene [Presence] in Unspecified specimen by Nucleic acid amplification using primer-probe set N1', 
					'2019 Novel Coronavirus N gene [Presence] in Unspecified specimen by Nucleic acid amplification using primer-probe set N2', 
					'2019 Novel Coronavirus RdRp gene [Presence] in Unspecified specimen by NAA with probe detection', 
					'2019 Novel Coronavirus RNA panel - Unspecified specimen by NAA with probe detection', 
					'2019-nCoV RNA XXX NAA+probe-Imp: Imp: Pt: XXX: Ord: Probe.amp.tar', 
					'COVID-19 N gene in Saliva by NAA with probe detection',
					'COVID-19 RNA in Saliva by NAA with probe detection', 
					'COVID-19 RNA in Saliva by sequencing', 
					'FLUABV + SARS-CoV-2 Resp NAA+probe', 
					'FLUABV+SARS-CoV-2+RSV Pnl Resp NAA+probe',
					'HCoV 229E RNA Nph Ql NAA+non-probe: PrThr: Pt: Nph: Ord: Non-probe.amp.tar',
					'HCoV HKU1 RNA Nph Ql NAA+non-probe: PrThr: Pt: Nph: Ord: Non-probe.amp.tar', 
					'HCoV NL63 RNA Nph Ql NAA+non-probe: PrThr: Pt: Nph: Ord: Non-probe.amp.tar', 
					'HCoV OC43 RNA Nph Ql NAA+non-probe: PrThr: Pt: Nph: Ord: Non-probe.amp.tar', 
					'Influenza virus A and B and SARS-CoV-2 (COVID-19) RNA panel NAA+probe (Resp)', 
					'SARS coronavirus 2 ORF1ab region [Presence] in Respiratory specimen by NAA with probe detection', 
					'SARS coronavirus 2 RNA PrThr Pt Nph', 
					'SARS coronavirus 2 RNA PrThr Pt Respiratory Ord Probe.amp.tar', 
					'SARS coronavirus RNA [Presence] in Unspecified specimen by NAA with probe detection', 
					'SARS-CoV-2 (COVID-19) and SARS-related CoV RNA panel and Influenza virus A and B – Respiratory specimen by NAA with probe detectionActive', 
					'SARS-CoV-2 (COVID-19) N gene [Presence] in Respiratory specimen by NAA with probe detection', 
					'SARS-CoV-2 (COVID-19) RdRp gene NAA+probe (Unsp spec) [ThreshNum]', 
					'SARS-CoV-2 (COVID-19) RNA NAA+probe Ql (Nph)', 
					'SARS-CoV-2 N gene Nose Ql NAA+non-probe', 
					'SARS-CoV-2 N gene Nose Ql NAA+probe', 
					'SARS-CoV-2 N gene Nph Ql NAA+probe', 
					'SARS-CoV-2 N gene Resp Ql NAA N1', 
					'SARS-CoV-2 N gene Resp Ql NAA N2', 
					'SARS-CoV-2 N gene Sal Ql NAA N1', 
					'SARS-CoV-2 ORF1ab Sal Ql NAA+probe', 
					'SARS-CoV-2 RdRp Resp Ql NAA+probe', 
					'SARS-CoV-2 RNA Nose Ql NAA+probe', 
					'SARS-CoV-2 RNA Resp Ql NAA+non-probe', 
					'SARS-CoV-2 RNA Resp Ql NAA+probe', 
					'SARS-CoV-2 RNA Resp Ql Seq', 
					'SARS-like Coronavirus N gene [Cycle Threshold #] in Unspecified specimen by NAA with probe detection', 
					'SARS-like Coronavirus N gene [Presence] in Unspecified specimen by NAA with probe detection', 
					'SARS-related CoV RNA NAA+probe Ql (Resp)',
					'SARS-CoV-2 (COVID-19) and SARS-related CoV RNA panel and Influenza virus A and B - Respiratory specimen by NAA with probe detection')
			)
		
		AND 
				
			(
		(RESULT = 'Positive' OR RESULT = 'Detected' OR RESULT = 'Above Threshold' OR RESULT = 'Reactive'
			OR RESULT_VALUE LIKE '%pos%' OR RESULT_VALUE LIKE 'Det%' OR RESULT_VALUE LIKE 'COVID det%' OR RESULT_VALUE LIKE 'SARS_CoV-2 Det%' 
			OR RESULT_VALUE LIKE 'SARS-CoV-2 Det%' OR RESULT_VALUE LIKE '%(+)%') 
			)
		)
		
		THEN 1
		ELSE 0
		END AS 'RNA Positive Count',
		
		/*Tally Ag Positive Results*/
		CASE WHEN
		(  
		   TEST IN ('SARS-CoV+SARS-CoV-2 Ag Resp Ql IA.rapid',
					'SARS-CoV-2 (COVID-19) Ag [Presence] in Respiratory specimen by Rapid immunoassay',
					'SARS-CoV-2 (COVID-19) Ag [Presence] in Upper respiratory specimen by Immunoassay',
					'SARS-CoV-2 (COVID-19) Ag panel and Influenza virus A and B - Upper respiratory specimen by Rapid immunoassay',
					'SARS-CoV-2 Ag Upper resp Ql IA.rapid')
		
		AND
		
			(
		(RESULT = 'Positive' OR RESULT = 'Detected' OR RESULT = 'Above Threshold' OR RESULT = 'Reactive' OR RESULT LIKE '%antigen%'
			OR RESULT_VALUE LIKE '%pos%' OR RESULT_VALUE LIKE '%above thresh%' OR RESULT_VALUE LIKE 'React%' 
			OR RESULT_VALUE LIKE 'Det%' OR RESULT_VALUE LIKE 'COVID det%' OR RESULT_VALUE LIKE 'SARS_CoV-2 Det%' 
			OR RESULT_VALUE LIKE 'SARS-CoV-2 Det%' OR RESULT_VALUE LIKE '%(+)%') OR RESULT_VALUE LIKE '%antigen%'
			)
		)
		
		THEN 1
		ELSE 0
		END AS 'Ag Positive Count',
		
		/*Tally Ab Positive Results*/
		CASE WHEN
		(
			TEST IN ('SARS coronavirus 2 Ab [Presence] in Serum or Plasma by Immunoassay',
					 'SARS coronavirus 2 Ab PrThr Pt Ser/Plas Ord IA',
					 'SARS CORONAVIRUS 2 AB.IGA',
					 'SARS coronavirus 2 Ab.IgA PrThr Pt Ser/Plas  Ord IA',
					 'SARS CORONAVIRUS 2 AB.IGG',
					 'SARS CORONAVIRUS 2 AB.IGM',
					 'SARS-CoV IgG Ser Ql IA',
					 'SARS-CoV-2 (COVID-19) IgG Ab [Presence] in DBS by Immunoassay',
					 'SARS-CoV-2 (COVID-19) S protein RBD neut ab IA Ql',
					 'SARS-CoV-2 (COVID-19) specific TCRB gene rearrangements PrThr Pt Bld Ord Sequencing (Ab)',
					 'SARS-CoV-2 Ab DBS Ql IA',
					 'SARS-CoV-2 Ab Pnl SerPl IA',
					 'SARS-CoV-2 Ab Pnl SerPlBld IA.rapid',
					 'SARS-CoV-2 Ab SerPl IA-aCnc',
					 'SARS-CoV-2 Ab SerPlBld Ql IA.rapid',
					 'SARS-CoV-2 Ab SerPl-Imp',
					 'SARS-CoV-2 IgG',
					 'SARS-CoV-2 IgG + IgM',
					 'SARS-CoV-2 IgG SerPl IA-aCnc',
					 'SARS-CoV-2 IgM',
					 'SARS-CoV-2 IgM SerPl IA-aCnc',
					 'SARS-CoV-2 NAb Titr Ser pVNT')
				
		AND
			(
		(RESULT = 'Positive' OR RESULT = 'Detected' OR RESULT = 'Above Threshold' OR RESULT = 'Reactive'
			OR RESULT_VALUE LIKE 'Pos%' OR RESULT_VALUE LIKE 'Det%' OR RESULT_VALUE LIKE '%above thresh%' OR RESULT_VALUE LIKE 'React%') 
			)
		)
		
		THEN 1
		ELSE 0
		END AS 'Ab Positive Count',

		/*Inconclusive or Indeterminate results that will pass qc*/
		CASE WHEN
		(
			(
			(TEST IN ('SARS-CoV+SARS-CoV-2 Ag Resp Ql IA.rapid',
					  'SARS-CoV-2 (COVID-19) Ag [Presence] in Respiratory specimen by Rapid immunoassay',
					  'SARS-CoV-2 (COVID-19) Ag [Presence] in Upper respiratory specimen by Immunoassay',
					  'SARS-CoV-2 (COVID-19) Ag panel and Influenza virus A and B - Upper respiratory specimen by Rapid immunoassay',
					  'SARS-CoV-2 Ag Upper resp Ql IA.rapid',
				      'SARS coronavirus 2 Ab [Presence] in Serum or Plasma by Immunoassay',
					  'SARS coronavirus 2 Ab PrThr Pt Ser/Plas Ord IA',
				      'SARS CORONAVIRUS 2 AB.IGA',
					  'SARS coronavirus 2 Ab.IgA PrThr Pt Ser/Plas  Ord IA',
					  'SARS CORONAVIRUS 2 AB.IGG',
					  'SARS CORONAVIRUS 2 AB.IGM',
					  'SARS-CoV IgG Ser Ql IA',
					  'SARS-CoV-2 (COVID-19) IgG Ab [Presence] in DBS by Immunoassay',
					  'SARS-CoV-2 (COVID-19) S protein RBD neut ab IA Ql',
					  'SARS-CoV-2 (COVID-19) specific TCRB gene rearrangements PrThr Pt Bld Ord Sequencing (Ab)',
					  'SARS-CoV-2 Ab DBS Ql IA',
					  'SARS-CoV-2 Ab Pnl SerPl IA',
					  'SARS-CoV-2 Ab Pnl SerPlBld IA.rapid',
					  'SARS-CoV-2 Ab SerPl IA-aCnc',
					  'SARS-CoV-2 Ab SerPlBld Ql IA.rapid',
					  'SARS-CoV-2 Ab SerPl-Imp',
					  'SARS-CoV-2 IgG',
					  'SARS-CoV-2 IgG + IgM',
					  'SARS-CoV-2 IgG SerPl IA-aCnc',
					  'SARS-CoV-2 IgM',
					  'SARS-CoV-2 IgM SerPl IA-aCnc',
					  'SARS-CoV-2 NAb Titr Ser pVNT')
		
		AND
		
		(
		(RESULT IN ('Inconclusive', 'Indeterminate', 'Invalid', 'Equivocal', 'No RNA detected', 'Specimen unsatisfactory for evaluation') OR RESULT IS NULL
		)
		OR 
		(
		RESULT_VALUE IN ('Inconclusive', 'Indeterminate', 'Invalid', 'Equivocal', 'No RNA detected', 'Specimen unsatisfactory for evaluation') OR RESULT_VALUE IS NULL
		)
		)
		
		AND 
		
				(
		(labs.DISEASE_STATUS IN ('Not a case', 'Suspect', 'Probable'))
					)
				)
			)
		)
		THEN 'Pass'

		/*Contact cases that will pass QC*/
			WHEN
			(
				(labs.DISEASE_STATUS = 'Contact')
				AND
				(TEST IS NULL)
			)
			THEN 'Pass'


		ELSE 'Fail'
	END AS 'Classification QC'

FROM 
	COVID_LABS_TABLE AS labs
	JOIN COVID_CASES_TABLE AS cases
		ON	labs.PATIENT_ID = cases.PATIENT_ID

WHERE
	(labs.DATE BETWEEN '{previous_month_start}' AND '{previous_month_end}')
	
	AND
	
	(cases.STATE = 'SC')
	
	AND
	
	(TEST IN	(/*All tests except WGS which have a separate QC report*/
				'2019 Novel Coronavirus E gene [Presence] in Unspecified specimen by NAA with probe detection', 
				'2019 Novel Coronavirus N gene [Presence] in Unspecified specimen by NAA with probe detection', 
				'2019 Novel Coronavirus N gene [Presence] in Unspecified specimen by Nucleic acid amplification using primer-probe set N1', 
				'2019 Novel Coronavirus N gene [Presence] in Unspecified specimen by Nucleic acid amplification using primer-probe set N2', 
				'2019 Novel Coronavirus RdRp gene [Presence] in Unspecified specimen by NAA with probe detection', 
				'2019 Novel Coronavirus RNA panel - Unspecified specimen by NAA with probe detection', 
				'2019-nCoV RNA XXX NAA+probe-Imp: Imp: Pt: XXX: Ord: Probe.amp.tar', 
				'COVID-19 N gene in Saliva by NAA with probe detection',
				'COVID-19 RNA in Saliva by NAA with probe detection', 
				'COVID-19 RNA in Saliva by sequencing', 
				'FLUABV + SARS-CoV-2 Resp NAA+probe', 
				'FLUABV+SARS-CoV-2+RSV Pnl Resp NAA+probe',
				'HCoV 229E RNA Nph Ql NAA+non-probe: PrThr: Pt: Nph: Ord: Non-probe.amp.tar',
				'HCoV HKU1 RNA Nph Ql NAA+non-probe: PrThr: Pt: Nph: Ord: Non-probe.amp.tar', 
				'HCoV NL63 RNA Nph Ql NAA+non-probe: PrThr: Pt: Nph: Ord: Non-probe.amp.tar', 
				'HCoV OC43 RNA Nph Ql NAA+non-probe: PrThr: Pt: Nph: Ord: Non-probe.amp.tar', 
				'Influenza virus A and B and SARS-CoV-2 (COVID-19) RNA panel NAA+probe (Resp)', 
				'SARS coronavirus 2 ORF1ab region [Presence] in Respiratory specimen by NAA with probe detection', 
				'SARS coronavirus 2 RNA PrThr Pt Nph', 
				'SARS coronavirus 2 RNA PrThr Pt Respiratory Ord Probe.amp.tar', 
				'SARS coronavirus RNA [Presence] in Unspecified specimen by NAA with probe detection', 
				'SARS-CoV-2 (COVID-19) and SARS-related CoV RNA panel and Influenza virus A and B – Respiratory specimen by NAA with probe detectionActive', 
				'SARS-CoV-2 (COVID-19) N gene [Presence] in Respiratory specimen by NAA with probe detection', 
				'SARS-CoV-2 (COVID-19) RdRp gene NAA+probe (Unsp spec) [ThreshNum]', 
				'SARS-CoV-2 (COVID-19) RNA NAA+probe Ql (Nph)', 
				'SARS-CoV-2 N gene Nose Ql NAA+non-probe', 
				'SARS-CoV-2 N gene Nose Ql NAA+probe', 
				'SARS-CoV-2 N gene Nph Ql NAA+probe', 
				'SARS-CoV-2 N gene Resp Ql NAA N1', 
				'SARS-CoV-2 N gene Resp Ql NAA N2', 
				'SARS-CoV-2 N gene Sal Ql NAA N1', 
				'SARS-CoV-2 ORF1ab Sal Ql NAA+probe', 
				'SARS-CoV-2 RdRp Resp Ql NAA+probe', 
				'SARS-CoV-2 RNA Nose Ql NAA+probe', 
				'SARS-CoV-2 RNA Resp Ql NAA+non-probe', 
				'SARS-CoV-2 RNA Resp Ql NAA+probe', 
				'SARS-CoV-2 RNA Resp Ql Seq', 
				'SARS-like Coronavirus N gene [Cycle Threshold #] in Unspecified specimen by NAA with probe detection', 
				'SARS-like Coronavirus N gene [Presence] in Unspecified specimen by NAA with probe detection', 
				'SARS-related CoV RNA NAA+probe Ql (Resp)',
				'SARS-CoV+SARS-CoV-2 Ag Resp Ql IA.rapid',
				'SARS-CoV-2 (COVID-19) Ag [Presence] in Respiratory specimen by Rapid immunoassay',
				'SARS-CoV-2 (COVID-19) Ag [Presence] in Upper respiratory specimen by Immunoassay',
				'SARS-CoV-2 (COVID-19) Ag panel and Influenza virus A and B - Upper respiratory specimen by Rapid immunoassay',
				'SARS-CoV-2 Ag Upper resp Ql IA.rapid',
				'SARS coronavirus 2 Ab [Presence] in Serum or Plasma by Immunoassay',
				'SARS coronavirus 2 Ab PrThr Pt Ser/Plas Ord IA',
				'SARS CORONAVIRUS 2 AB.IGA',
				'SARS coronavirus 2 Ab.IgA PrThr Pt Ser/Plas  Ord IA',
				'SARS CORONAVIRUS 2 AB.IGG',
				'SARS CORONAVIRUS 2 AB.IGM',
				'SARS-CoV IgG Ser Ql IA',
				'SARS-CoV-2 (COVID-19) IgG Ab [Presence] in DBS by Immunoassay',
				'SARS-CoV-2 (COVID-19) S protein RBD neut ab IA Ql',
				'SARS-CoV-2 (COVID-19) specific TCRB gene rearrangements PrThr Pt Bld Ord Sequencing (Ab)',
				'SARS-CoV-2 Ab DBS Ql IA',
				'SARS-CoV-2 Ab Pnl SerPl IA',
				'SARS-CoV-2 Ab Pnl SerPlBld IA.rapid',
				'SARS-CoV-2 Ab SerPl IA-aCnc',
				'SARS-CoV-2 Ab SerPlBld Ql IA.rapid',
				'SARS-CoV-2 Ab SerPl-Imp',
				'SARS-CoV-2 IgG',
				'SARS-CoV-2 IgG + IgM',
				'SARS-CoV-2 IgG SerPl IA-aCnc',
				'SARS-CoV-2 IgM',
				'SARS-CoV-2 IgM SerPl IA-aCnc',
				'SARS-CoV-2 NAb Titr Ser pVNT'
				)	
			)
		")
  )

# START CLEANING #####

print(paste("Cleaning started:", date()))
print(paste("Step 1: Computing total tests by CASE ID:", date()))

### tally total number of tests for each CASE ID, save as a separate data frame #####
QC_totaltests <- 
  QC_dirty %>%
  group_by(PATIENT_ID) %>%
  summarize('Total tests' = length(TEST)) %>%
  as_tibble()

### join the new data frame with Total Tests variable to the original data frame #####
print(paste("Step 2: Creating clean dataframe:", date()))
QC_clean <- 
  left_join(QC_dirty, QC_totaltests, by = "PATIENT_ID")

### tally total number of RNA tests for each CASE ID, save as separate data frame #####
print(paste("Step 3: Computing total RNA tests by CASE ID:", date()))
QC_totalRNA <- 
  QC_dirty %>%
  group_by(PATIENT_ID) %>%
  summarize('Total RNA Positives' = sum(`RNA Positive Count`)) %>%
  as_tibble()

### join the new RNA variable to the clean table #####
print(paste("Step 4: Joining RNA tests to clean dataframe:", date()))
QC_clean %<>%
  left_join(QC_totalRNA, by = "PATIENT_ID")

### tally total number of Ag tests for each CASE ID, save as separate data frame #####
print(paste("Step 5: Computing total Ag tests by CASE ID:", date()))
QC_totalAg <- 
  QC_dirty %>%
  group_by(PATIENT_ID) %>%
  summarize('Total Ag Positives' = sum(`Ag Positive Count`)) %>%
  as_tibble()

### join the new Ag variable to the clean table #####
print(paste("Step 6: Joining Ag tests to clean dataframe:", date()))
QC_clean %<>%
  left_join(QC_totalAg, by = "PATIENT_ID")

### tally total number of Ab tests for each CASE ID, save as separate data frame #####
print(paste("Step 7: Computing total Ab tests by CASE ID:", date()))
QC_totalAb <- 
  QC_dirty %>%
  group_by(PATIENT_ID) %>%
  summarize('Total Ab Positives' = sum(`Ab Positive Count`)) %>%
  as_tibble()

### join the new Ab variable to the clean table #####
print(paste("Step 8: Joining Ab tests to clean dataframe:", date()))
QC_clean %<>%
  left_join(QC_totalAb, by = "PATIENT_ID")

### clean up weird character values in RESULT_VALUE field containing <, >, or = #####
print(paste("Step 9 Cleaning RESULT_VALUE field (<, >, = values to Pass):", date()))
QC_clean %<>%
  mutate(
    `Classification QC` = case_when(
      grepl("^=[1-9]+ | ^>[1-9]+", RESULT_VALUE, ignore.case = TRUE) ~ "Pass", TRUE ~ NA
    )
  )

# FINISH CLEANING #####
print(paste("Cleaning finished:", date()))

# compute final QC #####
print(paste("Step 10: Computing Final QC variable:", date()))
QC_clean %<>%
  mutate(
    `Classification QC` = if_else(
      (`Total RNA Positives` >= 1 & DISEASE_STATUS == 'Confirmed') |
      (`Total RNA Positives` == 0 & `Total Ag Positives` >= 1 & `Total Ab Positives` == 0 & DISEASE_STATUS %in% c('Probable')) |
      (`Total RNA Positives` == 0 & `Total Ag Positives` >= 1 & `Total Ab Positives` >= 1 & DISEASE_STATUS %in% c('Probable')) |
      (`Total RNA Positives` == 0 & `Total Ag Positives` == 0 & `Total Ab Positives` >= 1 & DISEASE_STATUS %in% c('Suspect', 'Not a case')) |
      (`Total RNA Positives` == 0 & `Total Ag Positives` == 0 & `Total Ab Positives` == 0 & DISEASE_STATUS %in% c('Suspect', 'Not a case'))
      , "Pass", "Fail")
  ) %>% 
  relocate(
    `Classification QC`, .after = last_col()
  )

# save results #####
print(paste("Saving results:", date()))
write_xlsx(QC_clean, save_path)

# FINISH #####
print(paste("Finished:", date()))
