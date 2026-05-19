/*==============================================================================
  KDHS 2014 DATA CLEANING DO-FILE
  Thesis: ML vs REML Estimation in Multilevel Models of Maternal Healthcare
          Underutilization in Kenya
  Author:   [Your Name]
  Date:     April 2026
  Stata:    18+ (uses melogit, margins, mixed syntax)
  
  STRUCTURE OF THIS DO-FILE:
    PART 0 — Global settings and file paths
    PART 1 — Births Recode (BR): Outcome construction (transition probabilities)
    PART 2 — Individual Recode (IR): Behavioral & individual-level determinants
    PART 3 — Household Member Recode (PR): Community education & literacy
    PART 4 — Household Recode (HR): Community poverty & infrastructure
    PART 5 — Merge all files into analysis dataset
    PART 6 — Community-level aggregate construction (cluster level)
    PART 7 — Final variable labelling, checks, and save
    
  KEY DESIGN DECISIONS:
    - Base unit: Last birth in past 5 years (BR file, index birth)
    - Outcome: 4 sequential transition probabilities (continuum of care)
    - Community variables: Constructed from cluster aggregates (HR + PR + IR)
    - Cluster ID: v001 (primary sampling unit = community level 2)
    - Survey weights applied at analysis stage, NOT cleaning stage
==============================================================================*/


/*==============================================================================
  PART 0 — GLOBAL SETTINGS AND FILE PATHS
==============================================================================*/

version 18
clear all
set more off
set linesize 120
capture log close

* ── Edit these paths to match your folder structure ──────────────────────────
global root  "C:\Users\norah.ogonji_evidenc\Desktop\Norah_Learning\Strathmore_Semester Study\MY MAsters Project\DHS Data\2014"
global br    "$root\KEBR72DT"
global ir    "$root\KEIR72DT"
global hr    "$root\KEHR72DT"
global pr    "$root\KEPR72DT"
global clean "C:\Users\norah.ogonji_evidenc\Desktop\Norah_Learning\Strathmore_Semester Study\MY MAsters Project\DHS Data\2-Clean Data"
global logs  "$root\logs"

*cap mkdir "$clean"
cap mkdir "$logs"

log using "$logs\\kdhs2014_prep_`c(current_date)'.log", replace text
*log using "$log/KDHS2014_cleaning_log.txt", text replace

/*
  FILE NAMES — standard DHS Kenya 2014 naming convention:
    BR: KEBR70FL.dta   (Births Recode)
    IR: KEIR70FL.dta   (Individual/Women's Recode)
    PR: KEPR70FL.dta   (Household Member Recode)
    HR: KEHR70FL.dta   (Household Recode)
    
  NOTE: If your files use a different prefix (e.g., KE_2014_DHS), 
        update the filenames below accordingly.
*/

global br_file  "KEBR72FL.DTA"
global ir_file  "KEIR72FL.DTA"
global pr_file  "KEPR72FL.DTA"
global hr_file  "KEHR72FL.DTA"

di "======================================================================"
di "  KDHS 2014 CLEANING BEGINS — $(c(current_date)) $(c(current_time))"
di "======================================================================"


/*==============================================================================
  PART 1 — BIRTHS RECODE (BR FILE)
  
  Purpose:  Construct the 4-stage continuum of care transition variables.
            One row per birth. We restrict to:
              (a) Last/index birth in past 5 years (b5=1 & b19<=59)
              (b) Non-missing outcome components
              
  Transitions modeled:
    T1: ANC initiation        — did she attend ≥1 ANC visit?
    T2: ANC adequacy          — conditional on T1=1, did she get 4+ visits?
    T3: Skilled/facility delivery — conditional on T2=1, skilled attendant?
    T4: Postnatal care (PNC)  — conditional on T3=1, PNC within 48 hours?
    
  Key BR variables used:
    b19    — child age in months (use <60 for last 5 years)
    b5     — child alive (1=yes; use all births, not just alive)
    m14    — number of antenatal visits
    m15    — place of delivery
    m3a–m3n — who assisted at delivery
    m62    — timing of postnatal check for mother (hours/days)
    m66    — postnatal check within 2 days (constructed)
    v001   — cluster number (community identifier)
    v002   — household number
    v003   — respondent line number
    caseid — unique case identifier (links to IR)
==============================================================================*/

di _newline "--- PART 1: BIRTHS RECODE ---"
use "$br/$br_file", clear

di "BR file loaded. Observations: `c(N)', Variables: `c(k)'"

/* ---- 1.1 RESTRICT TO INDEX BIRTH (last birth, past 5 years) -------------- */
/*
  DHS standard: bord = birth order; bidx = birth index within interview
  bidx==1 identifies the most recent birth within the woman's record.
  b19 = months since birth; restrict to <60 (5 years) per DHS convention.
*/

keep if bidx == 1            // most recent birth only
*keep if b19 < 60             // born in last 5 years
keep if b5 == 1 | b5 == 0    // keep all (alive and dead; ANC happens before death)

di "After restricting to index birth (last 5 yrs): N = `c(N)'"

/* ---- 1.2 ANC VISITS — m14 ------------------------------------------------ */
/*
  m14: number of antenatal visits during pregnancy
  Values: 0=none, 1–90=visits, 98=don't know, 99=missing
  
  NOTE for 2014: DHS codes 98 as "don't know" — treat as missing.
  WHO 2002 standard (in use for 2014 data): 4+ visits = adequate
  WHO 2016 standard: 8+ visits — we use 4+ for 2014, 8+ for 2022.
  This divergence will be documented in your methodology chapter.
*/

/* Flag don't know / missing */
replace m14 = .  if m14 == 98 | m14 == 99

/* T1: ANC initiation — attended at least 1 visit */
gen     anc_any = 0  if m14 != .
replace anc_any = 1  if m14 >= 1 & m14 != .
label variable anc_any "T1: Any ANC visit (>=1)"
label define yesno 0 "No" 1 "Yes"
label values anc_any yesno

/* T2: ANC adequacy — 4+ visits (conditional; set to missing if no ANC) */
gen     anc4plus = .
replace anc4plus = 0  if anc_any == 1 & m14 < 4
replace anc4plus = 1  if anc_any == 1 & m14 >= 4 & m14 != .
label variable anc4plus "T2: 4+ ANC visits (conditional on any ANC)"
label values anc4plus yesno

/* ---- 1.3 SKILLED BIRTH ATTENDANT — m3a to m3n ---------------------------- */
/*
  DHS records who assisted at delivery using a series of binary flags:
    m3a = doctor
    m3b = nurse/midwife
    m3c = auxiliary nurse/midwife
    m3d = auxiliary midwife
    m3e = traditional birth attendant (TBA)
    m3f = other person
    m3n = no one
    
  Skilled birth attendant (SBA) = doctor OR nurse/midwife OR auxiliary nurse/midwife
  WHO classification for Kenya 2014 DHS.
*/

/* Check variables exist */
foreach v in m3a m3b m3c m3d m3e m3f m3n {
    capture confirm variable `v'
    if _rc != 0 {
        di "WARNING: Variable `v' not found in BR file"
    }
}

/* Recode to numeric if needed (some DHS files store as byte with labels) */
foreach v in m3a m3b m3c m3d m3e m3f m3n {
    capture replace `v' = 0 if `v' == .
}

/* Skilled birth attendant indicator */
gen     sba = 0
replace sba = 1  if m3a == 1 | m3b == 1 | m3c == 1  // doctor or nurse/midwife
label variable sba "T3: Skilled birth attendant (doctor/nurse/midwife)"
label values sba yesno

/* ---- 1.4 FACILITY DELIVERY — m15 ----------------------------------------- */
/*
  m15: place of delivery
  Values vary slightly between DHS rounds. Common Kenya 2014 codes:
    10–19 = home
    20–29 = government hospital/health centre/clinic/dispensary
    30–39 = private hospital/clinic
    96    = other
    
  Facility delivery = any institutional delivery (public or private)
  We construct BOTH sba and facility delivery; use facility delivery for T3
  as it is more directly policy-relevant and avoids recall bias on attendant.
*/

gen     facdelivery = .
replace facdelivery = 0  if m15 >= 10 & m15 <= 19  // home
replace facdelivery = 1  if m15 >= 20 & m15 <= 39  // any health facility
replace facdelivery = 0  if m15 == 96              // other (non-facility)
replace facdelivery = .  if m15 == 99 | m15 == .   // missing
label variable facdelivery "T3-alt: Facility delivery (any institution)"
label values facdelivery yesno

/*
  THESIS NOTE: You will use facdelivery as primary T3, sba as sensitivity check.
  In your methods section, state:
  "Facility delivery was used as the primary T3 indicator; skilled birth
   attendance was used in sensitivity analyses."
*/

/* T3 primary: facility delivery, conditional on anc4plus */
gen     t3_delivery = .
replace t3_delivery = facdelivery  if anc4plus != .
label variable t3_delivery "T3: Facility delivery (conditional on 4+ ANC)"
label values t3_delivery yesno

/* ---- 1.5 POSTNATAL CARE — m62, m66, m67 ----------------------------------- */
/*
  DHS 2014 Kenya: postnatal care variables
    m62: timing of first postnatal check for MOTHER (not baby)
         Codes: 0=during delivery stay, 100–148 = hours after birth
                200–248 = days after birth, 300–348 = weeks after birth
                993=not with health provider, 994=never, 995=still in facility
                998=don't know, 999=missing
                
    WHO standard: PNC within 48 hours = adequate
    48 hours = within first 2 days post-delivery
    
  NOTE: Some DHS files also have m66 (check within 2 days, direct binary).
        Check if m66 exists in your file; if so, use it directly.
        If not, construct from m62.
*/

capture confirm variable m66
if _rc == 0 {
    /* m66 exists — use directly */
    gen     pnc48 = .
    replace pnc48 = 0  if m66 == 0
    replace pnc48 = 1  if m66 == 1
    di "  PNC: Using m66 (direct 2-day check variable)"
}
else {
    /* Construct from m62 */
    gen     pnc48 = .
    /* During delivery stay or within 48 hours (<= 148 in hour coding) */
    replace pnc48 = 1  if m62 == 0                         // during delivery stay
    replace pnc48 = 1  if m62 >= 100 & m62 <= 148         // within 48 hours
    replace pnc48 = 0  if m62 >= 149 & m62 <= 248         // 49+ hours to 48 days
    replace pnc48 = 0  if m62 >= 249 & m62 <= 395         // >48 days (weeks coding)
    replace pnc48 = 0  if m62 == 993 | m62 == 994         // no check / never
    replace pnc48 = .  if m62 == 995                       // still in facility — ambiguous
    replace pnc48 = .  if m62 == 998 | m62 == 999 | m62 == . // missing/DK
    di "  PNC: Constructed from m62 (hours/days since birth)"
}

label variable pnc48 "T4: Postnatal care within 48 hours (mother)"
label values pnc48 yesno

/* T4: PNC conditional on facility delivery */
gen     t4_pnc = .
replace t4_pnc = pnc48  if t3_delivery != .
label variable t4_pnc "T4: PNC within 48hrs (conditional on facility delivery)"
label values t4_pnc yesno

/* ---- 1.6 CONTINUUM OF CARE SCORE (unconditional, for descriptives) -------- */
/*
  Sum of 4 components: anc_any + anc4plus + facdelivery + pnc48
  Range 0–4. Used for descriptive purposes only.
  Transition probability models use conditional versions above.
*/

gen coc_score = 0
replace coc_score = coc_score + 1  if anc_any     == 1
replace coc_score = coc_score + 1  if anc4plus    == 1
replace coc_score = coc_score + 1  if facdelivery == 1
replace coc_score = coc_score + 1  if pnc48       == 1
replace coc_score = .  if anc_any == . & anc4plus == . & facdelivery == . & pnc48 == .
label variable coc_score "Continuum of care score (0-4)"

/* Binary underutilization composite (for sensitivity / robustness models) */
gen     underutilized = .
replace underutilized = 0  if coc_score >= 3 & coc_score != .  // adequate
replace underutilized = 1  if coc_score < 3  & coc_score != .  // underutilized
label variable underutilized "Underutilized (CoC score <3): 1=underutilized"
label values underutilized yesno

/*
  THESIS NOTE:
  The binary 'underutilized' variable is used for:
    (a) Prevalence estimation (Objective 1)
    (b) Sensitivity model comparing to transition-probability approach
  Primary analytical framework remains the 4-stage transitions (T1–T4).
*/

/* ---- 1.7 DROP-OFF INDICATORS (descriptive) -------------------------------- */
/*
  For Table: "Where women drop out of the continuum"
  Useful for policy narrative in your results chapter.
*/

gen dropoff_stage = .
replace dropoff_stage = 0  if anc_any     == 1 & anc4plus == 1 & facdelivery == 1 & pnc48 == 1
replace dropoff_stage = 1  if anc_any     == 0                   // dropped at T1
replace dropoff_stage = 2  if anc_any     == 1 & anc4plus == 0   // dropped at T2
replace dropoff_stage = 3  if anc4plus    == 1 & facdelivery == 0 // dropped at T3
replace dropoff_stage = 4  if facdelivery == 1 & pnc48 == 0      // dropped at T4

label variable dropoff_stage "Stage where CoC dropped (0=complete)"
label define dropstage 0 "Complete CoC" 1 "Dropped: No ANC" ///
    2 "Dropped: ANC <4" 3 "Dropped: No facility delivery" 4 "Dropped: No PNC"
label values dropoff_stage dropstage

/* ---- 1.8 KEY LINKING VARIABLES FROM BR ------------------------------------ */
/*
  We keep only what we need for the merge with IR.
  caseid links BR to IR (one-to-one; both have one row per woman per birth).
  v001 = cluster number (critical for multilevel structure).
*/

/* Additional birth-level controls from BR */

/* Parity — total children ever born (proxy behavioral fatigue) */
/* NOTE: In IR we have v220 (living children) — more appropriate */
/* In BR, bord gives birth order of this specific birth */
gen parity_br = bord
label variable parity_br "Birth order of this birth (from BR)"

/* Child sex (for covariate balance checks) */
gen     child_male = .
replace child_male = 0  if b4 == 2   // female
replace child_male = 1  if b4 == 1   // male
label variable child_male "Child is male (1=yes)"
label values child_male yesno

/* Preceding birth interval (months) */
gen     prev_binterval = b11
replace prev_binterval = .  if b11 == 9996 | b11 == 9998 | b11 == 9999
label variable prev_binterval "Preceding birth interval (months)"

/* ---- 1.9 KEEP VARIABLES FROM BR FOR MERGE --------------------------------- */

keep caseid v001 v002 v003 bidx  b5 b4 ///
     anc_any anc4plus facdelivery sba pnc48 ///
     t3_delivery t4_pnc coc_score underutilized dropoff_stage ///
     parity_br child_male prev_binterval m14

/* Quick check */
di _newline "BR outcomes summary:"
tab anc_any,      missing
tab anc4plus,     missing
tab facdelivery,  missing
tab pnc48,        missing
tab underutilized, missing
tab dropoff_stage, missing

/* Save BR cleaned */
save "$clean/br_2014_outcomes.dta", replace
di "BR 2014 cleaned and saved. N = `c(N)'"


/*==============================================================================
  PART 2 — INDIVIDUAL RECODE (IR FILE)
  Behavioral and individual-level determinants
  
  Variables constructed:
    BEHAVIORAL DETERMINANTS (framed as determinants, not covariates):
      edu_level      — woman's education level (v106)
      media_exposure — composite media exposure index (v157, v158, v159)
      autonomy_index — decision-making autonomy composite (v743a–v743f)
      parity         — number of living children (v220) — behavioral fatigue
      marital_status — current union status (v501)
      
    INDIVIDUAL SES (confounders):
      wealth_index   — household wealth quintile (v190)
      residence      — urban/rural (v025)
      region         — administrative region (v024)
      age_group      — woman's age in 5-year groups (v013)
      
    DESIGN/LINKING:
      v001           — cluster (PSU) — Level 2 in multilevel model
      v005           — survey weight
      v021           — primary sampling unit
      v022           — sampling stratum
      v023           — domain
==============================================================================*/

di _newline "--- PART 2: INDIVIDUAL RECODE ---"
use "$ir/$ir_file", clear
di "IR file loaded. N = `c(N)'"

/* ---- 2.1 RESTRICT TO WOMEN WITH A BIRTH IN LAST 5 YEARS ------------------ */
/*
  We only keep women who appear in the BR file (had a birth in last 5 years).
  This is done via the merge in Part 5 — here we keep all IR records
  but flag eligibility.
  
  v208 = number of births in last 5 years
*/

gen ir_eligible = (v208 >= 1 & v208 != .)
label variable ir_eligible "Has birth in last 5 years (1=eligible for analysis)"

/* ---- 2.2 EDUCATION — v106 ------------------------------------------------- */
/*
  v106: highest educational level
    0 = no education
    1 = primary
    2 = secondary
    3 = higher
    
  Framing in thesis: education as BEHAVIORAL DETERMINANT
  Mechanism: information-seeking behaviour, health literacy, empowerment
  Interaction tested: Education × community poverty
*/

gen     edu_level = v106
replace edu_level = .  if v106 == 9 | v106 == .
label variable edu_level "Education level (0=none, 1=primary, 2=secondary, 3=higher)"
label define edu_lbl 0 "No education" 1 "Primary" 2 "Secondary" 3 "Higher"
label values edu_level edu_lbl

/* Education binary (for interaction models) */
gen     edu_secondary = .
replace edu_secondary = 0  if edu_level <= 1  // none or primary
replace edu_secondary = 1  if edu_level >= 2 & edu_level != .  // secondary+
label variable edu_secondary "Secondary+ education (1=yes)"
label values edu_secondary yesno

/* ---- 2.3 MEDIA EXPOSURE — v157, v158, v159 -------------------------------- */
/*
  v157: reads newspaper (0=not at all, 1=<once/wk, 2=≥once/wk, 3=almost daily)
  v158: watches television (same codes)
  v159: listens to radio (same codes)
  
  Framing: media exposure as behavioral determinant
  Mechanism: health information access, demand creation, awareness of services
  
  We construct:
    (a) media_any    — any regular media exposure (≥1 item ≥weekly)
    (b) media_index  — count of media types accessed ≥weekly (0–3)
*/

foreach v in v157 v158 v159 {
    replace `v' = .  if `v' == 9 | `v' == .
}

/* Any regular media: at least one source accessed ≥ weekly */
gen     media_any = 0
replace media_any = 1  if v157 >= 2 | v158 >= 2 | v159 >= 2
replace media_any = .  if v157 == . & v158 == . & v159 == .
label variable media_any "Any media exposure ≥weekly (1=yes)"
label values media_any yesno

/* Media index: number of sources accessed at least weekly */
gen media_index = 0

replace media_index = media_index + 1  if v157 >= 2 & v157 != .
replace media_index = media_index + 1  if v158 >= 2 & v158 != .
replace media_index = media_index + 1  if v159 >= 2 & v159 != .
replace media_index = .  if v157 == . & v158 == . & v159 == .
label variable media_index "Media exposure index (0-3 sources >=weekly)"

/* ---- 2.4 AUTONOMY — v743a to v743f --------------------------------------- */
/*
  Decision-making autonomy variables (who has final say):
    v743a: own health care
    v743b: large household purchases
    v743c: visits to family/relatives
    v743d: food to be cooked daily (not always in all DHS versions)
    v743e: own earnings (if working)
    v743f: husband's earnings
    
  Codes: 1=respondent alone, 2=respondent and husband, 3=husband alone,
         4=someone else, 5=other, 6=not applicable
  
  Framing: autonomy as behavioral determinant (not just covariate)
  Mechanism: ability to decide to seek care independently
  
  We focus on v743a (own healthcare) as PRIMARY autonomy indicator
  — most directly linked to maternal healthcare seeking.
  We also construct a composite autonomy index.
*/

foreach v in v743a v743b v743c {
    capture confirm variable `v'
    if _rc == 0 {
        replace `v' = .  if `v' == 9 | `v' == 6 | `v' == .
    }
    else {
        di "WARNING: `v' not found in IR file"
    }
}

/* Healthcare autonomy: respondent has say in own healthcare */
gen     auto_health = .
capture {
    replace auto_health = 0  if v743a == 3 | v743a == 4 | v743a == 5  // no/little say
    replace auto_health = 1  if v743a == 1 | v743a == 2                // has say (alone or joint)
}
label variable auto_health "Autonomy: say in own healthcare (1=has say)"
label values auto_health yesno

/* Composite autonomy index (0–3 domains) */
gen autonomy_index = 0
capture replace autonomy_index = autonomy_index + 1  if v743a == 1 | v743a == 2
capture replace autonomy_index = autonomy_index + 1  if v743b == 1 | v743b == 2
capture replace autonomy_index = autonomy_index + 1  if v743c == 1 | v743c == 2
replace autonomy_index = .  if auto_health == . & v743b == . & v743c == .
label variable autonomy_index "Autonomy index (0-3 domains with say)"

/* ---- 2.5 PARITY — v220 --------------------------------------------------- */
/*
  v220: number of living children
  Framing: BEHAVIORAL FATIGUE — higher parity → reduced perceived need for care
  Mechanism: familiarity with pregnancy, reduced anxiety, opportunity cost
  
  Recode into groups to allow for non-linear effects:
    0–1 children (primipara / low parity)
    2–3 children (moderate parity)
    4+ children (high parity, behavioral fatigue most likely)
*/

gen     parity_grp = .
replace parity_grp = 0  if v220 <= 1
replace parity_grp = 1  if v220 >= 2 & v220 <= 3
replace parity_grp = 2  if v220 >= 4 & v220 != .
replace parity_grp = .  if v220 == 9 | v220 == .
label variable parity_grp "Parity group (0=0-1, 1=2-3, 2=4+)"
label define parity_lbl 0 "0-1 children" 1 "2-3 children" 2 "4+ children (high parity)"
label values parity_grp parity_lbl

/* Continuous parity (for interaction models and for mixed models) */
gen parity_cont = v220
replace parity_cont = .  if v220 == 9 | v220 == .
label variable parity_cont "Number of living children (continuous)"

/* ---- 2.6 MARITAL STATUS — v501 ------------------------------------------- */
/*
  v501: current marital status
    0 = never married
    1 = currently married
    2 = living together
    3 = widowed
    4 = divorced
    5 = not living together / separated
    
  Framing: BEHAVIORAL DETERMINANT
  Mechanism: partner support, social norms, financial access
*/

gen     marital_status = .
replace marital_status = 0  if v501 == 0                 // never married
replace marital_status = 1  if v501 == 1                 // married
replace marital_status = 2  if v501 == 2                 // cohabiting
replace marital_status = 3  if v501 >= 3 & v501 <= 5    // previously married
replace marital_status = .  if v501 == 9 | v501 == .
label variable marital_status "Marital status (0=never, 1=married, 2=cohabit, 3=prev)"
label define marital_lbl 0 "Never married" 1 "Currently married" ///
    2 "Living together" 3 "Previously married/separated"
label values marital_status marital_lbl

/* Binary: currently in union (married or cohabiting) */
gen     in_union = .
replace in_union = 0  if marital_status == 0 | marital_status == 3
replace in_union = 1  if marital_status == 1 | marital_status == 2
label variable in_union "Currently in union (married or cohabiting)"
label values in_union yesno

/* ---- 2.7 INDIVIDUAL SES CONTROLS ----------------------------------------- */

/* Age group — v013 (5-year groups, already coded in DHS) */
gen     age_group = v013
replace age_group = .  if v013 == 7 | v013 == .  // DHS uses 7 for 45-49 sometimes
label variable age_group "Age group (5-year categories, v013)"
label define age_lbl 1 "15-19" 2 "20-24" 3 "25-29" 4 "30-34" ///
    5 "35-39" 6 "40-44" 7 "45-49"
label values age_group age_lbl

/* Wealth index — v190 */
gen     wealth_q = v190
replace wealth_q = .  if v190 == 9 | v190 == .
label variable wealth_q "Wealth quintile (1=poorest, 5=richest)"
label define wealth_lbl 1 "Poorest" 2 "Poorer" 3 "Middle" 4 "Richer" 5 "Richest"
label values wealth_q wealth_lbl

/* Wealth binary (poor vs non-poor) for interaction models */
gen     wealth_poor = .
replace wealth_poor = 1  if wealth_q <= 2
replace wealth_poor = 0  if wealth_q >= 3 & wealth_q != .
label variable wealth_poor "Poor (bottom 2 wealth quintiles, 1=yes)"
label values wealth_poor yesno

/* Residence — v025 */
gen     residence = v025
replace residence = .  if v025 == 9 | v025 == .
label variable residence "Residence (1=urban, 2=rural)"
label define res_lbl 1 "Urban" 2 "Rural"
label values residence res_lbl

gen     rural = .
replace rural = 0  if v025 == 1  // urban
replace rural = 1  if v025 == 2  // rural
label variable rural "Rural residence (1=rural)"
label values rural yesno

/* Region — v024 */
gen     region = v024
replace region = .  if v024 == .
label variable region "Region (v024)"
/* Labels from DHS codebook — Kenya 2014 has 8 provinces/regions */

/* ---- 2.8 SURVEY DESIGN VARIABLES ----------------------------------------- */
/*
  These are NOT recoded — used as-is for svyset and for multilevel models.
  v001 = cluster number (PSU and Level-2 community identifier)
  v005 = individual sampling weight (divide by 1,000,000 for Stata)
  v021 = primary sampling unit
  v022 = sampling stratum
*/

gen wt = v005 / 1000000
label variable wt "Survey weight (v005/1,000,000)"
label variable v001 "Cluster number (PSU = community)"
label variable v021 "Primary sampling unit"
label variable v022 "Sampling stratum"

/* ---- 2.9 INTERACTION TERM INDICATORS FOR ANALYSIS STAGE ------------------ */
/*
  Pre-construct interaction terms to facilitate model building.
  Following your thesis framework:
    edu × community poverty  → edu_secondary × comm_poverty (merged later)
    wealth × rural context   → wealth_poor × rural
*/

/* Wealth × Rural interaction */
gen     wealth_rural = wealth_poor * rural
replace wealth_rural = .  if wealth_poor == . | rural == .
label variable wealth_rural "Interaction: Poor × Rural"

/* ---- 2.10 KEEP IR VARIABLES ----------------------------------------------- */

keep caseid v001 v002 v003 v005 v021 v022 wt ///
     edu_level edu_secondary media_any media_index ///
     auto_health autonomy_index parity_grp parity_cont ///
     marital_status in_union age_group wealth_q wealth_poor ///
     residence rural region wealth_rural ir_eligible ///
     v208 v024 v025 v190 v106 v013 v012

/* Quick check */
di _newline "IR behavioral determinants summary:"
foreach v in edu_level media_any auto_health parity_grp marital_status wealth_q residence {
    di "  `v':"
    tab `v', missing
}

save "$clean/ir_2014_determinants.dta", replace
di "IR 2014 cleaned and saved. N = `c(N)'"


/*==============================================================================
  PART 3 — HOUSEHOLD MEMBER RECODE (PR FILE)
  
  Purpose: Construct cluster-level community education and literacy variables.
           Restrict to women 15–49 to match target population.
           Aggregate to cluster level in Part 6.
           
  Key PR variables:
    hv001  = cluster number (matches v001 in IR/BR)
    hv104  = sex of household member
    hv105  = age of household member
    hv106  = highest education level
    hv110  = educational attainment (for alternative measure)
    hv121  = member currently attending school
    hv108  = years of education
==============================================================================*/

di _newline "--- PART 3: HOUSEHOLD MEMBER RECODE (PR) ---"
use "$pr/$pr_file", clear
di "PR file loaded. N = `c(N)'"

/* ---- 3.1 RESTRICT TO WOMEN 15–49 ----------------------------------------- */
keep if hv104 == 2                          // female
keep if hv105 >= 15 & hv105 <= 49          // age 15-49

di "After restricting to women 15-49: N = `c(N)'"

/* ---- 3.2 EDUCATION VARIABLES --------------------------------------------- */

/* Education level */
gen     edu_pr = hv106
replace edu_pr = .  if hv106 == 9 | hv106 == .

/* Secondary or higher (binary, for community proportion) */
gen     edu_secondary_pr = .
replace edu_secondary_pr = 0  if edu_pr <= 1
replace edu_secondary_pr = 1  if edu_pr >= 2 & edu_pr != .

/* Years of education */
gen     edu_years = hv108
replace edu_years = .  if hv108 == 97 | hv108 == 98 | hv108 == 99

/* ---- 3.3 LITERACY --------------------------------------------------------- */
/*
  Literacy proxy: secondary+ education
  Some PR files have hv110 = literacy (can read/write)
  Check availability:
*/

capture confirm variable hv110
if _rc == 0 {
    gen     literate = .
    replace literate = 0  if hv110 == 0 | hv110 == 3   // cannot read
    replace literate = 1  if hv110 == 1 | hv110 == 2   // can read
    replace literate = .  if hv110 == 9 | hv110 == .
    di "  Literacy: Using hv110 (direct literacy variable)"
}
else {

    gen literate = edu_secondary_pr
	
    di "  Literacy: Proxied by secondary+ education (hv110 not found)"
}

replace literate = edu_secondary_pr

/* ---- 3.4 KEEP FOR AGGREGATION -------------------------------------------- */

gen cluster_id = hv001    // PSU/cluster identifier
label variable cluster_id "Cluster ID (from PR = hv001)"

keep cluster_id edu_pr edu_secondary_pr edu_years literate

/* Save for aggregation in Part 6 */
save "$clean/pr_2014_education.dta", replace
di "PR 2014 (women 15-49) saved for aggregation. N = `c(N)'"


/*==============================================================================
  PART 4 — HOUSEHOLD RECODE (HR FILE)
  
  Purpose: Construct cluster-level community poverty and infrastructure variables.
           One row per household. Aggregate to cluster level in Part 6.
           
  Key HR variables:
    hv001  = cluster number
    hv270  = wealth index (quintile)
    hv025  = type of place of residence (urban/rural)
    hv206  = has electricity
    hv207  = has radio
    hv208  = has television
    hv221  = has mobile telephone
    hv204  = time to water source (minutes)
    hv226  = type of cooking fuel
    hv205  = toilet facility type
    hv210  = has bicycle
    hv211  = has motorcycle/scooter
    hv212  = has car/truck
==============================================================================*/

di _newline "--- PART 4: HOUSEHOLD RECODE (HR) ---"
use "$hr/$hr_file", clear
di "HR file loaded. N = `c(N)'"

/* ---- 4.1 POVERTY INDICATOR ----------------------------------------------- */

gen     wealth_hh = hv270
replace wealth_hh = .  if hv270 == 9 | hv270 == .

/* Poor = bottom 2 quintiles */
gen     poor_hh = .
replace poor_hh = 1  if wealth_hh <= 2
replace poor_hh = 0  if wealth_hh >= 3 & wealth_hh != .
label variable poor_hh "Household in poorest 2 wealth quintiles"

/* ---- 4.2 INFRASTRUCTURE INDICATORS --------------------------------------- */

/* Electricity */
gen     has_electricity = .
capture replace has_electricity = 0  if hv206 == 0
capture replace has_electricity = 1  if hv206 == 1

/* Improved toilet */
/*
  hv205 values: 10–15 = flush toilet, 21–23 = pit latrine improved
  Non-improved: 23 (open pit), 31 (no facility)
  JMP classification used here
*/
gen     improved_toilet = .
capture {
    replace improved_toilet = 0  if hv205 == 31 | hv205 == 30  // no facility/field
    replace improved_toilet = 1  if hv205 >= 10 & hv205 <= 29  // any facility
    replace improved_toilet = .  if hv205 == 99 | hv205 == .
}

/* Improved water source */
/*
  hv201: main source of drinking water
  Improved: piped (10–14), tube well/borehole (20–21), protected well (31),
            protected spring (41), rainwater (51), bottled (71)
*/
gen     improved_water = .
capture {
    replace improved_water = 1  if hv201 <= 14                    // piped
    replace improved_water = 1  if hv201 >= 20 & hv201 <= 21     // tube well
    replace improved_water = 1  if hv201 == 31 | hv201 == 41     // protected
    replace improved_water = 1  if hv201 == 51 | hv201 == 71     // rain/bottled
    replace improved_water = 0  if improved_water == . & hv201 != . & hv201 != 99
    replace improved_water = .  if hv201 == 99 | hv201 == .
}

/* ---- 4.3 RESIDENCE -------------------------------------------------------- */
gen     urban_hh = .
replace urban_hh = 1  if hv025 == 1  // urban
replace urban_hh = 0  if hv025 == 2  // rural

/* ---- 4.4 DISTANCE TO FACILITY PROXY ------------------------------------- */
/*
  hv204: time to water source in minutes (imperfect proxy for access)
  For health facility access: some DHS files have sh18a/sh18b (distance)
  Check availability — otherwise use electricity/road as infrastructure proxy.
*/

gen     water_time = hv204
replace water_time = .  if hv204 == 998 | hv204 == 999  // on premises coded 0

/* ---- 4.5 KEEP FOR AGGREGATION -------------------------------------------- */

gen cluster_id = hv001
label variable cluster_id "Cluster ID (from HR = hv001)"

keep cluster_id poor_hh wealth_hh has_electricity improved_toilet ///
     improved_water urban_hh water_time

save "$clean/hr_2014_household.dta", replace
di "HR 2014 saved for aggregation. N = `c(N)'"


/*==============================================================================
  PART 5 — MERGE BR + IR INTO INDIVIDUAL ANALYSIS FILE
  
  Strategy:
    1. Start with BR cleaned (outcomes) — one row per birth
    2. Merge IR determinants using caseid (one-to-one)
    3. Result: one row per woman per last birth, with outcomes + determinants
    4. Community variables added in Part 6 after cluster-level aggregation
==============================================================================*/

di _newline "--- PART 5: MERGING BR AND IR ---"
use "$clean/br_2014_outcomes.dta", clear

/* ---- 5.1 MERGE BR WITH IR ------------------------------------------------- */
/*
  caseid is the unique woman identifier in both BR and IR.
  Sort both files by caseid before merging.
  Expect 1:1 merge (one record per woman for most recent birth).
*/

sort caseid
merge 1:1 caseid using "$clean/ir_2014_determinants.dta", ///
    keep(match master) nogen

di "After BR-IR merge: N = `c(N)'"

/* ---- 5.2 CONSISTENCY CHECKS ---------------------------------------------- */

/* Check cluster IDs are identical across files */
count if v001 == .
di "  Obs with missing cluster ID: `r(N)'"

/* Check: women restricted to those with birth in last 5 years */
count if ir_eligible == 0
di "  IR records without birth (should be 0 after merge): `r(N)'"

/* Check parity consistency between BR and IR */
gen parity_diff = abs(parity_br - parity_cont)
count if parity_diff > 3 & parity_diff != .
di "  Large parity discrepancy between BR and IR: `r(N)' obs"
drop parity_diff

/* ---- 5.3 RENAME CLUSTER FOR MULTILEVEL MODELLING ------------------------- */
/*
  In multilevel models, we need a clean Level-2 identifier.
  v001 = cluster (PSU). This is our community level.
  Some models may also use v022 (stratum) as Level-3 if needed.
*/

rename v001 cluster_id
label variable cluster_id "Cluster ID (community, Level 2 in MLM)"

/* ---- 5.4 SAVE INDIVIDUAL-LEVEL MERGED FILE ------------------------------- */

save "$clean/kdhs2014_individual.dta", replace
di "Individual merged file saved. N = `c(N)'"


/*==============================================================================
  PART 6 — COMMUNITY-LEVEL AGGREGATE CONSTRUCTION
  
  We construct 5 cluster-level (community) variables:
  
    comm_poverty      — % households in poorest 2 wealth quintiles per cluster
    comm_edu          — % women 15-49 with secondary+ education per cluster
    comm_literacy     — % literate women 15-49 per cluster
    comm_urban        — % urban households per cluster
    comm_anc_coverage — % women in cluster with 4+ ANC (leave-one-out method)
    comm_infra        — infrastructure index (electricity + water + toilet)
    
  NOTE on leave-one-out (LOO) for comm_anc_coverage:
    Standard cluster aggregation would include the woman herself, creating
    mechanical correlation between individual ANC and community ANC.
    LOO removes each woman's own observation when computing her cluster mean.
    We implement LOO after the merge.
    
  THESIS NOTE:
    These variables operationalise your theoretical community constructs:
      Community poverty → social norms, financial barriers
      Community education → information environment, norm formation
      Community ANC coverage → local health system norm / supply proxy
      Infrastructure → physical access to services
==============================================================================*/

di _newline "--- PART 6: COMMUNITY VARIABLE CONSTRUCTION ---"

/* ---- 6.1 POVERTY AGGREGATE (from HR) -------------------------------------- */
use "$clean/hr_2014_household.dta", clear

bysort cluster_id: egen comm_poverty   = mean(poor_hh)
bysort cluster_id: egen comm_urban     = mean(urban_hh)
bysort cluster_id: egen comm_elec      = mean(has_electricity)
bysort cluster_id: egen comm_water     = mean(improved_water)
bysort cluster_id: egen comm_toilet    = mean(improved_toilet)

/* Infrastructure composite index (mean of 3 binary indicators) */
gen comm_infra = (comm_elec + comm_water + comm_toilet) / 3
replace comm_infra = . if comm_elec == . & comm_water == . & comm_toilet == .

/* Collapse to one row per cluster */
collapse (mean) comm_poverty comm_urban comm_elec comm_water ///
    comm_toilet comm_infra, by(cluster_id)

label variable comm_poverty "Community poverty (% poor HHs in cluster)"
label variable comm_urban   "Community urbanicity (% urban HHs)"
label variable comm_infra   "Community infrastructure index (0-1)"
label variable comm_elec    "% households with electricity in cluster"

save "$clean/community_hr_2014.dta", replace
di "Community HR aggregates saved. N clusters = `c(N)'"

/* ---- 6.2 EDUCATION AGGREGATE (from PR) ------------------------------------ */
use "$clean/pr_2014_education.dta", clear

bysort cluster_id: egen comm_edu      = mean(edu_secondary_pr)
bysort cluster_id: egen comm_literacy = mean(literate)
bysort cluster_id: egen comm_edu_yrs  = mean(edu_years)

collapse (mean) comm_edu comm_literacy comm_edu_yrs, by(cluster_id)

label variable comm_edu      "Community education (% women w/ secondary+)"
label variable comm_literacy "Community literacy (% literate women 15-49)"
label variable comm_edu_yrs  "Community mean years of education"

save "$clean/community_pr_2014.dta", replace
di "Community PR aggregates saved. N clusters = `c(N)'"

/* ---- 6.3 ANC COVERAGE AGGREGATE (from individual file — LOO method) ------ */
/*
  Leave-one-out community ANC coverage:
    For each woman i in cluster j:
    comm_anc_cov_LOO[i] = mean(anc4plus) among all women in cluster j EXCEPT i
    
  This avoids endogeneity: a woman's own ANC does not predict itself.
  Standard formula: LOO_mean = (cluster_sum - individual_value) / (cluster_n - 1)
*/

use "$clean/kdhs2014_individual.dta", clear

/* Step 1: Compute cluster sum and count of anc4plus */
bysort cluster_id: egen cluster_anc_sum = total(anc4plus), missing
bysort cluster_id: egen cluster_anc_n   = count(anc4plus)

/* Step 2: LOO mean */
gen comm_anc_loo = (cluster_anc_sum - anc4plus) / (cluster_anc_n - 1)
replace comm_anc_loo = . if cluster_anc_n <= 1  // can't compute LOO with 1 obs
replace comm_anc_loo = . if anc4plus == .        // missing individual, keep as missing

label variable comm_anc_loo "Community ANC coverage (LOO, % with 4+ ANC)"

/* Collapse to cluster level (LOO already individual-level, keep as is for merge) */
/* We keep at individual level — merge back after community aggregates */

drop cluster_anc_sum cluster_anc_n

save "$clean/kdhs2014_individual_loo.dta", replace

/* ---- 6.4 MERGE COMMUNITY AGGREGATES INTO INDIVIDUAL FILE ------------------ */

use "$clean/kdhs2014_individual_loo.dta", clear

/* Merge poverty and infrastructure (from HR) */
merge m:1 cluster_id using "$clean/community_hr_2014.dta", ///
    keep(master match) nogen keepusing(comm_poverty comm_urban comm_infra comm_elec)

di "After merging community HR aggregates: N = `c(N)'"

/* Merge education and literacy (from PR) */
merge m:1 cluster_id using "$clean/community_pr_2014.dta", ///
    keep(master match) nogen keepusing(comm_edu comm_literacy comm_edu_yrs)

di "After merging community PR aggregates: N = `c(N)'"

di "Community variables successfully merged."

/* ---- 6.5 INTERACTION TERMS WITH COMMUNITY VARIABLES ----------------------- */
/*
  Your thesis specifies two key interactions:
    1. Education × community poverty
    2. Wealth × rural context (already created as wealth_rural in Part 2)
    
  Construct these now as clean variables.
*/

/* Interaction 1: Individual education × Community poverty */
gen     edu_comm_poverty = edu_secondary * comm_poverty
replace edu_comm_poverty = .  if edu_secondary == . | comm_poverty == .
label variable edu_comm_poverty "Interaction: Edu_secondary × Comm_poverty"

/* Interaction 2: Wealth × Rural (already have wealth_rural from IR, but restate) */
/* wealth_rural was constructed in Part 2 — confirm it is present */
capture confirm variable wealth_rural
if _rc != 0 {
    gen wealth_rural = wealth_poor * rural
    replace wealth_rural = . if wealth_poor == . | rural == .
}
label variable wealth_rural "Interaction: Poor × Rural"

/* Interaction 3: Education × Community education (cross-level) */
gen     edu_comm_edu = edu_secondary * comm_edu
replace edu_comm_edu = .  if edu_secondary == . | comm_edu == .
label variable edu_comm_edu "Interaction: Edu_secondary × Comm_edu (cross-level)"


/*==============================================================================
  PART 7 — FINAL CLEANING, LABELLING, CHECKS, AND SAVE
==============================================================================*/

di _newline "--- PART 7: FINAL CHECKS AND SAVE ---"

/* ---- 7.1 MISSING DATA SUMMARY -------------------------------------------- */

di _newline "=== MISSING DATA SUMMARY ==="
di "OUTCOMES:"
foreach v in anc_any anc4plus facdelivery pnc48 underutilized {
    quietly count if `v' == .
    di "  `v': `r(N)' missing (`=string(round(`r(N)'/`c(N)'*100,.1))'%)"
}

di _newline "BEHAVIORAL DETERMINANTS:"
foreach v in edu_level media_any auto_health parity_grp marital_status {
    quietly count if `v' == .
    di "  `v': `r(N)' missing (`=string(round(`r(N)'/`c(N)'*100,.1))'%)"
}

di _newline "INDIVIDUAL SES:"
foreach v in wealth_q residence region age_group {
    quietly count if `v' == .
    di "  `v': `r(N)' missing (`=string(round(`r(N)'/`c(N)'*100,.1))'%)"
}

di _newline "COMMUNITY VARIABLES:"
foreach v in comm_poverty comm_edu comm_literacy comm_anc_loo comm_infra comm_urban {
    quietly count if `v' == .
    di "  `v': `r(N)' missing (`=string(round(`r(N)'/`c(N)'*100,.1))'%)"
}

/* ---- 7.2 COMPLETE CASE INDICATOR ----------------------------------------- */
/*
  Flag observations with complete data on all key analysis variables.
  You will use this in sensitivity analyses (complete case vs. multiple imputation).
*/

gen complete_case = 1
foreach v in anc_any anc4plus facdelivery pnc48 ///
             edu_level media_any auto_health parity_grp marital_status ///
             wealth_q residence age_group ///
             comm_poverty comm_edu comm_anc_loo {
    replace complete_case = 0  if `v' == .
}
label variable complete_case "Complete case on all analysis variables (1=complete)"
label values complete_case yesno

quietly count if complete_case == 1
di _newline "Complete cases: `r(N)' of `c(N)' total (`=string(round(`r(N)'/`c(N)'*100,.1))'%)"

/* ---- 7.3 ANALYTICAL SAMPLE RESTRICTION ----------------------------------- */
/*
  PRIMARY ANALYTICAL SAMPLE:
    - Women with a birth in the last 5 years (already restricted via BR merge)
    - Age 15–49 (DHS eligibility criterion — already in IR)
    - Non-missing on at least T1 (ANC any) — minimum outcome requirement
    
  We create an analysis flag but do NOT drop observations at this stage.
  Apply restriction at analysis stage for transparency.
*/

gen analysis_sample = 1
replace analysis_sample = 0  if anc_any == .      // missing on T1
replace analysis_sample = 0  if cluster_id == .   // missing cluster (can't use in MLM)
replace analysis_sample = 0  if age_group == .    // missing age
label variable analysis_sample "Primary analysis sample (1=included)"
label values analysis_sample yesno

quietly count if analysis_sample == 1
di "Primary analysis sample: `r(N)' women"

/* ---- 7.4 SURVEY DESIGN DECLARATION (for descriptive prevalence estimation) */
/*
  Use svyset for Objective 1 (prevalence estimation with correct SEs).
  For multilevel models (Objectives 3-4), use melogit with pw(wt).
  
  Kenya 2014 DHS uses two-stage stratified cluster sampling:
    Strata = v022 (sampling strata)
    PSU    = v021 (primary sampling unit, often same as cluster)
    Weight = wt   (survey weight = v005/1,000,000)
*/

di _newline "Survey design:"
di "  svyset v021 [pw=wt], strata(v022) singleunit(centered)"

/* ---- 7.5 VARIABLE LABELS FINALISATION ------------------------------------ */

label variable caseid      "Case ID (unique woman identifier)"
label variable cluster_id  "Cluster ID (PSU = community, Level 2)"
label variable v021        "Primary sampling unit"
label variable v022        "Sampling stratum"
label variable wt          "Analytic weight (v005/1,000,000)"

/* ---- 7.6 CREATE VARIABLE CODEBOOK ---------------------------------------- */

di _newline "=== FINAL VARIABLE LIST ==="
describe

/* ---- 7.7 BASIC DESCRIPTIVE CHECK ----------------------------------------- */
/*
  Quick tabulations to confirm recodes are sensible before saving.
  These are for YOUR review — do not include raw output in thesis.
*/

di _newline "=== QUICK VALIDATION TABULATIONS ==="
di "T1 ANC initiation:"
tab anc_any if analysis_sample == 1, missing
di "T2 ANC 4+:"
tab anc4plus if analysis_sample == 1, missing
di "T3 Facility delivery:"
tab facdelivery if analysis_sample == 1, missing
di "T4 PNC 48hrs:"
tab pnc48 if analysis_sample == 1, missing

di _newline "Education:"
tab edu_level if analysis_sample == 1, missing
di "Wealth:"
tab wealth_q if analysis_sample == 1, missing
di "Residence:"
tab residence if analysis_sample == 1, missing
di "Parity group:"
tab parity_grp if analysis_sample == 1, missing

/* ---- 7.8 CLUSTER COUNT CHECK FOR MLM -------------------------------------- */
/*
  Multilevel models require adequate cluster representation.
  Rule of thumb: ≥30 clusters, ≥5 obs per cluster preferred.
  Document in your methods section.
*/

preserve
keep if analysis_sample == 1
bysort cluster_id: gen n_in_cluster = _N
quietly su n_in_cluster
di _newline "=== CLUSTER DIAGNOSTICS ==="
*di "Number of clusters: " quietly distinct cluster_id; di `r(ndistinct)'
di "Mean women per cluster: `r(mean)'"
di "Min women per cluster:  `r(min)'"
di "Max women per cluster:  `r(max)'"
restore

/* ---- 7.9 SAVE FINAL ANALYSIS FILE ---------------------------------------- */

save "$clean/KDHS2014_analysis_final.dta", replace

di _newline "======================================================================"
di "  KDHS 2014 CLEANING COMPLETE"
di "  Final dataset: $clean/KDHS2014_analysis_final.dta"
di "  N (total): `c(N)'"
quietly count if analysis_sample == 1
di "  N (analysis sample): `r(N)'"
quietly distinct cluster_id
di "  Clusters: `r(ndistinct)'"
di "  Completed: $(c(current_date)) $(c(current_time))"
di "======================================================================"

log close


/*==============================================================================
  APPENDIX — VARIABLE CROSSWALK TABLE FOR 2022 COMPARISON
  
  The following variables will CHANGE between 2014 and 2022.
  Document these in your methodology chapter.
  
  VARIABLE         | 2014 (KDHS 6)       | 2022 (KDHS 8)       | ACTION
  -----------------+---------------------+---------------------+------------------
  ANC visits       | m14 (max ~90)       | m14 (similar)       | Same recode
  ANC standard     | WHO 2002: 4+        | WHO 2016: 8+        | USE 4+ for BOTH
                   |                     |                     | (sensitivity: 8+)
  PNC timing       | m62 (hours/days)    | m62 / m66 direct    | Check availability
  Delivery place   | m15 (codes 10-39)   | m15 (expanded)      | Map codes carefully
  Autonomy         | v743a-v743f         | v743a-v743f         | Same
  Wealth index     | v190 (quintile)     | v190 (quintile)     | Same
  Media            | v157/158/159        | v157/158/159        | Same
  Region           | v024 (8 regions)    | v024 (47 counties)  | RECODE 2022 to 8
  Cluster          | v001                | v001                | Same
  Weight           | v005/1000000        | v005/1000000        | Same
  
  NOTE: 2022 KDHS (KDHS 8) uses county-level strata (47 counties).
  You will need to harmonise region coding for cross-time comparison.
  The cleaning do-file for 2022 will address this explicitly.
==============================================================================*/

/*  END OF KDHS 2014 CLEANING DO-FILE  */
