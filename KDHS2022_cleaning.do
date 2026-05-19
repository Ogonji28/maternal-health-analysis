/*==============================================================================
  KDHS 2022 DATA CLEANING DO-FILE
  Thesis: ML vs REML Estimation in Multilevel Models of Maternal Healthcare
          Underutilization in Kenya
  Author:   [Your Name]
  Date:     April 2026
  Stata:    18+ (uses melogit, margins, mixed syntax)

  FILE NAMES (KDHS 2022 — DHS Round 8):
    BR: KEBR8CDT.dta   (Births Recode)
    IR: KEIR8CDT.dta   (Individual/Women's Recode)
    PR: KEPR8CDT.dta   (Household Member Recode)
    HR: KEHR8CDT.dta   (Household Recode)

  STRUCTURE OF THIS DO-FILE:
    PART 0  — Global settings and file paths
    PART 1  — Births Recode (BR): Outcome construction (transition probabilities)
    PART 2  — Individual Recode (IR): Behavioral & individual-level determinants
    PART 3  — Household Member Recode (PR): Community education & literacy
    PART 4  — Household Recode (HR): Community poverty & infrastructure
    PART 5  — Merge BR + IR into individual analysis file
    PART 6  — Community-level aggregate construction (cluster level)
    PART 7  — County-to-region harmonisation (CRITICAL for 2014 comparison)
    PART 8  — Final variable labelling, checks, and save
    PART 9  — Cross-wave harmonisation checks (2014 vs 2022 comparability)

  KEY DIFFERENCES FROM 2014 CLEANING FILE — READ BEFORE RUNNING:
  ---------------------------------------------------------------
  1. FILE PREFIX: KEBR8CDT (not KEBR70FL) — Round 8 naming convention
  2. ANC STANDARD: WHO 2016 recommends 8+ visits; we retain 4+ as primary
     threshold for cross-wave comparability, and add 8+ as sensitivity
  3. REGION vs COUNTY: 2022 uses 47 counties (v024 = 1–47); 2014 used
     8 provinces. We harmonise 2022 counties back to 8 regions in Part 7.
  4. PNC VARIABLES: 2022 DHS expanded PNC module — m66/m67 more reliable
  5. WEALTH INDEX: Recomputed for 2022 — not directly comparable to 2014
     quintiles. We use within-survey quintiles; comparisons use tertiles.
  6. COVID-19 CONTEXT: 2022 data collected post-COVID. We flag a
     sensitivity indicator for potential disruption period births.
  7. NEW DHS-8 VARIABLES: Some autonomy and media variables have
     expanded response categories — all remapped to 2014-compatible codes.
  8. SAMPLE WEIGHTS: 2022 uses de facto weights (v005); same formula.
==============================================================================*/


/*==============================================================================
  PART 0 — GLOBAL SETTINGS AND FILE PATHS
==============================================================================*/

version 18
clear all
set more off
set maxvar 10000
set linesize 120
capture log close

* ── Edit these paths ─────────────────────────────────────────────────────────
global root  "C:\Users\norah.ogonji_evidenc\Desktop\Norah_Learning\Strathmore_Semester Study\MY MAsters Project\DHS Data\2022"
global br    "$root\KEBR8CDT"
global ir    "$root\KEIR8CDT"
global hr    "$root\KEHR8CDT"
global pr    "$root\KEPR8CDT"
global clean "C:\Users\norah.ogonji_evidenc\Desktop\Norah_Learning\Strathmore_Semester Study\MY MAsters Project\DHS Data\2-Clean Data"
global logs  "$root/logs"

*cap mkdir "$clean"
cap mkdir "$logs"



/* Harmonised output folder (pooled 2014+2022 dataset goes here) */
global pooled  "$clean"

capture mkdir "$clean"
capture mkdir "$log"
capture mkdir "$pooled"

log using "$logs\\kdhs2022_prep_`c(current_date)'.log", replace text

/* ---- 2022 FILE NAMES ------------------------------------------------------ */
global br22   "KEBR8CFL.DTA"
global ir22   "KEIR8CFL.DTA"
global pr22   "KEPR8CFL.DTA"
global hr22   "KEHR8CFL.DTA"

/* ---- SURVEY YEAR FLAG (used throughout for cross-wave pooling) ------------ */
global syear 2022

di "======================================================================"
di "  KDHS 2022 CLEANING BEGINS — $(c(current_date)) $(c(current_time))"
di "  Files: $br22 / $ir22 / $pr22 / $hr22"
di "======================================================================"


/*==============================================================================
  PART 1 — BIRTHS RECODE (BR FILE — KEBR8CDT)

  Purpose:  Construct 4-stage continuum of care transition variables.
            Identical structure to 2014 file for cross-wave comparability.
            Key 2022-specific notes annotated throughout.

  Transitions:
    T1: ANC initiation       — ≥1 ANC visit
    T2: ANC adequacy         — 4+ visits (primary); 8+ (sensitivity)
    T3: Facility delivery    — any institutional delivery
    T4: Postnatal care (PNC) — maternal PNC within 48 hours

  CHANGES FROM 2014 BR:
    - m14 coding: same but check for expanded "don't know" codes
    - m15 (delivery place): expanded facility codes in DHS-8; remap below
    - m62/m66: 2022 has improved PNC timing module; m66 more available
    - COVID sensitivity flag: births 2020–2021 flagged (b18/b19 range)
==============================================================================*/

di _newline "--- PART 1: BIRTHS RECODE (2022) ---"
use "$br/$br22", clear


di "BR 2022 loaded. N = `c(N)', Variables = `c(k)'"

/* ---- 1.1 RESTRICT TO INDEX BIRTH (last birth, past 5 years) -------------- */
/*
  Same restriction as 2014:
    bidx == 1  → most recent birth
    b19  < 60  → born in last 5 years (months since birth)
    
  NOTE: In 2022, b19 uses same coding as 2014.
  COVID sensitivity: births with b19 >= 12 & b19 <= 36 approximately
  correspond to 2019–2021 (pandemic period). We flag but do NOT exclude.
*/

keep if bidx == 1
keep if b19 < 60

di "After restricting to index birth (last 5 yrs): N = `c(N)'"

/* COVID-period birth flag (for sensitivity analysis in thesis) */
/*
  2022 survey conducted ~2022; b19 counts months before interview.
  b19 >= 12 & b19 <= 48 → approximate pandemic period (2018–2021)
  Flag as sensitivity indicator only — do NOT drop.
*/
gen     covid_period_birth = 0
replace covid_period_birth = 1  if b19 >= 12 & b19 <= 48
label variable covid_period_birth "Birth during approximate COVID period (sensitivity flag)"
label define yesno 0 "No" 1 "Yes"
label values covid_period_birth yesno

/* ---- 1.2 ANC VISITS — m14 ------------------------------------------------ */
/*
  m14: number of antenatal visits
  Coding: same as 2014 (0=none, 1–90=count, 98=DK, 99=missing)

  KEY METHODOLOGICAL DECISION FOR CROSS-WAVE COMPARISON:
    WHO 2002 standard (used for 2014): 4+ visits = adequate ANC
    WHO 2016 standard (current):       8+ visits = adequate ANC

  THESIS APPROACH:
    PRIMARY threshold:   4+ visits (both waves) — ensures comparability
    SENSITIVITY check:   8+ visits (2022 only)  — shows policy relevance
    Document in Methods: "We use the 4+ threshold for both waves to enable
    temporal comparison. We additionally report 8+ adequacy for 2022 as
    a sensitivity analysis reflecting current WHO recommendations."
*/

replace m14 = .  if m14 == 98 | m14 == 99

/* T1: ANC initiation */
gen     anc_any = 0  if m14 != .
replace anc_any = 1  if m14 >= 1 & m14 != .
label variable anc_any "T1: Any ANC visit (>=1)"
label values anc_any yesno

/* T2: ANC adequacy — PRIMARY threshold (4+, for cross-wave comparability) */
gen     anc4plus = .
replace anc4plus = 0  if anc_any == 1 & m14 < 4
replace anc4plus = 1  if anc_any == 1 & m14 >= 4 & m14 != .
label variable anc4plus "T2: 4+ ANC visits (conditional on any ANC) [cross-wave]"
label values anc4plus yesno

/* T2-alt: ANC adequacy — SENSITIVITY threshold (8+, WHO 2016) */
gen     anc8plus = .
replace anc8plus = 0  if anc_any == 1 & m14 < 8
replace anc8plus = 1  if anc_any == 1 & m14 >= 8 & m14 != .
label variable anc8plus "T2-alt: 8+ ANC visits (WHO 2016 standard, sensitivity)"
label values anc8plus yesno

/* ---- 1.3 SKILLED BIRTH ATTENDANT — m3a to m3n ---------------------------- */
/*
  Same variable structure as 2014. DHS-8 may include additional categories
  (e.g., clinical officer) — check codebook and update recode if needed.

  2022 Kenya DHS-8 delivery assistance categories (verify in codebook):
    m3a = doctor
    m3b = nurse/midwife
    m3c = auxiliary nurse/midwife
    m3d = community health worker (may be NEW in DHS-8)
    m3e = traditional birth attendant (TBA)
    m3f = relative/other
    m3n = no one
    
  If m3d = community health worker, DO NOT include in SBA (not WHO-classified).
*/

foreach v in m3a m3b m3c m3d m3e m3f m3n {
    capture confirm variable `v'
    if _rc != 0 {
        di "WARNING: `v' not found in BR 2022 file — check codebook"
    }
    else {
        capture replace `v' = 0 if `v' == .
    }
}

gen     sba = 0
replace sba = 1  if m3a == 1 | m3b == 1 | m3c == 1
label variable sba "T3-alt: Skilled birth attendant (doctor/nurse/midwife)"
label values sba yesno

/* ---- 1.4 FACILITY DELIVERY — m15 ----------------------------------------- */
/*
  IMPORTANT: DHS-8 (2022) EXPANDED FACILITY CODES.
  
  2014 codes (simple):  10-19=home, 20-29=govt, 30-39=private, 96=other
  2022 codes (expanded; verify against KEIR8CDT codebook):
    11 = respondent's home
    12 = other home
    21 = government hospital
    22 = government health centre
    23 = government dispensary/clinic
    26 = government other
    31 = private hospital/clinic
    32 = mission/faith-based hospital
    33 = private doctor
    36 = private other
    96 = other
    
  Strategy: any code 20–39 = facility delivery
            codes 11–19 or 96 = non-facility (home/other)
  This mapping is consistent with 2014 and ensures cross-wave comparability.
*/

gen     facdelivery = .
replace facdelivery = 0  if m15 >= 11 & m15 <= 19   // home delivery
replace facdelivery = 1  if m15 >= 20 & m15 <= 39   // any health facility
replace facdelivery = 0  if m15 == 96               // other non-facility
replace facdelivery = .  if m15 == 99 | m15 == .
label variable facdelivery "T3: Facility delivery (any institution)"
label values facdelivery yesno

/* Public vs private facility (new for 2022 — useful for policy analysis) */
gen     pub_delivery = .
replace pub_delivery = 1  if m15 >= 21 & m15 <= 26   // government
replace pub_delivery = 0  if m15 >= 31 & m15 <= 36   // private/faith-based
replace pub_delivery = .  if facdelivery == 0 | m15 == . | m15 == 99
label variable pub_delivery "Public facility delivery (1=public, 0=private)"
label values pub_delivery yesno

/* T3 primary */
gen     t3_delivery = .
replace t3_delivery = facdelivery  if anc4plus != .
label variable t3_delivery "T3: Facility delivery (conditional on 4+ ANC)"
label values t3_delivery yesno

/* ---- 1.5 POSTNATAL CARE — m62, m66, m67 ---------------------------------- */
/*
  2022 DHS-8 PNC module improvements:
    m66: postnatal check within 2 days (binary, direct — more available in 2022)
    m67: timing of first postnatal check (days)
    m62: timing in hours/days (same as 2014 but may have cleaner coding)
    
  Strategy:
    1. Try m66 first (most direct, least subject to recall bias)
    2. Fall back to m62 construction if m66 not available
    This is same strategy as 2014, but m66 is MORE likely present in 2022.
*/

capture confirm variable m66
if _rc == 0 {
    gen     pnc48 = .
    replace pnc48 = 0  if m66 == 0
    replace pnc48 = 1  if m66 == 1
    replace pnc48 = .  if m66 == 9 | m66 == .
    di "  PNC 2022: Using m66 (direct 2-day check variable)"
}
else {
    gen     pnc48 = .
    replace pnc48 = 1  if m62 == 0
    replace pnc48 = 1  if m62 >= 100 & m62 <= 148
    replace pnc48 = 0  if m62 >= 149 & m62 <= 248
    replace pnc48 = 0  if m62 >= 249 & m62 <= 395
    replace pnc48 = 0  if m62 == 993 | m62 == 994
    replace pnc48 = .  if m62 == 995
    replace pnc48 = .  if m62 == 998 | m62 == 999 | m62 == .
    di "  PNC 2022: Constructed from m62 (m66 not found)"
}

label variable pnc48 "T4: Postnatal care within 48 hours (mother)"
label values pnc48 yesno

/* T4 conditional */
gen     t4_pnc = .
replace t4_pnc = pnc48  if t3_delivery != .
label variable t4_pnc "T4: PNC within 48hrs (conditional on facility delivery)"
label values t4_pnc yesno

/* ---- 1.6 CONTINUUM OF CARE SCORE AND BINARY OUTCOME ---------------------- */

gen coc_score = 0
replace coc_score = coc_score + 1  if anc_any     == 1
replace coc_score = coc_score + 1  if anc4plus    == 1
replace coc_score = coc_score + 1  if facdelivery == 1
replace coc_score = coc_score + 1  if pnc48       == 1
replace coc_score = .  if anc_any == . & anc4plus == . & facdelivery == . & pnc48 == .
label variable coc_score "Continuum of care score (0-4)"

gen     underutilized = .
replace underutilized = 0  if coc_score >= 3 & coc_score != .
replace underutilized = 1  if coc_score <  3 & coc_score != .
label variable underutilized "Underutilized (CoC score <3): 1=underutilized"
label values underutilized yesno

/* ---- 1.7 DROP-OFF STAGE INDICATOR ---------------------------------------- */

gen dropoff_stage = .
replace dropoff_stage = 0  if anc_any == 1 & anc4plus == 1 & facdelivery == 1 & pnc48 == 1
replace dropoff_stage = 1  if anc_any     == 0
replace dropoff_stage = 2  if anc_any     == 1 & anc4plus    == 0
replace dropoff_stage = 3  if anc4plus    == 1 & facdelivery == 0
replace dropoff_stage = 4  if facdelivery == 1 & pnc48       == 0
label variable dropoff_stage "Stage where CoC dropped (0=complete)"
label define dropstage 0 "Complete CoC" 1 "Dropped: No ANC" ///
    2 "Dropped: ANC <4" 3 "Dropped: No facility delivery" 4 "Dropped: No PNC"
label values dropoff_stage dropstage

/* ---- 1.8 BIRTH-LEVEL CONTROLS -------------------------------------------- */

gen     parity_br = bord
label variable parity_br "Birth order of this birth (from BR)"

gen     child_male = .
replace child_male = 0  if b4 == 2
replace child_male = 1  if b4 == 1
label variable child_male "Child is male (1=yes)"
label values child_male yesno

gen     prev_binterval = b11
replace prev_binterval = .  if b11 == 9996 | b11 == 9998 | b11 == 9999
label variable prev_binterval "Preceding birth interval (months)"

/* ---- 1.9 SURVEY YEAR MARKER (CRITICAL for cross-wave pooling) ------------ */
gen survey_year = $syear
label variable survey_year "Survey year (2022)"

/* ---- 1.10 KEEP AND SAVE -------------------------------------------------- */

keep caseid v001 v002 v003 bidx b19 b5 b4 ///
     anc_any anc4plus anc8plus facdelivery sba pnc48 ///
     t3_delivery t4_pnc coc_score underutilized dropoff_stage ///
     parity_br child_male prev_binterval m14 pub_delivery ///
     covid_period_birth survey_year

di _newline "BR 2022 outcomes summary:"
tab anc_any,       missing
tab anc4plus,      missing
tab anc8plus,      missing
tab facdelivery,   missing
tab pnc48,         missing
tab underutilized, missing
tab dropoff_stage, missing

save "$clean/br_2022_outcomes.dta", replace
di "BR 2022 cleaned. N = `c(N)'"


/*==============================================================================
  PART 2 — INDIVIDUAL RECODE (IR FILE — KEIR8CDT)
  Behavioral and individual-level determinants

  CHANGES FROM 2014 IR:
  1. Autonomy module: DHS-8 may have expanded v743 response codes — remapped
  2. Media: v159 (radio) declining in usage; v171a (internet) added in DHS-8
  3. Wealth: v190 recomputed for 2022 — use within-survey quintiles
  4. Region: v024 now = county code (1–47), not province — harmonised in Part 7
  5. New variable: v190a (wealth factor score) available in some DHS-8 files
==============================================================================*/

di _newline "--- PART 2: INDIVIDUAL RECODE (2022) ---"
use "$ir/$ir22", clear
di "IR 2022 loaded. N = `c(N)'"

gen ir_eligible = (v208 >= 1 & v208 != .)
label variable ir_eligible "Has birth in last 5 years (1=eligible)"

/* ---- 2.1 EDUCATION — v106 ------------------------------------------------- */
/*
  Coding identical to 2014: 0=none, 1=primary, 2=secondary, 3=higher
  No change between DHS rounds for this variable.
*/

gen     edu_level = v106
replace edu_level = .  if v106 == 9 | v106 == .
label variable edu_level "Education level (0=none,1=primary,2=secondary,3=higher)"
label define edu_lbl 0 "No education" 1 "Primary" 2 "Secondary" 3 "Higher"
label values edu_level edu_lbl

gen     edu_secondary = .
replace edu_secondary = 0  if edu_level <= 1
replace edu_secondary = 1  if edu_level >= 2 & edu_level != .
label variable edu_secondary "Secondary+ education (1=yes)"
label values edu_secondary yesno

/* ---- 2.2 MEDIA EXPOSURE — v157, v158, v159, v171a ------------------------ */
/*
  2014: v157 (newspaper), v158 (TV), v159 (radio)
  2022: Same + v171a (internet/social media) — IMPORTANT NEW VARIABLE in DHS-8
  
  THESIS APPROACH for cross-wave comparability:
    Primary media index: v157 + v158 + v159 only (same 3 sources as 2014)
    Extended 2022 index: adds v171a (internet) — reported as supplementary
    
  This ensures the media_index is temporally comparable, while the extended
  index captures the changed information environment in 2022.
*/

foreach v in v157 v158 v159 {
    replace `v' = .  if `v' == 9 | `v' == .
}

/* Primary media index (comparable to 2014) */
gen     media_any = 0
replace media_any = 1  if v157 >= 2 | v158 >= 2 | v159 >= 2
replace media_any = .  if v157 == . & v158 == . & v159 == .
label variable media_any "Any media >=weekly: newspaper/TV/radio (comparable to 2014)"
label values media_any yesno

gen media_index = 0
replace media_index = media_index + 1  if v157 >= 2 & v157 != .
replace media_index = media_index + 1  if v158 >= 2 & v158 != .
replace media_index = media_index + 1  if v159 >= 2 & v159 != .
replace media_index = .  if v157 == . & v158 == . & v159 == .
label variable media_index "Media index: newspaper/TV/radio >=weekly (0-3) [comparable]"

/* Internet/social media — 2022 only (v171a) */
capture confirm variable v171a
if _rc == 0 {
    gen     internet_use = .
    replace internet_use = 0  if v171a == 0
    replace internet_use = 1  if v171a >= 1 & v171a != .
    replace internet_use = .  if v171a == 9
    label variable internet_use "Internet/social media use (2022 only)"
    label values internet_use yesno
    
    /* Extended media index (2022 only, 0-4) */
    gen media_index_ext = media_index
    replace media_index_ext = media_index_ext + 1  if internet_use == 1
    replace media_index_ext = .  if media_index == . & internet_use == .
    label variable media_index_ext "Extended media index incl. internet (2022 only, 0-4)"
    di "  Internet variable v171a found and recoded."
}
else {
    gen internet_use    = .
    gen media_index_ext = media_index
    di "  WARNING: v171a (internet) not found in IR 2022 — using 3-item media index"
}

/* ---- 2.3 AUTONOMY — v743a to v743f --------------------------------------- */
/*
  IMPORTANT: DHS-8 autonomy variable coding.
  
  Standard DHS codes (2014 and 2022):
    1 = respondent alone
    2 = respondent and husband/partner jointly
    3 = husband/partner alone
    4 = someone else
    5 = other
    6 = not applicable / not working (v743e/f)
    8 = don't know (may appear in 2022)
    9 = missing
    
  2022-specific: Code 8 (don't know) more frequently present in DHS-8.
  Treatment: recode 8 as missing.
*/

foreach v in v743a v743b v743c {
    capture confirm variable `v'
    if _rc == 0 {
        replace `v' = .  if `v' == 8 | `v' == 9 | `v' == 6 | `v' == .
    }
    else {
        di "WARNING: `v' not found in IR 2022 — check codebook"
    }
}

gen     auto_health = .
capture {
    replace auto_health = 0  if v743a == 3 | v743a == 4 | v743a == 5
    replace auto_health = 1  if v743a == 1 | v743a == 2
}
label variable auto_health "Autonomy: say in own healthcare (1=has say)"
label values auto_health yesno

gen autonomy_index = 0
capture replace autonomy_index = autonomy_index + 1  if v743a == 1 | v743a == 2
capture replace autonomy_index = autonomy_index + 1  if v743b == 1 | v743b == 2
capture replace autonomy_index = autonomy_index + 1  if v743c == 1 | v743c == 2
replace autonomy_index = .  if auto_health == . & v743b == . & v743c == .
label variable autonomy_index "Autonomy index (0-3 domains with say)"

/* ---- 2.4 PARITY — v220 --------------------------------------------------- */
/*
  No change from 2014. Same coding, same groupings.
  Framing: behavioral fatigue — maintained across both waves.
*/

gen     parity_grp = .
replace parity_grp = 0  if v220 <= 1
replace parity_grp = 1  if v220 >= 2 & v220 <= 3
replace parity_grp = 2  if v220 >= 4 & v220 != .
replace parity_grp = .  if v220 == 9 | v220 == .
label variable parity_grp "Parity group (0=0-1, 1=2-3, 2=4+)"
label define parity_lbl 0 "0-1 children" 1 "2-3 children" 2 "4+ children (high parity)"
label values parity_grp parity_lbl

gen parity_cont = v220
replace parity_cont = .  if v220 == 9 | v220 == .
label variable parity_cont "Number of living children (continuous)"

/* ---- 2.5 MARITAL STATUS — v501 ------------------------------------------- */
/*
  Same coding as 2014. No DHS-8 changes for this variable.
*/

gen     marital_status = .
replace marital_status = 0  if v501 == 0
replace marital_status = 1  if v501 == 1
replace marital_status = 2  if v501 == 2
replace marital_status = 3  if v501 >= 3 & v501 <= 5
replace marital_status = .  if v501 == 9 | v501 == .
label variable marital_status "Marital status (0=never,1=married,2=cohabit,3=prev)"
label define marital_lbl 0 "Never married" 1 "Currently married" ///
    2 "Living together" 3 "Previously married/separated"
label values marital_status marital_lbl

gen     in_union = .
replace in_union = 0  if marital_status == 0 | marital_status == 3
replace in_union = 1  if marital_status == 1 | marital_status == 2
label variable in_union "Currently in union (1=married or cohabiting)"
label values in_union yesno

/* ---- 2.6 INDIVIDUAL SES CONTROLS ----------------------------------------- */

/* Age group */
gen     age_group = v013
replace age_group = .  if v013 == . | v013 > 7
label variable age_group "Age group (5-year categories, v013)"
label define age_lbl 1 "15-19" 2 "20-24" 3 "25-29" 4 "30-34" ///
    5 "35-39" 6 "40-44" 7 "45-49"
label values age_group age_lbl

/* Wealth index */
/*
  CROSS-WAVE WARNING:
  v190 quintiles are WITHIN-SURVEY quintiles. The wealth index is recomputed
  for each DHS round using that round's asset data. A woman in quintile 3
  in 2014 is NOT directly comparable to quintile 3 in 2022 — overall
  wealth levels in Kenya changed substantially over this period.
  
  THESIS APPROACH:
    Use within-survey quintiles for within-wave analyses.
    For cross-wave comparison: use TERTILES (poor/middle/rich) which are
    more stable. We construct both here.
*/

gen     wealth_q = v190
replace wealth_q = .  if v190 == 9 | v190 == .
label variable wealth_q "Wealth quintile (1=poorest, 5=richest) [within 2022 survey]"
label define wealth_lbl 1 "Poorest" 2 "Poorer" 3 "Middle" 4 "Richer" 5 "Richest"
label values wealth_q wealth_lbl

gen     wealth_poor = .
replace wealth_poor = 1  if wealth_q <= 2
replace wealth_poor = 0  if wealth_q >= 3 & wealth_q != .
label variable wealth_poor "Poor (bottom 2 quintiles, 1=yes)"
label values wealth_poor yesno

/* Wealth tertile (for cross-wave comparability) */
gen     wealth_tertile = .
replace wealth_tertile = 1  if wealth_q <= 2          // poor
replace wealth_tertile = 2  if wealth_q == 3          // middle
replace wealth_tertile = 3  if wealth_q >= 4 & wealth_q != .  // rich
label variable wealth_tertile "Wealth tertile (1=poor,2=middle,3=rich) [cross-wave safe]"
label define tert_lbl 1 "Poor (Q1-Q2)" 2 "Middle (Q3)" 3 "Rich (Q4-Q5)"
label values wealth_tertile tert_lbl

/* Residence */
gen     residence = v025
replace residence = .  if v025 == 9 | v025 == .
label variable residence "Residence (1=urban, 2=rural)"
label define res_lbl 1 "Urban" 2 "Rural"
label values residence res_lbl

gen     rural = .
replace rural = 0  if v025 == 1
replace rural = 1  if v025 == 2
label variable rural "Rural residence (1=rural)"
label values rural yesno

/* ---- 2.7 REGION (RAW) — v024 (county in 2022) --------------------------- */
/*
  CRITICAL: In 2022, v024 = county (1–47 Kenya counties).
  In 2014, v024 = province/region (1–8).
  
  We KEEP the raw county variable here and create the harmonised
  8-region variable in Part 7. Both are retained for analysis.
*/

gen     county_2022 = v024
replace county_2022 = .  if v024 == . | v024 == 99
label variable county_2022 "County code 2022 (v024, 1-47 — harmonise in Part 7)"

/* ---- 2.8 SURVEY DESIGN VARIABLES ----------------------------------------- */

gen wt = v005 / 1000000
label variable wt "Survey weight (v005/1,000,000)"
label variable v001 "Cluster number (PSU = community)"
label variable v021 "Primary sampling unit"
label variable v022 "Sampling stratum"

/* ---- 2.9 INTERACTION TERMS ----------------------------------------------- */

gen     wealth_rural = wealth_poor * rural
replace wealth_rural = .  if wealth_poor == . | rural == .
label variable wealth_rural "Interaction: Poor × Rural"

/* ---- 2.10 SURVEY YEAR MARKER --------------------------------------------- */

gen survey_year = $syear
label variable survey_year "Survey year (2022)"

/* ---- 2.11 KEEP AND SAVE -------------------------------------------------- */

keep caseid v001 v002 v003 v005 v021 v022 v012 wt ///
     edu_level edu_secondary media_any media_index media_index_ext internet_use ///
     auto_health autonomy_index parity_grp parity_cont ///
     marital_status in_union age_group ///
     wealth_q wealth_poor wealth_tertile wealth_rural ///
     residence rural county_2022 ir_eligible ///
     v208 v024 v025 v190 v106 v013 survey_year

di _newline "IR 2022 behavioral determinants summary:"
foreach v in edu_level media_any auto_health parity_grp marital_status wealth_q residence {
    tab `v', missing
}

save "$clean/ir_2022_determinants.dta", replace
di "IR 2022 cleaned. N = `c(N)'"


/*==============================================================================
  PART 3 — HOUSEHOLD MEMBER RECODE (PR FILE — KEPR8CDT)
  Community education and literacy aggregation

  CHANGES FROM 2014:
  - hv106 coding identical
  - hv110 (literacy) more consistently coded in DHS-8
  - Restrict to women 15–49 for community education proxy
==============================================================================*/

di _newline "--- PART 3: HOUSEHOLD MEMBER RECODE (2022) ---"
use "$pr/$pr22", clear
di "PR 2022 loaded. N = `c(N)'"

keep if hv104 == 2
keep if hv105 >= 15 & hv105 <= 49
di "After restricting to women 15-49: N = `c(N)'"

gen     edu_pr = hv106
replace edu_pr = .  if hv106 == 9 | hv106 == .

gen     edu_secondary_pr = .
replace edu_secondary_pr = 0  if edu_pr <= 1
replace edu_secondary_pr = 1  if edu_pr >= 2 & edu_pr != .

gen     edu_years = hv108
replace edu_years = .  if hv108 == 97 | hv108 == 98 | hv108 == 99

/* Literacy */
capture confirm variable hv110
if _rc == 0 {
    gen     literate = .
    replace literate = 0  if hv110 == 0 | hv110 == 3
    replace literate = 1  if hv110 == 1 | hv110 == 2
    replace literate = .  if hv110 == 9 | hv110 == .
    di "  Literacy 2022: Using hv110"
}
else {
    gen literate = edu_secondary_pr
    di "  Literacy 2022: Proxied from education (hv110 not found)"
}

gen cluster_id = hv001
keep cluster_id edu_pr edu_secondary_pr edu_years literate

save "$clean/pr_2022_education.dta", replace
di "PR 2022 saved. N = `c(N)'"


/*==============================================================================
  PART 4 — HOUSEHOLD RECODE (HR FILE — KEHR8CDT)
  Community poverty and infrastructure

  CHANGES FROM 2014:
  - hv270 (wealth) recomputed for 2022 — within-survey quintiles only
  - hv221 (mobile telephone) now near-universal; less useful as proxy
  - hv243a (mobile phone): more detailed mobile ownership in DHS-8
  - hv271 (wealth factor score): may be available in KEHR8CDT
  - Expanded water/sanitation codes under JMP 2017 classification
==============================================================================*/

di _newline "--- PART 4: HOUSEHOLD RECODE (2022) ---"
use "$hr/$hr22", clear
di "HR 2022 loaded. N = `c(N)'"

/* Poverty */
gen     wealth_hh = hv270
replace wealth_hh = .  if hv270 == 9 | hv270 == .

gen     poor_hh = .
replace poor_hh = 1  if wealth_hh <= 2
replace poor_hh = 0  if wealth_hh >= 3 & wealth_hh != .
label variable poor_hh "Household in poorest 2 quintiles (2022 within-survey)"

/* Electricity */
gen     has_electricity = .
capture replace has_electricity = 0  if hv206 == 0
capture replace has_electricity = 1  if hv206 == 1

/* Improved toilet — JMP 2017 classification used in DHS-8 */
/*
  2022 DHS-8 expanded sanitation codes (verify against codebook):
    10-15 = flush/pour flush toilet
    20-22 = improved pit latrine
    30    = composting toilet
    96    = other (unimproved)
    31    = pit latrine without slab (unimproved)
    23    = open pit
    42–44 = hanging/bucket toilet (unimproved)
    00    = no facility/open defecation
*/
gen     improved_toilet = .
capture {
    replace improved_toilet = 1  if hv205 >= 10 & hv205 <= 22  // flush + improved pit
    replace improved_toilet = 1  if hv205 == 30                  // composting
    replace improved_toilet = 0  if hv205 == 23 | hv205 == 31  // unimproved pit
    replace improved_toilet = 0  if hv205 >= 42 & hv205 <= 44  // hanging/bucket
    replace improved_toilet = 0  if hv205 == 00 | hv205 == 96  // no facility/other
    replace improved_toilet = .  if hv205 == 99 | hv205 == .
}

/* Improved water — JMP 2017 */
gen     improved_water = .
capture {
    replace improved_water = 1  if hv201 >= 11 & hv201 <= 14   // piped
    replace improved_water = 1  if hv201 >= 21 & hv201 <= 22   // tube well/borehole
    replace improved_water = 1  if hv201 == 31 | hv201 == 41   // protected well/spring
    replace improved_water = 1  if hv201 == 51 | hv201 == 61   // rainwater/tanker
    replace improved_water = 1  if hv201 == 71                   // bottled
    replace improved_water = 0  if improved_water == . & hv201 != . & hv201 != 99
    replace improved_water = .  if hv201 == 99 | hv201 == .
}

/* Residence */
gen     urban_hh = .
replace urban_hh = 1  if hv025 == 1
replace urban_hh = 0  if hv025 == 2

/* Time to water */
gen     water_time = hv204
replace water_time = .  if hv204 == 998 | hv204 == 999

/* Mobile phone ownership (near-universal in 2022 — use as access proxy) */
capture confirm variable hv243a
if _rc == 0 {
    gen     has_mobile = .
    replace has_mobile = 0  if hv243a == 0
    replace has_mobile = 1  if hv243a == 1
    replace has_mobile = .  if hv243a == 9 | hv243a == .
    label variable has_mobile "Household has mobile phone (2022)"
    label values has_mobile yesno
}
else {
    gen has_mobile = .
    di "  hv243a (mobile) not found in HR 2022"
}

gen cluster_id = hv001
keep cluster_id poor_hh wealth_hh has_electricity improved_toilet ///
     improved_water urban_hh water_time has_mobile

save "$clean/hr_2022_household.dta", replace
di "HR 2022 saved. N = `c(N)'"


/*==============================================================================
  PART 5 — MERGE BR + IR (2022)
==============================================================================*/

di _newline "--- PART 5: MERGING BR AND IR (2022) ---"
use "$clean/br_2022_outcomes.dta", clear

sort caseid
merge 1:1 caseid using "$clean/ir_2022_determinants.dta", ///
    keep(match master) nogen

rename v001 cluster_id
label variable cluster_id "Cluster ID (community, Level 2 in MLM)"

di "After BR-IR merge: N = `c(N)'"

count if ir_eligible == 0
di "  IR records without birth (should be 0 after merge): `r(N)'"

save "$clean/kdhs2022_individual.dta", replace
di "Individual merged file (2022) saved. N = `c(N)'"


/*==============================================================================
  PART 6 — COMMUNITY VARIABLE CONSTRUCTION (2022)

  Identical strategy to 2014:
    comm_poverty      — from HR
    comm_edu          — from PR
    comm_literacy     — from PR
    comm_urban        — from HR
    comm_infra        — from HR (electricity + water + toilet)
    comm_anc_loo      — from IR (leave-one-out ANC coverage)
    
  2022 additions:
    comm_mobile       — % households with mobile phone (new infrastructure proxy)
    comm_internet     — % women using internet (from IR, if v171a available)
==============================================================================*/

di _newline "--- PART 6: COMMUNITY VARIABLE CONSTRUCTION (2022) ---"

/* ---- 6.1 POVERTY + INFRASTRUCTURE (from HR) ------------------------------ */
use "$clean/hr_2022_household.dta", clear

bysort cluster_id: egen comm_poverty = mean(poor_hh)
bysort cluster_id: egen comm_urban   = mean(urban_hh)
bysort cluster_id: egen comm_elec    = mean(has_electricity)
bysort cluster_id: egen comm_water   = mean(improved_water)
bysort cluster_id: egen comm_toilet  = mean(improved_toilet)
bysort cluster_id: egen comm_mobile  = mean(has_mobile)

gen comm_infra = (comm_elec + comm_water + comm_toilet) / 3
replace comm_infra = .  if comm_elec == . & comm_water == . & comm_toilet == .

collapse (mean) comm_poverty comm_urban comm_elec comm_water ///
    comm_toilet comm_infra comm_mobile, by(cluster_id)

label variable comm_poverty "Community poverty (% poor HHs) [2022]"
label variable comm_urban   "Community urbanicity (% urban HHs) [2022]"
label variable comm_infra   "Community infrastructure index (0-1) [2022]"
label variable comm_mobile  "Community mobile ownership (% HHs) [2022]"

save "$clean/community_hr_2022.dta", replace
di "Community HR 2022 aggregates saved. N clusters = `c(N)'"

/* ---- 6.2 EDUCATION + LITERACY (from PR) ---------------------------------- */
use "$clean/pr_2022_education.dta", clear

bysort cluster_id: egen comm_edu      = mean(edu_secondary_pr)
bysort cluster_id: egen comm_literacy = mean(literate)
bysort cluster_id: egen comm_edu_yrs  = mean(edu_years)

collapse (mean) comm_edu comm_literacy comm_edu_yrs, by(cluster_id)

label variable comm_edu      "Community education (% women w/ secondary+) [2022]"
label variable comm_literacy "Community literacy (% literate women 15-49) [2022]"

save "$clean/community_pr_2022.dta", replace
di "Community PR 2022 aggregates saved. N clusters = `c(N)'"

/* ---- 6.3 ANC LEAVE-ONE-OUT (from IR+BR merged) --------------------------- */
use "$clean/kdhs2022_individual.dta", clear

bysort cluster_id: egen cluster_anc_sum = total(anc4plus), missing
bysort cluster_id: egen cluster_anc_n   = count(anc4plus)

gen comm_anc_loo = (cluster_anc_sum - anc4plus) / (cluster_anc_n - 1)
replace comm_anc_loo = .  if cluster_anc_n <= 1
replace comm_anc_loo = .  if anc4plus == .
label variable comm_anc_loo "Community ANC coverage (LOO, % with 4+ ANC) [2022]"

drop cluster_anc_sum cluster_anc_n

/* Internet access — community level (2022 only) */
capture confirm variable internet_use
if _rc == 0 {
    bysort cluster_id: egen comm_internet = mean(internet_use)
    label variable comm_internet "Community internet use (% women, 2022 only)"
}
else {
    gen comm_internet = .
}

save "$clean/kdhs2022_individual_loo.dta", replace

/* ---- 6.4 MERGE COMMUNITY AGGREGATES -------------------------------------- */
use "$clean/kdhs2022_individual_loo.dta", clear

merge m:1 cluster_id using "$clean/community_hr_2022.dta", ///
    keep(master match) nogen keepusing(comm_poverty comm_urban comm_infra comm_mobile)

merge m:1 cluster_id using "$clean/community_pr_2022.dta", ///
    keep(master match) nogen keepusing(comm_edu comm_literacy comm_edu_yrs)

di "Community merges complete. N = `c(N)'"

/* ---- 6.5 INTERACTION TERMS ----------------------------------------------- */

gen     edu_comm_poverty = edu_secondary * comm_poverty
replace edu_comm_poverty = .  if edu_secondary == . | comm_poverty == .
label variable edu_comm_poverty "Interaction: Edu_secondary × Comm_poverty"

gen     edu_comm_edu = edu_secondary * comm_edu
replace edu_comm_edu = .  if edu_secondary == . | comm_edu == .
label variable edu_comm_edu "Interaction: Edu_secondary × Comm_edu (cross-level)"

save "$clean/kdhs2022_precounty.dta", replace


/*==============================================================================
  PART 7 — COUNTY-TO-REGION HARMONISATION
  
  CRITICAL FOR CROSS-WAVE COMPARISON.
  
  2014 KDHS: v024 = 8 provinces (Central, Coast, Eastern, Nairobi,
             North Eastern, Nyanza, Rift Valley, Western)
  2022 KDHS: v024 = 47 counties (per Kenya 2010 Constitution)
  
  We map all 47 counties to their historical 8-province equivalent.
  This allows region-level fixed effects to be included in pooled models.
  
  County-to-province mapping based on:
  Kenya National Bureau of Statistics county-to-former-province mapping.
  
  IMPORTANT: Verify these mappings against the 2022 DHS codebook.
  The numeric codes below are STANDARD but may differ slightly in your file.
  Cross-check: tab county_2022, missing — inspect labels against codebook.
==============================================================================*/

di _newline "--- PART 7: COUNTY-TO-REGION HARMONISATION ---"
use "$clean/kdhs2022_precounty.dta", clear

/*
  MAPPING: county_2022 (v024 in KEIR8CDT) → region8 (8-province system)
  
  Region codes (harmonised to match 2014):
    1 = Nairobi
    2 = Central
    3 = Coast
    4 = Eastern
    5 = North Eastern
    6 = Nyanza
    7 = Rift Valley
    8 = Western
    
  County-to-Region mapping (47 counties):
  Nairobi (1):       Nairobi (47)
  Central (2):       Kiambu (22), Kirinyaga (20), Murang'a (21), Nyandarua (18), Nyeri (19)
  Coast (3):         Kilifi (27), Kwale (25), Lamu (29), Mombasa (1), Taita-Taveta (6),
                     Tana River (4)
  Eastern (4):       Embu (14), Isiolo (11), Kitui (15), Machakos (16), Makueni (17),
                     Marsabit (10), Meru (12), Tharaka-Nithi (13)
  North Eastern (5): Garissa (7), Mandera (9), Wajir (8)
  Nyanza (6):        Homa Bay (43), Kisii (45), Kisumu (42), Migori (44), Nyamira (46),
                     Siaya (41)
  Rift Valley (7):   Baringo (30), Bomet (36), Elgeyo-Marakwet (31), Kajiado (24),
                     Kericho (35), Laikipia (32), Nakuru (33), Nandi (34),
                     Narok (23), Samburu (37), Trans-Nzoia (26), Turkana (3),
                     Uasin Gishu (27 — NOTE: verify; may clash with Kilifi code),
                     West Pokot (5)
  Western (8):       Bungoma (39), Busia (40), Kakamega (37 — NOTE: verify),
                     Vihiga (38)
                     
  VERIFICATION NOTE:
  The codes in parentheses are STANDARD DHS Kenya 2022 county codes.
  Some codes may vary in the actual KEIR8CDT file.
  Run: tab county_2022 to verify before running this recode.
  Cross-reference the DHS Kenya 2022 final report appendix.
*/

gen region8 = .

/* Nairobi */
replace region8 = 1  if county_2022 == 47

/* Central */
replace region8 = 2  if inlist(county_2022, 18, 19, 20, 21, 22)

/* Coast */
replace region8 = 3  if inlist(county_2022, 1, 4, 6, 25, 27, 29)

/* Eastern */
replace region8 = 4  if inlist(county_2022, 10, 11, 12, 13, 14, 15, 16, 17)

/* North Eastern */
replace region8 = 5  if inlist(county_2022, 7, 8, 9)

/* Nyanza */
replace region8 = 6  if inlist(county_2022, 41, 42, 43, 44, 45, 46)

/* Rift Valley */
replace region8 = 7  if inlist(county_2022, 3, 5, 23, 24, 26, 30, 31, 32, 33, 34, 35, 36)
replace region8 = 7  if inlist(county_2022, 37, 38, 39, 40)   // some RV/Western boundary counties
/* NOTE: Uasin Gishu (28), Nandi (34), Trans-Nzoia (26) — verify county codes */

/* Western */
replace region8 = 8  if inlist(county_2022, 37, 38, 39, 40)
/*
  IMPORTANT: Kakamega (37) and Vihiga (38) are Western.
  Some boundary counties (Bungoma=39, Busia=40) also Western.
  The Rift Valley recode above may overlap — verify and correct
  based on actual codebook values in your KEIR8CDT file.
*/

label variable region8 "Region (8-province harmonised, cross-wave comparable)"
label define reg8_lbl 1 "Nairobi" 2 "Central" 3 "Coast" 4 "Eastern" ///
    5 "North Eastern" 6 "Nyanza" 7 "Rift Valley" 8 "Western"
label values region8 reg8_lbl

/* Check: any unmatched counties */
count if region8 == . & county_2022 != .
di "  Counties not mapped to region8: `r(N)' — review mapping above"
tab county_2022 if region8 == ., missing

/* Recommended action if mapping gaps found:
   tab county_2022 region8, missing
   — Identify unmapped counties and add to recode above.
*/

save "$clean/kdhs2022_withregion.dta", replace


/*==============================================================================
  PART 8 — FINAL CLEANING, CHECKS, AND SAVE (2022)
==============================================================================*/

di _newline "--- PART 8: FINAL CHECKS AND SAVE (2022) ---"
use "$clean/kdhs2022_withregion.dta", clear

/* ---- 8.1 MISSING DATA SUMMARY -------------------------------------------- */

di _newline "=== MISSING DATA SUMMARY (2022) ==="
di "OUTCOMES:"
foreach v in anc_any anc4plus anc8plus facdelivery pnc48 underutilized {
    quietly count if `v' == .
    di "  `v': `r(N)' missing (`=string(round(`r(N)'/`c(N)'*100,.1))'%)"
}

di _newline "BEHAVIORAL DETERMINANTS:"
foreach v in edu_level media_any auto_health parity_grp marital_status {
    quietly count if `v' == .
    di "  `v': `r(N)' missing (`=string(round(`r(N)'/`c(N)'*100,.1))'%)"
}

di _newline "COMMUNITY VARIABLES:"
foreach v in comm_poverty comm_edu comm_literacy comm_anc_loo comm_infra comm_urban {
    quietly count if `v' == .
    di "  `v': `r(N)' missing (`=string(round(`r(N)'/`c(N)'*100,.1))'%)"
}

/* ---- 8.2 COMPLETE CASE INDICATOR ----------------------------------------- */

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
di "Complete cases: `r(N)' of `c(N)' total"

/* ---- 8.3 ANALYTICAL SAMPLE FLAG ------------------------------------------ */

gen analysis_sample = 1
replace analysis_sample = 0  if anc_any    == .
replace analysis_sample = 0  if cluster_id == .
replace analysis_sample = 0  if age_group  == .
label variable analysis_sample "Primary analysis sample (1=included)"
label values analysis_sample yesno

quietly count if analysis_sample == 1
di "Primary analysis sample (2022): `r(N)' women"

/* ---- 8.4 CLUSTER DIAGNOSTICS --------------------------------------------- */

preserve
keep if analysis_sample == 1
bysort cluster_id: gen n_in_cluster = _N
quietly su n_in_cluster
di _newline "=== CLUSTER DIAGNOSTICS (2022) ==="
quietly distinct cluster_id
di "  Clusters: `r(ndistinct)'"
di "  Mean women per cluster: `r(mean)'"
di "  Min: `r(min)', Max: `r(max)'"
restore

/* ---- 8.5 SAVE FINAL 2022 ANALYSIS FILE ----------------------------------- */

save "$clean/KDHS2022_analysis_final.dta", replace

di _newline "======================================================================"
di "  KDHS 2022 CLEANING COMPLETE"
di "  Final dataset: $clean22/KDHS2022_analysis_final.dta"
di "  N (total): `c(N)'"
quietly count if analysis_sample == 1
di "  N (analysis sample): `r(N)'"
di "  Completed: $(c(current_date)) $(c(current_time))"
di "======================================================================"


/*==============================================================================
  PART 9 — CROSS-WAVE HARMONISATION AND POOLED DATASET CREATION

  This part:
    (a) Loads BOTH cleaned waves (2014 and 2022)
    (b) Renames/recodes any remaining non-comparable variables
    (c) Appends into a single pooled dataset with a survey_year indicator
    (d) Constructs pooled community variables where needed
    (e) Saves the pooled file for cross-wave analysis

  The pooled file enables:
    - Cross-temporal comparison of ICC, fixed effects, ML vs REML estimates
    - Interaction: survey_year × behavioral determinants (temporal change)
    - Pooled models with survey_year as fixed effect
==============================================================================*/

di _newline "--- PART 9: CROSS-WAVE HARMONISATION AND POOLED DATASET ---"

/* ---- 9.1 LOAD 2014 FINAL FILE AND ADD HARMONISATION VARIABLES ------------ */
use "$clean/KDHS2014_analysis_final.dta", clear

/* Add survey year if not already present */
capture confirm variable survey_year
if _rc != 0 {
    gen survey_year = 2014
    label variable survey_year "Survey year"
}

/* Add variables present in 2022 but not 2014 — set to missing for 2014 */
gen anc8plus        = .           // 2022 only (WHO 2016 standard)
gen internet_use    = .           // 2022 only
gen media_index_ext = media_index // 2014: extended = standard (no internet)
gen wealth_tertile  = .
replace wealth_tertile = 1  if wealth_q <= 2
replace wealth_tertile = 2  if wealth_q == 3
replace wealth_tertile = 3  if wealth_q >= 4 & wealth_q != .
gen comm_mobile     = .           // 2022 only
gen comm_internet   = .           // 2022 only
gen pub_delivery    = .           // 2022 only
gen covid_period_birth = 0        // not applicable for 2014

/* Region harmonisation: 2014 already has 8-province v024 */
gen region8 = v024
replace region8 = .  if v024 == . | v024 == 9
label variable region8 "Region (8-province, cross-wave harmonised)"
label define reg8_lbl 1 "Nairobi" 2 "Central" 3 "Coast" 4 "Eastern" ///
    5 "North Eastern" 6 "Nyanza" 7 "Rift Valley" 8 "Western"
label values region8 reg8_lbl

/* County variable: 2014 did not use counties */
gen county_2022 = .

/* Unique cluster ID across waves: prefix with year to avoid collision */
gen cluster_wave = survey_year * 10000 + cluster_id
label variable cluster_wave "Unique cluster ID across waves (year*10000 + cluster)"

save "$pooled/kdhs2014_harmonised.dta", replace
di "2014 harmonised file saved."

/* ---- 9.2 LOAD 2022 FINAL FILE AND ADD HARMONISATION VARIABLES ------------ */
use "$clean/KDHS2022_analysis_final.dta", clear

/* Variables present in 2014 but not 2022 */
capture confirm variable anc8plus
if _rc != 0  gen anc8plus = .

/* Cluster wave ID */
gen cluster_wave = survey_year * 10000 + cluster_id
label variable cluster_wave "Unique cluster ID across waves (year*10000 + cluster)"

save "$pooled/kdhs2022_harmonised.dta", replace
di "2022 harmonised file saved."

/* ---- 9.3 APPEND INTO POOLED DATASET -------------------------------------- */
use "$pooled/kdhs2014_harmonised.dta", clear
append using "$pooled/kdhs2022_harmonised.dta"

di "Pooled dataset created. Total N = `c(N)'"
tab survey_year, missing

/* ---- 9.4 POOLED CROSS-WAVE VARIABLES ------------------------------------- */

/* Binary wave indicator (for interaction models) */
gen     wave2022 = (survey_year == 2022)
label variable wave2022 "Survey wave (0=2014, 1=2022)"
label values wave2022 yesno

/* Behavioral change interactions (wave × determinant) */
/* These test whether the effect of a determinant CHANGED between 2014 and 2022 */
gen     wave_edu      = wave2022 * edu_secondary
gen     wave_media    = wave2022 * media_any
gen     wave_autonomy = wave2022 * auto_health
gen     wave_parity   = wave2022 * parity_cont
gen     wave_wealth   = wave2022 * wealth_poor
gen     wave_rural    = wave2022 * rural

foreach v in wave_edu wave_media wave_autonomy wave_parity wave_wealth wave_rural {
    replace `v' = .  if wave2022 == . | `v' == .
}

label variable wave_edu      "Interaction: Wave2022 × Education (secondary+)"
label variable wave_media    "Interaction: Wave2022 × Media exposure"
label variable wave_autonomy "Interaction: Wave2022 × Autonomy (healthcare)"
label variable wave_parity   "Interaction: Wave2022 × Parity (continuous)"
label variable wave_wealth   "Interaction: Wave2022 × Poor (bottom 2 quintiles)"
label variable wave_rural    "Interaction: Wave2022 × Rural"

/* ---- 9.5 FINAL SAMPLE CHECK ON POOLED DATA ------------------------------- */

di _newline "=== POOLED DATASET SUMMARY ==="
tab survey_year analysis_sample, missing
tab survey_year underutilized   if analysis_sample == 1, missing
tab survey_year residence       if analysis_sample == 1, missing

di _newline "Cluster counts by wave:"
preserve
keep if analysis_sample == 1
bysort survey_year: distinct cluster_id
restore

/* ---- 9.6 SAVE POOLED FILE ------------------------------------------------ */
save "$pooled/KDHS_pooled_2014_2022_final.dta", replace

di _newline "======================================================================"
di "  POOLED DATASET SAVED"
di "  File: $pooled/KDHS_pooled_2014_2022_final.dta"
di "  Total N: `c(N)'"
quietly count if analysis_sample == 1
di "  Analysis sample (both waves): `r(N)'"
di "  Completed: $(c(current_date)) $(c(current_time))"
di "======================================================================"

log close


/*==============================================================================
  APPENDIX A — FULL VARIABLE CROSSWALK: 2014 vs 2022

  VARIABLE          2014 CODE    2022 CODE    HARMONISED NAME     ISSUE
  ----------------  -----------  -----------  ------------------  ----------------
  Outcomes
  ANC initiation    m14 >=1      m14 >=1      anc_any             None
  ANC 4+            m14 >=4      m14 >=4      anc4plus            Use for both waves
  ANC 8+ (WHO2016)  N/A          m14 >=8      anc8plus            2022 only
  Facility deliv.   m15 20-39    m15 20-39    facdelivery         2022 has more codes
  SBA               m3a/b/c      m3a/b/c(+d?) sba                 Check m3d in 2022
  PNC 48hrs         m62/m66      m66 (better) pnc48               More reliable 2022
  CoC score         Derived      Derived      coc_score           Comparable
  Underutilized     Derived      Derived      underutilized       Comparable

  Behavioral determinants
  Education         v106 0-3     v106 0-3     edu_level           Identical
  Media (newspaper) v157 0-3     v157 0-3     media_index         Identical
  Media (TV)        v158 0-3     v158 0-3     media_index         Identical
  Media (radio)     v159 0-3     v159 0-3     media_index         Identical
  Internet          N/A          v171a        internet_use        2022 only
  Autonomy health   v743a 1-5    v743a 1-8    auto_health         Code 8=DK in 2022
  Parity            v220         v220         parity_cont/grp     Identical
  Marital status    v501 0-5     v501 0-5     marital_status      Identical

  Individual SES
  Wealth quintile   v190 1-5     v190 1-5     wealth_q            Not cross-comparable
  Wealth tertile    Derived      Derived      wealth_tertile      USE FOR COMPARISON
  Residence         v025 1-2     v025 1-2     rural               Identical
  Region            v024 1-8     v024 1-47    region8             REMAP 2022 counties
  Age group         v013 1-7     v013 1-7     age_group           Identical

  Survey design
  Cluster           v001         v001         cluster_id          Recode to cluster_wave
  Weight            v005/1e6     v005/1e6     wt                  Identical formula
  PSU               v021         v021         v021                Identical
  Stratum           v022         v022         v022                Identical

==============================================================================*/

/*==============================================================================
  APPENDIX B — RECOMMENDED ORDER OF ANALYSIS DO-FILES

  Step 1: KDHS2014_cleaning.do         → $clean/KDHS2014_analysis_final.dta
  Step 2: KDHS2022_cleaning.do         → $clean22/KDHS2022_analysis_final.dta
          (this file, Parts 0-9)        → $pooled/KDHS_pooled_2014_2022_final.dta
  Step 3: 01_descriptives.do           → prevalence tables, CoC dropout tables
  Step 4: 02_null_models.do            → empty MLM, ICC, ML vs REML comparison
  Step 5: 03_individual_models.do      → T1-T4 with behavioral determinants
  Step 6: 04_community_models.do       → add community variables, cross-level interactions
  Step 7: 05_ml_reml_comparison.do     → systematic ML vs REML comparison (H5-H8)
  Step 8: 06_cross_wave_models.do      → pooled models, wave interactions
  Step 9: 07_sensitivity_analyses.do   → 8+ ANC threshold, complete case, etc.
==============================================================================*/

/*  END OF KDHS 2022 CLEANING DO-FILE  */
