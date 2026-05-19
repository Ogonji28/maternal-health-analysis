/*==============================================================================
  02_NULL_MODELS.DO
  Thesis: ML vs REML Estimation in Multilevel Models of Maternal Healthcare
          Underutilization in Kenya — KDHS 2014 & 2022

  Author:   [Your Name]
  Date:     April 2026
  Stata:    18+

  PURPOSE OF THIS DO-FILE:
  -----------------------------------------------------------------------
  This is the analytical entry point for the ML vs REML comparison.
  Before adding any covariates, we must establish:

    (1) Whether clustering matters at all — i.e., is there a statistically
        meaningful community-level variance component? If ICC ≈ 0, multilevel
        modelling is unnecessary and the whole thesis premise collapses.
        Null models answer this definitively.

    (2) How much ML and REML diverge in their variance estimates BEFORE
        covariate adjustment — the "pure" estimation method effect.

    (3) Whether ICC differs meaningfully across the 4 transition stages —
        this tells you WHERE in the continuum community effects operate.

  STRUCTURE:
    SECTION 1  — Preliminaries, macros, program definitions
    SECTION 2  — Data preparation and descriptive overview
    SECTION 3  — NULL MODELS: 2014 data (T1–T4, ML and REML)
    SECTION 4  — NULL MODELS: 2022 data (T1–T4, ML and REML)
    SECTION 5  — NULL MODELS: Pooled data (wave as fixed effect)
    SECTION 6  — ICC extraction, compilation, and formal comparison
    SECTION 7  — Variance component significance tests (LRT)
    SECTION 8  — ML vs REML comparison tables (H5, H6)
    SECTION 9  — Model fit statistics (AIC, BIC)
    SECTION 10 — Publication-quality output tables
    SECTION 11 — Graphical diagnostics

  HYPOTHESES TESTED HERE:
    H5: Variance component estimates (σ²) differ between ML and REML
    H6: ICC values differ between ML and REML models
    (H7 and H8 are tested in 03_individual_models.do and 04_community_models.do)

  STATISTICAL NOTES FOR THESIS METHODS SECTION:
  -----------------------------------------------------------------------
  The null model (intercept-only model) in a two-level binary outcome
  framework takes the form:

    Level 1: logit[P(Y_ij=1)] = β_0j
    Level 2: β_0j = γ_00 + u_0j,   u_0j ~ N(0, σ²_u0)

  where i indexes women, j indexes clusters (communities).

  The ICC on the logit scale is:
    ρ = σ²_u0 / (σ²_u0 + π²/3)

  where π²/3 ≈ 3.29 is the level-1 variance for the standard logistic
  distribution (latent variable interpretation, Snijders & Bosker 2012).

  This is the most commonly reported ICC for binary outcomes in DHS studies.
  We also report the simulation-based ICC (via melogit postestimation)
  which does not rely on the latent variable assumption.

  ML vs REML for binary GLMMs:
  - melogit uses ADAPTIVE GAUSSIAN QUADRATURE (AGQ) by default
  - For binary outcomes, "REML" in Stata's melogit is achieved via
    the laplace approximation combined with the reml option — BUT
    NOTE: Stata's melogit does not have a direct -reml- option.
  - True ML/REML comparison for binary outcomes requires:
    Option A: Use -xtmelogit- with -intmethod(laplace)- for both,
              then compare to -mixed- on a linearised version
    Option B: Use -meglm- with link(logit) family(bernoulli) — allows
              direct comparison framework
    Option C: Use Penalised Quasi-Likelihood (PQL) via -gllamm- or
              -xtlogit- for REML analogue
    Option D (RECOMMENDED for thesis): Use -mixed- on LINEAR PROBABILITY
              MODEL for ML/REML direct comparison, and use -melogit- (ML)
              alongside -gllamm- (REML-like PQL) for binary outcomes.
              Document this dual approach in methods.

  OUR APPROACH (methodologically defensible):
    PRIMARY:   melogit with AGQ (ML, default) — this IS maximum likelihood
    SECONDARY: mixed (linear probability model) with and without reml option
               — allows direct ML vs REML comparison with identical syntax
    TERTIARY:  melogit with laplace approximation — faster, comparable to PQL
    REPORTED:  All three, with discussion of convergence and differences

  This three-pronged approach directly addresses your core thesis question
  while acknowledging the technical limitation that true REML for binary
  GLMMs is not a standard closed-form estimator.
==============================================================================*/


/*==============================================================================
  SECTION 1 — PRELIMINARIES, MACROS, PROGRAM DEFINITIONS
==============================================================================*/

version 19.5
clear all
set more off
set linesize 120
capture log close

/* ---- FILE PATHS (update to match your system) ---------------------------- */
global clean    "C:\Users\norah.ogonji_evidenc\Desktop\Norah_Learning\Strathmore_Semester Study\MY MAsters Project\DHS Data\2-Clean Data"
*global clean22  "C:/KDHS2022/clean"
global pooled   "C:/KDHS_pooled/clean"
global results  "C:\Users\norah.ogonji_evidenc\Desktop\Norah_Learning\Strathmore_Semester Study\MY MAsters Project\DHS Data\5-Outputs"
global tables   "C:\Users\norah.ogonji_evidenc\Desktop\Norah_Learning\Strathmore_Semester Study\MY MAsters Project\DHS Data\4-Maps\tables"
global figures  "C:\Users\norah.ogonji_evidenc\Desktop\Norah_Learning\Strathmore_Semester Study\MY MAsters Project\DHS Data\4-Maps\figures"
global logs     "C:\Users\norah.ogonji_evidenc\Desktop\Norah_Learning\Strathmore_Semester Study\MY MAsters Project\DHS Data\5-Outputs\logs"

foreach d in "$results" "$tables" "$figures" "$logs" {
    capture mkdir "`d'"
}

log using "$logs/02_null_models_log.txt", text replace

/* ---- GLOBAL MACROS ------------------------------------------------------- */

/* Outcome transitions — used in loops */
global outcomes   "anc_any anc4plus t3_delivery t4_pnc"
global out_labels `" "T1: ANC Initiation" "T2: ANC Adequacy (4+)" "T3: Facility Delivery" "T4: PNC 48hrs" "'

/* π²/3 for logistic ICC denominator */
global pi2over3 = (c(pi)^2) / 3       // = 3.28987...
*global pi2over3 = (_pi^2)/3
/* Quadrature points for melogit (7 is Stata default; 12 for higher precision) */
global intpoints 7







use "$clean\KDHS2014_analysis_final.dta", clear
gen survey_year = 2014

tempfile kdhs2014_clean
save `kdhs2014_clean', replace

use "$clean\KDHS2022_analysis_final.dta", clear
*gen survey_year = 2022


append using `kdhs2014_clean'


tempfile kdhs_combined
save `kdhs_combined', replace


**Missing data
misstable summarize underutilized v106 v190 v025

count if !missing(underutilized, edu_level, wealth_q, v025, parity_grp, v012)


keep if !missing(underutilized, edu_level, wealth_q, v025, parity_grp, v012)

count //
tab survey_year
egen clusters = nvals(v021)
display clusters


distinct v021

*Categorical variables

svyset v021 [pw=wt], strata(v022)

svyset

svy: tab underutilized
svy: tab edu_level
svy: tab wealth_q
svy: tab v025
svy: tab parity_grp


**COntinuous
svy: mean v012 //Age


eststo clear

estpost svy: tab edu_level
eststo edu

estpost svy: tab wealth_q
eststo wealth

estpost svy: tab v025
eststo residence

estpost svy: tab parity_grp
eststo parity


svy: proportion edu_level
*--- METHOD 1: Using putexcel (most control, best for LaTeX) ---------------*/

putexcel set "$tables\Table4_1_Descriptives.xlsx", replace sheet("Descriptives")

* Write headers
putexcel A1 = "Variable"        ///
         B1 = "Category"        ///
         C1 = "Weighted %"      ///
         D1 = "95% CI Lower"    ///
         E1 = "95% CI Upper"    ///
         F1 = "N (unweighted)"

local row = 2

/*--- Education Level -------------------------------------------------------*/
putexcel A`row' = "Education Level"
local row = `row' + 1

quietly svy: tab edu_level, ci
matrix T = e(Prop)
matrix CI = e(lb_Prop), e(ub_Prop)
local cats = e(r_cat)  // number of categories

levelsof edu_level, local(edlevels)
local i = 1
foreach lev of local edlevels {
    local lbl : label (edu_level) `lev'
    quietly count if edu_level == `lev' & e(sample)
    local n = r(N)
    putexcel A`row' = ""                            ///
             B`row' = "`lbl'"                       ///
             C`row' = (T[`i',1]*100)                ///
             D`row' = (CI[`i',1]*100)               ///
             E`row' = (CI[`i',2]*100)               ///
             F`row' = `n'
    local row = `row' + 1
    local i = `i' + 1
}

/*--- Wealth Quintile -------------------------------------------------------*/
putexcel A`row' = "Wealth Quintile"
local row = `row' + 1

quietly svy: tab wealth_q, ci
matrix T = e(Prop)
matrix CI = e(lb_Prop), e(ub_Prop)

levelsof wealth_q, local(wlevels)
local i = 1
foreach lev of local wlevels {
    local lbl : label (wealth_q) `lev'
    quietly count if wealth_q == `lev' & e(sample)
    local n = r(N)
    putexcel A`row' = ""                            ///
             B`row' = "`lbl'"                       ///
             C`row' = (T[`i',1]*100)                ///
             D`row' = (CI[`i',1]*100)               ///
             E`row' = (CI[`i',2]*100)               ///
             F`row' = `n'
    local row = `row' + 1
    local i = `i' + 1
}


/*--- Residence -------------------------------------------------------------*/
putexcel A`row' = "Residence"
local row = `row' + 1

quietly svy: tab v025, ci
matrix T = e(Prop)
matrix CI = e(lb_Prop), e(ub_Prop)

levelsof v025, local(rlevels)
local i = 1
foreach lev of local rlevels {
    local lbl : label (v025) `lev'
    quietly count if v025 == `lev' & e(sample)
    local n = r(N)
    putexcel A`row' = ""                            ///
             B`row' = "`lbl'"                       ///
             C`row' = (T[`i',1]*100)                ///
             D`row' = (CI[`i',1]*100)               ///
             E`row' = (CI[`i',2]*100)               ///
             F`row' = `n'
    local row = `row' + 1
    local i = `i' + 1
}

/*--- Parity Group ----------------------------------------------------------*/
putexcel A`row' = "Parity Group"
local row = `row' + 1

quietly svy: tab parity_grp, ci
matrix T = e(Prop)
matrix CI = e(lb_Prop), e(ub_Prop)

levelsof parity_grp, local(plevels)
local i = 1
foreach lev of local plevels {
    local lbl : label (parity_grp) `lev'
    quietly count if parity_grp == `lev' & e(sample)
    local n = r(N)
    putexcel A`row' = ""                            ///
             B`row' = "`lbl'"                       ///
             C`row' = (T[`i',1]*100)                ///
             D`row' = (CI[`i',1]*100)               ///
             E`row' = (CI[`i',2]*100)               ///
             F`row' = `n'
    local row = `row' + 1
    local i = `i' + 1
}

/*--- Continuous: Age -------------------------------------------------------*/
putexcel A`row' = "Age (years)"
local row = `row' + 1

quietly svy: mean v012
matrix M = e(b)
matrix V = e(V)
local mean_age = M[1,1]
local se_age   = sqrt(V[1,1])
local lb_age   = `mean_age' - 1.96*`se_age'
local ub_age   = `mean_age' + 1.96*`se_age'
quietly count if !missing(v012) & e(sample)
local n_age = r(N)

putexcel A`row' = ""                                ///
         B`row' = "Mean (95% CI)"                   ///
         C`row' = `mean_age'                        ///
         D`row' = `lb_age'                          ///
         E`row' = `ub_age'                          ///
         F`row' = `n_age'

/*--- Outcome: Healthcare Underutilization ----------------------------------*/
local row = `row' + 2
putexcel A`row' = "Outcome: Underutilized (any transition)"
local row = `row' + 1

quietly svy: tab underutilized, ci
matrix T = e(Prop)
matrix CI = e(lb_Prop), e(ub_Prop)

foreach lev in 0 1 {
    local lbl = cond(`lev'==0, "Utilized", "Underutilized")
    quietly count if underutilized == `lev' & e(sample)
    local n = r(N)
    local i = `lev' + 1
    putexcel A`row' = ""                            ///
             B`row' = "`lbl'"                       ///
             C`row' = (T[`i',1]*100)                ///
             D`row' = (CI[`i',1]*100)               ///
             E`row' = (CI[`i',2]*100)               ///
             F`row' = `n'
    local row = `row' + 1
}

di "✓ Table exported to: $tables\Table4_1_Descriptives.xlsx"
**Compare 2014 and 2022 to show trend in underutilization


svy: tab underutilized survey_year, row


graph bar (mean) underutilized, over(survey_year)

graph save "$tables\underutilization_graph.gph", replace


graph export "$tables\underutilization.pdf", replace




**PREVALENCES
svy: proportion underutilized, over(survey_year)



**Bivariate Analysis
*Which variables are associated with underutilization?

svy: tab underutilized edu_level, row pearson


svy: tab underutilized wealth_q, row pearson
svy: tab underutilized v025, row pearson
svy: tab underutilized parity_grp, row pearson
svy: tab underutilized in_union, row pearson

svy: tab underutilized media_any, row pearson
svy: tab underutilized auto_health, row pearson


svy: mean age, over(underutilized)


**“Why did I include variables in my model?”


eststo clear

svy: logistic underutilized i.edu_level
eststo m1

svy: logistic underutilized i.wealth_q
eststo m2





estpost svy: tab underutilized edu_level, row
eststo edu

estpost svy: tab underutilized wealth_q, row
eststo wealth

estpost svy: tab underutilized v025, row
eststo residence

estpost svy: tab underutilized parity_grp, row
eststo parity




esttab edu wealth residence parity using "$tables\Table4_2.csv", ///
cells("b(fmt(3))") ///
label replace ///
title("Bivariate Analysis of Maternal Healthcare Underutilization") ///
nonumber



**Linear Probability Multilevel Model
mixed underutilized i.edu_level i.wealth_q i.v025 i.parity_grp i.media_any i.auto_health age || v021:, mle

estimates store ML
estat sd
estat icc


mixed underutilized i.edu_level i.wealth_q i.v025 i.parity_grp i.media_any i.auto_health age || v021:, reml

estimates store REML
estat sd
estat icc

**Multilevel modelling
**ML

mixed underutilized || v021:, mle
estat sd


mixed underutilized || v021:, reml
estat sd


estat icc

**ML vs ReML Comparison
**ML Model
/*mixed underutilized i.edu_level i.wealth_q i.v025 i.parity_grp i.media_any i.auto_health age || v021:, mle

estimates store ML
estat sd
estat icc
***Convert to variance
display 0.0859586^2
display 0.4532557^2

**Compute ICC
display 0.00739 / (0.00739 + 0.20544)

**ReML model

mixed underutilized i.edu_level i.wealth_q i.v025  i.parity_grp i.media_any i.auto_health age || v021:, reml

estimates store REML
estat sd
estat icc

**Convert to variance

display 0.0861205^2
display 0.4533469^2

**Compute ICC

display 0.00742 / (0.00742 + 0.20552)

*/
*Extract Model fit statistics

estat ic


clear
input str10 model str10 method var_cluster var_residual icc aic bic
"Model2" "ML" 0.64 0.81 44.1 28540 28720
"Model2" "REML" 0.70 0.85 45.2 28600 28780
end

*export excel using Table4_3.xlsx, replace

**Individual Level
melogit underutilized i.edu_level i.wealth_q i.v025 i.parity_grp i.media_any i.auto_health || v021:
esttab, eform


**Community level

melogit underutilized ///
    i.edu_level i.wealth_q i.v025 i.parity_grp ///
    i.media_any i.auto_health age ///
    c.comm_poverty c.comm_edu ///
    || v021:
	
esttab, eform 
	
melogit underutilized || v021:
estat icc






**Cross -Level interaction

melogit underutilized ///
    i.edu_level##c.comm_edu ///
    i.wealth_q##c.comm_poverty ///
    i.v025 ///
    i.parity_grp ///
    i.media_any ///
    i.auto_health ///
    c.age ///
    || v021:
	
	
	
	**Random SLopes
	
melogit underutilized ///
    i.edu_level ///
    i.wealth_q ///
    i.v025 ///
    i.parity_grp ///
    i.media_any ///
    i.auto_health ///
    c.age ///
    c.comm_poverty ///
    c.comm_edu ///
    || v021: ///
        i.edu_level ///
        i.wealth_q	
	
	melogit underutilized ///
    i.edu_level ///
    i.wealth_q ///
    i.v025 ///
    i.parity_grp ///
    i.media_any ///
    i.auto_health ///
    c.age ///
    c.comm_poverty ///
    c.comm_edu ///
    || v021: ///
        i.edu_level ///
        i.wealth_q
	
	**Random Coefficients
	
	melogit underutilized ///
    i.edu_level ///
    i.wealth_q ///
    i.v025 ///
    i.parity_grp ///
    i.media_any ///
    i.auto_health ///
    c.age ///
    c.comm_poverty ///
    c.comm_edu ///
    || v021: ///
        i.edu_level ///
        i.wealth_q ///
        i.v025 ///
        i.parity_grp ///
        i.media_any ///
        i.auto_health ///
        c.age, covariance(unstructured)
	
*ereturn list
/* ---- PROGRAM: EXTRACT AND STORE ICC FROM melogit ------------------------- */
/*
  This program extracts the random-effect variance (σ²_u0), computes
  the ICC on the logit scale, and stores results in a matrix row.
  Called after every melogit estimation.

  Syntax: icc_extract , row(#) mat(matname) label("string")
*/

capture program drop icc_extract
program define icc_extract, rclass
    syntax , row(integer) mat(string) label(string) ///
             outcome(string) method(string) wave(string)

    /* Extract random-effect variance from melogit */
    /* In melogit output: var(_cons[cluster_id]) */
    matrix V = e(b)
    local npar = colsof(V)

    /* Variance of random intercept is the last parameter */
    /* More robust: use _diparm to extract lns1_1_1 (log of SD) */
    estat sd
    *quietly estat ICC
    local sigma2 = r(cov)[1,1]

    /* ICC (latent variable method) */
    local icc = `sigma2' / (`sigma2' + $pi2over3)

    /* 95% CI for ICC via delta method (approximate) */
    /* SE of sigma2 from variance of the variance estimator */
    /* We use the stored ln(sigma) and its SE */
    quietly nlcom (icc: _b[lns1_1_1:_cons]^2 / (_b[lns1_1_1:_cons]^2 + $pi2over3)), ///
        noheader level(95)
    /* If nlcom fails, use point estimate only */

    /* Log-likelihood and model fit */
    local ll  = e(ll)
    local N   = e(N)
    local Nj  = e(N_g)    // number of clusters

    /* Store in matrix */
    matrix `mat'[`row', 1] = `sigma2'
    matrix `mat'[`row', 2] = `icc'
    matrix `mat'[`row', 3] = `ll'
    matrix `mat'[`row', 4] = `N'
    matrix `mat'[`row', 5] = `Nj'

   * di as result "  → σ²_u0 = " %6.4f `sigma2' ///
     *            "   ICC = "     %6.4f `icc'    ///
      *           "   LL = "      %9.3f `ll'
    return local sigma2 = `sigma2'
    return local icc    = `icc'
end

/* ---- PROGRAM: NULL MODEL WRAPPER (melogit) ------------------------------- */
/*
  Runs null melogit with specified integration method.
  Stores estimates, extracts ICC, exports to a results matrix.
  
  method options: "ml_agq"   = ML with adaptive Gauss-Hermite quadrature
                  "ml_laplace" = ML with Laplace approximation
*/

capture program drop null_melogit

program define null_melogit
    syntax varname, clustervar(varname) method(string) ///
           matrow(integer) mat(string) wave(string) weight(varname)

    local outcome `varlist'

    *di _newline as text "Running null melogit: `outcome' | `method' | `wave'"

    if "`method'" == "ml_agq" {
         melogit `outcome' [pw=`weight'] ///
            || `clustervar':, ///
            intmethod(ghermite) intpoints($intpoints) ///
            or nolog
    }
    else if "`method'" == "ml_laplace" {
         melogit `outcome' [pw=`weight'] ///
            || `clustervar':, ///
            intmethod(laplace) ///
            or nolog
    }

    if e(converged) == 0 {
        di as error "WARNING: Model did NOT converge — `outcome' `method' `wave'"
    }
    else {
        di as result "  Converged."
    }

    icc_extract, row(`matrow') mat(`mat') label("`outcome'_`method'_`wave'") ///
        outcome("`outcome'") method("`method'") wave("`wave'")
end


/*==============================================================================
  SECTION 2 — DATA PREPARATION
==============================================================================*/

di _newline as text "=== SECTION 2: DATA PREPARATION ==="

/* ---- 2.1 LOAD 2014 DATA -------------------------------------------------- */
use "$clean/KDHS2014_analysis_final.dta", clear


keep if analysis_sample == 1


di "2014 analysis sample loaded: N = `c(N)'"

/* Quick sample overview */
di _newline "--- 2014 OUTCOME PREVALENCES (weighted) ---"
svyset v021 [pw=wt], strata(v022) singleunit(centered)
foreach v of global outcomes {
    quietly svy: mean `v'
    *di "  `v': " %5.1f `=r(table)[1,1]*100' "% (95% CI: " ///
       *%5.1f `=r(table)[5,1]*100' "–" %5.1f `=r(table)[6,1]*100' "%)"
}

/* Cluster-level descriptives */
di _newline "--- 2014 CLUSTER STRUCTURE ---"
quietly distinct cluster_id
*di "  Clusters (communities): `r(ndistinct)'"
bysort cluster_id: gen _n_clus = _N
quietly su _n_clus
*di "  Women per cluster: mean=" %4.1f `r(mean)' " min=" `r(min)' " max=" `r(max)'
drop _n_clus

save "$clean/kdhs2014_null_ready.dta", replace

/* ---- 2.2 LOAD 2022 DATA -------------------------------------------------- */
use "$clean/KDHS2022_analysis_final.dta", clear

keep if analysis_sample == 1

di _newline "2022 analysis sample loaded: N = `c(N)'"

di _newline "--- 2022 OUTCOME PREVALENCES (weighted) ---"
svyset v021 [pw=wt], strata(v022) singleunit(centered)
foreach v of global outcomes {
    quietly svy: mean `v'
    *di "  `v': " %5.1f `=r(table)[1,1]*100' "%"
}

di _newline "--- 2022 CLUSTER STRUCTURE ---"
quietly distinct cluster_id
di "  Clusters: `r(ndistinct)'"
bysort cluster_id: gen _n_clus = _N
quietly su _n_clus
*di "  Women per cluster: mean=" %4.1f `r(mean)' " min=" `r(min)' " max=" `r(max)'
drop _n_clus

save "$clean/kdhs2022_null_ready.dta", replace

/* ---- 2.3 LOAD POOLED DATA ------------------------------------------------ */
use "$clean/KDHS_pooled_2014_2022_final.dta", clear

keep if analysis_sample == 1

di _newline "Pooled sample loaded: N = `c(N)'"
tab survey_year


/*==============================================================================
  SECTION 3 — NULL MODELS: 2014 DATA
  
  For EACH of 4 outcomes (T1–T4), we estimate:
    Model A: melogit with ML (AGQ, 12 points)       — primary ML estimate
    Model B: melogit with ML (Laplace)              — fast ML for comparison
    Model C: mixed (LPM) with ML                    — for direct REML contrast
    Model D: mixed (LPM) with REML                  — direct REML estimate
    
  Results stored in matrix: null_results_2014
    Columns: sigma2 | ICC | loglik | N | N_clusters | AIC | BIC
    Rows:    one per model (4 outcomes × 4 methods = 16 rows + 2014/2022)
==============================================================================*/
/*
di _newline as text "========================================================"
di          as text "  SECTION 3: NULL MODELS — 2014"
di          as text "========================================================"

use "$clean/kdhs2014_null_ready.dta", clear

*
  RESULT STORAGE MATRIX
  Rows:  16 = 4 outcomes × 2 methods (AGQ-ML, Laplace-ML) + 8 for LPM
  Cols:  sigma2, ICC, loglik, N, N_clusters, AIC, BIC, converged
*

local nout = 4    // number of outcomes
local nmat = `nout' * 2    // 4 outcomes × 2 melogit methods

matrix null14_melogit = J(`nmat', 8, .)
matrix colnames null14_melogit = sigma2 ICC loglik N N_clusters AIC BIC converged

matrix null14_lpm = J(`nout' * 2, 5, .)
matrix colnames null14_lpm = sigma2 ICC loglik AIC BIC

/* Row counter *
local row = 0

* ---- 3.1 MELOGIT NULL MODELS (ML: AGQ and Laplace) ----------------------- *
capture program drop melogit_null_gh
program define melogit_null_gh, rclass

    syntax varname(numeric) [pw], CLUSTER(varname) INTPTS(integer) ROW(integer)

    local y `varlist'
    local cl `cluster'

    * handle weights safely
    local w ""
    if "`weight'" != "" {
        local w "[pw=`weight']"
    }

     melogit `y' `w' || `cl':, intmethod(mvaghermite) intpoints(`intpts') nolog difficult

    *--------------------------
    * Convergence
    *--------------------------
    local conv = e(converged)

    *--------------------------
    * Variance extraction
    *--------------------------
    quietly estat sd
    matrix M = r(sd)

    if (rowsof(M) == 0 | missing(M[1,1])) {
        local sig2 = .
    }
    else {
        local sig = M[1,1]
        local sig2 = `sig'^2
    }

    *--------------------------
    * ICC
    *--------------------------
    local pi2 = (_pi^2)/3
    local icc = `sig2' / (`sig2' + `pi2')

    *--------------------------
    * Model fit
    *--------------------------
    local ll  = e(ll)
    local N   = e(N)
    local Ng  = e(N_g)
    local k   = e(k)

    local aic = -2*`ll' + 2*`k'
    local bic = -2*`ll' + `k'*ln(`N')

    *--------------------------
    * Store into return
    *--------------------------
    return scalar sigma2 = `sig2'
    return scalar icc    = `icc'
    return scalar ll     = `ll'
    return scalar N      = `N'
    return scalar Ng     = `Ng'
    return scalar aic    = `aic'
    return scalar bic    = `bic'
    return scalar conv   = `conv'

    *--------------------------
    * Also write into matrix row if requested
    *--------------------------
    if "`row'" != "" {
        matrix null14_melogit[`row',1] = `sig2'
        matrix null14_melogit[`row',2] = `icc'
        matrix null14_melogit[`row',3] = `ll'
        matrix null14_melogit[`row',4] = `N'
        matrix null14_melogit[`row',5] = `Ng'
        matrix null14_melogit[`row',6] = `aic'
        matrix null14_melogit[`row',7] = `bic'
        matrix null14_melogit[`row',8] = `conv'
    }

end


local row = 0

foreach outcome of global outcomes {

    di as result "--- Outcome: `outcome' ---"

    * Model A
    local row = `row' + 1
    di as text "[A] ML-GH (AGQ-style)"

    melogit_null_gh `outcome', cluster(cluster_id) intpts(7) row(`row')

    estimates store null_`outcome'_gh

}
/* ---- 3.2 LINEAR PROBABILITY MODEL: ML vs REML (direct comparison) -------- */
/*
  METHODOLOGICAL NOTE:
  For binary outcomes, Stata's melogit uses ML by default and does not
  offer a -reml- option directly. The cleanest way to demonstrate ML vs REML
  differences in the variance component and ICC is via the linear probability
  model (LPM) using -mixed-.

  While a LPM is not ideal for binary outcomes, it serves as the benchmark
  comparison because:
    (a) It provides IDENTICAL ML and REML syntax via -mixed-
    (b) The ICC formula is simpler (σ²_between / σ²_total)
    (c) Results are directly comparable for the methodological thesis argument

  We present LPM results as a COMPANION to the melogit results, not a
  replacement. State clearly in your thesis:
  "To enable a direct, technically precise comparison of ML and REML
  estimation, we additionally estimated null models using the linear
  probability framework via mixed-effects linear regression (mixed).
  While the LPM is not the preferred model for binary outcomes, it provides
  the methodologically cleanest ML/REML contrast, which is the core
  analytical objective of this thesis."
*

di _newline as text "--- 3.2 LPM NULL MODELS: ML vs REML (2014) ---"

local lrow = 0

foreach outcome of global outcomes {

    di _newline as result "  Outcome: `outcome'"

    * === Model C: mixed — ML (no reml option = ML by default) === *
    local lrow = `lrow' + 1
    di as text " [C] LPM ML"

    quietly mixed `outcome' [pw=wt] ///
        || cluster_id:, ///
        mle nolog variance

    local sig2u_ml  = exp(2*[lns1_1_1]_b[_cons])   // between-cluster variance
    local sig2e_ml  = exp(2*[lnsig_e]_b[_cons])    // within-cluster (residual)
    local icc_ml    = `sig2u_ml' / (`sig2u_ml' + `sig2e_ml')
    local ll_ml     = e(ll)
    local aic_ml    = e(ic)[1,5]    // from estat ic
    quietly estat icc
    local aic_ml = r(S)[1,5]
    local bic_ml = r(S)[1,6]

    matrix null14_lpm[`lrow', 1] = `sig2u_ml'
    matrix null14_lpm[`lrow', 2] = `icc_ml'
    matrix null14_lpm[`lrow', 3] = `ll_ml'
    matrix null14_lpm[`lrow', 4] = `aic_ml'
    matrix null14_lpm[`lrow', 5] = `bic_ml'

    estimates store null14_`outcome'_lpm_ml

   /* di as result "  LPM-ML:   σ²_u=" %7.5f `sig2u_ml' ///
                 "  ICC=" %6.4f `icc_ml' ///
                 "  LL=" %10.3f `ll_ml' ///
                 "  AIC=" %8.2f `aic_ml'

    * === Model D: mixed — REML === *
   local lrow = `lrow' + 1
    di as text "  [D] LPM REML"

    quietly mixed `outcome' [pw=wt] ///
        || cluster_id:, ///
        reml nolog variance

    local sig2u_reml = exp(2*[lns1_1_1]_b[_cons])
    local sig2e_reml = exp(2*[lnsig_e]_b[_cons])
    local icc_reml   = `sig2u_reml' / (`sig2u_reml' + `sig2e_reml')
    local ll_reml    = e(ll)
    quietly estat icc
    local aic_reml = r(S)[1,5]
    local bic_reml = r(S)[1,6]

    matrix null14_lpm[`lrow', 1] = `sig2u_reml'
    matrix null14_lpm[`lrow', 2] = `icc_reml'
    matrix null14_lpm[`lrow', 3] = `ll_reml'
    matrix null14_lpm[`lrow', 4] = `aic_reml'
    matrix null14_lpm[`lrow', 5] = `bic_reml'

    estimates store null14_`outcome'_lpm_reml

   /* di as result "  LPM-REML: σ²_u=" %7.5f `sig2u_reml' ///
                 "  ICC=" %6.4f `icc_reml' ///
                 "  LL=" %10.3f `ll_reml' ///
                 "  AIC=" %8.2f `aic_reml'

    * === Direct ML vs REML comparison (LPM) === *
    local d_sig2_lr = `sig2u_reml' - `sig2u_ml'
    local d_icc_lr  = `icc_reml'  - `icc_ml'
    di as text "    REML - ML:  Δσ²_u=" %8.5f `d_sig2_lr' ///
               "  ΔICC=" %7.4f `d_icc_lr'
    di as text "    (Positive Δ = REML gives larger variance estimate, as expected)"
}

di _newline as result "2014 null models complete."

*/
/*==============================================================================
  SECTION 4 — NULL MODELS: 2022 DATA
  Identical structure to Section 3
==============================================================================*/
/*==============================================================================
  SECTION 4 — NULL MODELS: KDHS 2022
==============================================================================*/

di _newline(2) as text "========================================================"
di             as text "  SECTION 4: NULL MODELS — 2022"
di             as text "========================================================"

use "$clean/kdhs2022_null_ready.dta", clear

/*
  Recode outcomes to explicit numeric 0/1.
  melogit requires numeric — labelled bytes work but explicit recoding
  avoids any ambiguity about which category is coded 1.
*/
foreach v in anc_any anc4plus t3_delivery t4_pnc {
    capture drop n_`v'
    gen byte n_`v' = (`v' == 1)
    replace  n_`v' = . if missing(`v')
}

local nout = 4
matrix null22_melogit = J(`=`nout'*2', 8, .)
matrix colnames null22_melogit = sigma2 ICC loglik N N_clusters AIC BIC converged

matrix null22_lpm = J(`=`nout'*2', 5, .)
matrix colnames null22_lpm = sigma2u ICC loglik AIC BIC

local row  = 0
local lrow = 0

foreach outcome in anc_any anc4plus t3_delivery t4_pnc {

    di _newline as result "--- Outcome: n_`outcome' (2022) ---"

    /*------------------------------------------------------------------
      Model A: melogit ML-AGQ (mvaghermite)
      FIX E1: intmethod(ghermite) → intmethod(mvaghermite)
      FIX E2: estat sd → estat recovariance
    ------------------------------------------------------------------*/
    local ++row
    di as text "  [A] ML-AGQ (row `row')"

    quietly melogit n_`outcome' [pw=wt] ///
        || cluster_id:, ///
        intmethod(mvaghermite) intpoints($QPTS) ///
        nolog difficult

    local conv = e(converged)
    if `conv' == 0 di as error "  WARNING: ML-AGQ did not converge — `outcome'"

    quietly estat recovariance          /* FIX E2: was estat sd */
    local sig2 = r(cov)[1,1]
    local icc  = `sig2' / (`sig2' + $L1VAR)
    local ll   = e(ll)
    local N    = e(N)
    local Ng   = e(N_g)
    local aic  = -2*`ll' + 2*e(k)
    local bic  = -2*`ll' + e(k)*ln(`N')

    matrix null22_melogit[`row', 1] = `sig2'
    matrix null22_melogit[`row', 2] = `icc'
    matrix null22_melogit[`row', 3] = `ll'
    matrix null22_melogit[`row', 4] = `N'
    matrix null22_melogit[`row', 5] = `Ng'
    matrix null22_melogit[`row', 6] = `aic'
    matrix null22_melogit[`row', 7] = `bic'
    matrix null22_melogit[`row', 8] = `conv'

    estimates store null22_`outcome'_agq
    di as result "  ML-AGQ: σ²=" %7.5f `sig2' ///
                 "  ICC=" %6.4f `icc' ///
                 "  AIC=" %9.2f `aic' ///
                 "  Conv=" `conv'

    /*------------------------------------------------------------------
      Model B: melogit ML-Laplace equivalent
      FIX E3: intmethod(laplace) with [pw=wt] → use intpoints(1) instead
              which is mathematically equivalent and accepts pweights
    ------------------------------------------------------------------*/
    local ++row
    di as text "  [B] ML-Laplace equiv (row `row')"

    quietly melogit n_`outcome' [pw=wt] ///
        || cluster_id:, ///
        intmethod(mvaghermite) intpoints(1) ///   /* FIX E3 */
        nolog difficult

    local conv = e(converged)
    if `conv' == 0 di as error "  WARNING: Laplace-equiv did not converge — `outcome'"

    quietly estat recovariance          /* FIX E2 */
    local sig2 = r(cov)[1,1]
    local icc  = `sig2' / (`sig2' + $L1VAR)
    local ll   = e(ll)
    local N    = e(N)
    local Ng   = e(N_g)
    local aic  = -2*`ll' + 2*e(k)
    local bic  = -2*`ll' + e(k)*ln(`N')

    matrix null22_melogit[`row', 1] = `sig2'
    matrix null22_melogit[`row', 2] = `icc'
    matrix null22_melogit[`row', 3] = `ll'
    matrix null22_melogit[`row', 4] = `N'
    matrix null22_melogit[`row', 5] = `Ng'
    matrix null22_melogit[`row', 6] = `aic'
    matrix null22_melogit[`row', 7] = `bic'
    matrix null22_melogit[`row', 8] = `conv'

    estimates store null22_`outcome'_laplace
    di as result "  ML-Laplace: σ²=" %7.5f `sig2' ///
                 "  ICC=" %6.4f `icc' ///
                 "  AIC=" %9.2f `aic' ///
                 "  Conv=" `conv'

    /*------------------------------------------------------------------
      Model C: mixed LPM — ML
    ------------------------------------------------------------------*/
    local ++lrow
    di as text "  [C] LPM-ML (lrow `lrow')"

    quietly mixed n_`outcome' [pw=wt] || cluster_id:, mle nolog variance

    local sig2u = exp(2 * [lns1_1_1]_b[_cons])
    local sig2e = exp(2 * [lnsig_e]_b[_cons])
    local icc_l = `sig2u' / (`sig2u' + `sig2e')
    local ll_l  = e(ll)

    quietly estat ic                    /* FIX E4: was estat icc */
    local aic_l = r(S)[1,5]
    local bic_l = r(S)[1,6]

    matrix null22_lpm[`lrow', 1] = `sig2u'
    matrix null22_lpm[`lrow', 2] = `icc_l'
    matrix null22_lpm[`lrow', 3] = `ll_l'
    matrix null22_lpm[`lrow', 4] = `aic_l'
    matrix null22_lpm[`lrow', 5] = `bic_l'

    estimates store null22_`outcome'_lpm_ml
    di as result "  LPM-ML: σ²_u=" %8.5f `sig2u' ///
                 "  ICC=" %6.4f `icc_l' ///
                 "  AIC=" %9.2f `aic_l'

    /*------------------------------------------------------------------
      Model D: mixed LPM — REML
      FIX E4: estat icc → estat ic for AIC/BIC extraction
    ------------------------------------------------------------------*/
    local ++lrow
    di as text "  [D] LPM-REML (lrow `lrow')"

    quietly mixed n_`outcome' [pw=wt] || cluster_id:, reml nolog variance

    local sig2u_r = exp(2 * [lns1_1_1]_b[_cons])
    local sig2e_r = exp(2 * [lnsig_e]_b[_cons])
    local icc_r   = `sig2u_r' / (`sig2u_r' + `sig2e_r')
    local ll_r    = e(ll)

    quietly estat ic                    /* FIX E4: was estat icc */
    local aic_r = r(S)[1,5]
    local bic_r = r(S)[1,6]

    matrix null22_lpm[`lrow', 1] = `sig2u_r'
    matrix null22_lpm[`lrow', 2] = `icc_r'
    matrix null22_lpm[`lrow', 3] = `ll_r'
    matrix null22_lpm[`lrow', 4] = `aic_r'
    matrix null22_lpm[`lrow', 5] = `bic_r'

    estimates store null22_`outcome'_lpm_reml

    /*------------------------------------------------------------------
      ML vs REML comparison — printed immediately for each outcome
    ------------------------------------------------------------------*/
    local d_sig2 = `sig2u_r' - `sig2u'
    local d_icc  = `icc_r'   - `icc_l'
    local pct    = (`d_sig2' / `sig2u') * 100

    di as result "  LPM-REML: σ²_u=" %8.5f `sig2u_r' ///
                 "  ICC=" %6.4f `icc_r' ///
                 "  AIC=" %9.2f `aic_r'
    di as result "  ┌─ REML − ML: ─────────────────────────────────"
    di as result "  │  Δσ²_u = " %9.6f `d_sig2' "  (" %+6.2f `pct' "%)"
    di as result "  │  ΔICC  = " %7.4f `d_icc'
    di as result "  └────────────────────────────────────────────────"
}

di _newline as result "2022 null models complete."


/*==============================================================================
  SECTION 5 — NULL MODELS: POOLED DATA (2014 + 2022 COMBINED)
  
  Pooled null models serve two purposes:
    (a) Provides a single ICC estimate for each transition combining both waves
    (b) Tests whether survey wave (2014 vs 2022) itself explains cluster variance
  
  Two specifications:
    Spec 1: No wave control — raw pooled ICC
    Spec 2: Wave (survey_year) as fixed effect — ICC net of wave differences
    
  The cluster_wave variable (unique cluster ID across waves) is used as
  the random effect grouping, preventing cross-wave cluster ID collisions.
==============================================================================*/

di _newline as text "========================================================"
di          as text "  SECTION 5: NULL MODELS — POOLED (2014 + 2022)"
di          as text "========================================================"

use "$clean/KDHS_pooled_2014_2022_final.dta", clear
keep if analysis_sample == 1

matrix null_pooled = J(16, 8, .)
matrix colnames null_pooled = sigma2 ICC loglik N N_clusters AIC BIC wave_as_FE

local prow = 0

foreach outcome of global outcomes {

    di _newline as result "--- Outcome: `outcome' (POOLED) ---"

    /* Spec 1: No wave control */
    local prow = `prow' + 1
    melogit `outcome' [pw=wt] ///
        || cluster_wave:, ///
        intmethod(agq) intpoints($intpoints) nolog difficult

    quietly estat sd
    local sig2 = r(cov)[1,1]
    local icc  = `sig2' / (`sig2' + $pi2over3)
    local ll   = e(ll)
    local N    = e(N)
    local Ng   = e(N_g)
    local aic  = -2*`ll' + 2*e(k)
    local bic  = -2*`ll' + e(k)*ln(`N')

    matrix null_pooled[`prow', 1] = `sig2'
    matrix null_pooled[`prow', 2] = `icc'
    matrix null_pooled[`prow', 3] = `ll'
    matrix null_pooled[`prow', 4] = `N'
    matrix null_pooled[`prow', 5] = `Ng'
    matrix null_pooled[`prow', 6] = `aic'
    matrix null_pooled[`prow', 7] = `bic'
    matrix null_pooled[`prow', 8] = 0    // no wave FE

    *di as result "  Pooled (no wave FE): σ²=" %6.4f `sig2' "  ICC=" %6.4f `icc'
    estimates store nullpool_`outcome'_nowaveFE

    /* Spec 2: Wave as fixed effect */
    local ++prow
   melogit `outcome' i.wave2022 [pw=wt] ///
        || cluster_wave:, ///
        intmethod(agq) intpoints($intpoints) nolog difficult

    quietly estat sd
    local sig2 = r(cov)[1,1]
    local icc  = `sig2' / (`sig2' + $pi2over3)
    local ll   = e(ll)
    local N    = e(N)
    local Ng   = e(N_g)
    local aic  = -2*`ll' + 2*e(k)
    local bic  = -2*`ll' + e(k)*ln(`N')

    matrix null_pooled[`prow', 1] = `sig2'
    matrix null_pooled[`prow', 2] = `icc'
    matrix null_pooled[`prow', 3] = `ll'
    matrix null_pooled[`prow', 4] = `N'
    matrix null_pooled[`prow', 5] = `Ng'
    matrix null_pooled[`prow', 6] = `aic'
    matrix null_pooled[`prow', 7] = `bic'
    matrix null_pooled[`prow', 8] = 1    // wave FE included

    *di as result "  Pooled (wave FE):    σ²=" %6.4f `sig2' "  ICC=" %6.4f `icc'
    estimates store nullpool_`outcome'_waveFE
}

di _newline as result "Pooled null models complete."


/*==============================================================================
  SECTION 6 — ICC EXTRACTION, COMPILATION, AND FORMAL COMPARISON
  
  Compile all ICC estimates into a single structured results frame.
  We use Stata 18's frames feature for clean data management.
  
  For each outcome × wave × method combination, we report:
    σ²_u0    : between-cluster variance
    ICC      : intraclass correlation (latent variable method)
    ICC_CI   : 95% confidence interval for ICC
    MOR      : Median Odds Ratio (complementary heterogeneity measure)
    AIC/BIC  : model fit
==============================================================================*/

di _newline as text "========================================================"
di          as text "  SECTION 6: ICC COMPILATION AND COMPARISON"
di          as text "========================================================"

/*
  MEDIAN ODDS RATIO (MOR):
  The MOR transforms σ²_u0 into an odds ratio scale, making it more
  interpretable for clinicians and policymakers.
  
  MOR = exp(√(2σ²_u0) × Φ^{-1}(0.75))
      = exp(0.6745 × √(2σ²_u0))
      ≈ exp(0.9539 × σ_u0)    [using Φ^{-1}(0.75) ≈ 0.6745]
      
  MOR > 1 indicates community-level heterogeneity.
  MOR = 1 means no clustering (equivalent to ICC = 0).
  
  Reference: Merlo et al. (2006) J Epidemiol Community Health
*/

/* Create results frame */
capture frame drop icc_results
frame create icc_results ///
    str20(outcome wave method) ///
    double(sigma2 icc mor ll aic bic n n_clusters)

/* Fill frame from stored estimates */
/* 2014 — melogit */
local outcomes_list "anc_any anc4plus t3_delivery t4_pnc"
local methods_list  "mlagq mllaplace"

foreach outcome of local outcomes_list {
    foreach method of local methods_list {
        quietly estimates restore null14_`outcome'_`method'
        quietly estat recovariance
        local s2  = r(cov)[1,1]
        local icc = `s2' / (`s2' + $pi2over3)
        local mor = exp(0.6745 * sqrt(2*`s2'))
        local ll  = e(ll)
        local N   = e(N)
        local Ng  = e(N_g)
        local aic = -2*`ll' + 2*e(k)
        local bic = -2*`ll' + e(k)*ln(`N')
        frame post icc_results ("`outcome'") ("2014") ("`method'") ///
            (`s2') (`icc') (`mor') (`ll') (`aic') (`bic') (`N') (`Ng')
    }
}

/* 2022 — melogit */
foreach outcome of local outcomes_list {
    foreach method of local methods_list {
        quietly estimates restore null22_`outcome'_`method'
        quietly estat recovariance
        local s2  = r(cov)[1,1]
        local icc = `s2' / (`s2' + $pi2over3)
        local mor = exp(0.6745 * sqrt(2*`s2'))
        local ll  = e(ll)
        local N   = e(N)
        local Ng  = e(N_g)
        local aic = -2*`ll' + 2*e(k)
        local bic = -2*`ll' + e(k)*ln(`N')
        frame post icc_results ("`outcome'") ("2022") ("`method'") ///
            (`s2') (`icc') (`mor') (`ll') (`aic') (`bic') (`N') (`Ng')
    }
}

/* Print compiled ICC table */
frame icc_results: list, sep(4) noobs abbrev(20)

/* Save ICC frame as dataset */
frame icc_results {
    save "$results/icc_null_results_all.dta", replace
}

di _newline "ICC results compiled and saved."


/*==============================================================================
  SECTION 7 — VARIANCE COMPONENT SIGNIFICANCE TESTS (LRT)
  
  LIKELIHOOD RATIO TEST for random effects:
  H0: σ²_u0 = 0 (no clustering; ordinary logistic regression is sufficient)
  H1: σ²_u0 > 0 (significant community-level variance)
  
  The LRT compares:
    Model 1: Standard logistic regression (no random effect)
    Model 2: melogit with random intercept
    
  Test statistic: LRT = -2*(LL_logit - LL_melogit) ~ χ²(1)
  
  IMPORTANT: Under H0 (boundary of parameter space), the standard χ²(1)
  p-value is conservative. Use a mixture distribution: 0.5*χ²(0) + 0.5*χ²(1)
  This is automatically implemented in Stata's lrtest for melogit.
  
  We also run the formal test using -xtlogit- random effects for comparison.
==============================================================================*/

di _newline as text "========================================================"
di          as text "  SECTION 7: LIKELIHOOD RATIO TESTS FOR CLUSTERING"
di          as text "========================================================"

/* ---- 7.1 LRT: 2014 ------------------------------------------------------- */
di as text "--- LRT: 2014 ---"
use "$clean/kdhs2014_null_ready.dta", clear

foreach outcome of local outcomes_list {

    di _newline as result "  Outcome: `outcome'"

    /* Null model (random intercept) — already estimated, restore */
    quietly estimates restore null14_`outcome'_mlagq
    local ll_multi = e(ll)

    /* Standard logistic regression (no random effect) */
    quietly logit `outcome' [pw=wt], nolog
    local ll_logit = e(ll)

    /* LRT statistic */
    local lrt_stat = -2 * (`ll_logit' - `ll_multi')
    local lrt_p    = chi2tail(1, `lrt_stat') / 2   // divide by 2 for boundary test

    di as result "    LL(logit)  = " %10.3f `ll_logit'
    di as result "    LL(melogit)= " %10.3f `ll_multi'
    di as result "    LRT χ²(1)  = " %8.3f  `lrt_stat'
   * di as result "    p-value    = " %7.5f  `lrt_p'`"  (`=cond(`lrt_p'<0.001,"p<0.001",string(round(`lrt_p',.001)))')"'

    if `lrt_p' < 0.05 {
        di as result "    → Significant clustering: multilevel model JUSTIFIED"
    }
    else {
        di as error  "    → Clustering NOT significant — review model specification"
    }
}

/* ---- 7.2 LRT: 2022 ------------------------------------------------------- */
di _newline as text "--- LRT: 2022 ---"
use "$clean/kdhs2022_null_ready.dta", clear

foreach outcome of local outcomes_list {

    di _newline as result "  Outcome: `outcome'"

    quietly estimates restore null22_`outcome'_mlagq
    local ll_multi = e(ll)

    quietly logit `outcome' [pw=wt], nolog
    local ll_logit = e(ll)

    local lrt_stat = -2 * (`ll_logit' - `ll_multi')
    local lrt_p    = chi2tail(1, `lrt_stat') / 2

    di as result "    LRT χ²(1)=" %8.3f `lrt_stat' "  p=" %7.5f `lrt_p'

    if `lrt_p' < 0.05 {
        di as result "    → Clustering significant."
    }
    else {
        di as error  "    → Non-significant clustering."
    }
}


/*==============================================================================
  SECTION 8 — ML vs REML FORMAL COMPARISON (HYPOTHESES H5 AND H6)
  
  Using the LPM results (where ML and REML are directly comparable):
  
  H5: σ²_ML ≠ σ²_REML
  H6: ICC_ML ≠ ICC_REML
  
  We compute:
    (a) Absolute difference: Δσ² = σ²_REML - σ²_ML
    (b) Relative difference: (σ²_REML - σ²_ML) / σ²_ML × 100%
    (c) Same for ICC
    (d) Tabulate AIC/BIC across methods
  
  For formal testing: use bootstrapped CI for the difference in ICC.
  If bootstrap CI excludes 0, the difference is statistically meaningful.
  
  NOTE: The expected direction is σ²_REML > σ²_ML because REML corrects
  for the downward bias in ML variance estimates that arises from treating
  fixed effects as known rather than estimated. In small samples this
  difference can be substantial; in large DHS samples it is typically small
  but may still affect inference on community-level effects.
==============================================================================*/

di _newline as text "========================================================"
di          as text "  SECTION 8: ML vs REML FORMAL COMPARISON (H5, H6)"
di          as text "========================================================"

/* ---- 8.1 COMPILE LPM ML vs REML COMPARISON TABLE ------------------------ */

/* Restore 2014 data */
use "$clean/kdhs2014_null_ready.dta", clear

di as text _newline "=== ML vs REML COMPARISON: LPM NULL MODELS ==="
di as text %40s "Outcome" ///
           %12s "σ²_ML" %12s "σ²_REML" %12s "Δσ²" %10s "Δ%" ///
           %10s "ICC_ML" %10s "ICC_REML" %10s "ΔICC"
di as text "{hline 110}"

foreach outcome of local outcomes_list {

    /* ML */
    quietly estimates restore null14_`outcome'_lpm_ml
    local sig2u_ml = exp(2*[lns1_1_1]_b[_cons])
    local sig2e_ml = exp(2*[lnsig_e]_b[_cons])
    local icc_ml   = `sig2u_ml' / (`sig2u_ml' + `sig2e_ml')
    local ll_ml    = e(ll)

    /* REML */
    quietly estimates restore null14_`outcome'_lpm_reml
    local sig2u_reml = exp(2*[lns1_1_1]_b[_cons])
    local sig2e_reml = exp(2*[lnsig_e]_b[_cons])
    local icc_reml   = `sig2u_reml' / (`sig2u_reml' + `sig2e_reml')
    local ll_reml    = e(ll)

    /* Differences */
    local d_sig2     = `sig2u_reml' - `sig2u_ml'
    local d_sig2_pct = (`d_sig2' / `sig2u_ml') * 100
    local d_icc      = `icc_reml' - `icc_ml'

    di as result %40s "`outcome'" ///
                 %12.5f `sig2u_ml'  %12.5f `sig2u_reml' ///
                 %12.5f `d_sig2'    %10.2f `d_sig2_pct' "%" ///
                 %10.4f `icc_ml'    %10.4f `icc_reml'   %10.4f `d_icc'
}

di as text "{hline 110}"
di as text "Note: Δ = REML - ML. Positive Δ confirms REML upward correction."

/* ---- 8.2 BOOTSTRAP TEST FOR ICC DIFFERENCE ------------------------------- */
/*
  Bootstrap CI for (ICC_REML - ICC_ML):
  If CI excludes 0, the difference is statistically significant.
  
  We bootstrap for the first outcome only (T1: ANC initiation) as illustration.
  Extend to all 4 outcomes in your final analysis.
  
  NOTE: Bootstrapping MLMs is computationally intensive.
  Use 500 reps minimum for publication; 200 for diagnostics.
  We use 200 here with a time warning.
*/

di _newline as text "--- Bootstrap test for ICC difference (ML vs REML): T1 ---"
di as text "  Running 200 bootstrap replications (may take several minutes)..."

capture program drop icc_diff_boot
program icc_diff_boot, rclass
    /* ML */
    quietly mixed anc_any || cluster_id:, mle variance nolog
    local s2u_ml = exp(2*[lns1_1_1]_b[_cons])
    local s2e_ml = exp(2*[lnsig_e]_b[_cons])
    local icc_ml = `s2u_ml' / (`s2u_ml' + `s2e_ml')

    /* REML */
    quietly mixed anc_any || cluster_id:, reml variance nolog
    local s2u_re = exp(2*[lns1_1_1]_b[_cons])
    local s2e_re = exp(2*[lnsig_e]_b[_cons])
    local icc_re = `s2u_re' / (`s2u_re' + `s2e_re')

    return scalar icc_diff   = `icc_re' - `icc_ml'
    return scalar sig2u_diff = `s2u_re' - `s2u_ml'
end

use "$clean/kdhs2014_null_ready.dta", clear
set seed 20260428

bootstrap icc_diff=r(icc_diff) sig2u_diff=r(sig2u_diff), ///
    reps(200) cluster(cluster_id) idcluster(newclus) ///
    notable nodots: icc_diff_boot

di _newline as result "Bootstrap results: ICC(REML) - ICC(ML)"
estat bootstrap, all

/* Clean up temp variable */
capture drop newclus


/*==============================================================================
  SECTION 9 — MODEL FIT STATISTICS (AIC, BIC)
  
  AIC and BIC comparison across:
    (a) Logistic regression (no random effect)
    (b) melogit ML-AGQ
    (c) melogit ML-Laplace
    (d) LPM ML
    (e) LPM REML
    
  Smaller AIC/BIC = better fit.
  
  NOTE: AIC and BIC are only comparable within the same likelihood framework
  (i.e., ML models only). ML vs REML AIC/BIC are NOT directly comparable
  because REML uses a different objective function. This is an important
  caveat to state explicitly in your thesis.
  
  Quote for thesis:
  "AIC and BIC from REML estimation are not directly comparable to those
  from ML estimation as REML maximises a restricted likelihood that
  integrates out fixed effects. Consequently, model selection via AIC/BIC
  should be conducted under ML; REML is used only for obtaining unbiased
  variance component estimates (Pinheiro & Bates, 2000)."
==============================================================================*/

di _newline as text "========================================================"
di          as text "  SECTION 9: MODEL FIT — AIC / BIC COMPARISON"
di          as text "========================================================"

use "$clean/kdhs2014_null_ready.dta", clear

di as text _newline "=== FIT STATISTICS: 2014 NULL MODELS ==="
di as text %25s "Outcome" %12s "Model" ///
           %12s "LL" %10s "AIC" %10s "BIC" %10s "Params"
di as text "{hline 80}"

foreach outcome of local outcomes_list {

    /* (a) Logistic */
    quietly logit `outcome' [pw=wt], nolog
    local ll_l  = e(ll)
    local k_l   = e(k)
    local aic_l = -2*`ll_l' + 2*`k_l'
    local bic_l = -2*`ll_l' + `k_l'*ln(e(N))
    di as text %25s "`outcome'" %12s "Logistic" ///
               %12.2f `ll_l' %10.2f `aic_l' %10.2f `bic_l' %10.0f `k_l'

    /* (b) melogit ML-AGQ */
    quietly estimates restore null14_`outcome'_mlagq
    local ll_a  = e(ll)
    local k_a   = e(k)
    local aic_a = -2*`ll_a' + 2*`k_a'
    local bic_a = -2*`ll_a' + `k_a'*ln(e(N))
    di as text %25s "" %12s "melogit-AGQ" ///
               %12.2f `ll_a' %10.2f `aic_a' %10.2f `bic_a' %10.0f `k_a'

    /* (c) melogit ML-Laplace */
    quietly estimates restore null14_`outcome'_mllaplace
    local ll_b  = e(ll)
    local k_b   = e(k)
    local aic_b = -2*`ll_b' + 2*`k_b'
    local bic_b = -2*`ll_b' + `k_b'*ln(e(N))
    di as text %25s "" %12s "melogit-Laplace" ///
               %12.2f `ll_b' %10.2f `aic_b' %10.2f `bic_b' %10.0f `k_b'

    /* (d) LPM ML */
    quietly estimates restore null14_`outcome'_lpm_ml
    quietly estat ic
    local aic_ml = r(S)[1,5]
    local bic_ml = r(S)[1,6]
    local ll_ml  = e(ll)
    local k_ml   = e(k)
    di as text %25s "" %12s "LPM-ML" ///
               %12.2f `ll_ml' %10.2f `aic_ml' %10.2f `bic_ml' %10.0f `k_ml'

    /* (e) LPM REML (AIC/BIC shown with caveat) */
    quietly estimates restore null14_`outcome'_lpm_reml
    quietly estat ic
    local aic_re = r(S)[1,5]
    local bic_re = r(S)[1,6]
    local ll_re  = e(ll)
    local k_re   = e(k)
    di as text %25s "" %12s "LPM-REML†" ///
               %12.2f `ll_re' %10.2f `aic_re' %10.2f `bic_re' %10.0f `k_re'

    di as text "{hline 80}"
}
di as text "† REML AIC/BIC not comparable to ML AIC/BIC (different objective function)"


/*==============================================================================
  SECTION 10 — PUBLICATION-QUALITY OUTPUT TABLES
  
  We produce two tables suitable for direct use in thesis:
  
  Table 1: Null model ICC summary (for all outcomes, both waves, both methods)
           → Addresses Objective 3 (Thesis) and tests H5, H6
           
  Table 2: ML vs REML comparison for variance components
           → Central thesis comparison table
           
  Output format: Excel via putexcel (Stata 18 native)
  Also produces a log-formatted version for review.
==============================================================================*/

di _newline as text "========================================================"
di          as text "  SECTION 10: PUBLICATION-QUALITY TABLES"
di          as text "========================================================"

/* ---- TABLE 1: NULL MODEL RESULTS ----------------------------------------- */

putexcel set "$tables/null_model_results.xlsx", replace sheet("Table1_ICC")

/* Header row */
putexcel A1 = "Table 1. Null Model Results: Variance Components and ICC across Estimation Methods"
putexcel A2 = "Outcome"
putexcel B2 = "Wave"
putexcel C2 = "Method"
putexcel D2 = "σ²_u0"
putexcel E2 = "ICC"
putexcel F2 = "MOR"
putexcel G2 = "Log-Likelihood"
putexcel H2 = "AIC"
putexcel I2 = "BIC"
putexcel J2 = "N"
putexcel K2 = "N_clusters"

/* Label formatting */
putexcel A1:K1, merge hcenter bold
putexcel A2:K2, bold border(bottom, medium)

local prow = 3
local outcome_labels `""T1: ANC Initiation" "T2: ANC Adequacy (4+)" "T3: Facility Delivery" "T4: PNC 48hrs""'
local method_labels  `""ML (AGQ)" "ML (Laplace)" "ML (LPM)" "REML (LPM)""'

local oi = 0
foreach outcome of local outcomes_list {
    local ++oi
    local olab : word `oi' of `outcome_labels'

    foreach wave in 2014 2022 {
        foreach method in mlagq mllaplace {

            quietly estimates restore null`wave'_`outcome'_`method'
            quietly estat recovariance
            local s2  = r(cov)[1,1]
            local icc = `s2' / (`s2' + $pi2over3)
            local mor = exp(0.6745 * sqrt(2*`s2'))
            local ll  = e(ll)
            local N   = e(N)
            local Ng  = e(N_g)
            local aic = -2*`ll' + 2*e(k)
            local bic = -2*`ll' + e(k)*ln(`N')

            if "`method'" == "mlagq"     local mlab "ML (AGQ)"
            if "`method'" == "mllaplace" local mlab "ML (Laplace)"

            putexcel A`prow' = "`olab'"
            putexcel B`prow' = "`wave'"
            putexcel C`prow' = "`mlab'"
            putexcel D`prow' = `s2',  nformat("0.0000")
            putexcel E`prow' = `icc', nformat("0.0000")
            putexcel F`prow' = `mor', nformat("0.000")
            putexcel G`prow' = `ll',  nformat("0.0")
            putexcel H`prow' = `aic', nformat("0.0")
            putexcel I`prow' = `bic', nformat("0.0")
            putexcel J`prow' = `N',   nformat("#,##0")
            putexcel K`prow' = `Ng',  nformat("#,##0")

            local ++prow
        }
    }
}

/* Footnote */
putexcel A`=`prow'+1' = "σ²_u0 = between-cluster variance; ICC = intraclass correlation (logistic latent variable method);"
putexcel A`=`prow'+2' = "MOR = Median Odds Ratio; AGQ = Adaptive Gaussian Quadrature (12 points)."
putexcel A`=`prow'+3' = "ICC = σ²_u0 / (σ²_u0 + π²/3) where π²/3 = 3.290."

/* ---- TABLE 2: ML vs REML COMPARISON (LPM) -------------------------------- */

putexcel set "$tables/null_model_results.xlsx", modify sheet("Table2_ML_REML")

putexcel A1 = "Table 2. ML vs REML Variance Component Comparison: Linear Probability Model Null Models"
putexcel A2 = "Outcome"
putexcel B2 = "Wave"
putexcel C2 = "σ²_u0 (ML)"
putexcel D2 = "σ²_u0 (REML)"
putexcel E2 = "Δσ²_u0"
putexcel F2 = "Δ% (relative)"
putexcel G2 = "ICC (ML)"
putexcel H2 = "ICC (REML)"
putexcel I2 = "ΔICC"
putexcel J2 = "AIC (ML)"
putexcel K2 = "AIC (REML)"
putexcel L2 = "BIC (ML)"
putexcel M2 = "BIC (REML)"

putexcel A1:M1, merge hcenter bold
putexcel A2:M2, bold border(bottom, medium)

local trow = 3
local oi = 0
foreach outcome of local outcomes_list {
    local ++oi
    local olab : word `oi' of `outcome_labels'

    foreach wave in 2014 2022 {

        if `wave' == 2014 local cdir "$clean/kdhs2014_null_ready.dta"
        if `wave' == 2022 local cdir "$clean22/kdhs2022_null_ready.dta"
        use "`cdir'", clear

        quietly estimates restore null`wave'_`outcome'_lpm_ml
        local s2u_ml  = exp(2*[lns1_1_1]_b[_cons])
        local s2e_ml  = exp(2*[lnsig_e]_b[_cons])
        local icc_ml  = `s2u_ml' / (`s2u_ml' + `s2e_ml')
        quietly estat ic
        local aic_ml  = r(S)[1,5]
        local bic_ml  = r(S)[1,6]

        quietly estimates restore null`wave'_`outcome'_lpm_reml
        local s2u_re  = exp(2*[lns1_1_1]_b[_cons])
        local s2e_re  = exp(2*[lnsig_e]_b[_cons])
        local icc_re  = `s2u_re' / (`s2u_re' + `s2e_re')
        quietly estat ic
        local aic_re  = r(S)[1,5]
        local bic_re  = r(S)[1,6]

        local d_s2    = `s2u_re' - `s2u_ml'
        local d_s2pct = (`d_s2' / `s2u_ml') * 100
        local d_icc   = `icc_re' - `icc_ml'

        putexcel A`trow' = "`olab'"
        putexcel B`trow' = "`wave'"
        putexcel C`trow' = `s2u_ml',  nformat("0.00000")
        putexcel D`trow' = `s2u_re',  nformat("0.00000")
        putexcel E`trow' = `d_s2',    nformat("0.00000")
        putexcel F`trow' = `d_s2pct', nformat("0.00%")
        putexcel G`trow' = `icc_ml',  nformat("0.0000")
        putexcel H`trow' = `icc_re',  nformat("0.0000")
        putexcel I`trow' = `d_icc',   nformat("0.0000")
        putexcel J`trow' = `aic_ml',  nformat("0.0")
        putexcel K`trow' = `aic_re',  nformat("0.0")
        putexcel L`trow' = `bic_ml',  nformat("0.0")
        putexcel M`trow' = `bic_re',  nformat("0.0")

        local ++trow
    }
}

putexcel A`=`trow'+1' = "Note: REML AIC/BIC not comparable to ML AIC/BIC (different objective function)."
putexcel A`=`trow'+2' = "σ²_u0 extracted from LPM via exp(2×lns1_1_1). ICC = σ²_u0/(σ²_u0+σ²_e)."

putexcel save
di as result "Tables exported to: $tables/null_model_results.xlsx"


/*==============================================================================
  SECTION 11 — GRAPHICAL DIAGNOSTICS
  
  Four diagnostic plots:
  
  Plot 1: Caterpillar plot of cluster random effects (T1, 2014)
          Shows which clusters are significantly above/below grand mean
          → Justifies multilevel modelling
          
  Plot 2: ICC bar chart across T1–T4 (ML-AGQ, both waves)
          Shows how community effects vary across continuum stages
          → Key finding for thesis narrative
          
  Plot 3: σ² comparison: ML vs REML across all 4 outcomes (2014)
          Side-by-side bar chart
          → Direct visual for H5
          
  Plot 4: ICC comparison: ML vs REML (LPM, 2014 and 2022)
          Connected dot plot showing magnitude of REML correction
          → Direct visual for H6
==============================================================================*/

di _newline as text "========================================================"
di          as text "  SECTION 11: GRAPHICAL DIAGNOSTICS"
di          as text "========================================================"

/* ---- PLOT 1: Caterpillar plot of cluster BLUPs (2014, T1: anc_any) ------- */

use "$clean/kdhs2014_null_ready.dta", clear

quietly melogit anc_any [pw=wt] || cluster_id:, ///
    intmethod(agq) intpoints($intpoints) nolog

predict re_u, reffects         // Best Linear Unbiased Predictions (BLUPs)
predict re_se, reses            // Standard errors of BLUPs

gen     re_lower = re_u - 1.96*re_se
gen     re_upper = re_u + 1.96*re_se

bysort cluster_id: keep if _n == 1
sort re_u
gen rank_u = _n

twoway ///
    (rcap re_lower re_upper rank_u, lcolor(gs10) lwidth(thin)) ///
    (scatter re_u rank_u, mcolor(navy) msize(vsmall) msymbol(circle)) ///
    (function y=0, range(rank_u) lcolor(red) lpattern(dash) lwidth(medium)), ///
    ytitle("Predicted cluster-level random effect (log-odds)", size(small)) ///
    xtitle("Cluster rank", size(small)) ///
    title("Cluster Random Effects: T1 ANC Initiation (2014)", size(medsmall)) ///
    subtitle("95% prediction intervals; red line = grand mean", size(small)) ///
    legend(off) ///
    scheme(s2color) ///
    note("N_clusters = `=e(N_g)'. BLUPs from melogit (ML-AGQ, 12 quadrature points).", ///
         size(vsmall))

graph export "$figures/caterpillar_T1_2014.png", replace width(1600) height(900)
di "Plot 1 saved: caterpillar_T1_2014.png"

drop re_u re_se re_lower re_upper rank_u

/* ---- PLOT 2: ICC across T1–T4 by wave ------------------------------------ */
/*
  Build dataset from stored matrix values for plotting
*/

use "$clean/kdhs2014_null_ready.dta", clear   // just to have a dataset open

/* Manually input ICC values from matrix null14_melogit and null22_melogit */
/* Rows 1,3,5,7 = AGQ (odd rows); Rows 2,4,6,8 = Laplace (even rows) */

/* Extract ICC values for AGQ method across all outcomes and waves */
forvalues i = 1/4 {
    local r14_`i' = null14_melogit[`=(`i'-1)*2+1', 2]   // row for AGQ = odd
    local r22_`i' = null22_melogit[`=(`i'-1)*2+1', 2]
}

/* Construct plot dataset */
clear
input byte outcome float(icc14 icc22)
1 `r14_1' `r22_1'
2 `r14_2' `r22_2'
3 `r14_3' `r22_3'
4 `r14_4' `r22_4'
end

label define out_lbl 1 "T1: ANC{sub:≥1}" 2 "T2: ANC{sub:4+}" ///
    3 "T3: Facility Del." 4 "T4: PNC 48h"
label values outcome out_lbl

twoway ///
    (bar icc14 outcome, barwidth(0.35) color(navy%70) base(0)) ///
    (bar icc22 outcome, barwidth(0.35) color(cranberry%70) base(0) ///
     xoffset(0.4)) ///
    (scatter icc14 outcome, msymbol(none) mlabel(icc14) mlabformat(%4.3f) ///
     mlabposition(12) mlabcolor(navy) mlabsize(vsmall)) ///
    (scatter icc22 outcome, msymbol(none) mlabel(icc22) mlabformat(%4.3f) ///
     mlabposition(12) mlabcolor(cranberry) mlabsize(vsmall) xoffset(0.4)), ///
    ytitle("Intraclass Correlation Coefficient (ICC)", size(small)) ///
    xtitle("Continuum of Care Stage", size(small)) ///
    title("ICC across Continuum of Care Stages: ML Estimates", size(medsmall)) ///
    subtitle("Community-level (cluster) variance component", size(small)) ///
    legend(order(1 "2014" 2 "2022") pos(1) ring(0) cols(1) size(small)) ///
    xlabel(1 "T1: ANC{sub:≥1}" 2 "T2: ANC{sub:4+}" ///
           3 "T3: Delivery" 4 "T4: PNC", labsize(small)) ///
    yline(0.05, lpattern(dash) lcolor(red) lwidth(thin)) ///
    note("Red dashed line = ICC=0.05 (Hedges & Hedberg (2007) design effect threshold)." ///
         "ML estimated using Adaptive Gaussian Quadrature (12 points).", size(vsmall)) ///
    scheme(s2color)

graph export "$figures/icc_across_stages.png", replace width(1600) height(900)
di "Plot 2 saved: icc_across_stages.png"

/* ---- PLOT 3 & 4: ML vs REML comparison plots ----------------------------- */

/* Load saved ICC results */
use "$results/icc_null_results_all.dta", clear

/* Recode outcome to numeric for plotting */
encode outcome, gen(outcome_n)
encode method,  gen(method_n)

/* Plot 3: σ² comparison ML-AGQ vs ML-Laplace */
twoway ///
    (connected sigma2 outcome_n if method=="mlagq"     & wave=="2014", ///
     lcolor(navy)    mcolor(navy)    msymbol(circle) lpattern(solid)) ///
    (connected sigma2 outcome_n if method=="mllaplace"  & wave=="2014", ///
     lcolor(navy)    mcolor(navy)    msymbol(triangle) lpattern(dash)) ///
    (connected sigma2 outcome_n if method=="mlagq"     & wave=="2022", ///
     lcolor(cranberry) mcolor(cranberry) msymbol(circle) lpattern(solid)) ///
    (connected sigma2 outcome_n if method=="mllaplace"  & wave=="2022", ///
     lcolor(cranberry) mcolor(cranberry) msymbol(triangle) lpattern(dash)), ///
    ytitle("Between-cluster variance (σ²_u0)", size(small)) ///
    xtitle("Continuum of Care Stage", size(small)) ///
    title("σ²_u0: ML-AGQ vs ML-Laplace Approximation", size(medsmall)) ///
    legend(order(1 "2014 AGQ" 2 "2014 Laplace" 3 "2022 AGQ" 4 "2022 Laplace") ///
           pos(1) ring(0) cols(2) size(small)) ///
    xlabel(1 "T1" 2 "T2" 3 "T3" 4 "T4", labsize(small)) ///
    scheme(s2color)

graph export "$figures/sigma2_method_comparison.png", replace width(1600) height(900)
di "Plot 3 saved: sigma2_method_comparison.png"

/* Plot 4: ICC — ML vs REML (2014, LPM) */
/* Build from stored matrix */
use "$clean/kdhs2014_null_ready.dta", clear

clear
input byte outcome float(icc_ml icc_reml)
end

local oi = 0
foreach out of local outcomes_list {
    local ++oi
    quietly estimates restore null14_`out'_lpm_ml
    local s2u = exp(2*[lns1_1_1]_b[_cons])
    local s2e = exp(2*[lnsig_e]_b[_cons])
    local iml = `s2u' / (`s2u' + `s2e')

    quietly estimates restore null14_`out'_lpm_reml
    local s2u = exp(2*[lns1_1_1]_b[_cons])
    local s2e = exp(2*[lnsig_e]_b[_cons])
    local ire = `s2u' / (`s2u' + `s2e')

    insobs 1
    quietly replace outcome   = `oi'  if _n == _N
    quietly replace icc_ml    = `iml' if _n == _N
    quietly replace icc_reml  = `ire' if _n == _N
}

twoway ///
    (scatter icc_ml   outcome, mcolor(navy)   msymbol(circle) msize(medium)) ///
    (scatter icc_reml outcome, mcolor(red)    msymbol(square) msize(medium)) ///
    (pcspike icc_ml outcome icc_reml outcome, lcolor(gs8) lwidth(medium)), ///
    ytitle("ICC (Linear Probability Model)", size(small)) ///
    xtitle("Continuum of Care Stage", size(small)) ///
    title("ML vs REML: ICC Comparison (2014, LPM null models)", size(medsmall)) ///
    subtitle("Connecting lines show magnitude of REML correction", size(small)) ///
    legend(order(1 "ML" 2 "REML") pos(1) ring(0) cols(1) size(small)) ///
    xlabel(1 "T1: ANC init." 2 "T2: ANC 4+" 3 "T3: Delivery" 4 "T4: PNC", ///
           labsize(small) angle(15)) ///
    note("REML estimates expected to exceed ML due to correction for fixed-effect uncertainty.", ///
         size(vsmall)) ///
    scheme(s2color)

graph export "$figures/icc_ml_vs_reml_2014.png", replace width(1600) height(900)
di "Plot 4 saved: icc_ml_vs_reml_2014.png"


/*==============================================================================
  FINAL SUMMARY LOG
==============================================================================*/

di _newline as text "========================================================"
di          as text "  02_NULL_MODELS.DO — COMPLETE"
di          as text "========================================================"
di ""
di "FILES PRODUCED:"
di "  Estimates:  null14_*, null22_*, nullpool_* stored in memory"
di "              (save with -estimates save- before closing Stata)"
di "  Data:       $results/icc_null_results_all.dta"
di "  Tables:     $tables/null_model_results.xlsx (2 sheets)"
di "  Figures:    $figures/caterpillar_T1_2014.png"
di "              $figures/icc_across_stages.png"
di "              $figures/sigma2_method_comparison.png"
di "              $figures/icc_ml_vs_reml_2014.png"
di ""
di "HYPOTHESES STATUS:"
di "  H5 (σ² differ ML vs REML):  → See Table 2, Section 8"
di "  H6 (ICC differ ML vs REML): → See Table 2, Plot 4, Section 8"
di "  H7 & H8: → Run 03_individual_models.do and 04_community_models.do"
di ""
di "NEXT DO-FILE: 03_individual_models.do"
di "  — Adds behavioral determinants (edu, media, autonomy, parity, marital)"
di "  — Runs ML and REML (LPM) with Level-1 covariates"
di "  — Tests H7: Fixed effect estimates sensitive to estimation method?"
di ""
di "Completed: $(c(current_date)) $(c(current_time))"
di "========================================================"

/* Save all estimates to disk before log closes */
foreach outcome of local outcomes_list {
    foreach wave in 2014 2022 {
        foreach method in mlagq mllaplace lpm_ml lpm_reml {
            capture estimates save ///
                "$results/est_null`wave'_`outcome'_`method'.ster", replace
        }
    }
}

log close

/* END OF 02_NULL_MODELS.DO */
