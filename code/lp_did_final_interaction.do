/*
    LP-DiD: Interaction with Rule of Law — prosecution outcome
    Estimators: TWFE, LP-DiD, Reweighted LP-DiD, PMD LP-DiD, Reweighted PMD LP-DiD
    Interaction: high_rol_pre (above-median Rule of Law, pre-treatment)
    Outputs: pooled tables (LaTeX), event-study graph, event-study tables (xlsx/tex/docx)
*/

clear all
set more off
macro drop _all
capture log close

use "C:\Users\MSAM01\Documents\stata19\dataset\ds_small_incl_richard_unp_adapt_not_nan.dta", clear

set scheme s2color
set trace off

***************************************************
* Rename & restrict
***************************************************
rename prosecution    Y
rename new_ccode      unit
rename approxTrafYear time
rename unp_adapt      treat

keep if time <= 2013
keep if time >= 1990
keep if Y != .

levelsof unit, local(units_list)
local units : word count `units_list'

levelsof time, local(time_list)
local time_periods : word count `time_list'

***************************************************
* Panel & event-time setup
***************************************************
local post_window 5
local pre_window  5

xtset unit time

gen dtreat = D.treat

bysort unit (time): egen treat_year = min(cond(dtreat == 1, time, .))
gen t = time - treat_year
replace t = . if t < -`pre_window' | t > `post_window'

gen first_treat = time if dtreat == 1
bysort unit: egen cohort = min(first_treat)

***************************************************
* TWFE EVENT-STUDY ESTIMATES (with interaction)
***************************************************
gen b_twfe      = .
gen se_twfe     = .
gen N_twfe      = .
gen b_twfe_low  = .
gen se_twfe_low = .
gen b_twfe_high = .
gen se_twfe_high= .
gen b_twfe_diff = .
gen se_twfe_diff= .

xtset unit time
reghdfe Y c.L(-5/5).D.treat##i.high_rol_pre, absorb(unit time) cluster(unit)

* Pre-periods (normalised to t=-1)
forval j = 2/`pre_window' {
    lincom F`j'D.treat - FD.treat
    replace b_twfe  = r(estimate) if t == -`j'
    replace se_twfe = r(se)       if t == -`j'
    replace N_twfe  = e(N)        if t == -`j'

    replace b_twfe_low  = _b[F`j'D.treat]               if t == -`j'
    replace se_twfe_low = _se[F`j'D.treat]              if t == -`j'
    replace b_twfe_diff = _b[1.high_rol_pre#cF`j'D.treat]  if t == -`j'
    replace se_twfe_diff= _se[1.high_rol_pre#cF`j'D.treat] if t == -`j'

    lincom F`j'D.treat + 1.high_rol_pre#cF`j'D.treat
    replace b_twfe_high  = r(estimate) if t == -`j'
    replace se_twfe_high = r(se)       if t == -`j'
}

replace b_twfe      = 0 if t == -1
replace b_twfe_low  = 0 if t == -1
replace b_twfe_high = 0 if t == -1
replace b_twfe_diff = 0 if t == -1

* t=0
lincom D.treat - FD.treat
replace b_twfe  = r(estimate) if t == 0
replace se_twfe = r(se)       if t == 0
replace N_twfe  = e(N)        if t == 0
replace b_twfe_low  = _b[D.treat]                   if t == 0
replace se_twfe_low = _se[D.treat]                  if t == 0
replace b_twfe_diff = _b[1.high_rol_pre#cD.treat]   if t == 0
replace se_twfe_diff= _se[1.high_rol_pre#cD.treat]  if t == 0
lincom D.treat + 1.high_rol_pre#cD.treat
replace b_twfe_high  = r(estimate) if t == 0
replace se_twfe_high = r(se)       if t == 0

* t=1
lincom LD.treat - FD.treat
replace b_twfe  = r(estimate) if t == 1
replace se_twfe = r(se)       if t == 1
replace N_twfe  = e(N)        if t == 1
replace b_twfe_low  = _b[LD.treat]                  if t == 1
replace se_twfe_low = _se[LD.treat]                 if t == 1
replace b_twfe_diff = _b[1.high_rol_pre#cLD.treat]  if t == 1
replace se_twfe_diff= _se[1.high_rol_pre#cLD.treat] if t == 1
lincom LD.treat + 1.high_rol_pre#cLD.treat
replace b_twfe_high  = r(estimate) if t == 1
replace se_twfe_high = r(se)       if t == 1

* Post-periods
forval j = 2/`post_window' {
    lincom L`j'D.treat - FD.treat
    replace b_twfe  = r(estimate) if t == `j'
    replace se_twfe = r(se)       if t == `j'
    replace N_twfe  = e(N)        if t == `j'
    replace b_twfe_low  = _b[L`j'D.treat]               if t == `j'
    replace se_twfe_low = _se[L`j'D.treat]              if t == `j'
    replace b_twfe_diff = _b[1.high_rol_pre#cL`j'D.treat]  if t == `j'
    replace se_twfe_diff= _se[1.high_rol_pre#cL`j'D.treat] if t == `j'
    lincom L`j'D.treat + 1.high_rol_pre#cL`j'D.treat
    replace b_twfe_high  = r(estimate) if t == `j'
    replace se_twfe_high = r(se)       if t == `j'
}

gen t_twfe_diff = b_twfe_diff / se_twfe_diff
gen p_twfe_diff = 2*ttail(N_twfe, abs(t_twfe_diff))
replace se_twfe_diff = 0 if t == -1

gen tstat_twfe = b_twfe / se_twfe
gen pval_twfe  = 2*ttail(N_twfe, abs(tstat_twfe))

***************************************************
* LP-DiD ESTIMATES (with interaction)
***************************************************
* Long differences
forval j = 0/`post_window' {
    gen D`j'y = F`j'.Y - L.Y
}
forval j = 2/`pre_window' {
    gen Dm`j'y = L`j'.Y - L.Y
}

gen b_lpdid      = .
gen se_lpdid     = .
gen N_lpdid      = .
gen N_units      = .
gen b_lpdid_low  = .
gen se_lpdid_low = .
gen b_lpdid_high = .
gen se_lpdid_high= .
gen b_lpdid_diff = .
gen se_lpdid_diff= .

forval j = 0/`post_window' {
    reghdfe D`j'y c.D.treat##i.high_rol_pre ///
        if D.treat==1 | F`j'.treat==0, ///
        absorb(time) vce(cluster unit)

    preserve
        keep if e(sample)
        levelsof unit, local(u)
        local nunits : word count `u'
    restore

    replace b_lpdid      = _b[D.treat]                    if t == `j'
    replace se_lpdid     = _se[D.treat]                   if t == `j'
    replace N_lpdid      = e(N)                            if t == `j'
    replace N_units      = `nunits'                        if t == `j'
    replace b_lpdid_low  = _b[D.treat]                    if t == `j'
    replace se_lpdid_low = _se[D.treat]                   if t == `j'
    replace b_lpdid_diff = _b[1.high_rol_pre#c.D.treat]   if t == `j'
    replace se_lpdid_diff= _se[1.high_rol_pre#c.D.treat]  if t == `j'
    lincom D.treat + 1.high_rol_pre#c.D.treat
    replace b_lpdid_high  = r(estimate) if t == `j'
    replace se_lpdid_high = r(se)       if t == `j'

    if `j' > 1 & `j' <= `pre_window' {
        reghdfe Dm`j'y c.D.treat##i.high_rol_pre ///
            if D.treat==1 | treat==0, ///
            absorb(time) vce(cluster unit)

        preserve
            keep if e(sample)
            levelsof unit, local(u)
            local nunits : word count `u'
        restore

        replace b_lpdid      = _b[D.treat]                    if t == -`j'
        replace se_lpdid     = _se[D.treat]                   if t == -`j'
        replace N_lpdid      = e(N)                            if t == -`j'
        replace N_units      = `nunits'                        if t == -`j'
        replace b_lpdid_low  = _b[D.treat]                    if t == -`j'
        replace se_lpdid_low = _se[D.treat]                   if t == -`j'
        replace b_lpdid_diff = _b[1.high_rol_pre#c.D.treat]   if t == -`j'
        replace se_lpdid_diff= _se[1.high_rol_pre#c.D.treat]  if t == -`j'
        lincom D.treat + 1.high_rol_pre#c.D.treat
        replace b_lpdid_high  = r(estimate) if t == -`j'
        replace se_lpdid_high = r(se)       if t == -`j'
    }
}

gen t_lpdid_diff = b_lpdid_diff / se_lpdid_diff
gen p_lpdid_diff = 2*ttail(N_lpdid, abs(t_lpdid_diff))
replace b_lpdid_diff  = 0 if t == -1
replace se_lpdid_diff = 0 if t == -1

replace b_lpdid  = 0 if t == -1
replace se_lpdid = 0 if t == -1

gen tstat_lpdid = b_lpdid / se_lpdid
gen pval_lpdid  = 2*ttail(N_lpdid, abs(tstat_lpdid))

summ N_units
local min_units = r(min)
local max_units = r(max)

summ N_lpdid
local Nmin = r(min)
local Nmax = r(max)

***************************************************
* REWEIGHTING WEIGHTS
***************************************************
forval j = 0/`post_window' {
    qui gen group_h`j' = .
    qui replace group_h`j' = time if (D.treat==1 | F`j'.treat==0)
    qui reghdfe D.treat if (D.treat==1 | F`j'.treat==0), absorb(time) residuals(num_weights_`j')
    qui replace num_weights_`j' = . if D.treat != 1
    qui egen den_weights_`j' = total(num_weights_`j')
    qui gen weight_`j' = num_weights_`j' / den_weights_`j'
    qui bysort group_h`j': egen gweight_`j' = max(weight_`j')
    qui replace weight_`j' = gweight_`j' if weight_`j' == .
    qui replace weight_`j' = round(weight_`j', 0.00000001)
    qui gen reweight_`j' = 1 / weight_`j'
    qui sort unit time
}

***************************************************
* REWEIGHTED LP-DiD (with interaction)
***************************************************
gen b_rw_lpdid      = .
gen se_rw_lpdid     = .
gen N_rw_lpdid      = .
gen b_rw_lpdid_low  = .
gen se_rw_lpdid_low = .
gen b_rw_lpdid_high = .
gen se_rw_lpdid_high= .
gen b_rw_lpdid_diff = .
gen se_rw_lpdid_diff= .

forval j = 0/`post_window' {
    reghdfe D`j'y c.D.treat##i.high_rol_pre ///
        if (D.treat==1 | F`j'.treat==0) ///
        [pweight = reweight_`j'], ///
        absorb(time) vce(cluster unit)

    replace b_rw_lpdid      = _b[D.treat]                    if t == `j'
    replace se_rw_lpdid     = _se[D.treat]                   if t == `j'
    replace N_rw_lpdid      = e(N)                            if t == `j'
    replace b_rw_lpdid_low  = _b[D.treat]                    if t == `j'
    replace se_rw_lpdid_low = _se[D.treat]                   if t == `j'
    replace b_rw_lpdid_diff = _b[1.high_rol_pre#c.D.treat]   if t == `j'
    replace se_rw_lpdid_diff= _se[1.high_rol_pre#c.D.treat]  if t == `j'
    lincom D.treat + 1.high_rol_pre#c.D.treat
    replace b_rw_lpdid_high  = r(estimate) if t == `j'
    replace se_rw_lpdid_high = r(se)       if t == `j'

    if `j' > 1 & `j' <= `pre_window' {
        reghdfe Dm`j'y c.D.treat##i.high_rol_pre ///
            if (D.treat==1 | treat==0) ///
            [pweight = reweight_0], ///
            absorb(time) vce(cluster unit)

        replace b_rw_lpdid      = _b[D.treat]                    if t == -`j'
        replace se_rw_lpdid     = _se[D.treat]                   if t == -`j'
        replace N_rw_lpdid      = e(N)                            if t == -`j'
        replace b_rw_lpdid_low  = _b[D.treat]                    if t == -`j'
        replace se_rw_lpdid_low = _se[D.treat]                   if t == -`j'
        replace b_rw_lpdid_diff = _b[1.high_rol_pre#c.D.treat]   if t == -`j'
        replace se_rw_lpdid_diff= _se[1.high_rol_pre#c.D.treat]  if t == -`j'
        lincom D.treat + 1.high_rol_pre#c.D.treat
        replace b_rw_lpdid_high  = r(estimate) if t == -`j'
        replace se_rw_lpdid_high = r(se)       if t == -`j'
    }
}

gen t_rw_lpdid_diff = b_rw_lpdid_diff / se_rw_lpdid_diff
gen p_rw_lpdid_diff = 2*ttail(N_rw_lpdid, abs(t_rw_lpdid_diff))
replace b_rw_lpdid_diff  = 0 if t == -1
replace se_rw_lpdid_diff = 0 if t == -1

gen tstat_rw_lpdid = b_rw_lpdid / se_rw_lpdid
gen pval_rw_lpdid  = 2*ttail(N_rw_lpdid, abs(tstat_rw_lpdid))

***************************************************
* POOLED LP-DiD + TABLE
***************************************************
xtset unit time

gen LPDID_pool = 0
forval j = 0/`post_window' {
    replace LPDID_pool = LPDID_pool + D`j'y
}
replace LPDID_pool = LPDID_pool / (`post_window' + 1)

gen pooled_sample = (dtreat==1 | F`post_window'.treat==0)

* Reweighted (equally-weighted ATT)
reghdfe LPDID_pool c.dtreat##i.high_rol_pre ///
    if (dtreat==1 | F`post_window'.treat==0) ///
    [pweight = reweight_`post_window'], ///
    absorb(time) vce(cluster unit)

scalar tbl_N_rw   = e(N)
scalar tbl_r2w_rw = e(r2_within)
scalar tbl_df_rw  = e(df_r)
scalar b_rw_low   = _b[dtreat]
scalar se_rw_low  = _se[dtreat]
scalar t_rw_low   = _b[dtreat] / _se[dtreat]
scalar p_rw_low   = 2*ttail(e(df_r), abs(t_rw_low))
scalar df_rw      = e(df_r)

scalar b_rw_diff  = _b[1.high_rol_pre#c.dtreat]
scalar se_rw_diff = _se[1.high_rol_pre#c.dtreat]
scalar t_rw_diff  = b_rw_diff / se_rw_diff
scalar p_rw_diff  = 2*ttail(df_rw, abs(t_rw_diff))

lincom dtreat + 1.high_rol_pre#c.dtreat
scalar b_rw_high  = r(estimate)
scalar se_rw_high = r(se)
scalar t_rw_high  = b_rw_high / se_rw_high
scalar p_rw_high  = 2*ttail(df_rw, abs(t_rw_high))

* Variance-weighted ATT
reghdfe LPDID_pool c.D.treat##i.high_rol_pre ///
    if pooled_sample, ///
    absorb(time) vce(cluster unit)

scalar tbl_N_vw   = e(N)
scalar tbl_r2w_vw = e(r2_within)
scalar tbl_df_vw  = e(df_r)
scalar b_vw_low   = _b[D.treat]
scalar se_vw_low  = _se[D.treat]
scalar t_vw_low   = _b[D.treat] / _se[D.treat]
scalar p_vw_low   = 2*ttail(e(df_r), abs(t_vw_low))
scalar df_vw      = e(df_r)

scalar b_vw_diff  = _b[1.high_rol_pre#c.D.treat]
scalar se_vw_diff = _se[1.high_rol_pre#c.D.treat]
scalar t_vw_diff  = b_vw_diff / se_vw_diff
scalar p_vw_diff  = 2*ttail(df_vw, abs(t_vw_diff))

lincom D.treat + 1.high_rol_pre#c.D.treat
scalar b_vw_high  = r(estimate)
scalar se_vw_high = r(se)
scalar t_vw_high  = b_vw_high / se_vw_high
scalar p_vw_high  = 2*ttail(df_vw, abs(t_vw_high))

capture program drop get_stars
program define get_stars, rclass
    args pval
    if      `pval' < 0.01 return local stars "***"
    else if `pval' < 0.05 return local stars "**"
    else if `pval' < 0.10 return local stars "*"
    else                  return local stars ""
end

get_stars `=scalar(p_rw_low)'
local cell_rw_low  = "\makecell{" + string(scalar(b_rw_low),"%6.3f")  + "`r(stars)'" + "\\" + "(" + string(scalar(t_rw_low),"%5.2f")  + ")}"
get_stars `=scalar(p_rw_high)'
local cell_rw_high = "\makecell{" + string(scalar(b_rw_high),"%6.3f") + "`r(stars)'" + "\\" + "(" + string(scalar(t_rw_high),"%5.2f") + ")}"
get_stars `=scalar(p_rw_diff)'
local cell_rw_diff = "\makecell{" + string(scalar(b_rw_diff),"%6.3f") + "`r(stars)'" + "\\" + "(" + string(scalar(t_rw_diff),"%5.2f") + ")}"
get_stars `=scalar(p_vw_low)'
local cell_vw_low  = "\makecell{" + string(scalar(b_vw_low),"%6.3f")  + "`r(stars)'" + "\\" + "(" + string(scalar(t_vw_low),"%5.2f")  + ")}"
get_stars `=scalar(p_vw_high)'
local cell_vw_high = "\makecell{" + string(scalar(b_vw_high),"%6.3f") + "`r(stars)'" + "\\" + "(" + string(scalar(t_vw_high),"%5.2f") + ")}"
get_stars `=scalar(p_vw_diff)'
local cell_vw_diff = "\makecell{" + string(scalar(b_vw_diff),"%6.3f") + "`r(stars)'" + "\\" + "(" + string(scalar(t_vw_diff),"%5.2f") + ")}"

local tbl_N_rw   = string(scalar(tbl_N_rw))
local tbl_N_vw   = string(scalar(tbl_N_vw))
local tbl_r2w_rw = string(scalar(tbl_r2w_rw), "%6.3f")
local tbl_r2w_vw = string(scalar(tbl_r2w_vw), "%6.3f")

file open tex using pooled_lpdid_rol_prosecution_high_rol_pre.tex, write replace
file write tex "\begin{table}[htbp]" _n
file write tex "\centering" _n
file write tex "\caption{Pooled LP-DiD Estimates: Effect on Prosecution by Rule of Law}" _n
file write tex "\label{tab:pooled_lpdid_rol}" _n
file write tex "\begin{tabular}{lccc}" _n
file write tex "\toprule" _n
file write tex " & Low RoL & High RoL & High \$-\$ Low \\" _n
file write tex "\midrule" _n
file write tex "Reweighted LP-DiD        & `cell_rw_low' & `cell_rw_high' & `cell_rw_diff' \\" _n
file write tex "[6pt]" _n
file write tex "Variance-weighted LP-DiD & `cell_vw_low' & `cell_vw_high' & `cell_vw_diff' \\" _n
file write tex "\midrule" _n
file write tex "Observations (RW)   & \multicolumn{3}{c}{`tbl_N_rw'}   \\" _n
file write tex "Observations (VW)   & \multicolumn{3}{c}{`tbl_N_vw'}   \\" _n
file write tex "Within R\$^2\$ (RW) & \multicolumn{3}{c}{`tbl_r2w_rw'} \\" _n
file write tex "Within R\$^2\$ (VW) & \multicolumn{3}{c}{`tbl_r2w_vw'} \\" _n
file write tex "\bottomrule" _n
file write tex "\end{tabular}" _n
file write tex "\begin{tablenotes}" _n
file write tex "\small" _n
file write tex "\item \textit{Notes:} Pooled LP-DiD estimates. " _n
file write tex "Outcome: average of post-treatment long differences (h=0 to h=`post_window'). " _n
file write tex "High RoL = above-median Rule of Law (pre-treatment). " _n
file write tex "t-statistics in parentheses, clustered at country level. " _n
file write tex "* p\$<\$0.10, ** p\$<\$0.05, *** p\$<\$0.01." _n
file write tex "\end{tablenotes}" _n
file write tex "\end{table}" _n
file close tex
di "Saved: pooled_lpdid_rol_prosecution_high_rol_pre.tex"

***************************************************
* PMD LP-DiD (with interaction)
***************************************************
gen b_pmd_lpdid  = .
gen se_pmd_lpdid = .
gen N_pmd_lpdid  = .
gen b_pmd_low    = .
gen se_pmd_low   = .
gen b_pmd_high   = .
gen se_pmd_high  = .
gen b_pmd_diff   = .
gen se_pmd_diff  = .

xtset unit time

gen Y_pre = Y if time < treat_year
bysort unit (time): gen cumsum_pre = sum(Y_pre)
bysort unit (time): gen n_pre      = sum(Y_pre < .)
xtset unit time
gen aveLY = L.cumsum_pre / L.n_pre

forval j = 0/`post_window' {
    gen PMD`j'y = F`j'.Y - aveLY
}
forval j = 2/`pre_window' {
    gen PMDm`j'y = L`j'.Y - aveLY
}

forval j = 0/`post_window' {
    reghdfe PMD`j'y c.D.treat##i.high_rol_pre ///
        if D.treat==1 | F`j'.treat==0, ///
        absorb(time) vce(cluster unit)

    replace b_pmd_lpdid  = _b[D.treat]                   if t == `j'
    replace se_pmd_lpdid = _se[D.treat]                  if t == `j'
    replace N_pmd_lpdid  = e(N)                           if t == `j'
    replace b_pmd_low    = _b[D.treat]                   if t == `j'
    replace se_pmd_low   = _se[D.treat]                  if t == `j'
    replace b_pmd_diff   = _b[1.high_rol_pre#c.D.treat]  if t == `j'
    replace se_pmd_diff  = _se[1.high_rol_pre#c.D.treat] if t == `j'
    lincom D.treat + 1.high_rol_pre#c.D.treat
    replace b_pmd_high   = r(estimate) if t == `j'
    replace se_pmd_high  = r(se)       if t == `j'

    if `j' > 1 & `j' <= `pre_window' {
        reghdfe PMDm`j'y c.D.treat##i.high_rol_pre ///
            if D.treat==1 | treat==0, ///
            absorb(time) vce(cluster unit)

        replace b_pmd_lpdid  = _b[D.treat]                   if t == -`j'
        replace se_pmd_lpdid = _se[D.treat]                  if t == -`j'
        replace N_pmd_lpdid  = e(N)                           if t == -`j'
        replace b_pmd_low    = _b[D.treat]                   if t == -`j'
        replace se_pmd_low   = _se[D.treat]                  if t == -`j'
        replace b_pmd_diff   = _b[1.high_rol_pre#c.D.treat]  if t == -`j'
        replace se_pmd_diff  = _se[1.high_rol_pre#c.D.treat] if t == -`j'
        lincom D.treat + 1.high_rol_pre#c.D.treat
        replace b_pmd_high   = r(estimate) if t == -`j'
        replace se_pmd_high  = r(se)       if t == -`j'
    }
}

gen t_pmd_diff = b_pmd_diff / se_pmd_diff
gen p_pmd_diff = 2*ttail(N_pmd_lpdid, abs(t_pmd_diff))

replace b_pmd_lpdid  = 0 if t == -1
replace se_pmd_lpdid = 0 if t == -1
replace b_pmd_diff   = 0 if t == -1
replace se_pmd_diff  = 0 if t == -1

***************************************************
* REWEIGHTED PMD LP-DiD (with interaction)
***************************************************
gen b_rwpmd_lpdid  = .
gen se_rwpmd_lpdid = .
gen N_rwpmd_lpdid  = .
gen N_high         = .
gen N_low          = .
gen Nunits_high    = .
gen Nunits_low     = .
gen b_rwpmd_low    = .
gen se_rwpmd_low   = .
gen b_rwpmd_high   = .
gen se_rwpmd_high  = .
gen b_rwpmd_diff   = .
gen se_rwpmd_diff  = .

forval j = 0/`post_window' {

    quietly count if (D.treat==1 | F`j'.treat==0) & high_rol_pre==1
    replace N_high = r(N) if t == `j'
    quietly count if (D.treat==1 | F`j'.treat==0) & high_rol_pre==0
    replace N_low  = r(N) if t == `j'

    tempvar tagh tagl
    quietly egen `tagh' = tag(unit) if (D.treat==1 | F`j'.treat==0) & high_rol_pre==1
    quietly summarize `tagh'
    replace Nunits_high = r(sum) if t == `j'
    quietly egen `tagl' = tag(unit) if (D.treat==1 | F`j'.treat==0) & high_rol_pre==0
    quietly summarize `tagl'
    replace Nunits_low  = r(sum) if t == `j'

    reghdfe PMD`j'y c.D.treat##i.high_rol_pre ///
        if (D.treat==1 | F`j'.treat==0) ///
        [pweight = reweight_`j'], ///
        absorb(time) vce(cluster unit)

    replace N_rwpmd_lpdid  = e(N)                            if t == `j'
    replace b_rwpmd_lpdid  = _b[D.treat]                    if t == `j'
    replace se_rwpmd_lpdid = _se[D.treat]                   if t == `j'
    replace b_rwpmd_low    = _b[D.treat]                    if t == `j'
    replace se_rwpmd_low   = _se[D.treat]                   if t == `j'
    replace b_rwpmd_diff   = _b[1.high_rol_pre#c.D.treat]   if t == `j'
    replace se_rwpmd_diff  = _se[1.high_rol_pre#c.D.treat]  if t == `j'
    lincom D.treat + 1.high_rol_pre#c.D.treat
    replace b_rwpmd_high   = r(estimate) if t == `j'
    replace se_rwpmd_high  = r(se)       if t == `j'

    if `j' > 1 & `j' <= `pre_window' {

        quietly count if (D.treat==1 | F`j'.treat==0) & high_rol_pre==1
        replace N_high = r(N) if t == -`j'
        quietly count if (D.treat==1 | F`j'.treat==0) & high_rol_pre==0
        replace N_low  = r(N) if t == -`j'

        tempvar tagh tagl
        quietly egen `tagh' = tag(unit) if (D.treat==1 | F`j'.treat==0) & high_rol_pre==1
        quietly summarize `tagh'
        replace Nunits_high = r(sum) if t == -`j'
        quietly egen `tagl' = tag(unit) if (D.treat==1 | F`j'.treat==0) & high_rol_pre==0
        quietly summarize `tagl'
        replace Nunits_low  = r(sum) if t == -`j'

        reghdfe PMDm`j'y c.D.treat##i.high_rol_pre ///
            if (D.treat==1 | F`j'.treat==0) ///
            [pweight = reweight_0], ///
            absorb(time) vce(cluster unit)

        replace N_rwpmd_lpdid  = e(N)                            if t == -`j'
        replace b_rwpmd_lpdid  = _b[D.treat]                    if t == -`j'
        replace se_rwpmd_lpdid = _se[D.treat]                   if t == -`j'
        replace b_rwpmd_low    = _b[D.treat]                    if t == -`j'
        replace se_rwpmd_low   = _se[D.treat]                   if t == -`j'
        replace b_rwpmd_diff   = _b[1.high_rol_pre#c.D.treat]   if t == -`j'
        replace se_rwpmd_diff  = _se[1.high_rol_pre#c.D.treat]  if t == -`j'
        lincom D.treat + 1.high_rol_pre#c.D.treat
        replace b_rwpmd_high   = r(estimate) if t == -`j'
        replace se_rwpmd_high  = r(se)       if t == -`j'
    }
}

replace b_rwpmd_diff  = 0 if t == -1

gen t_rwpmd_diff = b_rwpmd_diff / se_rwpmd_diff
gen p_rwpmd_diff = 2*ttail(N_rwpmd_lpdid, abs(t_rwpmd_diff))

replace se_rwpmd_diff = 0 if t == -1

***************************************************
* POOLED PMD LP-DiD + TABLE
***************************************************
xtset unit time

gen PMD_pool = 0
forval j = 0/`post_window' {
    replace PMD_pool = PMD_pool + PMD`j'y
}
replace PMD_pool = PMD_pool / (`post_window' + 1)

* Reweighted (equally-weighted ATT)
reghdfe PMD_pool c.dtreat##i.high_rol_pre ///
    if (dtreat==1 | F`post_window'.treat==0) ///
    [pweight = reweight_`post_window'], ///
    absorb(time) vce(cluster unit)

scalar tbl_N_pmd_rw   = e(N)
scalar tbl_r2w_pmd_rw = e(r2_within)
scalar tbl_df_pmd_rw  = e(df_r)
scalar b_rw_low       = _b[dtreat]
scalar se_rw_low      = _se[dtreat]
scalar t_rw_low       = _b[dtreat] / _se[dtreat]
scalar p_rw_low       = 2*ttail(e(df_r), abs(t_rw_low))
scalar df_rw          = e(df_r)

scalar b_rw_diff  = _b[1.high_rol_pre#c.dtreat]
scalar se_rw_diff = _se[1.high_rol_pre#c.dtreat]
scalar t_rw_diff  = b_rw_diff / se_rw_diff
scalar p_rw_diff  = 2*ttail(df_rw, abs(t_rw_diff))

lincom dtreat + 1.high_rol_pre#c.dtreat
scalar b_rw_high  = r(estimate)
scalar se_rw_high = r(se)
scalar t_rw_high  = b_rw_high / se_rw_high
scalar p_rw_high  = 2*ttail(df_rw, abs(t_rw_high))

* Variance-weighted ATT
reghdfe PMD_pool c.D.treat##i.high_rol_pre ///
    if pooled_sample, ///
    absorb(time) vce(cluster unit)

scalar tbl_N_pmd_vw   = e(N)
scalar tbl_r2w_pmd_vw = e(r2_within)
scalar tbl_df_pmd_vw  = e(df_r)
scalar b_vw_low       = _b[D.treat]
scalar se_vw_low      = _se[D.treat]
scalar t_vw_low       = _b[D.treat] / _se[D.treat]
scalar p_vw_low       = 2*ttail(e(df_r), abs(t_vw_low))
scalar df_vw          = e(df_r)

scalar b_vw_diff  = _b[1.high_rol_pre#c.D.treat]
scalar se_vw_diff = _se[1.high_rol_pre#c.D.treat]
scalar t_vw_diff  = b_vw_diff / se_vw_diff
scalar p_vw_diff  = 2*ttail(df_vw, abs(t_vw_diff))

lincom D.treat + 1.high_rol_pre#c.D.treat
scalar b_vw_high  = r(estimate)
scalar se_vw_high = r(se)
scalar t_vw_high  = b_vw_high / se_vw_high
scalar p_vw_high  = 2*ttail(df_vw, abs(t_vw_high))

get_stars `=scalar(p_rw_low)'
local cell_rw_low  = "\makecell{" + string(scalar(b_rw_low),"%6.3f")  + "`r(stars)'" + "\\" + "(" + string(scalar(t_rw_low),"%5.2f")  + ")}"
get_stars `=scalar(p_rw_high)'
local cell_rw_high = "\makecell{" + string(scalar(b_rw_high),"%6.3f") + "`r(stars)'" + "\\" + "(" + string(scalar(t_rw_high),"%5.2f") + ")}"
get_stars `=scalar(p_rw_diff)'
local cell_rw_diff = "\makecell{" + string(scalar(b_rw_diff),"%6.3f") + "`r(stars)'" + "\\" + "(" + string(scalar(t_rw_diff),"%5.2f") + ")}"
get_stars `=scalar(p_vw_low)'
local cell_vw_low  = "\makecell{" + string(scalar(b_vw_low),"%6.3f")  + "`r(stars)'" + "\\" + "(" + string(scalar(t_vw_low),"%5.2f")  + ")}"
get_stars `=scalar(p_vw_high)'
local cell_vw_high = "\makecell{" + string(scalar(b_vw_high),"%6.3f") + "`r(stars)'" + "\\" + "(" + string(scalar(t_vw_high),"%5.2f") + ")}"
get_stars `=scalar(p_vw_diff)'
local cell_vw_diff = "\makecell{" + string(scalar(b_vw_diff),"%6.3f") + "`r(stars)'" + "\\" + "(" + string(scalar(t_vw_diff),"%5.2f") + ")}"

local tbl_N_pmd_rw   = string(scalar(tbl_N_pmd_rw))
local tbl_N_pmd_vw   = string(scalar(tbl_N_pmd_vw))
local tbl_r2w_pmd_rw = string(scalar(tbl_r2w_pmd_rw), "%6.3f")
local tbl_r2w_pmd_vw = string(scalar(tbl_r2w_pmd_vw), "%6.3f")

file open tex using pooled_pmd_lpdid_rol_prosecution_high_rol_pre.tex, write replace
file write tex "\begin{table}[htbp]" _n
file write tex "\centering" _n
file write tex "\caption{Pooled PMD LP-DiD Estimates: Effect on Prosecution by Rule of Law}" _n
file write tex "\label{tab:pooled_pmd_lpdid_rol}" _n
file write tex "\begin{tabular}{lccc}" _n
file write tex "\toprule" _n
file write tex " & Low RoL & High RoL & High \$-\$ Low \\" _n
file write tex "\midrule" _n
file write tex "Reweighted PMD LP-DiD        & `cell_rw_low' & `cell_rw_high' & `cell_rw_diff' \\" _n
file write tex "[6pt]" _n
file write tex "Variance-weighted PMD LP-DiD & `cell_vw_low' & `cell_vw_high' & `cell_vw_diff' \\" _n
file write tex "\midrule" _n
file write tex "Observations (RW)   & \multicolumn{3}{c}{`tbl_N_pmd_rw'}   \\" _n
file write tex "Observations (VW)   & \multicolumn{3}{c}{`tbl_N_pmd_vw'}   \\" _n
file write tex "Within R\$^2\$ (RW) & \multicolumn{3}{c}{`tbl_r2w_pmd_rw'} \\" _n
file write tex "Within R\$^2\$ (VW) & \multicolumn{3}{c}{`tbl_r2w_pmd_vw'} \\" _n
file write tex "\bottomrule" _n
file write tex "\end{tabular}" _n
file write tex "\begin{tablenotes}" _n
file write tex "\small" _n
file write tex "\item \textit{Notes:} Pooled PMD LP-DiD estimates. " _n
file write tex "Outcome: average of post-treatment PMD long differences (h=0 to h=`post_window'). " _n
file write tex "High RoL = above-median Rule of Law (pre-treatment). " _n
file write tex "t-statistics in parentheses, clustered at country level. " _n
file write tex "* p\$<\$0.10, ** p\$<\$0.05, *** p\$<\$0.01." _n
file write tex "\end{tablenotes}" _n
file write tex "\end{table}" _n
file close tex
di "Saved: pooled_pmd_lpdid_rol_prosecution_high_rol_pre.tex"

***************************************************
* COUNTS TABLE (N by RoL group and event time)
***************************************************
preserve
keep t N_high N_low Nunits_high Nunits_low
bysort t: keep if _n == 1
sort t

export excel using counts_table_gov_eff_diff_table_rwpmd_prosecution.xlsx, firstrow(variables) replace

table t, ///
    statistic(mean N_high) ///
    statistic(mean N_low)  ///
    statistic(mean Nunits_high) ///
    statistic(mean Nunits_low)

collect export counts_table_gov_eff_diff_table_rwpmd_prosecution.tex, replace
restore

***************************************************
* RWPMD DIFFERENTIAL TABLE (xlsx + tex)
***************************************************
gen stars     = cond(p_rwpmd_diff<0.01,"***", cond(p_rwpmd_diff<0.05,"**", cond(p_rwpmd_diff<0.10,"*","")))
gen coef_star = string(b_rwpmd_diff, "%9.3f") + stars
gen t_paren   = "(" + string(t_rwpmd_diff, "%9.2f") + ")"

preserve
keep t coef_star t_paren
bysort t: keep if _n == 1
sort t

export excel t coef_star t_paren using diff_table_rwpmd_prosecution.xlsx, ///
    firstrow(variables) replace

file open tex using diff_table_rwpmd_prosecution.tex, write replace
file write tex "\begin{tabular}{ccc}" _n
file write tex "Event time & Coef. & t-stat \\" _n
file write tex "\hline" _n
local N = _N
quietly forvalues i = 1/`N' {
    local et = t[`i']
    local b  = coef_star[`i']
    local tt = t_paren[`i']
    file write tex "`et' & `b' & `tt' \\" _n
}
file write tex "\end{tabular}" _n
file close tex
restore

***************************************************
* ESTIMATES GRAPH
***************************************************
preserve
bysort t: keep if _n == 1
sort t
keep if t >= -`pre_window' & t <= `post_window'

tw (line b_twfe_diff      t, lcolor(red)    lpattern(dash)      lwidth(medium)) ///
   (line b_lpdid_diff     t, lcolor(green)  lpattern(longdash)  lwidth(medium)) ///
   (line b_rw_lpdid_diff  t, lcolor(orange) lpattern(solid)     lwidth(thin))   ///
   (line b_pmd_diff       t, lcolor(purple) lpattern(longdash)  lwidth(thin))   ///
   (line b_rwpmd_diff     t, lcolor(blue)   lpattern(longdash)  lwidth(thin)),  ///
    legend(order(1 "TWFE Distributed Lags" 2 "LP-DiD" 3 "Reweighted LP-DiD" 4 "PMD LP-DiD" 5 "Reweighted PMD LP-DiD") pos(6) col(2)) ///
    xtitle("Event Time") ytitle("Differential Treatment Effect (High - Low RoL)") ///
    note("Units: `units', Periods: `time_periods'" "Obs range: `Nmin'–`Nmax', Countries: `min_units'–`max_units'", ring(0) pos(10) size(*0.8) box) ///
    title("Heterogeneous Treatment by Rule of Law on Prosecution") ///
    name(effects_6, replace)

restore

***************************************************
* EVENT-STUDY DIFFERENTIAL TABLE
***************************************************
preserve
bysort t: keep if _n == 1
sort t
keep if t >= -`pre_window' & t <= `post_window'

gen star_twfe  = cond(p_twfe_diff<0.01,"***",    cond(p_twfe_diff<0.05,"**",    cond(p_twfe_diff<0.10,"*","")))
gen star_lpdid = cond(p_lpdid_diff<0.01,"***",   cond(p_lpdid_diff<0.05,"**",   cond(p_lpdid_diff<0.10,"*","")))
gen star_rw    = cond(p_rw_lpdid_diff<0.01,"***",cond(p_rw_lpdid_diff<0.05,"**",cond(p_rw_lpdid_diff<0.10,"*","")))
gen star_pmd   = cond(p_pmd_diff<0.01,"***",     cond(p_pmd_diff<0.05,"**",     cond(p_pmd_diff<0.10,"*","")))
gen star_rwpmd = cond(p_rwpmd_diff<0.01,"***",   cond(p_rwpmd_diff<0.05,"**",   cond(p_rwpmd_diff<0.10,"*","")))

gen cell_twfe  = string(b_twfe_diff,"%9.3f")     + star_twfe  + " (" + string(t_twfe_diff,"%9.2f")     + ")"
gen cell_lpdid = string(b_lpdid_diff,"%9.3f")    + star_lpdid + " (" + string(t_lpdid_diff,"%9.2f")    + ")"
gen cell_rw    = string(b_rw_lpdid_diff,"%9.3f") + star_rw    + " (" + string(t_rw_lpdid_diff,"%9.2f") + ")"
gen cell_pmd   = string(b_pmd_diff,"%9.3f")      + star_pmd   + " (" + string(t_pmd_diff,"%9.2f")      + ")"
gen cell_rwpmd = string(b_rwpmd_diff,"%9.3f")    + star_rwpmd + " (" + string(t_rwpmd_diff,"%9.2f")    + ")"

keep t cell_twfe cell_lpdid cell_rw cell_pmd cell_rwpmd
sort t

export excel using lp_did_diff_methods_interaction_high_rol_pre_diff_prosecution.xlsx, firstrow(variables) replace

outfile t cell_twfe cell_lpdid cell_rw cell_pmd cell_rwpmd ///
    using lp_did_diff_methods_interaction_high_rol_pre.txt, replace wide

putdocx begin
putdocx paragraph, halign(center)
putdocx text ("Event Study Results — Differential by Rule of Law"), bold

count
local n = r(N)
local rows = `n' + 1

putdocx table tbl = (`rows', 6)
putdocx table tbl(1,1) = ("t")
putdocx table tbl(1,2) = ("TWFE")
putdocx table tbl(1,3) = ("LP-DiD")
putdocx table tbl(1,4) = ("RW LP-DiD")
putdocx table tbl(1,5) = ("PMD")
putdocx table tbl(1,6) = ("RWPMD")

forvalues i = 1/`n' {
    local r = `i' + 1
    putdocx table tbl(`r',1) = (t[`i'])
    putdocx table tbl(`r',2) = (cell_twfe[`i'])
    putdocx table tbl(`r',3) = (cell_lpdid[`i'])
    putdocx table tbl(`r',4) = (cell_rw[`i'])
    putdocx table tbl(`r',5) = (cell_pmd[`i'])
    putdocx table tbl(`r',6) = (cell_rwpmd[`i'])
}
putdocx save "lp_did_diff_methods_interaction_high_rol_pre_diff_prosecution.docx", replace

listtex using high_rol_pre_prosecution.tex, replace ///
    rstyle(tabular) ///
    head("\begin{tabular}{lccccc}" ///
         "\toprule" ///
         "t & TWFE & LP-DiD & RW LP-DiD & PMD & RWPMD \\" ///
         "\midrule") ///
    foot("\bottomrule" "\end{tabular}")

restore
