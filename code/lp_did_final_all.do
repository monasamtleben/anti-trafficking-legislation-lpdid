/*
    LP-DiD: All estimators — prosecution outcome
    Estimators: TWFE, LP-DiD, Reweighted LP-DiD, PMD LP-DiD, Reweighted PMD LP-DiD
    Outputs: pooled tables (LaTeX), event-study graph, event-study table (xlsx/tex/docx)
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
* TWFE EVENT-STUDY ESTIMATES
***************************************************
gen b_twfe  = .
gen se_twfe = .
gen N_twfe  = .

xtset unit time
reghdfe Y L(-5/5).D.treat, absorb(unit time) cluster(unit)

* Pre-periods (normalised to t=-1)
forval j = 2/`pre_window' {
    lincom F`j'D.treat - FD.treat
    replace b_twfe  = r(estimate) if t == -`j'
    replace se_twfe = r(se)       if t == -`j'
    replace N_twfe  = e(N)        if t == -`j'
}

replace b_twfe  = 0 if t == -1
replace se_twfe = 0 if t == -1

* t=0
lincom D.treat - FD.treat
replace b_twfe  = r(estimate) if t == 0
replace se_twfe = r(se)       if t == 0
replace N_twfe  = e(N)        if t == 0

* t=1
lincom LD.treat - FD.treat
replace b_twfe  = r(estimate) if t == 1
replace se_twfe = r(se)       if t == 1
replace N_twfe  = e(N)        if t == 1

* Post-periods
forval j = 2/`post_window' {
    lincom L`j'D.treat - FD.treat
    replace b_twfe  = r(estimate) if t == `j'
    replace se_twfe = r(se)       if t == `j'
    replace N_twfe  = e(N)        if t == `j'
}

gen tstat_twfe = b_twfe / se_twfe
gen pval_twfe  = 2*ttail(N_twfe, abs(tstat_twfe))

***************************************************
* LP-DiD ESTIMATES
***************************************************
gen b_lpdid  = .
gen se_lpdid = .
gen df_lpdid = .
gen N_lpdid  = .
gen N_units  = .

* Long differences
forval j = 0/`post_window' {
    gen D`j'y = F`j'.Y - L.Y
}
forval j = 2/`pre_window' {
    gen Dm`j'y = L`j'.Y - L.Y
}

forval j = 0/`post_window' {

    reghdfe D`j'y D.treat ///
        if D.treat==1 | F`j'.treat==0, ///
        absorb(time) vce(cluster unit)

    preserve
        keep if e(sample)
        levelsof unit, local(u)
        local nunits : word count `u'
    restore

    replace b_lpdid  = _b[D.treat] if t == `j'
    replace se_lpdid = _se[D.treat] if t == `j'
    replace df_lpdid = e(df_r)      if t == `j'
    replace N_lpdid  = e(N)         if t == `j'
    replace N_units  = `nunits'     if t == `j'

    if `j' > 1 & `j' <= `pre_window' {
        reghdfe Dm`j'y D.treat ///
            if D.treat==1 | treat==0, ///
            absorb(time) vce(cluster unit)

        replace b_lpdid  = _b[D.treat] if t == -`j'
        replace se_lpdid = _se[D.treat] if t == -`j'
        replace df_lpdid = e(df_r)      if t == -`j'
        replace N_lpdid  = e(N)         if t == -`j'

        preserve
            keep if e(sample)
            levelsof unit, local(u)
            local nunits : word count `u'
        restore
        replace N_units = `nunits' if t == -`j'
    }
}

replace b_lpdid  = 0 if t == -1
replace se_lpdid = 0 if t == -1

gen tstat_lpdid = b_lpdid / se_lpdid
gen pval_lpdid  = 2*ttail(df_lpdid, abs(tstat_lpdid))

summ N_units
local min_units = r(min)
local max_units = r(max)

summ N_lpdid
local Nmin = r(min)
local Nmax = r(max)

***************************************************
* REWEIGHTED LP-DiD
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

gen b_rw_lpdid  = .
gen se_rw_lpdid = .
gen df_rw_lpdid = .

forval j = 0/`post_window' {
    reghdfe D`j'y D.treat ///
        if (D.treat==1 | F`j'.treat==0) ///
        [pweight = reweight_`j'], ///
        absorb(time) vce(cluster unit)

    replace b_rw_lpdid  = _b[D.treat] if t == `j'
    replace se_rw_lpdid = _se[D.treat] if t == `j'
    replace df_rw_lpdid = e(df_r)      if t == `j'

    if `j' > 1 & `j' <= `pre_window' {
        reghdfe Dm`j'y D.treat ///
            if (D.treat==1 | treat==0) ///
            [pweight = reweight_0], ///
            absorb(time) vce(cluster unit)

        replace b_rw_lpdid  = _b[D.treat] if t == -`j'
        replace se_rw_lpdid = _se[D.treat] if t == -`j'
        replace df_rw_lpdid = e(df_r)      if t == -`j'
    }
}

replace b_rw_lpdid  = 0 if t == -1
replace se_rw_lpdid = 0 if t == -1

gen tstat_rw_lpdid = b_rw_lpdid / se_rw_lpdid
gen pval_rw_lpdid  = 2*ttail(df_rw_lpdid, abs(tstat_rw_lpdid))

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
reghdfe LPDID_pool D.treat ///
    if (D.treat==1 | F`post_window'.treat==0) ///
    [pweight = reweight_`post_window'], ///
    absorb(time) vce(cluster unit)

scalar pool_lpdid_N_rw   = e(N)
scalar pool_lpdid_r2w_rw = e(r2_within)
scalar pool_lpdid_df_rw  = e(df_r)
scalar pool_lpdid_b_rw   = _b[D.treat]
scalar pool_lpdid_se_rw  = _se[D.treat]
scalar pool_lpdid_t_rw   = _b[D.treat] / _se[D.treat]
scalar pool_lpdid_p_rw   = 2 * ttail(e(df_r), abs(pool_lpdid_t_rw))

* Variance-weighted ATT
reghdfe LPDID_pool D.treat ///
    if pooled_sample, ///
    absorb(time) vce(cluster unit)

scalar pool_lpdid_N_vw   = e(N)
scalar pool_lpdid_r2w_vw = e(r2_within)
scalar pool_lpdid_df_vw  = e(df_r)
scalar pool_lpdid_b_vw   = _b[D.treat]
scalar pool_lpdid_se_vw  = _se[D.treat]
scalar pool_lpdid_t_vw   = _b[D.treat] / _se[D.treat]
scalar pool_lpdid_p_vw   = 2 * ttail(e(df_r), abs(pool_lpdid_t_vw))

capture program drop get_stars
program define get_stars, rclass
    args pval
    if      `pval' < 0.01 return local stars "***"
    else if `pval' < 0.05 return local stars "**"
    else if `pval' < 0.10 return local stars "*"
    else                  return local stars ""
end

get_stars `=scalar(pool_lpdid_p_rw)'
local cell_pool_lpdid_rw = "\makecell{" + string(scalar(pool_lpdid_b_rw), "%6.3f") + "`r(stars)'" + "\\" + "(" + string(scalar(pool_lpdid_t_rw), "%5.2f") + ")}"

get_stars `=scalar(pool_lpdid_p_vw)'
local cell_pool_lpdid_vw = "\makecell{" + string(scalar(pool_lpdid_b_vw), "%6.3f") + "`r(stars)'" + "\\" + "(" + string(scalar(pool_lpdid_t_vw), "%5.2f") + ")}"

local pool_lpdid_N_rw   = string(scalar(pool_lpdid_N_rw))
local pool_lpdid_N_vw   = string(scalar(pool_lpdid_N_vw))
local pool_lpdid_r2w_rw = string(scalar(pool_lpdid_r2w_rw), "%6.3f")
local pool_lpdid_r2w_vw = string(scalar(pool_lpdid_r2w_vw), "%6.3f")

file open tex using pooled_lpdid_overall3p.tex, write replace
file write tex "\begin{table}[htbp]" _n
file write tex "\centering" _n
file write tex "\caption{Pooled LP-DiD Estimates: Effect on Overall 3P Index}" _n
file write tex "\label{tab:pooled_lpdid_overall3p}" _n
file write tex "\begin{tabular}{lc}" _n
file write tex "\toprule" _n
file write tex " & Overall ATT \\" _n
file write tex "\midrule" _n
file write tex "Reweighted LP-DiD         & `cell_pool_lpdid_rw' \\" _n
file write tex "[6pt]" _n
file write tex "Variance-weighted LP-DiD  & `cell_pool_lpdid_vw' \\" _n
file write tex "\midrule" _n
file write tex "Observations (RW)   & `pool_lpdid_N_rw'   \\" _n
file write tex "Observations (VW)   & `pool_lpdid_N_vw'   \\" _n
file write tex "Within R\$^2\$ (RW) & `pool_lpdid_r2w_rw' \\" _n
file write tex "Within R\$^2\$ (VW) & `pool_lpdid_r2w_vw' \\" _n
file write tex "\bottomrule" _n
file write tex "\end{tabular}" _n
file write tex "\begin{tablenotes}" _n
file write tex "\small" _n
file write tex "\item \textit{Notes:} Pooled LP-DiD estimates. " _n
file write tex "Outcome: average of post-treatment long differences in Overall 3P Index (h=0 to h=`post_window'). " _n
file write tex "t-statistics in parentheses, clustered at country level. " _n
file write tex "* p\$<\$0.10, ** p\$<\$0.05, *** p\$<\$0.01." _n
file write tex "\end{tablenotes}" _n
file write tex "\end{table}" _n
file close tex
di "Saved: pooled_lpdid_overall3p.tex"

***************************************************
* PMD LP-DiD
***************************************************
gen b_pmd_lpdid  = .
gen se_pmd_lpdid = .
gen df_pmd       = .

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
    reghdfe PMD`j'y D.treat ///
        if D.treat==1 | F`j'.treat==0, ///
        absorb(time) vce(cluster unit)

    replace b_pmd_lpdid  = _b[D.treat] if t == `j'
    replace se_pmd_lpdid = _se[D.treat] if t == `j'
    replace df_pmd       = e(df_r)      if t == `j'

    if `j' > 1 & `j' <= `pre_window' {
        reghdfe PMDm`j'y D.treat ///
            if (D.treat==1 | treat==0) ///
            [pweight = reweight_`j'], ///
            absorb(time) vce(cluster unit)

        replace b_pmd_lpdid  = _b[D.treat] if t == -`j'
        replace se_pmd_lpdid = _se[D.treat] if t == -`j'
        replace df_pmd       = e(df_r)      if t == -`j'
    }
}

replace b_pmd_lpdid  = 0 if t == -1
replace se_pmd_lpdid = 0 if t == -1

gen tstat_pmd_lpdid = b_pmd_lpdid / se_pmd_lpdid
gen pval_pmd_lpdid  = 2*ttail(df_pmd, abs(tstat_pmd_lpdid))

***************************************************
* REWEIGHTED PMD LP-DiD
***************************************************
gen b_rwpmd_lpdid  = .
gen se_rwpmd_lpdid = .
gen df_rwpmd       = .

forval j = 0/`post_window' {
    reghdfe PMD`j'y D.treat ///
        if (D.treat==1 | F`j'.treat==0) ///
        [pweight = reweight_`j'], ///
        absorb(time) vce(cluster unit)

    replace b_rwpmd_lpdid  = _b[D.treat] if t == `j'
    replace se_rwpmd_lpdid = _se[D.treat] if t == `j'
    replace df_rwpmd       = e(df_r)      if t == `j'

    if `j' > 1 & `j' <= `pre_window' {
        reghdfe PMDm`j'y D.treat ///
            if (D.treat==1 | F`j'.treat==0) ///
            [pweight = reweight_0], ///
            absorb(time) vce(cluster unit)

        replace b_rwpmd_lpdid  = _b[D.treat] if t == -`j'
        replace se_rwpmd_lpdid = _se[D.treat] if t == -`j'
        replace df_rwpmd       = e(df_r)      if t == -`j'
    }
}

replace b_rwpmd_lpdid  = 0 if t == -1
replace se_rwpmd_lpdid = 0 if t == -1

gen tstat_rwpmd_lpdid = b_rwpmd_lpdid / se_rwpmd_lpdid
gen pval_rwpmd_lpdid  = 2*ttail(df_rwpmd, abs(tstat_rwpmd_lpdid))

***************************************************
* POOLED PMD LP-DiD + TABLE
***************************************************
gen PMD_pool = 0
forval j = 0/`post_window' {
    replace PMD_pool = PMD_pool + PMD`j'y
}
replace PMD_pool = PMD_pool / (`post_window' + 1)

* Reweighted (equally-weighted ATT)
reghdfe PMD_pool c.D.treat ///
    if (D.treat==1 | F`post_window'.treat==0) ///
    [pweight = reweight_`post_window'], ///
    absorb(time) vce(cluster unit)

scalar pool_pmd_N_rw   = e(N)
scalar pool_pmd_r2w_rw = e(r2_within)
scalar pool_pmd_df_rw  = e(df_r)
scalar pool_pmd_b_rw   = _b[D.treat]
scalar pool_pmd_se_rw  = _se[D.treat]
scalar pool_pmd_t_rw   = _b[D.treat] / _se[D.treat]
scalar pool_pmd_p_rw   = 2 * ttail(e(df_r), abs(pool_pmd_t_rw))

* Variance-weighted ATT
reghdfe PMD_pool c.D.treat ///
    if pooled_sample, ///
    absorb(time) vce(cluster unit)

scalar pool_pmd_N_vw   = e(N)
scalar pool_pmd_r2w_vw = e(r2_within)
scalar pool_pmd_df_vw  = e(df_r)
scalar pool_pmd_b_vw   = _b[D.treat]
scalar pool_pmd_se_vw  = _se[D.treat]
scalar pool_pmd_t_vw   = _b[D.treat] / _se[D.treat]
scalar pool_pmd_p_vw   = 2 * ttail(e(df_r), abs(pool_pmd_t_vw))

get_stars `=scalar(pool_pmd_p_rw)'
local cell_pool_pmd_rw = "\makecell{" + string(scalar(pool_pmd_b_rw), "%6.3f") + "`r(stars)'" + "\\" + "(" + string(scalar(pool_pmd_t_rw), "%5.2f") + ")}"

get_stars `=scalar(pool_pmd_p_vw)'
local cell_pool_pmd_vw = "\makecell{" + string(scalar(pool_pmd_b_vw), "%6.3f") + "`r(stars)'" + "\\" + "(" + string(scalar(pool_pmd_t_vw), "%5.2f") + ")}"

local pool_pmd_N_rw   = string(scalar(pool_pmd_N_rw))
local pool_pmd_N_vw   = string(scalar(pool_pmd_N_vw))
local pool_pmd_r2w_rw = string(scalar(pool_pmd_r2w_rw), "%6.3f")
local pool_pmd_r2w_vw = string(scalar(pool_pmd_r2w_vw), "%6.3f")

file open tex using pooled_pmd_overall3p.tex, write replace
file write tex "\begin{table}[htbp]" _n
file write tex "\centering" _n
file write tex "\caption{Pooled PMD LP-DiD Estimates: Effect on Overall 3P Index}" _n
file write tex "\label{tab:pooled_pmd_overall3p}" _n
file write tex "\begin{tabular}{lc}" _n
file write tex "\toprule" _n
file write tex " & Overall ATT \\" _n
file write tex "\midrule" _n
file write tex "Reweighted PMD LP-DiD         & `cell_pool_pmd_rw' \\" _n
file write tex "[6pt]" _n
file write tex "Variance-weighted PMD LP-DiD  & `cell_pool_pmd_vw' \\" _n
file write tex "\midrule" _n
file write tex "Observations (RW)   & `pool_pmd_N_rw'   \\" _n
file write tex "Observations (VW)   & `pool_pmd_N_vw'   \\" _n
file write tex "Within R\$^2\$ (RW) & `pool_pmd_r2w_rw' \\" _n
file write tex "Within R\$^2\$ (VW) & `pool_pmd_r2w_vw' \\" _n
file write tex "\bottomrule" _n
file write tex "\end{tabular}" _n
file write tex "\begin{tablenotes}" _n
file write tex "\small" _n
file write tex "\item \textit{Notes:} Pooled PMD LP-DiD estimates. " _n
file write tex "Outcome: average of post-treatment PMD long differences in Overall 3P Index (h=0 to h=`post_window'). " _n
file write tex "t-statistics in parentheses, clustered at country level. " _n
file write tex "* p\$<\$0.10, ** p\$<\$0.05, *** p\$<\$0.01." _n
file write tex "\end{tablenotes}" _n
file write tex "\end{table}" _n
file close tex
di "Saved: pooled_pmd_overall3p.tex"

***************************************************
* ESTIMATES GRAPH
***************************************************
preserve
bysort t: keep if _n == 1
sort t
keep if t >= -`pre_window' & t <= `post_window'

tw (line b_twfe        t, lcolor(red)    lpattern(dash)      lwidth(medium)) ///
   (line b_lpdid       t, lcolor(green)  lpattern(longdash)  lwidth(medium)) ///
   (line b_rw_lpdid    t, lcolor(orange) lpattern(solid)     lwidth(thin))   ///
   (line b_pmd_lpdid   t, lcolor(purple) lpattern(longdash)  lwidth(thin))   ///
   (line b_rwpmd_lpdid t, lcolor(blue)   lpattern(longdash)  lwidth(thin)),  ///
    legend(order(1 "TWFE Distributed Lags" 2 "LP-DiD" 3 "Reweighted LP-DiD" 4 "PMD LP-DiD" 5 "Reweighted PMD LP-DiD") pos(6) col(2)) ///
    xtitle("Event Time") ytitle("Treatment Effect") ///
    note("Units: `units', Periods: `time_periods'" "Obs range: `Nmin'–`Nmax', Countries: `min_units'–`max_units'", ring(0) pos(10) size(*0.8) box) ///
    title("Treatment Effects for Protocol Ratification on Prosecution") ///
    name(effects_5, replace)

restore

***************************************************
* SAMPLE TABLE
***************************************************
preserve
bysort t: keep if _n == 1
sort t

export excel t N_lpdid N_units ///
    using sample_by_event_time.xlsx, ///
    firstrow(variables) replace

listtex t N_lpdid N_units using sample_by_event_time_prosecution.tex, replace ///
    rstyle(tabular) ///
    head("\begin{tabular}{lcc}" ///
         "\toprule" ///
         "Event time & Observations & Units \\" ///
         "\midrule") ///
    foot("\bottomrule" "\end{tabular}")

putdocx begin
putdocx table tbl = data(t N_lpdid N_units)
putdocx save sample_by_event_time_prosecution.docx, replace

restore

***************************************************
* EVENT-STUDY TABLE
***************************************************
preserve
bysort t: keep if _n == 1
sort t
keep if t >= -`pre_window' & t <= `post_window'

format b_*     %9.3f
format tstat_* %9.2f

gen star_twfe  = cond(pval_twfe<0.01,"***",  cond(pval_twfe<0.05,"**",  cond(pval_twfe<0.10,"*","")))
gen star_lpdid = cond(pval_lpdid<0.01,"***", cond(pval_lpdid<0.05,"**", cond(pval_lpdid<0.10,"*","")))
gen star_rw    = cond(pval_rw_lpdid<0.01,"***",   cond(pval_rw_lpdid<0.05,"**",   cond(pval_rw_lpdid<0.10,"*","")))
gen star_pmd   = cond(pval_pmd_lpdid<0.01,"***",  cond(pval_pmd_lpdid<0.05,"**",  cond(pval_pmd_lpdid<0.10,"*","")))
gen star_rwpmd = cond(pval_rwpmd_lpdid<0.01,"***",cond(pval_rwpmd_lpdid<0.05,"**",cond(pval_rwpmd_lpdid<0.10,"*","")))

gen cell_twfe  = string(b_twfe,"%9.3f")        + star_twfe  + " (" + string(tstat_twfe,"%9.2f")        + ")"
gen cell_lpdid = string(b_lpdid,"%9.3f")       + star_lpdid + " (" + string(tstat_lpdid,"%9.2f")       + ")"
gen cell_rw    = string(b_rw_lpdid,"%9.3f")    + star_rw    + " (" + string(tstat_rw_lpdid,"%9.2f")    + ")"
gen cell_pmd   = string(b_pmd_lpdid,"%9.3f")   + star_pmd   + " (" + string(tstat_pmd_lpdid,"%9.2f")   + ")"
gen cell_rwpmd = string(b_rwpmd_lpdid,"%9.3f") + star_rwpmd + " (" + string(tstat_rwpmd_lpdid,"%9.2f") + ")"

keep t cell_twfe cell_lpdid cell_rw cell_pmd cell_rwpmd
sort t

export excel using lpdid_all_prosecution.xlsx, firstrow(variables) replace

outfile t cell_twfe cell_lpdid cell_rw cell_pmd cell_rwpmd ///
    using lpdid_all_prosecution.txt, replace wide

putdocx begin
putdocx paragraph, halign(center)
putdocx text ("Event Study Results"), bold

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
putdocx save "lpdid_all_prosecution.docx", replace

listtex using lpdid_all_prosecution.tex, replace ///
    rstyle(tabular) ///
    head("\begin{tabular}{lccccc}" ///
         "\toprule" ///
         "t & TWFE & LP-DiD & RW LP-DiD & PMD & RWPMD \\" ///
         "\midrule") ///
    foot("\bottomrule" "\end{tabular}")

restore
