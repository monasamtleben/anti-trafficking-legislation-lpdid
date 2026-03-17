/*
    LP-DiD: Pre-Mean-Differenced (PMD) Interaction Analysis
    Four outcome variables: overall_3p, prevention, protection, prosecution
    Interaction: high_rol_pre (above-median Rule of Law, pre-treatment)
	
	Method: A Local Projections Approach To Difference-In-Differences (Dube, Girardi, Jorda' and Taylor, 2023) 
    Code based on: https://github.com/danielegirardi/lpdid/blob/main/STATA%20example%20files/LP_DiD_examplefile.do 
	Author: Daniele Girardi (University of Massachusetts Amherst), dgirardi@umass.edu (2/24/2023)

    Outputs:
    - pmd_eventstudy_overall_3p.tex / _prevention / _protection / _prosecution
        → per-outcome event-study table t=-5 to t=5 (Low / High / Diff columns)
    - pmd_eventstudy_overview.tex
        → combined overview table: High-Low differential for all 4 outcomes
    - pmd_pooled_all_outcomes.tex
        → combined pooled RWPMD table: treat / high / high×treat / N / R²
*/

clear all
set more off
macro drop _all
capture log close

use "C:\Users\dataset.dta", clear

set scheme s2color
set trace off

***************************************************
* Rename & restrict
***************************************************
rename new_ccode    unit
rename approxTrafYear time
rename unp_adapt    treat

keep if time <= 2013
keep if time >= 1990

***************************************************
* Panel & event-time setup
***************************************************
local post_window 5
local pre_window  5

xtset unit time

gen dtreat = D.treat    // create BEFORE any bysort

bysort unit (time): egen treat_year = min(cond(dtreat == 1, time, .))

gen t = time - treat_year
replace t = . if t < -`pre_window' | t > `post_window'

gen first_treat = time if dtreat == 1
bysort unit: egen cohort = min(first_treat)

***************************************************
* Reweighting weights (used for pooled RWPMD)
***************************************************
xtset unit time

forval j = 0/`post_window' {
    qui gen group_h`j' = .
    qui replace group_h`j' = time if (dtreat==1 | F`j'.treat==0)
    qui reghdfe dtreat if (dtreat==1 | F`j'.treat==0), absorb(time) residuals(num_weights_`j')
    qui replace num_weights_`j' = . if dtreat != 1
    qui egen den_weights_`j' = total(num_weights_`j')
    qui gen weight_`j' = num_weights_`j' / den_weights_`j'
    qui bysort group_h`j': egen gweight_`j' = max(weight_`j')
    qui replace weight_`j' = gweight_`j' if weight_`j' == .
    qui replace weight_`j' = round(weight_`j', 0.00000001)
    qui gen reweight_`j' = 1 / weight_`j'
    qui sort unit time
}

***************************************************
* Create outcome-specific long differences & PMD vars
***************************************************
local outcomes "overall_3p prevention protection prosecution"

xtset unit time

foreach ovar of local outcomes {

    * Long differences
    forval j = 0/`post_window' {
        gen D`j'y_`ovar' = F`j'.`ovar' - L.`ovar'
    }
    forval j = 2/`pre_window' {
        gen Dm`j'y_`ovar' = L`j'.`ovar' - L.`ovar'
    }

    * Pre-treatment mean (correct aveLY: cumsum / n of pre-treatment obs, lagged)
    gen Y_pre_`ovar' = `ovar' if time < treat_year
    bysort unit (time): gen cumsum_`ovar' = sum(Y_pre_`ovar')
    bysort unit (time): gen n_pre_`ovar'  = sum(Y_pre_`ovar' < .)
    xtset unit time
    gen aveLY_`ovar' = L.cumsum_`ovar' / L.n_pre_`ovar'

    * PMD long differences
    forval j = 0/`post_window' {
        gen PMD`j'y_`ovar' = F`j'.`ovar' - aveLY_`ovar'
    }
    forval j = 2/`pre_window' {
        gen PMDm`j'y_`ovar' = L`j'.`ovar' - aveLY_`ovar'
    }
}

***************************************************
* Storage variables for PMD event-study results
***************************************************
foreach ovar of local outcomes {
    gen b_pmd_low_`ovar'   = .
    gen se_pmd_low_`ovar'  = .
    gen b_pmd_high_`ovar'  = .
    gen se_pmd_high_`ovar' = .
    gen b_pmd_diff_`ovar'  = .
    gen se_pmd_diff_`ovar' = .
    gen N_pmd_`ovar'       = .
}

***************************************************
* PMD event-study regressions (per outcome, per horizon)
***************************************************
xtset unit time

foreach ovar of local outcomes {

    forval j = 0/`post_window' {

        * Post-period regression (also j=0 and j=1)
        reghdfe PMD`j'y_`ovar' c.D.treat##i.high_rol_pre ///
            if D.treat==1 | F`j'.treat==0, ///
            absorb(time) vce(cluster unit)

        replace b_pmd_low_`ovar'   = _b[D.treat]                   if t==`j'
        replace se_pmd_low_`ovar'  = _se[D.treat]                  if t==`j'
        replace b_pmd_diff_`ovar'  = _b[1.high_rol_pre#c.D.treat]  if t==`j'
        replace se_pmd_diff_`ovar' = _se[1.high_rol_pre#c.D.treat] if t==`j'
        replace N_pmd_`ovar'       = e(N)                          if t==`j'

        lincom D.treat + 1.high_rol_pre#c.D.treat
        replace b_pmd_high_`ovar'  = r(estimate) if t==`j'
        replace se_pmd_high_`ovar' = r(se)       if t==`j'

        * Pre-period regression (j>=2 only; t=-1 is the omitted baseline)
        if `j'>1 & `j'<=`pre_window' {
            reghdfe PMDm`j'y_`ovar' c.D.treat##i.high_rol_pre ///
                if D.treat==1 | treat==0, ///
                absorb(time) vce(cluster unit)

            replace b_pmd_low_`ovar'   = _b[D.treat]                   if t==-`j'
            replace se_pmd_low_`ovar'  = _se[D.treat]                  if t==-`j'
            replace b_pmd_diff_`ovar'  = _b[1.high_rol_pre#c.D.treat]  if t==-`j'
            replace se_pmd_diff_`ovar' = _se[1.high_rol_pre#c.D.treat] if t==-`j'
            replace N_pmd_`ovar'       = e(N)                          if t==-`j'

            lincom D.treat + 1.high_rol_pre#c.D.treat
            replace b_pmd_high_`ovar'  = r(estimate) if t==-`j'
            replace se_pmd_high_`ovar' = r(se)       if t==-`j'
        }
    }

    * t=-1 baseline: normalize to zero
    replace b_pmd_low_`ovar'   = 0 if t==-1
    replace b_pmd_high_`ovar'  = 0 if t==-1
    replace b_pmd_diff_`ovar'  = 0 if t==-1
    replace se_pmd_diff_`ovar' = 0 if t==-1

    * t-statistics and p-values
    gen t_pmd_low_`ovar'  = b_pmd_low_`ovar'  / se_pmd_low_`ovar'
    gen t_pmd_high_`ovar' = b_pmd_high_`ovar' / se_pmd_high_`ovar'
    gen t_pmd_diff_`ovar' = b_pmd_diff_`ovar' / se_pmd_diff_`ovar'

    gen p_pmd_low_`ovar'  = 2*ttail(N_pmd_`ovar', abs(t_pmd_low_`ovar'))
    gen p_pmd_high_`ovar' = 2*ttail(N_pmd_`ovar', abs(t_pmd_high_`ovar'))
    gen p_pmd_diff_`ovar' = 2*ttail(N_pmd_`ovar', abs(t_pmd_diff_`ovar'))
}

***************************************************
* Stars helper program
***************************************************
capture program drop get_stars
program define get_stars, rclass
    args pval
    if      `pval' < 0.01  return local stars "***"
    else if `pval' < 0.05  return local stars "**"
    else if `pval' < 0.10  return local stars "*"
    else                   return local stars ""
end

***************************************************
* Per-outcome event-study LaTeX tables
* Rows = event time; Columns = Low RoL / High RoL / High-Low
***************************************************
foreach ovar of local outcomes {

    preserve
    bysort t: keep if _n==1
    sort t
    keep if t >= -`pre_window' & t <= `post_window'

    file open tex using pmd_eventstudy_`ovar'.tex, write replace

    file write tex "\begin{table}[htbp]" _n
    file write tex "\centering" _n
    file write tex "\caption{PMD LP-DiD Event Study: \texttt{`ovar'} (Rule of Law Interaction)}" _n
    file write tex "\label{tab:pmd_es_`ovar'}" _n
    file write tex "\begin{threeparttable}" _n
    file write tex "\begin{tabular}{lccc}" _n
    file write tex "\toprule" _n
    file write tex "Event Time & Low RoL & High RoL & High \$-\$ Low \\" _n
    file write tex "\midrule" _n

    local Nrows = _N
    forvalues i = 1/`Nrows' {
        local tt = t[`i']

        if `tt' == -1 {
            file write tex "-1 & 0 (base) & 0 (base) & 0 (base) \\" _n
        }
        else {
            * Low group
            local bl = b_pmd_low_`ovar'[`i']
            local tl = t_pmd_low_`ovar'[`i']
            local pl = p_pmd_low_`ovar'[`i']
            get_stars `pl'
            local sl "`r(stars)'"

            * High group
            local bh = b_pmd_high_`ovar'[`i']
            local th = t_pmd_high_`ovar'[`i']
            local ph = p_pmd_high_`ovar'[`i']
            get_stars `ph'
            local sh "`r(stars)'"

            * Differential
            local bd = b_pmd_diff_`ovar'[`i']
            local td = t_pmd_diff_`ovar'[`i']
            local pd = p_pmd_diff_`ovar'[`i']
            get_stars `pd'
            local sd "`r(stars)'"

            local cl = "\makecell{" + string(`bl', "%6.3f") + "`sl'" + "\\" + "(" + string(`tl', "%5.2f") + ")}"
            local ch = "\makecell{" + string(`bh', "%6.3f") + "`sh'" + "\\" + "(" + string(`th', "%5.2f") + ")}"
            local cd = "\makecell{" + string(`bd', "%6.3f") + "`sd'" + "\\" + "(" + string(`td', "%5.2f") + ")}"

            file write tex "`tt' & `cl' & `ch' & `cd' \\" _n
        }
    }

    * Observation counts at two horizons
    quietly summ N_pmd_`ovar' if t==0
    local Nt0 = string(r(mean), "%5.0f")
    quietly summ N_pmd_`ovar' if t==`post_window'
    local Nt5 = string(r(mean), "%5.0f")

    file write tex "\midrule" _n
    file write tex "N (t=0) & \multicolumn{3}{c}{`Nt0'} \\" _n
    file write tex "N (t=`post_window') & \multicolumn{3}{c}{`Nt5'} \\" _n
    file write tex "\bottomrule" _n
    file write tex "\end{tabular}" _n
    file write tex "\begin{tablenotes}" _n
    file write tex "\small" _n
    file write tex "\item \textit{Notes:} PMD LP-DiD estimates, outcome: \texttt{`ovar'}." _n
    file write tex " Interaction: above-median Rule of Law (pre-treatment)." _n
    file write tex " t-statistics in parentheses, SE clustered at country level." _n
    file write tex " t=$-1$ is the omitted baseline. * p\$<\$0.10, ** p\$<\$0.05, *** p\$<\$0.01." _n
    file write tex "\end{tablenotes}" _n
    file write tex "\end{threeparttable}" _n
    file write tex "\end{table}" _n

    file close tex
    di "Saved: pmd_eventstudy_`ovar'.tex"
    restore
}

***************************************************
* Combined overview table: High-Low differential for all 4 outcomes
* Rows = event time; Columns = 4 outcomes
***************************************************
preserve
bysort t: keep if _n==1
sort t
keep if t >= -`pre_window' & t <= `post_window'

file open tex using pmd_eventstudy_overview_high_rol_pre.tex, write replace

file write tex "\begin{table}[htbp]" _n
file write tex "\centering" _n
file write tex "\caption{PMD LP-DiD: Differential Effect (High \$-\$ Low Rule of Law) by Outcome}" _n
file write tex "\label{tab:pmd_overview}" _n
file write tex "\begin{threeparttable}" _n
file write tex "\begin{tabular}{lcccc}" _n
file write tex "\toprule" _n
file write tex "Event Time & Overall 3P & Prevention & Protection & Prosecution \\" _n
file write tex "\midrule" _n

local Nrows = _N
forvalues i = 1/`Nrows' {
    local tt = t[`i']

    if `tt' == -1 {
        file write tex "-1 & 0 & 0 & 0 & 0 \\" _n
    }
    else {
        foreach ovar in overall_3p prevention protection prosecution {
            local bd = b_pmd_diff_`ovar'[`i']
            local td = t_pmd_diff_`ovar'[`i']
            local pd = p_pmd_diff_`ovar'[`i']
            get_stars `pd'
            local c_`ovar' = "\makecell{" + string(`bd', "%6.3f") + "`r(stars)'" + "\\" + "(" + string(`td', "%5.2f") + ")}"
        }
        file write tex "`tt' & `c_overall_3p' & `c_prevention' & `c_protection' & `c_prosecution' \\" _n
    }
}

* N row at t=0
foreach ovar in overall_3p prevention protection prosecution {
    quietly summ N_pmd_`ovar' if t==0
    local N0_`ovar' = string(r(mean), "%5.0f")
}
file write tex "\midrule" _n
file write tex "N (t=0) & `N0_overall_3p' & `N0_prevention' & `N0_protection' & `N0_prosecution' \\" _n
file write tex "\bottomrule" _n
file write tex "\end{tabular}" _n
file write tex "\begin{tablenotes}" _n
file write tex "\small" _n
file write tex "\item \textit{Notes:} PMD LP-DiD. Each cell reports the differential treatment effect" _n
file write tex " (High RoL minus Low RoL) with t-statistic (SE clustered at country) in parentheses." _n
file write tex " Interaction: above-median Rule of Law (pre-treatment). t=$-1$ is the omitted baseline." _n
file write tex " * p\$<\$0.10, ** p\$<\$0.05, *** p\$<\$0.01." _n
file write tex "\end{tablenotes}" _n
file write tex "\end{threeparttable}" _n
file write tex "\end{table}" _n

file close tex
di "Saved: pmd_eventstudy_overview_high_rol_pre.tex"
restore

***************************************************
* Pooled RWPMD: create averaged post-treatment PMD outcome
***************************************************
xtset unit time

foreach ovar of local outcomes {
    gen PMD_pool_`ovar' = 0
    forval j = 0/`post_window' {
        replace PMD_pool_`ovar' = PMD_pool_`ovar' + PMD`j'y_`ovar'
    }
    replace PMD_pool_`ovar' = PMD_pool_`ovar' / (`post_window' + 1)
}

***************************************************
* Pooled RWPMD regressions (per outcome)
* Coefficients: treat (Low RoL baseline), high (level diff), high x treat (differential)
***************************************************
foreach ovar of local outcomes {

    reghdfe PMD_pool_`ovar' ///
        c.dtreat##i.high_rol_pre ///
        if (dtreat==1 | F`post_window'.treat==0) ///
        [pweight = reweight_`post_window'], ///
        absorb(time) vce(cluster unit)

    scalar b_treat_`ovar'  = _b[dtreat]
    scalar se_treat_`ovar' = _se[dtreat]
    scalar t_treat_`ovar'  = _b[dtreat] / _se[dtreat]
    scalar p_treat_`ovar'  = 2*ttail(e(df_r), abs(t_treat_`ovar'))

    scalar b_high_`ovar'   = _b[1.high_rol_pre]
    scalar se_high_`ovar'  = _se[1.high_rol_pre]
    scalar t_high_`ovar'   = _b[1.high_rol_pre] / _se[1.high_rol_pre]
    scalar p_high_`ovar'   = 2*ttail(e(df_r), abs(t_high_`ovar'))

    scalar b_diff_`ovar'   = _b[1.high_rol_pre#c.dtreat]
    scalar se_diff_`ovar'  = _se[1.high_rol_pre#c.dtreat]
    scalar t_diff_`ovar'   = _b[1.high_rol_pre#c.dtreat] / _se[1.high_rol_pre#c.dtreat]
    scalar p_diff_`ovar'   = 2*ttail(e(df_r), abs(t_diff_`ovar'))

    scalar N_pool_`ovar'   = e(N)
    scalar r2w_pool_`ovar' = e(r2_within)
}

***************************************************
* Format pooled cells
***************************************************
foreach ovar in overall_3p prevention protection prosecution {

    get_stars `=scalar(p_treat_`ovar')'
    local cell_treat_`ovar' = "\makecell{" + string(scalar(b_treat_`ovar'), "%6.3f") + ///
        "`r(stars)'" + "\\" + "(" + string(scalar(t_treat_`ovar'), "%5.2f") + ")}"

    get_stars `=scalar(p_high_`ovar')'
    local cell_high_`ovar' = "\makecell{" + string(scalar(b_high_`ovar'), "%6.3f") + ///
        "`r(stars)'" + "\\" + "(" + string(scalar(t_high_`ovar'), "%5.2f") + ")}"

    get_stars `=scalar(p_diff_`ovar')'
    local cell_diff_`ovar' = "\makecell{" + string(scalar(b_diff_`ovar'), "%6.3f") + ///
        "`r(stars)'" + "\\" + "(" + string(scalar(t_diff_`ovar'), "%5.2f") + ")}"

    local N_`ovar'   = string(scalar(N_pool_`ovar'))
    local r2w_`ovar' = string(scalar(r2w_pool_`ovar'), "%6.3f")
}

***************************************************
* Write combined pooled RWPMD table
***************************************************
file open tex using pmd_pooled_all_outcomes_high_rol_pre.tex, write replace

file write tex "\begin{table}[htbp]" _n
file write tex "\centering" _n
file write tex "\caption{Pooled Reweighted PMD LP-DiD: Treatment Effects by Outcome (Rule of Law Interaction)}" _n
file write tex "\label{tab:pmd_pooled}" _n
file write tex "\begin{threeparttable}" _n
file write tex "\begin{tabular}{lcccc}" _n
file write tex "\toprule" _n
file write tex " & Overall 3P & Prevention & Protection & Prosecution \\" _n
file write tex "\midrule" _n
file write tex "Treat (Low RoL)          & `cell_treat_overall_3p' & `cell_treat_prevention' & `cell_treat_protection' & `cell_treat_prosecution' \\" _n
file write tex "[6pt]" _n
file write tex "High RoL                 & `cell_high_overall_3p'  & `cell_high_prevention'  & `cell_high_protection'  & `cell_high_prosecution'  \\" _n
file write tex "[6pt]" _n
file write tex "High RoL \$\times\$ Treat & `cell_diff_overall_3p'  & `cell_diff_prevention'  & `cell_diff_protection'  & `cell_diff_prosecution'  \\" _n
file write tex "\midrule" _n
file write tex "Observations             & `N_overall_3p'   & `N_prevention'   & `N_protection'   & `N_prosecution'   \\" _n
file write tex "Within R\$^2\$           & `r2w_overall_3p' & `r2w_prevention' & `r2w_protection' & `r2w_prosecution' \\" _n
file write tex "\bottomrule" _n
file write tex "\end{tabular}" _n
file write tex "\begin{tablenotes}" _n
file write tex "\small" _n
file write tex "\item \textit{Notes:} Reweighted PMD LP-DiD pooled estimates." _n
file write tex " Pooled outcome = average of post-treatment PMD long differences h=0 to h=`post_window'." _n
file write tex " Treat (Low RoL) = treatment effect for below-median Rule of Law countries." _n
file write tex " High RoL = level difference between High and Low RoL in pre-mean-differenced outcome." _n
file write tex " High RoL \$\times\$ Treat = differential treatment effect (High minus Low)." _n
file write tex " t-statistics in parentheses, SE clustered at country level." _n
file write tex " * p\$<\$0.10, ** p\$<\$0.05, *** p\$<\$0.01." _n
file write tex "\end{tablenotes}" _n
file write tex "\end{threeparttable}" _n
file write tex "\end{table}" _n

file close tex
di "Saved: pmd_pooled_all_outcomes_high_rol_pre.tex"

***************************************************
* ESTIMATES GRAPHS — Development over event time
***************************************************
preserve
bysort t: keep if _n == 1
sort t
keep if t >= -`pre_window' & t <= `post_window'

* 95% CI bounds
foreach ovar of local outcomes {
    gen ci_hi_low_`ovar'  = b_pmd_low_`ovar'  + 1.96 * se_pmd_low_`ovar'
    gen ci_lo_low_`ovar'  = b_pmd_low_`ovar'  - 1.96 * se_pmd_low_`ovar'
    gen ci_hi_high_`ovar' = b_pmd_high_`ovar' + 1.96 * se_pmd_high_`ovar'
    gen ci_lo_high_`ovar' = b_pmd_high_`ovar' - 1.96 * se_pmd_high_`ovar'
    gen ci_hi_diff_`ovar' = b_pmd_diff_`ovar' + 1.96 * se_pmd_diff_`ovar'
    gen ci_lo_diff_`ovar' = b_pmd_diff_`ovar' - 1.96 * se_pmd_diff_`ovar'
}

* --- Graph 1: Differential effect (High - Low RoL) for all 4 outcomes ---
tw (line b_pmd_diff_overall_3p  t, lcolor(blue)   lpattern(solid)     lwidth(medium)) ///
   (line b_pmd_diff_prevention  t, lcolor(red)    lpattern(dash)      lwidth(medium)) ///
   (line b_pmd_diff_protection  t, lcolor(green)  lpattern(longdash)  lwidth(medium)) ///
   (line b_pmd_diff_prosecution t, lcolor(orange) lpattern(shortdash) lwidth(medium)), ///
    yline(0, lcolor(black) lpattern(dot)) xline(-0.5, lcolor(gs10) lpattern(dot)) ///
    legend(order(1 "Overall 3P" 2 "Prevention" 3 "Protection" 4 "Prosecution") pos(6) col(2)) ///
    xtitle("Event Time") ytitle("Differential Treatment Effect (High - Low RoL)") ///
    title("PMD LP-DiD: Differential by Rule of Law") ///
    xlabel(-5(1)5) name(diff_all_outcomes, replace)

graph export "pmd_diff_all_outcomes_high_rol_pre.png", replace

* --- Graph 2: Low RoL group — all 4 outcomes ---
tw (line b_pmd_low_overall_3p  t, lcolor(blue)   lpattern(solid)     lwidth(medium)) ///
   (line b_pmd_low_prevention  t, lcolor(red)    lpattern(dash)      lwidth(medium)) ///
   (line b_pmd_low_protection  t, lcolor(green)  lpattern(longdash)  lwidth(medium)) ///
   (line b_pmd_low_prosecution t, lcolor(orange) lpattern(shortdash) lwidth(medium)), ///
    yline(0, lcolor(black) lpattern(dot)) xline(-0.5, lcolor(gs10) lpattern(dot)) ///
    legend(order(1 "Overall 3P" 2 "Prevention" 3 "Protection" 4 "Prosecution") pos(6) col(2)) ///
    xtitle("Event Time") ytitle("Treatment Effect") ///
    title("PMD LP-DiD: Low Rule of Law Group") ///
    xlabel(-5(1)5) name(low_all_outcomes, replace)

graph export "pmd_low_all_outcomes_high_rol_pre.png", replace

* --- Graph 3: High RoL group — all 4 outcomes ---
tw (line b_pmd_high_overall_3p  t, lcolor(blue)   lpattern(solid)     lwidth(medium)) ///
   (line b_pmd_high_prevention  t, lcolor(red)    lpattern(dash)      lwidth(medium)) ///
   (line b_pmd_high_protection  t, lcolor(green)  lpattern(longdash)  lwidth(medium)) ///
   (line b_pmd_high_prosecution t, lcolor(orange) lpattern(shortdash) lwidth(medium)), ///
    yline(0, lcolor(black) lpattern(dot)) xline(-0.5, lcolor(gs10) lpattern(dot)) ///
    legend(order(1 "Overall 3P" 2 "Prevention" 3 "Protection" 4 "Prosecution") pos(6) col(2)) ///
    xtitle("Event Time") ytitle("Treatment Effect") ///
    title("PMD LP-DiD: High Rule of Law Group") ///
    xlabel(-5(1)5) name(high_all_outcomes, replace)

graph export "pmd_high_all_outcomes_high_rol_pre.png", replace

* --- Graph 4: Per-outcome event study with 95% CI bands (combined 2×2) ---
foreach ovar in overall_3p prevention protection prosecution {
    tw (rarea ci_hi_low_`ovar'  ci_lo_low_`ovar'  t, color(blue%20)) ///
       (rarea ci_hi_high_`ovar' ci_lo_high_`ovar' t, color(red%20))  ///
       (line b_pmd_low_`ovar'   t, lcolor(blue) lwidth(medium))       ///
       (line b_pmd_high_`ovar'  t, lcolor(red)  lwidth(medium)),      ///
        yline(0, lcolor(black) lpattern(dot)) xline(-0.5, lcolor(gs10) lpattern(dot)) ///
        legend(order(3 "Low RoL" 4 "High RoL") pos(6) col(2)) ///
        xtitle("Event Time") ytitle("Treatment Effect") ///
        title("`ovar'") xlabel(-5(1)5) name(es_`ovar', replace)
}

graph combine es_overall_3p es_prevention es_protection es_prosecution, ///
    title("PMD LP-DiD Event Study by Outcome (Rule of Law Interaction)") ///
    cols(2) name(es_combined, replace)

graph export "pmd_eventstudy_combined_high_rol_pre.png", replace

* --- Graph 4b: Differential effect with 95% CI bands — 2×2 combined ---
foreach ovar in overall_3p prevention protection prosecution {
    tw (rarea ci_hi_diff_`ovar' ci_lo_diff_`ovar' t, color(blue%20)) ///
       (line b_pmd_diff_`ovar'  t, lcolor(blue) lwidth(medium)),      ///
        yline(0, lcolor(black) lpattern(dot)) xline(-0.5, lcolor(gs10) lpattern(dot)) ///
        legend(order(2 "Diff (High - Low RoL)") pos(6) col(1)) ///
        xtitle("Event Time") ytitle("Differential Treatment Effect") ///
        title("`ovar'") xlabel(-5(1)5) name(es_diff_`ovar', replace)
}

graph combine es_diff_overall_3p es_diff_prevention es_diff_protection es_diff_prosecution, ///
    title("PMD LP-DiD: Differential Effect (High - Low RoL) by Outcome") ///
    cols(2) name(es_diff_combined, replace)

graph export "pmd_diff_eventstudy_combined_high_rol_pre.png", replace

restore
