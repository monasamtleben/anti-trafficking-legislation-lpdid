/*
    LP-DiD: Pre-Mean-Differenced (PMD) — No Interaction
    Four outcome variables: overall_3p, prevention, protection, prosecution

	Method: A Local Projections Approach To Difference-In-Differences (Dube, Girardi, Jorda' and Taylor, 2023) 
    Code based on: https://github.com/danielegirardi/lpdid/blob/main/STATA%20example%20files/LP_DiD_examplefile.do 
	Author: Daniele Girardi (University of Massachusetts Amherst), dgirardi@umass.edu (2/24/2023)

    Outputs:
    - pmd_eventstudy_no_int.tex
        → combined event-study overview table t=-5 to t=5 (one column per outcome)
    - pmd_pooled_no_int.tex
        → combined pooled RWPMD table: treat / N / R² for all 4 outcomes

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
    gen b_pmd_`ovar'  = .
    gen se_pmd_`ovar' = .
    gen N_pmd_`ovar'  = .
}

***************************************************
* PMD event-study regressions (per outcome, per horizon)
***************************************************
xtset unit time

foreach ovar of local outcomes {

    forval j = 0/`post_window' {

        * Post-period regression
        reghdfe PMD`j'y_`ovar' D.treat ///
            if D.treat==1 | F`j'.treat==0, ///
            absorb(time) vce(cluster unit)

        replace b_pmd_`ovar'  = _b[D.treat] if t==`j'
        replace se_pmd_`ovar' = _se[D.treat] if t==`j'
        replace N_pmd_`ovar'  = e(N)         if t==`j'

        * Pre-period regression (j>=2; t=-1 is the omitted baseline)
        if `j'>1 & `j'<=`pre_window' {
            reghdfe PMDm`j'y_`ovar' D.treat ///
                if D.treat==1 | treat==0, ///
                absorb(time) vce(cluster unit)

            replace b_pmd_`ovar'  = _b[D.treat] if t==-`j'
            replace se_pmd_`ovar' = _se[D.treat] if t==-`j'
            replace N_pmd_`ovar'  = e(N)         if t==-`j'
        }
    }

    * t=-1 baseline: normalize to zero
    replace b_pmd_`ovar'  = 0 if t==-1
    replace se_pmd_`ovar' = 0 if t==-1

    * t-statistics and p-values
    gen t_pmd_`ovar' = b_pmd_`ovar' / se_pmd_`ovar'
    gen p_pmd_`ovar' = 2*ttail(N_pmd_`ovar', abs(t_pmd_`ovar'))
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
* Combined event-study overview table
* Rows = event time; Columns = 4 outcomes
***************************************************
preserve
bysort t: keep if _n==1
sort t
keep if t >= -`pre_window' & t <= `post_window'

file open tex using pmd_eventstudy_no_int.tex, write replace

file write tex "\begin{table}[htbp]" _n
file write tex "\centering" _n
file write tex "\caption{PMD LP-DiD Event Study: Treatment Effect by Outcome (t=$-5$ to $5$)}" _n
file write tex "\label{tab:pmd_es_no_int}" _n
file write tex "\begin{threeparttable}" _n
file write tex "\begin{tabular}{lcccc}" _n
file write tex "\toprule" _n
file write tex "Event Time & Overall 3P & Prevention & Protection & Prosecution \\" _n
file write tex "\midrule" _n

local Nrows = _N
forvalues i = 1/`Nrows' {
    local tt = t[`i']

    if `tt' == -1 {
        file write tex "-1 & \multicolumn{4}{c}{0 (baseline)} \\" _n
    }
    else {
        foreach ovar in overall_3p prevention protection prosecution {
            local b = b_pmd_`ovar'[`i']
            local tv = t_pmd_`ovar'[`i']
            local p = p_pmd_`ovar'[`i']
            get_stars `p'
            local c_`ovar' = "\makecell{" + string(`b', "%6.3f") + "`r(stars)'" + "\\" + "(" + string(`tv', "%5.2f") + ")}"
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
file write tex "\item \textit{Notes:} PMD LP-DiD estimates. Each cell: \makecell{coefficient\\(t-statistic)}." _n
file write tex " t-statistics based on SE clustered at country level." _n
file write tex " t=$-1$ is the omitted baseline. * p\$<\$0.10, ** p\$<\$0.05, *** p\$<\$0.01." _n
file write tex "\end{tablenotes}" _n
file write tex "\end{threeparttable}" _n
file write tex "\end{table}" _n

file close tex
di "Saved: pmd_eventstudy_no_int.tex"
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
* Pooled RWPMD regressions (per outcome, no interaction)
***************************************************
foreach ovar of local outcomes {

    reghdfe PMD_pool_`ovar' dtreat ///
        if (dtreat==1 | F`post_window'.treat==0) ///
        [pweight = reweight_`post_window'], ///
        absorb(time) vce(cluster unit)

    scalar b_pool_`ovar'   = _b[dtreat]
    scalar se_pool_`ovar'  = _se[dtreat]
    scalar t_pool_`ovar'   = _b[dtreat] / _se[dtreat]
    scalar p_pool_`ovar'   = 2*ttail(e(df_r), abs(t_pool_`ovar'))
    scalar N_pool_`ovar'   = e(N)
    scalar r2w_pool_`ovar' = e(r2_within)
}

***************************************************
* Format pooled cells
***************************************************
foreach ovar in overall_3p prevention protection prosecution {
    get_stars `=scalar(p_pool_`ovar')'
    local cell_pool_`ovar' = "\makecell{" + string(scalar(b_pool_`ovar'), "%6.3f") + ///
        "`r(stars)'" + "\\" + "(" + string(scalar(t_pool_`ovar'), "%5.2f") + ")}"
    local N_`ovar'   = string(scalar(N_pool_`ovar'))
    local r2w_`ovar' = string(scalar(r2w_pool_`ovar'), "%6.3f")
}

***************************************************
* Write combined pooled RWPMD table
***************************************************
file open tex using pmd_pooled_no_int.tex, write replace

file write tex "\begin{table}[htbp]" _n
file write tex "\centering" _n
file write tex "\caption{Pooled Reweighted PMD LP-DiD: Average Treatment Effect by Outcome}" _n
file write tex "\label{tab:pmd_pooled_no_int}" _n
file write tex "\begin{threeparttable}" _n
file write tex "\begin{tabular}{lcccc}" _n
file write tex "\toprule" _n
file write tex " & Overall 3P & Prevention & Protection & Prosecution \\" _n
file write tex "\midrule" _n
file write tex "Treat & `cell_pool_overall_3p' & `cell_pool_prevention' & `cell_pool_protection' & `cell_pool_prosecution' \\" _n
file write tex "\midrule" _n
file write tex "Observations  & `N_overall_3p'   & `N_prevention'   & `N_protection'   & `N_prosecution'   \\" _n
file write tex "Within R\$^2\$ & `r2w_overall_3p' & `r2w_prevention' & `r2w_protection' & `r2w_prosecution' \\" _n
file write tex "\bottomrule" _n
file write tex "\end{tabular}" _n
file write tex "\begin{tablenotes}" _n
file write tex "\small" _n
file write tex "\item \textit{Notes:} Reweighted PMD LP-DiD pooled estimates." _n
file write tex " Pooled outcome = average of post-treatment PMD long differences h=0 to h=`post_window'." _n
file write tex " Each cell: \makecell{coefficient\\(t-statistic)}, SE clustered at country level." _n
file write tex " * p\$<\$0.10, ** p\$<\$0.05, *** p\$<\$0.01." _n
file write tex "\end{tablenotes}" _n
file write tex "\end{threeparttable}" _n
file write tex "\end{table}" _n

file close tex
di "Saved: pmd_pooled_no_int.tex"

***************************************************
* ESTIMATES GRAPHS — Development over event time
***************************************************
preserve
bysort t: keep if _n == 1
sort t
keep if t >= -`pre_window' & t <= `post_window'

* 95% CI bounds
foreach ovar of local outcomes {
    gen ci_hi_`ovar' = b_pmd_`ovar' + 1.96 * se_pmd_`ovar'
    gen ci_lo_`ovar' = b_pmd_`ovar' - 1.96 * se_pmd_`ovar'
}

* --- Graph 1: All 4 outcomes on one plot ---
tw (line b_pmd_overall_3p  t, lcolor(blue)   lpattern(solid)     lwidth(medium)) ///
   (line b_pmd_prevention  t, lcolor(red)    lpattern(dash)      lwidth(medium)) ///
   (line b_pmd_protection  t, lcolor(green)  lpattern(longdash)  lwidth(medium)) ///
   (line b_pmd_prosecution t, lcolor(orange) lpattern(shortdash) lwidth(medium)), ///
    yline(0, lcolor(black) lpattern(dot)) xline(-0.5, lcolor(gs10) lpattern(dot)) ///
    legend(order(1 "Overall 3P" 2 "Prevention" 3 "Protection" 4 "Prosecution") pos(6) col(2)) ///
    xtitle("Event Time") ytitle("Treatment Effect") ///
    title("PMD LP-DiD: Treatment Effect by Outcome") ///
    xlabel(-5(1)5) name(pmd_all_outcomes, replace)

graph export "pmd_all_outcomes_no_int.png", replace

* --- Graph 2: Per-outcome event study with 95% CI bands (combined 2×2) ---
foreach ovar in overall_3p prevention protection prosecution {
    tw (rarea ci_hi_`ovar' ci_lo_`ovar' t, color(blue%20)) ///
       (line b_pmd_`ovar'  t, lcolor(blue) lwidth(medium)), ///
        yline(0, lcolor(black) lpattern(dot)) xline(-0.5, lcolor(gs10) lpattern(dot)) ///
        legend(off) ///
        xtitle("Event Time") ytitle("Treatment Effect") ///
        title("`ovar'") xlabel(-5(1)5) name(es_`ovar', replace)
}

graph combine es_overall_3p es_prevention es_protection es_prosecution, ///
    title("PMD LP-DiD Event Study by Outcome") ///
    cols(2) name(es_combined_no_int, replace)

graph export "pmd_eventstudy_combined_no_int.png", replace

restore
