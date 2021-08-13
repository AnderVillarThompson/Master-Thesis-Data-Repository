cd "ROOT DIRECTORY"
*set trace on
********* INPUT PARAMETER SECTION **********
// Declare max degree for p and q in ARIMA(p,d,q)
local deg = 2

// Declare max year for prediction
local max_year = 2050
// Declare Reference Year
local base_year= 2005
// Declare target year
local prob_year = 2030 
// Declare reduction percentage
local reduction = 0.39 //[0,1]

// Variable for forecasting
local var_for = "EUETS"

// Options (1=yes, 0=no)
local limiter = 1 // Limit size to that of ESD for internal consistency
local cte = 0 // Declare if test should be done with constant
// Declare coefficients. If one of these is "." then the program will simply run tests to find d or to find p & q
local d = 1 // If ".", the program will display the desired degree +1 so that the user can judge if the order is adequate
local p =0
local q =0

// The program requires innitialization before overall forecasts can be acheived, or if the dataset needs to be reset. If this is the case, set the following to 1
local reset = 0

********** REQUIRED PACKAGES **********
ssc install adftest
ssc install sxpose


********** DATA PREPARATION **********

// Index master data to loop through and convert all contents, aided by: https://www.statalist.org/forums/forum/general-stata-discussion/general/1303673-looping-over-files-in-a-folder
if `reset' == 1{
local files : dir "Master Data" files "*.tsv"
foreach file in `files' {
	// Import all data from file for cutting
	local dirname = "Master Data" + "\" + "`file'"
	import delimited using "`dirname'"
  
	drop if (strpos(v1, "TIME_PERIOD")==0 &strpos(v1, "time")==0) & strpos(v1, "ES")==0 //Keep Spain & year data
	replace v1 = subinstr(v1, ",ES", "", .) // Remove ES from data (is implicit)
	replace v1 = "year" if _n == 1 // Change first entry to year
	
	// Select needed categories
	keep if strpos(v1, "year")>0 | strpos(v1, "MIO_T,GHG,")>0 | (strpos(v1, "GIC,")>0 & strpos(v1, ",TJ")>0) | strpos(v1, "T2020_35T")>0 
	
	// Transpose data
	sxpose, force clear // need to download sxpose
	
	// Rename variables to top observation
	local j = 1
	foreach v of var * {
		capture rename `v' `=strtoname(`v'[1])'
	}
	macro drop j
	drop in 1
	quietly destring *, replace force
	
	// Save as .dta with same file name
	local savename = "Cut Data" + "\" + "`file'"
	local savename1 = subinstr("`savename'",".tsv","",.)
	save "`savename1'", replace
  
	clear
}
macro drop files

// Combine all datasets
gen year = .
local files : dir "Cut Data" files "*.dta"
foreach file in `files' {
	local dirname = "Cut Data" + "\" + "`file'"
	merge 1:1 year using "`dirname'"
	drop _merge
}

// Eliminate empty years
egen check = rowtotal(A_MIO_T_GHG* A_T2020_35T)
quietly sum year
local i = 1
while `i'<=`=_N'{
	if check[`i']==0{
		drop in `i'
	}
	else{
		local i = `i' + 1
	}
	quietly sum year
}
drop check

sort year 

rename A_T2020_35T ESD
rename A_MIO_T_GHG_CRF4 LULUCF
// Generate EU ETS data
gen EUETS = A_MIO_T_GHG_TOTX4_MEMO-A_MIO_T_GHG_CRF1A3A-ESD if ESD!=.

save "Combined Data\Combined Data", replace
clear
}
use "Combined Data\Combined Data"


********** ARIMA Modeling **********
// Find earliest variable (USES ESD AS LIMITING FACTOR)
local start = .
local stop = .

if `limiter' == 1{
	gen check = ESD
}
else{
	gen check = `var_for'
}

sum check

forval i = 1/`=_N'{
	if check[`i'] != . &`start' ==.{
		local start = `i'
	}
	else if check[`i'] == . & `start'!=. & `stop' == .{
		local stop = `i'-1
	}
	
}
local size = `stop' - `start'

drop check


// Expand data selection
quietly sum year
local T = _N // Generate end index for years
local m = `T' + `max_year' - year[`=_N']
set obs `m'
gen x = _N
quietly sum year
local s = r(N) + 1
quietly sum x
local endd = r(N)
forval i = `s'/`endd'{
    replace year = year[`i'-1] + 1 in `i'
}
drop x 


// Find differentiation order
tsset year
if `d' == .{
local check = 0 // Binary variable that addresses if stationarity has been acheived at a 95% confidence level
local d = 0 // The order of differentiation being checked
gen diff = `var_for' if year>=year[`start'] & year<=year[`stop'] 

while `check' == 0{
	
	quietly varsoc diff // Test for lags in ADF test
	
	// Find min AIC lag
	local min_aic=.
	local min_index=.
	local lag_aic = `r(mlag)'+1
	forval i = 1/`lag_aic' {
		if `min_aic' == . | `min_aic'>r(stats)[`i',7] {
			local min_aic = r(stats)[`i',7]
			local min_index = `i'
		}
	}
	
	local min_index = r(stats)[`min_index',1]-1
	quietly dfuller diff, lags(`min_index')
	if r(Zt) < r(cv_10){
		local check = 1
		quietly sum diff
		local lagg = r(N) - 1 - `d'
		ac diff, lags(`lagg') graphregion(color(white)) bgcolor(white) ytitle("ACF") mcolor(b) lcolor(b) ylabel(-1(0.5)1) title(Autocorrelations for {&nabla}{sup:`d'}y{sub:t})
		local graphname = "Graphs\ACF_d" + "`d'" + "`var_for'" + ".jpg"
		graph export "`graphname'", replace
		
		local lagg = floor((`lagg'-3)/2)-1
		pac diff, lags(`lagg') graphregion(color(white)) bgcolor(white) ytitle("PACF") mcolor(b) lcolor(b) ylabel(-1(0.5)1) title(Partial Autocorrelations for {&nabla}{sup:`d'}y{sub:t})
		local graphname = "Graphs\PACF_d" + "`d'" + "`var_for'" + ".jpg"
		graph export "`graphname'", replace
		
	}
	else{
		quietly sum diff
		local lagg = r(N) - 1 - `d'
		ac diff, lags(`lagg') graphregion(color(white)) bgcolor(white) ytitle("ACF") mcolor(b) lcolor(b) ylabel(-1(0.5)1) title(Autocorrelations for {&nabla}{sup:`d'}y{sub:t})
		local graphname = "Graphs\ACF_d" + "`d'" + "`var_for'" + ".jpg"
		graph export "`graphname'", replace
		
		local lagg = floor((`lagg'-3)/2)-1
		pac diff, lags(`lagg') graphregion(color(white)) bgcolor(white) ytitle("PACF") mcolor(b) lcolor(b) ylabel(-1(0.5)1) title(Partial Autocorrelations for {&nabla}{sup:`d'}y{sub:t})
		local graphname = "Graphs\PACF_d" + "`d'" + "`var_for'" + ".jpg"
		graph export "`graphname'", replace
		
		local d = `d' + 1
		gen ddiff = d.diff
		replace diff = ddiff
		drop ddiff
	}
	if `check' == 1{ // Display the ACF and PACF for d+1
		gen ddiff = d.diff
		quietly sum ddiff
		local lagg = r(N) - 1 - (`d' + 1)
		local dd = `d' + 1
		di `dd'
		ac ddiff, lags(`lagg') graphregion(color(white)) bgcolor(white) ytitle("ACF") mcolor(b) lcolor(b) ylabel(-1(0.5)1) title(Autocorrelations for {&nabla}{sup:`dd'}y{sub:t})
		local graphname = "Graphs\ACF_d" + "`dd'" + "`var_for'" + ".jpg"
		graph export "`graphname'", replace
		
		local lagg = floor((`lagg'-3)/2)-1
		pac ddiff, lags(`lagg') graphregion(color(white)) bgcolor(white) ytitle("PACF") mcolor(b) lcolor(b) ylabel(-1(0.5)1) title(Partial Autocorrelations for {&nabla}{sup:`dd'}y{sub:t})
		local graphname = "Graphs\PACF_d" + "`dd'" + "`var_for'" + ".jpg"
		graph export "`graphname'", replace
		
	}
}
drop diff
}


// Obtain performance values for tests
local aic = .
local aicc = .
local bic = .
local logl = .
local rmse = .
else if `p'==. | `q'==.{
forval p = 0/`deg'{
	forval q = 0/`deg'{
		capture quietly drop pred res
		if `cte' == 1{
			capture quietly arima `var_for' if year>=year[`start'] & year<=year[`stop'], arima(`p',`d',`q') 
		}
		else{
			capture quietly arima `var_for' if year>=year[`start'] & year<=year[`stop'], arima(`p',`d',`q') noconstant
		}
		
		if _rc==0 {

		
		quietly predict pred if year>=year[`start'], xb
		quietly predict res if year>=year[`start'], residuals
		quietly estat ic // Finds AIC and BIC
		local aic = r(S)[1,5]
		local bic = r(S)[1,6]
		local k = .
		if `cte' == 1{
			local k = 1
		}
		else{
			local k = 0
		}		
		local aicc = `aic' + 2*(`p'+`q'+`k'+1)*(`p'+`q'+`k'+2)/(e(N)-`p'-`q'-`k'-2)
		local logl = e(ll)
		local rmse = 0
		
		forval j = `start'/`stop'{
			if res[`j'] != .{
			local rmse = `rmse' + res[`j']^2
			}
		}
		
		local rmse = sqrt(`rmse'/(`size'))	
		
		
		
		local aic = round(`aic',0.01)
		local aicc = round(`aicc',0.01)
		local bic = round(`bic',0.01)
		local logl = round(`logl',0.01)
		local rmse = round(`rmse',0.01)
		if `p'+`q' != 0{
			local stable = 1 // Innitially assumes stability
			quietly estat aroots, nograph 
			if `p'>0{
				local ar_size = colsof(r(Modulus_ar))
				forval check = 1/`ar_size'{
					if r(Modulus_ar)[1,`check']>1{
						local stable = 0
					}
				}
			}
			if `q'>0{
				local ma_size = colsof(r(Modulus_ma))
				forval check = 1/`ma_size'{
					if r(Modulus_ma)[1,`check']>1{
						local stable = 0
					}
				}
			}	
			if `stable' == 1{
				di "(""`p'"",""`d'"",""`q'"") & ""`aic'""&""`aicc'""&""`bic'""&""`logl'""&""`rmse'""\\"
				di "\hline"
			}
		}
		else{
			di "(""`p'"",""`d'"",""`q'"") & ""`aic'""&""`aicc'""&""`bic'""&""`logl'""&""`rmse'""\\"
			di "\hline"
		}	
		quietly drop pred res
		}
	}
}
}
else{
// Perform ARIMA prediction based on calculated values
gen check = `var_for' if year>=year[`start'] & year<=year[`stop']
tsset year
di `cte'
if `cte' == 1{
	quietly arima check, arima(`p',`d',`q') 
}
else{
	quietly arima check, arima(`p',`d',`q') noconstant
}

predict pred, dynamic(year[`stop']+1)  y
predict res, residuals
drop check

gen UL95 =.
gen LL95 =.
gen UL80 =.
gen LL80 =.
quietly sum year
local limit = `max_year'-year[`stop']
forval i = 1/`limit'{
    local at = `stop'+`i'
    quietly sum res
	if `cte' == 1{
		replace UL95 = pred[`at']+1.96*(r(sd)*sqrt(`i'*(1+`i'/r(N)))) in `at'
		replace LL95 = pred[`at']-1.96*(r(sd)*sqrt(`i'*(1+`i'/r(N)))) in `at'
		replace UL80 = pred[`at']+1.28*(r(sd)*sqrt(`i'*(1+`i'/r(N)))) in `at'
		replace LL80 = pred[`at']-1.28*(r(sd)*sqrt(`i'*(1+`i'/r(N)))) in `at'		
	}
	else{
		replace UL95 = pred[`at']+1.96*(r(sd)*sqrt(`i')) in `at'
		replace LL95 = pred[`at']-1.96*(r(sd)*sqrt(`i')) in `at'
		replace UL80 = pred[`at']+1.28*(r(sd)*sqrt(`i')) in `at'
		replace LL80 = pred[`at']-1.28*(r(sd)*sqrt(`i')) in `at'
	}
}

twoway (rarea LL95 UL95 year, color(gs15)) (rarea LL80 UL80 year, color(gs14)) (line `var_for' year if year<=year[`stop'], lcolor(b)) (line pred year,lcolor(red)) if year>=year[`start'], xlabel(2005(5)2050) ytitle(MtCO{sub:2}eq) legend(order(1 "95% conf." 2 "80% conf." 3 "Historical Trend"  4 "Fitted Model")) graphregion(color(white)) bgcolor(white) title(Projected GHG Emissions (ARIMA(`p',`d',`q')) (`var_for'))
local graphname = "Graphs\FORECAST_ARIMA(" + "`p'"+","+ "`d'" + ","+"`q'"+ ")"+"`var_for'"+".jpg"
		graph export "`graphname'", replace
		

capture estat aroots, graphregion(color(white)) bgcolor(white) mcolor(b) lcolor(b) title(Inverse roots for ARIMA(`p',`d',`q'))
local graphname = "Graphs\ROOTS_ARIMA(" + "`p'"+","+ "`d'" +","+ "`q'"+ ")_TOTX4_MEMO.jpg"
capture graph export "`graphname'", replace

local prob_at = 100
quietly sum res

// Display probabilities depending on the existance of a constant
if `cte' ==1{
	local prob_at = round(`prob_at'*normal(((1-`reduction')*`var_for'[`stop'-(year[`stop']-`base_year')]-pred[`stop'+`prob_year'-year[`stop']])/(r(sd)*sqrt((`prob_year'-year[`stop'])*(1+(`prob_year'-year[`stop'])/(`size'+1))))),0.01)
}
else{
	local prob_at = round(`prob_at'*normal(((1-`reduction')*`var_for'[`stop'-(year[`stop']-`base_year')]-pred[`stop'+`prob_year'-year[`stop']])/(r(sd)*sqrt(`prob_year'-year[`stop']))),0.01)
}


		
capture gen `var_for'pred = pred
capture replace `var_for'pred = pred
capture gen `var_for'res = res
capture replace `var_for'res = res
drop pred res UL* LL*
save "Combined Data\Combined Data", replace

}
// Create graphs. The following is hard coded to the specific tests for this study - it is not dynamic.
// EUETS
quietly sum EUETSres
local target = (1-0.61)*EUETS[`stop'-(year[`stop']-`base_year')]
local p = round(100*normal((`target'-EUETSpred[`stop'+`prob_year'-year[`stop']])/(r(sd)*sqrt(`prob_year'-year[`stop']))),0.01)


twoway (line EUETS year if year<=year[`stop'], lcolor(b)) (line EUETSpred year,lcolor(red)) (scatteri `target' `prob_year' (3) "Approx. `p'%", mcolor(blue) mlabcolor(blue)) if year>=year[`start'], xlabel(2005(5)2050) ytitle(MtCO{sub:2}eq) legend(order(1 "Historical Trend" 2 "Fitted Model" 3 "Internal Goal")) graphregion(color(white)) bgcolor(white) title(Projected Goal Attainment (EUETS))
local graphname = "Graphs\PROB_EUETS.jpg"
graph export "`graphname'", replace

// ESD
quietly sum ESDres
local target1 = (1-0.26)*ESD[`stop'-(year[`stop']-`base_year')]
local p1 = round(100*normal((`target1'-ESDpred[`stop'+`prob_year'-year[`stop']])/(r(sd)*sqrt((`prob_year'-year[`stop'])*(1+(`prob_year'-year[`stop'])/(`size'+1))))),0.01)
local target2 = (1-0.39)*ESD[`stop'-(year[`stop']-`base_year')]
local p2 = round(100*normal((`target2'-ESDpred[`stop'+`prob_year'-year[`stop']])/(r(sd)*sqrt((`prob_year'-year[`stop'])*(1+(`prob_year'-year[`stop'])/(`size'+1))))),0.01)

twoway (line ESD year if year<=year[`stop'], lcolor(b)) (line ESDpred year,lcolor(red)) (scatteri `target1' `prob_year' (3) "Approx. `p1'%", mcolor(dkgreen) mlabcolor(dkgreen)) (scatteri `target2' 2030 (3) "Approx. `p2'%", mcolor(blue) mlabcolor(blue)) if year>=year[`start'], xlabel(2005(5)2050) ytitle(MtCO{sub:2}eq) legend(order(1 "Historical Trend" 2 "Fitted Model" 3 "Mandated Goal" 4 "Internal Goal")) graphregion(color(white)) bgcolor(white) title(Projected Goal Attainment (ESD))
local graphname = "Graphs\PROB_ESD.jpg"
graph export "`graphname'", replace

// LULUCF
local target = 0
quietly sum LULUCFres
local p = round(100*normal((`target'-LULUCFpred[`stop'+`prob_year'-year[`stop']])/(r(sd)*sqrt(`prob_year'-year[`stop']))),0.01)


twoway (line LULUCF year if year<=year[`stop'], lcolor(b)) (line LULUCFpred year,lcolor(red)) (scatteri `target' `prob_year' (3) "Approx. `p' %", mcolor(dkgreen) mlabcolor(dkgreen)) if year>=year[`start'], xlabel(2005(5)2050) ytitle(MtCO{sub:2}eq) legend(order(1 "Historical Trend" 2 "Fitted Model" 3 "Mandated Goal")) graphregion(color(white)) bgcolor(white) title(Projected Goal Attainment (LULUCF))
local graphname = "Graphs\PROB_LULUCF.jpg"
graph export "`graphname'", replace



gen TOT_hist = EUETS + ESD + LULUCF if year>=year[`start'] & year<=year[`stop']
gen TOT_histX4 = EUETS + ESD if year>=year[`start'] & year<=year[`stop']
egen TOT_pred = rowtotal(EUETSpred ESDpred)

// Check 2030 goals
local tot_var = 0
quietly sum EUETSres
local tot_var = `tot_var' + (r(sd)*sqrt(`prob_year'-year[`stop']))^2
quietly sum ESDres
local tot_var = `tot_var' + (r(sd)*sqrt((`prob_year'-year[`stop'])*(1+(`prob_year'-year[`stop'])/(`size'+1))))^2
local tot_sd = sqrt(`tot_var')


local 2030target1 = (1-0.2)*(A_MIO_T_GHG_TOTX4_MEMO[1]-A_MIO_T_GHG_CRF1A3A[1])
local 2030target2 = (1-0.23)*(A_MIO_T_GHG_TOTX4_MEMO[1]-A_MIO_T_GHG_CRF1A3A[1])

local 2030p1 =round(100*normal((`2030target1'-TOT_pred[`stop'+`prob_year'-year[`stop']])/`tot_sd'),0.01)
local 2030p2 =round(100*normal((`2030target2'-TOT_pred[`stop'+`prob_year'-year[`stop']])/`tot_sd'),0.01)

twoway (line TOT_histX4 year if year<=year[`stop'], lcolor(b)) (line TOT_pred year if year>=year[`start']+1,lcolor(red)) (scatteri `2030target1' 2030 (12) "Approx. `2030p1' %", mcolor(dkgreen) mlabcolor(dkgreen)) (scatteri `2030target2' 2030 (6) "Approx. `2030p2' %", mcolor(blue) mlabcolor(blue)) if year>=year[`start'], xlabel(2005(5)2050) ylabel(0(50)450) ytitle(MtCO{sub:2}eq) legend(order(1 "Historical Trend" 2 "Fitted Model" 3 "Mandated Goal" 4 "Internal Goal")) graphregion(color(white)) bgcolor(white) title(Projected Goal Attainment (TOTAL))
local graphname = "Graphs\PROB_TOT2030.jpg"
graph export "`graphname'", replace


// Check 2050 goals
drop TOT_pred
egen TOT_pred = rowtotal(EUETSpred ESDpred LULUCFpred)
local tot_var = 0
quietly sum EUETSres
local tot_var = `tot_var' + (r(sd)*sqrt(`max_year'-year[`stop']))^2
quietly sum ESDres
local tot_var = `tot_var' + (r(sd)*sqrt((`max_year'-year[`stop'])*(1+(`max_year'-year[`stop'])/(`size'+1))))^2
quietly sum LULUCFres
local tot_var = `tot_var' + (r(sd)*sqrt(`max_year'-year[`stop']))^2
local tot_sd = sqrt(`tot_var')

local 2050target = 0
local 2050p = round(100*normal((`2050target'-TOT_pred[`stop'+`max_year'-year[`stop']])/`tot_sd'),0.01)


twoway (line TOT_hist year if year<=year[`stop'], lcolor(b)) (line TOT_pred year if year>=year[`start']+1,lcolor(red)) (scatteri `2050target' 2050 (9) "Approx. `2050p' %", mcolor(dkgreen) mlabcolor(dkgreen)) if year>=year[`start'], xlabel(2005(5)2050) ytitle(MtCO{sub:2}eq) legend(order(1 "Historical Trend" 2 "Fitted Model" 3 "Mandated Goal" )) graphregion(color(white)) bgcolor(white) title(Projected Goal Attainment (TOTAL))
local graphname = "Graphs\PROB_TOT2050.jpg"
graph export "`graphname'", replace
clear