quietly {
* Program to simulate data sets with fixed firm effects
*	Estimates coefficients with OLS and Fama-MacBeth
*	Estimates standard errors with OLS, Rogers/Cluster (by firm), and Fama Macbeth 
*
* To run program type do simulation followed by three parameters. They are:
* `1' = Firm effect in X (firm_x) -- Magnitude of firm in the independent variable (0.00-1.00)
* `2' = Firm effect in r (firm_r) -- Magnitude of firm in the residual (0.00-1.00)
* `3' = T = Number of years per firm

#delimit ;
* ----- The following lines keep track of how long it takes to run the program -----;
local start = c(current_time);
local startd = c(current_date);
local startn = 24*60*60*date(c(current_date),"dmy")
		+ 3600 * real(substr(c(current_time),1,2))
		+ 60 * real(substr(c(current_time),4,2)) 
		+ real(substr(c(current_time),7,2));
dis "Start time is $_startd at $_start;  Start index is: $_startn  " ;

drop _all;
program drop _all;
set more 1;
set memory 400m;

* ----- Set the number of observations per data set, and the number of iterations (datasets);
local numobs = 5000;		* 5000 observations;
local numiter = 5000;		* 5000 iterations;
local nmax = max($_numobs,$_numiter);

set obs $_nmax;		* this sets the number of obs in a blank data set equal to # observations numobs;
local ntime = `3';
gen firmid =  1+ int((_n-1)/$_ntime);	* Sets firm id;
gen ncount = _n;

* ----- Variables to capture estimated coefficients and standard errors in each iteration -----;
gen b_ols=0; 	 	/* Standard OLS regression, with regular standard errors */
gen se_ols=0;
gen b_cluster=0;	/* Standard OLS regression, with rogers/clustered standard errors */
gen se_cluster=0;
gen b_fm=0; 	 	/* Fama-MacBeth regression -- coefficients and standard errors*/
gen se_fm=0;
label var b_ols 	"OLS Coefficient";
label var b_cluster	"OLS Coefficient with Rogers/Cluster SE";
label var b_fm		"Fama-MacBeth Coefficient";
label var se_ols 	"OLS SE";
label var se_cluster	"Cluster/Rogers SE";
label var se_fm		"Fama-MacBeth SE";

gen year=0;
gen r=0;
gen x=0;
gen y=0;
gen firm_temp=0;;

* ----- Construct parameters for data structure: firm effect -----;
local firm_x = `1';	/* % of independent variable variance which is due to firm fixed effect */
local firm_r = `2';  	/* % of residual variance which is due to firm fixed effect */

local j=1;	
while `j'<=$_numiter {;							/*  Number of times the simulation is run */

	sort firmid;
	by firmid: replace year=_n;
	replace x = sqrt(1-$_firm_x)*invnorm(uniform());		/* generates non-firm specific portion of x */
	replace firm_temp = sqrt($_firm_x) * invnorm(uniform());	/* generates firm specific portion of x */
	by firmid: replace x = x + firm_temp[_N];			/* Var(X) = 1 */
										
	replace r = sqrt(1-$_firm_r)*invnorm(uniform()); 		/* generates non-firm specific portion of r */
	replace firm_temp = sqrt($_firm_r) * invnorm(uniform());	/* generates firm specific portion of r */
	by firmid: replace r = 2* (r + firm_temp[_N]);			/* Var(r) = 4 */
	
	replace y = 1 * x + r;						/* generates dependent variable */
									/* R2 = scaling = [V(X)/(V(X)+V(5)] = 1/(1+4) */	
									
	* Estimate Standard OLS regression, with regular standard errors, and save coefficient & SE; 
	reg y x;
	replace b_ols  = _b[x]  if ncount==`j';
	replace se_ols = _se[x] if ncount==`j';

	* Estimate Standard OLS regression, with rogers/clustered standard errors, and save coefficient & SE;  
	reg y x, robust cluster(firmid);
	replace b_cluster  = _b[x]  if ncount==`j';
	replace se_cluster = _se[x] if ncount==`j';

	* Run FM regression;
	local k=1;
	while `k'<=$_ntime {;	
		regress y x if year==`k';
		replace b_fm = b_fm + _b[x]/$_ntime if ncount==`j';
		replace se_fm = se_fm + _b[x]^2 if ncount==`j';
		local k = `k' + 1;
		};
	replace se_fm = sqrt((se_fm/($_ntime-1))-((b_fm^2)*($_ntime/($_ntime-1))))/sqrt($_ntime) if ncount==`j';
	
	local j = `j' + 1;
	};


gen firm_x = $_firm_x;
gen firm_r = $_firm_r;
keep b_ols se_ols b_cluster se_cluster b_fm se_fm firm_x firm_r ncount;

sum b_ols if ncount <= $_numiter;
local abols = r(mean);
local sbols = r(sd);
sum se_ols if ncount <= $_numiter;
local asols = r(mean);
sum se_cluster if ncount <= $_numiter;
local ascls = r(mean);
sum b_fm if ncount <= $_numiter;
local abfmb = r(mean);
local sbfmb = r(sd);
sum se_fm if ncount <= $_numiter;
local asfmb = r(mean); 

};

log using simulation, replace;
* se_simulation.do;
* 	Written by Mitchell Petersen, 2005;
*	Estimating Standard Errors In Finance Panel Data Sets: Comparing Approches; 

dis " Firm effect in X (firm_x) is: " `1' " & in r (firm_r) is: " `2';
dis " Time effect in X (time_x) is: " `3' " & in r (time_r) is: " `4';
dis " # of obs = " %5.0f $_numobs "  & time periods per firm = " %3.0f $_ntime "  # of iterations = " %5.0f $_numiter;
dis _newline " Avg(b-ols)= " %6.4f $_abols "  Std(b-ols)= " %6.4f $_sbols "  Avg(s-ols)= " %6.4f $_asols "  Avg(s-cls)= " %6.4f $_ascls;
dis _newline " Avg(b-fmb)= " %6.4f $_abfmb "  Std(b-fmb)= " %6.4f $_sbfmb "  Avg(s-fmb)= " %6.4f $_asfmb;

	local stop = c(current_time);
	local stopn = 24*60*60*date(c(current_date),"dmy")
		+ 3600 * real(substr(c(current_time),1,2))
		+ 60 * real(substr(c(current_time),4,2)) 
		+ real(substr(c(current_time),7,2));
	dis "Program took " %7.1f int(($_stopn - $_startn)/3600) " hours and " %4.1f mod(($_stopn - $_startn),3600)/60 " minutes to run ";

log close;