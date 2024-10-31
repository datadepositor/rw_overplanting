********************************************************************************
*				program: simulation, main model, one iteration      		   *									
********************************************************************************
global root "..."

global data "$root/1_data"
global ds "$data/source"
global dd "$data/data"
global dp "$data/proc"
global dt "$data/temp"
global do "$data/out"


*====================obtain ECDF==============================================*

*****ECDF
use $dd\estim.dta, clear
keep if bt==0
keep if soilinsecticide==0
keep if seedtrtrate==0
rename year year_trial
gen year=year_trial-1
keep if year>=2014 & year<=2016
keep year east rw

forv i=0/1{
preserve

cumul rw if east==`i', gen(cdf`i') 

keep rw cdf`i' 
keep if cdf`i'!=.
sort cdf`i'
gen cdf`i'_per=cdf`i'*100
tostring cdf`i'_per, force replace
split cdf`i'_per, gen(part) parse(".")
sort cdf`i'_per
bysort part1 (cdf`i'): gen is_minimum = cdf`i' == cdf`i'[1] //quantile level at interval=1
keep if is_minimum==1
keep rw part1 cdf`i'
destring part1, force replace
drop if part1==.
sort part1
rename part1 percentile 
rename rw rw`i'
save $dt\temp`i'.dta, replace
restore
}


clear
set obs 100
gen percentile=_n
merge 1:1 percentile using $dt\temp0.dta, nogen
merge 1:1 percentile using $dt\temp1.dta, nogen
for var cdf*: replace X=percentile/100 
export excel using "$dt\insert_v6.xls", firstrow(variables) replace
	//less than 100 observations for the eastern states, so rw0 values are missing for some percentile values; 
	//approximate the missing values based on the average distance per percentile between two nearest percentile rw0 values; saved in insert_ed.
import excel "$dt\insert_ed_v6.xls", sheet("Sheet1") firstrow clear
save $dt\ecdf.dta, replace
	//this data file is provided


*=====================compile parameters from dataset===========================*
*****maize cash price, state level
use $dd\estim.dta, clear
keep state* east
duplicates drop state, force
expand 18
bys state: gen year=2003+_n-1 //because trial field rw in 2004 represents the area rw in 2003.
save $dt\sample_statefips.dta, replace

use D:\Dropbox\1_Research\8_climate_biotech\1_data\temp\cash.dta, clear
//evelator-level maize cash price data; see "Data and materials availability" for data availability

duplicates drop year statefips, force
keep year statefips cash_statefips cash_na
merge 1:1 year statefips using $dt\sample_statefips.dta, keep(match using) nogen
save $dt\par_pcash.dta, replace
//aggregated state-level maize cash price data; see "Data and materials availability" for data availability

*****obtain yield potential, seed premium, Bt rate, and m.
use $dd\estim.dta, clear
keep if soilinsecticide==0 &seedtrtrate==0 

gen premium_acre_statefips=ps_acre_statefips1-ps_acre_statefips0
keep year east state* countyfips  bt rw premium_acre* *rw_rate_scrd* def16
collapse (mean) premium_acre* *rw_rate* rw def16, by(year east state* countyfips bt)
reshape wide rw, i(year east stateabbr statefips countyfips *rw_rate_scrd premium_acre_statefips def16) j(bt)

xtset countyfips year
for var rw0 rw1: rename X X_tminus1
gen rw0=f.rw0_tminus1
gen rw1=f.rw1_tminus1
drop rw*_tminus1

keep if year>=2014 & year<=2016
drop if rw0==. & rw1==.

replace rw1=rw0-${a}+${b}*lag1_rw_rate_scrd/100 if rw1==. 
replace rw0=rw1+${a}-${b}*lag1_rw_rate_scrd/100 if rw0==. 
for var rw0 rw1: replace X=0 if X<0

save $dt\threeyear.dta, replace



//match with county level yield.
import delimited $ds\yield_NASS\yield_2010to16_county.csv, clear
//county-level maize yield data from NASS, USDA; see "Data and materials availability" for data availability

keep if period=="YEAR"
keep year stateansi countyansi value
rename value yield
gen countyfips=stateansi*1000+countyansi
drop if countyfips==.
drop *ansi

merge 1:1 countyfips year using $dt\threeyear.dta, keep(match using) nogen

//recover yield potential
gen rw_rate=rw_rate_scrd/100
gen rw_loss=(rw_rate*rw1+(100-rw_rate)*rw0)/100 //state-average root injury weighted by bt area.
gen yield_potential=yield/(1-0.45*rw_loss)

preserve
collapse (mean) mean1= yield mean2=yield_potential ///
	(sd) sd1= yield sd2=yield_potential ///
	(count) n1=yield n2=yield_potential ///
	, by(state stateabbr statefips)
save $do\yield_state.dta, replace
restore

preserve
collapse (mean) mean1= yield mean2=yield_potential ///
	(sd) sd1= yield sd2=yield_potential ///
	(count) n1=yield n2=yield_potential ///
	, by(east)
tostring east, force replace
replace east="West" if east=="0"
replace east="East" if east=="1"
rename east state
save $do\yield_east.dta, replace
restore


//calculate c/m
merge m:1 year statefips using $dt\par_pcash.dta, keep(match master) nogen

replace cash_statefips=cash_na if cash_statefips==.
gen m=3*0.15*cash_statefips*yield_potential/def16
gen c=premium_acre_statefips/def16
gen cbym=c/m
collapse (mean) c m cash_statefips yield_potential cbym rw_rate, by(east)
sort east


//store c/m for east and west.
gl cbym0=cbym[1] //+${b}-${r} - the two terms are added in the old version for benefit calculation
gl cbym1=cbym[2] //+${b}-${r}
di "$cbym0, $cbym1"

save $dt\par.dta, replace


*------produce main model 1 iteration results using value iteration (table S14)----*
gl size=1
gl tol=1e-2 //convergence criterion


//define cbym parameter
use $dt/par.dta, clear
gl cbym0=cbym[1] 
gl cbym1=cbym[2] 

*------define mata functions outside of loops to avoid break------*
capture mata: mata drop funpar()
mata:
function funpar() {
    M = st_matrix("Mpar")
	R = M[1,.] // Realisations  in row 1
	P = M[2,.] // Probabilities in row 2
	n = ${size}
	D = rdiscrete(1,n,P)  // n discrete random variates from 1 to k
	V = R[1, D]           // V[i] = R[D[i]], 1 <= i <= n
	st_matrix("V", V)
    }
end

capture mata: mata drop fun0()
mata:
function fun0() {
    M = st_matrix("M0")
	R = M[1,.] // Realisations  in row 1
	P = M[2,.] // Probabilities in row 2
	n = ${size}
	D = rdiscrete(1,n,P)  // n discrete random variates from 1 to k
	V = R[1, D]           // V[i] = R[D[i]], 1 <= i <= n
	st_matrix("V", V)
    }
end

capture mata: mata drop fun1()
mata:
function fun1() {
    M = st_matrix("M1")
	R = M[1,.] // Realisations  in row 1
	P = M[2,.] // Probabilities in row 2
	n = ${size}
	D = rdiscrete(1,n,P)  // n discrete random variates from 1 to k
	V = R[1, D]           // V[i] = R[D[i]], 1 <= i <= n
	st_matrix("V", V)
    }
end




*-------------------random draws from abr parameters---------------*
set seed 123 //re run the program for every u value, so that the simulation for each u starts with the same seed

quietly{
qui use $dd/bootstrap_abr_5000itr.dta, clear
gen p=0.0002

recast double p
egen double check=total(p)
qui replace p=p/check //to ensure that p adds up to 1 (otherwise mata would report error)
cap drop check
egen double check=total(p)

local parlist "a r b"
qui foreach par of local parlist{
preserve
mkmat `par' p, mat(Mpar) 
mat Mpar=Mpar'
mata funpar()

mat Vp=V'
clear
svmat Vp, names(col)
rename r1 `par'
gen drawid=_n
save $dt/draw_`par'.dta, replace
restore
}



*-------------------------random draws from z cdf---------------------*
//west
qui use $dt/ecdf.dta, clear
gen p=0.01

recast double p
egen double check=total(p)
replace p=p/check //to ensure that p adds up to 1 (otherwise mata would report error)
cap drop check
egen double check=total(p)

mkmat rw0 p, mat(M0) 
mat M0=M0'

mata fun0()

mat Vp=V'
clear
svmat Vp, names(col)
rename r1 z_r0

gen drawid=_n
save $dt/draw_z_r0.dta, replace

/*
cumul z_r0, gen(cdf_check) //cross-check the cdf
sort cdf_check
twoway connected cdf_check z_r0 
*/

//east
qui use $dt/ecdf.dta, clear
gen p=0.01

recast double p
egen double check=total(p)
replace p=p/check 
cap drop check
egen double check=total(p)

mkmat rw1 p, mat(M1) 
mat M1=M1'

mata fun1()

mat Vp=V'
clear
svmat Vp, names(col)
rename r1 z_r1

gen drawid=_n
save $dt/draw_z_r1.dta, replace

/*
cumul z_r1, gen(cdf_check) //cross-check the cdf
sort cdf_check
twoway connected cdf_check z_r1
*/
}


//compile randow draw datasets
use $dt/draw_a.dta, clear
merge 1:1 drawid using $dt/draw_r.dta, nogen
merge 1:1 drawid using $dt/draw_b.dta, nogen
merge 1:1 drawid using $dt/draw_z_r0.dta, nogen
merge 1:1 drawid using $dt/draw_z_r1.dta, nogen



//expand dataset to accommodate st-1 values (0 and 1)
expand 2
bysort drawid: gen stm1=_n-1 //st-1=0 or st-1=1
sort stm1 drawid 
order stm1 drawid 

expand 100
bysort drawid stm1: gen P100=_n
bysort drawid stm1: gen P=(_n)/100
save $dt/base_main.dta, replace

clear
local vlist "rid east interest tao100 u100 pr0 pr1 P k Pequi EV0 EV1"
foreach v of local vlist{
	if "`v'"=="interest"{
		gen `v'=""
	}
	else{
		gen `v'=.
	}
}
save $dt/equi_infinite_main_1iter.dta, replace

local intlist "private group"
foreach j of numlist  0 1 {
foreach int of local intlist{
foreach tao100 of numlist 85 90 95 {
foreach u100 of numlist   5 25 50 {

gl tao=`tao100'/100
gl u=`u100'/100
gl j=`j'

use $dt/base_main.dta, clear
//add spillover effects
gen rtilda=r*${u}/(1-${u})
gen btilda=b*${u}/(1-${u})

//initialize 
gen V=0 //initial value
gen max_value = -1e10  // Initialize with a large negative number
gen st_optimal=.
gen value0=.
gen value1=.

* Value iteration loop
local diff = 1
local iter = 1
gen diff=1

cd $dt
while `diff' > $tol {
qui{
	cap drop EV EV0 EV1 intflag
	bys stm1 P: egen EV=mean(V)
	
	preserve
	keep stm1 P100 P EV
	collapse (mean) EV, by(stm1 P100)
	reshape wide EV, i(P100) j(stm1)
	tempfile temp
	save `temp', replace 
	restore
	merge m:1 P100 using `temp', nogen
	
	gen intflag=("`int'"=="group")
	replace value0=min(a*0+r*stm1-b*stm1*0+rtilda*P-btilda*P*0, z_r${j}) -${cbym${j}}*0 - (btilda-rtilda)*stm1*intflag +${tao} * EV0  //updated value for st=0
	replace value1=min(a*1+r*stm1-b*stm1*1+rtilda*P-btilda*P*1, z_r${j}) -${cbym${j}}*1 - (btilda-rtilda)*stm1*intflag +${tao} * EV1  //updated value for st=1
	replace st_optimal=1 if value1>=value0
	replace st_optimal=0 if value1<value0
	replace max_value=value1 if  value1>=value0
	replace max_value=value0 if  value1<value0

	//compare the updated value with initial value
	replace diff=abs(V-max_value)
	replace V=max_value
	sum diff
	loc diff = r(max)
}
    //di "Iteration `iter': max diff = `diff'"
    local iter = `iter' + 1

}

//display results
qui{
collapse (mean) st_optimal EV0 EV1, by(stm1 P100 P)
reshape wide st_optimal, i(P100 P) j(stm1)
rename st_optimal0 pr0
rename st_optimal1 pr1
gen k=P*pr1+(1-P)*pr0
gen Pdiff=abs(P-k)
sort Pdiff
loc Pequi = P in 1

gen Pequi=`Pequi'
gen rid=1
gen east=`j'
gen interest="`int'"
gen tao100=`tao100'
gen u100=`u100'

keep rid east interest tao100 u100 pr0 pr1 P k Pequi EV0 EV1
order rid east interest tao100 u100 pr0 pr1  P k Pequi  EV0 EV1

append using equi_infinite_main_1iter.dta
save equi_infinite_main_1iter.dta, replace
}

di "at east=`j', interest=`int' tao100=`tao100', u100=`u100': the optimal Bt planting rate is `Pequi'"

}
}
}
}

use $dt/equi_infinite_main_1iter.dta, clear
save  $do/equi_infinite_main_1iter.dta, replace


***export tables (table S14)
use equi_infinite_main_1iter.dta, clear
collapse (mean) pr0 pr1 P k Pequi EV0 EV1, by(east interest tao100 u100)
sort east interest tao u Pequi 
asdoc list east interest tao u Pequi   , save(Pequi_main_1iter.doc) replace dec(2) tzok 

//the table produced from R=1 is used to make table S14



