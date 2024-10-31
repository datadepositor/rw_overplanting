********************************************************************************
*				program: simulation, crop mix two-part model		     	   *
********************************************************************************
global root "..."

global data "$root/1_data"
global ds "$data/source"
global dd "$data/data"
global dp "$data/proc"
global dt "$data/temp"
global do "$data/out"



*------crop mix two-part model results using value iteration----------------------------*
gl R=1
gl size=100
gl tol=1e-2 //convergence criterion


*=============step 1: obtain distributions============*
*(since the parameter distribution has been obtained in the main analysis)
*(here only need to obtain the mono-mix version of ecdf, i.e., divide ECDF into mixd and monoinuous corn sections)
*****calibrate the mix rate M for east and west that corresponds to the regional ecdf;
//weighted by the number of non-bt trial samples to reflect the sampling weight in the regional ecdf.
//mix rate M is the percent acres of corn under crop mix system in total corn acres
//M is calibrated using the rotation rate


use $dd/estim_v7_use, clear //the cdl data are at county level
keep if bt==0
keep if year>=2014 & year<=2016

collapse (mean) lag1_cornshare_corn (count) sampleno=bt, by(countyfips east stateabbr statefips)

collapse (mean) lag1_cornshare_corn [fw=sampleno], by(east)

sort east

gl M0=100-round(lag1_cornshare_corn[1],1) //mix rate in the west, calibrated using the rotation rate in the west
gl M1=100-round(lag1_cornshare_corn[2],1) //mix rate in the east, calibrated using the rotation rate in the east
di "$M0, $M1"

*****divide the regional ecdf into mixd and monoinous corn sections
forv j=0/1{

use $dt/ecdf.dta, clear
keep percentile *`j'

//obtain the ecdf for the mixd corn section
preserve
keep in 1/${M`j'}
gen cdf`j'_mix=100*cdf`j'/${M`j'}
keep rw`j' cdf`j'_mix
gen percentile_mix=round(cdf`j'_mix*100,1)
save $dt/temp_mix.dta, replace
restore

//obtain the ecdf for the monoinous corn section
preserve
loc rown=${M`j'}+1
keep in `rown'/100
gen cdf`j'_temp=cdf`j'-${M`j'}/100
gen cdf`j'_mono=100*cdf`j'_temp/(100-${M`j'})
keep rw`j' cdf`j'_mono
gen percentile_mono=round(cdf`j'_mono*100,1)
save $dt/temp_mono.dta, replace
restore


//approximate the missing rw value using the next percentile value
clear
set obs 100
gen percentile_mix=_n
merge 1:1 percentile_mix using $dt/temp_mix.dta, nogen
tsset percentile_mix
replace rw`j'=f.rw`j' if rw`j'==.
replace rw`j'=f.rw`j' if rw`j'==.
replace rw`j'=f.rw`j' if rw`j'==.

rename percentile_mix percentile
rename rw`j' rw`j'_mix
save $dt/ecdf`j'_mix.dta, replace

clear
set obs 100
gen percentile_mono=_n
merge 1:1 percentile_mono using $dt/temp_mono.dta, nogen
tsset percentile_mono

replace rw`j'=f.rw`j' if rw`j'==.
replace rw`j'=f.rw`j' if rw`j'==.
replace rw`j'=f.rw`j' if rw`j'==.

rename percentile_mono percentile
rename rw`j' rw`j'_mono
save $dt/ecdf`j'_mono.dta, replace
}

use  $dt/ecdf0_mix.dta, clear
merge 1:1 percentile using  $dt/ecdf0_mono.dta, nogen
merge 1:1 percentile using  $dt/ecdf1_mix.dta, nogen
merge 1:1 percentile using  $dt/ecdf1_mono.dta, nogen
save $do/ecdf_mix.dta, replace






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

capture mata: mata drop fun0_mix()
mata:
function fun0_mix() {
    M = st_matrix("M0_mix")
	R = M[1,.] // Realisations  in row 1
	P = M[2,.] // Probabilities in row 2
	n = ${size}
	D = rdiscrete(1,n,P)  // n discrete random variates from 1 to k
	V = R[1, D]           // V[i] = R[D[i]], 1 <= i <= n
	st_matrix("V", V)
    }
end

capture mata: mata drop fun0_mono()
mata:
function fun0_mono() {
    M = st_matrix("M0_mono")
	R = M[1,.] // Realisations  in row 1
	P = M[2,.] // Probabilities in row 2
	n = ${size}
	D = rdiscrete(1,n,P)  // n discrete random variates from 1 to k
	V = R[1, D]           // V[i] = R[D[i]], 1 <= i <= n
	st_matrix("V", V)
    }
end

capture mata: mata drop fun1_mix()
mata:
function fun1_mix() {
    M = st_matrix("M1_mix")
	R = M[1,.] // Realisations  in row 1
	P = M[2,.] // Probabilities in row 2
	n = ${size}
	D = rdiscrete(1,n,P)  // n discrete random variates from 1 to k
	V = R[1, D]           // V[i] = R[D[i]], 1 <= i <= n
	st_matrix("V", V)
    }
end

capture mata: mata drop fun1_mono()
mata:
function fun1_mono() {
    M = st_matrix("M1_mono")
	R = M[1,.] // Realisations  in row 1
	P = M[2,.] // Probabilities in row 2
	n = ${size}
	D = rdiscrete(1,n,P)  // n discrete random variates from 1 to k
	V = R[1, D]           // V[i] = R[D[i]], 1 <= i <= n
	st_matrix("V", V)
    }
end






set seed 123 //re run the program for every u value, so that the simulation for each u starts with the same seed

quietly{
*-------------------random draws from abr parameters---------------*
use $dd/bootstrap_abr_5000itr.dta, clear
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
*****define mata function to avoid the problem of "end mata" causes "break" of looping.
local clist "mix mono"
qui forv j=0/1{
	foreach  c of local clist{

	qui use $do/ecdf_mix.dta, clear
	gen p=0.01

	recast double p
	egen double check=total(p)
	qui replace p=p/check //to ensure that p adds up to 1 (otherwise mata would report error)
	cap drop check
	egen double check=total(p)

	mkmat rw`j'_`c' p, mat(M`j'_`c') 
	mat M`j'_`c'=M`j'_`c''

	mata fun`j'_`c'()

	mat Vp=V'
	clear
	svmat Vp, names(col)
	rename r1 z_r`j'_`c'

	gen drawid=_n
	qui save $dt/draw_z_r`j'_`c'.dta, replace
}
}
}


*------------simulate the probablities and equilibriums-------------------*
//define cbym parameter
qui use $dt/par.dta, clear
gl cbym0=cbym[1] 
gl cbym1=cbym[2]  


//compile randow draw datasets
qui use $dt/draw_a.dta, clear
qui merge 1:1 drawid using $dt/draw_r.dta, nogen
qui merge 1:1 drawid using $dt/draw_b.dta, nogen
qui merge 1:1 drawid using $dt/draw_z_r0_mono.dta, nogen
qui merge 1:1 drawid using $dt/draw_z_r0_mix.dta, nogen
qui merge 1:1 drawid using $dt/draw_z_r1_mono.dta, nogen
qui merge 1:1 drawid using $dt/draw_z_r1_mix.dta, nogen
order drawid


//expand dataset to accommodate st-1 values (0 and 1)
expand 2
bysort drawid: gen stm1=_n-1 //st-1=0 or st-1=1
sort stm1 drawid 
order stm1 drawid 

expand 100
bysort drawid stm1: gen P100_mono=_n
bysort drawid stm1: gen P_mono=(_n)/100

expand 100
bysort drawid stm1 P100_mono: gen P100_mix=_n
bysort drawid stm1 P100_mono: gen P_mix=(_n)/100

qui save $dt/base_mix.dta, replace

//empty dataset to store the equilibrium results
clear
local vlist "rid east interest tao u delta Pequi_mono Pequi_mix Pequi"
foreach v of local vlist{
	if "`v'"=="interest"{
		gen `v'=""
	}
	else{
		gen `v'=.
	}
}
save $dt/Pequi_infinite_mix, replace

//loop over all parameter combinations

loc intlist "private group"
foreach delta of numlist 0 0.25 0.5 {
forv j=0/1{
foreach tao100 of numlist 95{
foreach u100 of numlist 5 25 50 {
foreach int of local intlist{

gl delta=`delta'
gl tao=`tao100'/100
gl u=`u100'/100
gl j=`j'

use $dt/base_mix.dta, clear
gen P=(100-${M${j}})*(P_mono)/100+${M${j}}*P_mix/100

//add spillover effects
gen rtilda=r*${u}/(1-${u})
gen btilda=b*${u}/(1-${u})

//initialize 
loc type "mono mix"
foreach t of loc type{
gen V_`t'=0 //initial value
gen max_value_`t' = -1e10  // Initialize with a large negative number
gen st_optimal_`t'=.
gen value0_`t'=.
gen value1_`t'=.

* Value iteration loop
gen diff_`t'=1
}

cd $dt

local diff = 1
local iter = 1
while `diff' > $tol {
qui{
	cap drop EV* intflag diff
	
	bys stm1 P_mono P_mix: egen EV_mono=mean(V_mono)
	bys stm1 P_mono P_mix: egen EV_mix=mean(V_mix)
	
	preserve
	keep stm1 P100* P* EV*
	collapse (mean) EV*, by(stm1 P100*)
	reshape wide EV*, i(P100*) j(stm1)
	save temp_mix.dta, replace 
	//sleep 500
		//increase the sleep time is error occurs: file temp.dta cannot be modified or erased; likely cause is read-only directory or file
	restore
	merge m:1 P100* using temp_mix.dta, nogen
	
	//gen intflag=1
	gen intflag=("`int'"=="group")
	replace value0_mono=min(r*stm1+rtilda*P, z_r${j}_mono) - (btilda-rtilda)*stm1*intflag +${tao} * EV_mono0  //updated value for st=0
	replace value1_mono=min(a+r*stm1-b*stm1+rtilda*P-btilda*P, z_r${j}_mono) -${cbym${j}} - (btilda-rtilda)*stm1*intflag +${tao} * EV_mono1  //updated value for st=1
	replace value0_mix=min(r*stm1+${delta}*rtilda*P, z_r${j}_mix) - ${delta}*(btilda-rtilda)*stm1*intflag +${tao} * EV_mix0  //updated value for st=0
	replace value1_mix=min(a+r*stm1-b*stm1+${delta}*(rtilda*P-btilda*P) , z_r${j}_mix) -${cbym${j}}*1 - ${delta}*(btilda-rtilda)*stm1*intflag +${tao} * EV_mix1  //updated value for st=1
	
	replace st_optimal_mono=1 if value1_mono>=value0_mono
	replace st_optimal_mono=0 if value1_mono<value0_mono
	replace max_value_mono=value1_mono if  value1_mono>=value0_mono
	replace max_value_mono=value0_mono if  value1_mono<value0_mono

	replace st_optimal_mix=1 if value1_mix>=value0_mix
	replace st_optimal_mix=0 if value1_mix<value0_mix
	replace max_value_mix=value1_mix if  value1_mix>=value0_mix
	replace max_value_mix=value0_mix if  value1_mix<value0_mix


	//compare the updated value with initial value
	replace diff_mono=abs(V_mono-max_value_mono)
	replace diff_mix=abs(V_mix-max_value_mix)
	egen diff=rowmax(diff_mono diff_mix)
	sum diff
	loc diff = r(max)
	
	replace V_mono=max_value_mono
	replace V_mix=max_value_mix

	
}
    //di "Iteration `iter': max diff = `diff'"
    local iter = `iter' + 1

}

//display results
qui{
collapse (mean) st_optimal*, by(stm1 P100* P*)
reshape wide st_optimal*, i(P100_mono P100_mix P_mono P_mix P) j(stm1)

rename st_optimal_mono0 pr_mono0 //prob of st=1 conditional on stm1=0
rename st_optimal_mono1 pr_mono1 //prob of st=1 conditional on stm1=1
rename st_optimal_mix0 pr_mix0
rename st_optimal_mix1 pr_mix1

gen k_mono=P_mono*pr_mono1+(1-P_mono)*pr_mono0
gen k_mix=P_mix*pr_mix1+(1-P_mix)*pr_mix0
gen Pdiff_mono=abs(P_mono-k_mono)
gen Pdiff_mix=abs(P_mix-k_mix)
gen Pdiff=Pdiff_mono+Pdiff_mix
sort Pdiff

keep in 1
loc Pequi=P in 1

gen Pequi_mono=P_mono in 1
gen Pequi_mix=P_mix in 1
gen Pequi=P in 1
gen rid=1
gen east=${j}
gen interest="`int'"
gen tao=${tao}
gen u=${u}
gen delta=${delta}

keep rid east interest tao u delta Pequi_mono Pequi_mix Pequi
order rid east interest tao u delta Pequi_mono Pequi_mix Pequi

append using Pequi_infinite_mix.dta
save Pequi_infinite_mix.dta, replace


}

di "east=`j', int=`int', tao100=`tao100', u100=`u100': Pequi=`Pequi'"

}
}
}
}
}


use $dt/Pequi_infinite_mix.dta, clear
save $do/Pequi_infinite_mix.dta, replace


*---------------------------export tables (table S16)-----------------------------------*
*****tables
cd $do

use Pequi_infinite_mix.dta, clear
sort east   u interest delta Pequi* 
asdoc list east   u interest delta Pequi*  , save(Pequi_mix.doc) replace dec(2) tzok 








