********************************************************************************
*				program: simulation, crop rotation two-part model			   *
********************************************************************************
global root "..."

global data "$root/1_data"
global ds "$data/source"
global dd "$data/data"
global dp "$data/proc"
global dt "$data/temp"
global do "$data/out"


*------crop rotation two-part model results using value iteration----------------------------*
gl R=1
gl size=100
gl tol=1e-2 //convergence criterion


*=============step 1: obtain distributions============*
*(since the parameter distribution has been obtained in the main analysis)
*(here only need to obtain the rotation version of ecdf, i.e., divide ECDF into rotated and continuous corn sections)
*****calculate the rotation rate R for east and west that corresponds to the regional ecdf;
//weighted by the number of non-bt trial samples to reflect the sampling weight in the regional ecdf.
// rotation rate R is the percent of corn planted to non-corn in the previous year 
use $dd/estim_v7_use, clear //the cdl data are at county level
keep if bt==0
keep if year>=2014 & year<=2016

collapse (mean) lag1_cornshare_corn (count) sampleno=bt, by(countyfips east stateabbr statefips)

collapse (mean) lag1_cornshare_corn [fw=sampleno], by(east)

sort east
gl R0=100-round(lag1_cornshare_corn[1],1) //rotation rate in the west
gl R1=100-round(lag1_cornshare_corn[2],1) //rotation rate in the east
di "$R0, $R1"

*****divide the regional ecdf into rotated and continous corn sections
forv j=0/1{

use $dt/ecdf.dta, clear
keep percentile *`j'

//obtain the ecdf for the rotated corn section
preserve
keep in 1/${R`j'}
gen cdf`j'_rotate=100*cdf`j'/${R`j'}
keep rw`j' cdf`j'_rotate
gen percentile_rotate=round(cdf`j'_rotate*100,1)
save $dt/temp_rotate.dta, replace
restore

//obtain the ecdf for the continous corn section
preserve
loc rown=${R`j'}+1
keep in `rown'/100
gen cdf`j'_temp=cdf`j'-${R`j'}/100
gen cdf`j'_cont=100*cdf`j'_temp/(100-${R`j'})
keep rw`j' cdf`j'_cont
gen percentile_cont=round(cdf`j'_cont*100,1)
save $dt/temp_cont.dta, replace
restore


//approximate the missing rw value using the next percentile value
clear
set obs 100
gen percentile_rotate=_n
merge 1:1 percentile_rotate using $dt/temp_rotate.dta, nogen
tsset percentile_rotate
replace rw`j'=f.rw`j' if rw`j'==.
replace rw`j'=f.rw`j' if rw`j'==.
replace rw`j'=f.rw`j' if rw`j'==.

rename percentile_rotate percentile
rename rw`j' rw`j'_rotate
save $dt/ecdf`j'_rotate.dta, replace

clear
set obs 100
gen percentile_cont=_n
merge 1:1 percentile_cont using $dt/temp_cont.dta, nogen
tsset percentile_cont

replace rw`j'=f.rw`j' if rw`j'==.
replace rw`j'=f.rw`j' if rw`j'==.
replace rw`j'=f.rw`j' if rw`j'==.

rename percentile_cont percentile
rename rw`j' rw`j'_cont
save $dt/ecdf`j'_cont.dta, replace
}

use  $dt/ecdf0_rotate.dta, clear
merge 1:1 percentile using  $dt/ecdf0_cont.dta, nogen
merge 1:1 percentile using  $dt/ecdf1_rotate.dta, nogen
merge 1:1 percentile using  $dt/ecdf1_cont.dta, nogen
save $do/ecdf_rotation.dta, replace





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

capture mata: mata drop fun0_rotate()
mata:
function fun0_rotate() {
    M = st_matrix("M0_rotate")
	R = M[1,.] // Realisations  in row 1
	P = M[2,.] // Probabilities in row 2
	n = ${size}
	D = rdiscrete(1,n,P)  // n discrete random variates from 1 to k
	V = R[1, D]           // V[i] = R[D[i]], 1 <= i <= n
	st_matrix("V", V)
    }
end

capture mata: mata drop fun0_cont()
mata:
function fun0_cont() {
    M = st_matrix("M0_cont")
	R = M[1,.] // Realisations  in row 1
	P = M[2,.] // Probabilities in row 2
	n = ${size}
	D = rdiscrete(1,n,P)  // n discrete random variates from 1 to k
	V = R[1, D]           // V[i] = R[D[i]], 1 <= i <= n
	st_matrix("V", V)
    }
end

capture mata: mata drop fun1_rotate()
mata:
function fun1_rotate() {
    M = st_matrix("M1_rotate")
	R = M[1,.] // Realisations  in row 1
	P = M[2,.] // Probabilities in row 2
	n = ${size}
	D = rdiscrete(1,n,P)  // n discrete random variates from 1 to k
	V = R[1, D]           // V[i] = R[D[i]], 1 <= i <= n
	st_matrix("V", V)
    }
end

capture mata: mata drop fun1_cont()
mata:
function fun1_cont() {
    M = st_matrix("M1_cont")
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

/*
use $dt/draw_a.dta, clear
cumul a, gen(cdf) //cross-check the cdf
sort cdf
twoway connected cdf a 
*/

*-------------------------random draws from z cdf---------------------*
*****define mata function to avoid the problem of "end mata" causes "break" of looping.
local clist "rotate cont"
qui forv j=0/1{
	foreach  c of local clist{

	qui use $do/ecdf_rotation.dta, clear
	gen p=0.01

	recast double p
	egen double check=total(p)
	qui replace p=p/check //to ensure that p adds up to 1 (otherwise mata would report error)
	cap drop check
	egen double check=total(p)


	//rotate
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
qui merge 1:1 drawid using $dt/draw_z_r0_rotate.dta, nogen
qui merge 1:1 drawid using $dt/draw_z_r0_cont.dta, nogen
qui merge 1:1 drawid using $dt/draw_z_r1_rotate.dta, nogen
qui merge 1:1 drawid using $dt/draw_z_r1_cont.dta, nogen
order drawid


//expand dataset to accommodate st-1 values (0 and 1)
expand 2
bysort drawid: gen stm1=_n-1 //st-1=0 or st-1=1
sort stm1 drawid 
order stm1 drawid 

expand 101
bysort drawid stm1: gen P100_rotate=_n-1
bysort drawid stm1: gen P_rotate=(_n-1)/100

expand 101
bysort drawid stm1 P100_rotate: gen P100_cont=_n-1
bysort drawid stm1 P100_rotate: gen P_cont=(_n-1)/100

qui save $dt/base_rotation.dta, replace


//empty dataset to store the equilibrium results
clear
local vlist "rid east interest tao u delta Pequi_cont Pequi_rotate Pequi"
foreach v of local vlist{
	if "`v'"=="interest"{
		gen `v'=""
	}
	else{
		gen `v'=.
	}
}
save $dt/Pequi_infinite_rotation, replace

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

use $dt/base.dta, clear
gen P=(100-${R${j}})*(P_cont)/100+${R${j}}*P_rotate/100

//add spillover effects
gen rtilda=r*${u}/(1-${u})
gen btilda=b*${u}/(1-${u})

//initialize 
loc type "cont rotate"
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
	
	bys stm1 P_cont P_rotate: egen EV_cont=mean(V_cont)
	bys stm1 P_cont P_rotate: egen EV_rotate=mean(V_rotate)
	
	preserve
	keep stm1 P100* P* EV*
	collapse (mean) EV*, by(stm1 P100*)
	reshape wide EV*, i(P100*) j(stm1)
	save temp_rotation.dta, replace 
	//sleep 500
		//increase the sleep time is error occurs: file temp.dta cannot be modified or erased; likely cause is read-only directory or file
	restore
	merge m:1 P100* using temp_rotation.dta, nogen
	
	//gen intflag=1
	gen intflag=("`int'"=="group")
	replace value0_cont=min(r*stm1+rtilda*P, z_r${j}_cont) - (btilda-rtilda)*stm1*intflag +${tao} * EV_cont0  //updated value for st=0
	replace value1_cont=min(a+r*stm1-b*stm1+rtilda*P-btilda*P, z_r${j}_cont) -${cbym${j}} - (btilda-rtilda)*stm1*intflag +${tao} * EV_cont1  //updated value for st=1
	replace value0_rotate=min(${delta}*r*stm1+rtilda*P, z_r${j}_rotate) - (btilda-rtilda)*stm1*intflag +${tao} * EV_rotate0  //updated value for st=0
	replace value1_rotate=min(a+${delta}*(r*stm1-b*stm1)+(rtilda*P-btilda*P) , z_r${j}_rotate) -${cbym${j}} - (btilda-rtilda)*stm1*intflag +${tao} * EV_rotate1  //updated value for st=1
	
	replace st_optimal_cont=1 if value1_cont>=value0_cont
	replace st_optimal_cont=0 if value1_cont<value0_cont
	replace max_value_cont=value1_cont if  value1_cont>=value0_cont
	replace max_value_cont=value0_cont if  value1_cont<value0_cont

	replace st_optimal_rotate=1 if value1_rotate>=value0_rotate
	replace st_optimal_rotate=0 if value1_rotate<value0_rotate
	replace max_value_rotate=value1_rotate if  value1_rotate>=value0_rotate
	replace max_value_rotate=value0_rotate if  value1_rotate<value0_rotate


	//compare the updated value with initial value
	replace diff_cont=abs(V_cont-max_value_cont)
	replace diff_rotate=abs(V_rotate-max_value_rotate)
	egen diff=rowmax(diff_cont diff_rotate)
	sum diff
	loc diff = r(max)
	
	replace V_cont=max_value_cont
	replace V_rotate=max_value_rotate

	
}
    //di "Iteration `iter': max diff = `diff'"
    local iter = `iter' + 1

}

//display results
qui{
collapse (mean) st_optimal*, by(stm1 P100* P*)
reshape wide st_optimal*, i(P100_rotate P100_cont P_rotate P_cont P) j(stm1)

rename st_optimal_cont0 pr_cont0 //prob of st=1 conditional on stm1=0
rename st_optimal_cont1 pr_cont1 //prob of st=1 conditional on stm1=1
rename st_optimal_rotate0 pr_rotate0
rename st_optimal_rotate1 pr_rotate1

gen k_cont=P_cont*pr_cont1+(1-P_cont)*pr_cont0
gen k_rotate=P_rotate*pr_rotate1+(1-P_rotate)*pr_rotate0
gen Pdiff_cont=abs(P_cont-k_cont)
gen Pdiff_rotate=abs(P_rotate-k_rotate)
gen Pdiff=Pdiff_cont+Pdiff_rotate
sort Pdiff

keep in 1
loc Pequi=P in 1

gen Pequi_cont=P_cont in 1
gen Pequi_rotate=P_rotate in 1
gen Pequi=P in 1
gen rid=1
gen east=${j}
gen interest="`int'"
gen tao=${tao}
gen u=${u}
gen delta=${delta}

keep rid east interest tao u delta Pequi_cont Pequi_rotate Pequi
order rid east interest tao u delta Pequi_cont Pequi_rotate Pequi

append using Pequi_infinite_rotation.dta
save Pequi_infinite_rotation.dta, replace


}
di "east=`j', int=`int', tao100=`tao100', u100=`u100': Pequi=`Pequi'"

}
}
}
}
}


use $dt/Pequi_infinite_rotation.dta, clear
save $do/Pequi_infinite_rotation.dta, replace


*---------------------------export tables (table S15)-----------------------------------*
*****tables
cd $do

use Pequi_infinite_rotation.dta, clear
sort east   u interest delta Pequi*
asdoc list east   u interest delta Pequi*  , save(Pequi_rotation.doc) replace dec(2) tzok 


