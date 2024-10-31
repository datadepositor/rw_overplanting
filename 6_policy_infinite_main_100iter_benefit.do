********************************************************************************
*				program: simulation to obtain shifting benefit       		   *
********************************************************************************
global root "..."

global data "$root/1_data"
global ds "$data/source"
global dd "$data/data"
global dp "$data/proc"
global dt "$data/temp"
global do "$data/out"



*---simulate for status quo scenario by calibrating the seed premium adjustment parameter---------*
gl R=100
gl size=100
gl tol=1e-2 //convergence criterion
gl chypo=-26 //seed premium adjustment parameter

//define cbym parameter
use $dt/par.dta, clear
gen chypo=${chypo}
gl cbym0=(c[1]+chypo)/m[1]
gl cbym1=(c[2]+chypo)/m[2]

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
loc R=${R}
forv rid=1/`R'{
di "working on iteration `rid'"
quietly{
qui use $dd/bootstrap_abr_5000itr.dta, clear
gen p=0.0002

recast double p
egen double check=total(p)
qui replace p=p/check //to ensure that p adds up to 1 (otherwise mata would report error)
cap drop check
egen double check=total(p)

local parlist "a r b"
foreach par of local parlist{
preserve
mkmat `par' p, mat(Mpar) 
mat Mpar=Mpar'
mata funpar()

mat Vp=V'
clear
svmat Vp, names(col)
rename r1 `par'
gen drawid=_n

tempfile draw_`par'
save  `draw_`par'', replace



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
tempfile draw_z_r0
save  `draw_z_r0', replace


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
tempfile draw_z_r1
save  `draw_z_r1', replace


//compile randow draw datasets
use `draw_a', clear
merge 1:1 drawid using  `draw_r', nogen
merge 1:1 drawid using  `draw_b', nogen
merge 1:1 drawid using  `draw_z_r0', nogen
merge 1:1 drawid using `draw_z_r1', nogen

//di "it does not break here"


//expand dataset to accommodate st-1 values (0 and 1)
expand 2
bysort drawid: gen stm1=_n-1 //st-1=0 or st-1=1
sort stm1 drawid 
order stm1 drawid 

expand 100
bysort drawid stm1: gen P100=_n
bysort drawid stm1: gen P=(_n)/100
tempfile base_main
save `base_main', replace 


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

tempfile equi_infinite_main_r`rid'
save `equi_infinite_main_r`rid'', replace

local intlist "group "
foreach j of numlist  1  {
foreach int of local intlist{
foreach tao100 of numlist  95  {
foreach u100 of numlist   5  {

gl tao=`tao100'/100
gl u=`u100'/100
gl j=`j'

use `base_main', clear
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

    //di "Iteration `iter': max diff = `diff'"
    local iter = `iter' + 1

}


collapse (mean) st_optimal EV0 EV1, by(stm1 P100 P)
reshape wide st_optimal, i(P100 P) j(stm1)
rename st_optimal0 pr0
rename st_optimal1 pr1
gen k=P*pr1+(1-P)*pr0
gen Pdiff=abs(P-k)
sort Pdiff
loc Pequi = P in 1

gen Pequi=`Pequi'
gen rid=`rid'
gen east=`j'
gen interest="`int'"
gen tao100=`tao100'
gen u100=`u100'

keep rid east interest tao100 u100 pr0 pr1 P k Pequi EV0 EV1
order rid east interest tao100 u100 pr0 pr1  P k Pequi  EV0 EV1

append using `equi_infinite_main_r`rid''
save `equi_infinite_main_r`rid'', replace
}


}
}
}
}

}

loc R=${R}
cd $dt
use `equi_infinite_main_r1', clear
forv rid=2/`R'{
	append using `equi_infinite_main_r`rid''
}
save  $do/equi_infinite_main_100iter_benefit.dta, replace



*---------------------------calculate benefit-----------------------------------*
*****tables
cd $do
use equi_infinite_main_100iter_benefit.dta, clear
collapse (mean) Pequi EV0 EV1, by(east interest tao100 u100)
gen chypo=${chypo}

//results obtained here are used to produce table S13, "seed premium adjustment" and "benefit from shifting" columns








