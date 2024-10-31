********************************************************************************
*				program: simulation, farm size two-part model		     	   *
********************************************************************************
global root "..."

global data "$root/1_data"
global ds "$data/source"
global dd "$data/data"
global dp "$data/proc"
global dt "$data/temp"
global do "$data/out"


*------farm size two-part model results using value iteration----------------------------*
gl R=1
gl size=100
gl tol=1e-2 //convergence criterion
gl L=50 


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

}


*------------simulate the probablities and equilibriums-------------------*
//define cbym parameter
qui use $dt/par.dta, clear
gl cbym0=cbym[1] 
gl cbym1=cbym[2]  

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
bysort drawid stm1: gen P100_small=_n
bysort drawid stm1: gen P_small=(_n)/100

expand 100
bysort drawid stm1 P100_small: gen P100_large=_n
bysort drawid stm1 P100_small: gen P_large=(_n)/100

save $dt/base_size.dta, replace


//empty dataset to store the equilibrium results
clear
local vlist "rid east interest tao u_large u_small  Pequi_small Pequi_large Pequi"
foreach v of local vlist{
	if "`v'"=="interest"{
		gen `v'=""
	}
	else{
		gen `v'=.
	}
}
save $dt/Pequi_infinite_size.dta, replace

//loop over all parameter combinations

loc intlist "private group"
forv j=0/1{
foreach int of local intlist{
foreach tao100 of numlist 95{
foreach u_small of numlist 0.30 0.50{
foreach u_large of numlist 0.05 0.15{

gl u_small=`u_small'
gl u_large=`u_large'
gl tao=`tao100'/100
gl j=`j'

use $dt/base_size.dta, clear

//add spillover effects
gen btilda_small=b*${u_small}/(1-${u_small})
gen btilda_large=b*${u_large}/(1-${u_large})

gen rtilda_small=r*${u_small}/(1-${u_small})
gen rtilda_large=r*${u_large}/(1-${u_large})

gen P=(100-${L})*(P_large)/100+${L}*P_small/100

//initialize 
loc type "small large"
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
	
	bys stm1 P_small P_large: egen EV_small=mean(V_small)
	bys stm1 P_small P_large: egen EV_large=mean(V_large)
	
	preserve
	keep stm1 P100* P* EV*
	collapse (mean) EV*, by(stm1 P100*)
	reshape wide EV*, i(P100*) j(stm1)
	save temp_size.dta, replace 
	restore
	merge m:1 P100* using temp_size.dta, nogen
	
	//gen intflag=1
	gen intflag=("`int'"=="group")
	replace value0_small=min(r*stm1+rtilda_small*P, z_r${j}) - (btilda_small-rtilda_small)*stm1*intflag +${tao} * EV_small0  //updated value for st=0
	replace value1_small=min(a+r*stm1-b*stm1+rtilda_small*P-btilda_small*P, z_r${j}) -${cbym${j}} - (btilda_small-rtilda_small)*stm1*intflag +${tao} * EV_small1  //updated value for st=1
	replace value0_large=min(r*stm1+rtilda_large*P, z_r${j}) - (btilda_large-rtilda_large)*stm1*intflag +${tao} * EV_large0  //updated value for st=0
	replace value1_large=min(a+(r*stm1-b*stm1)+(rtilda_large*P-btilda_large*P) , z_r${j}) -${cbym${j}} - (btilda_large-rtilda_large)*stm1*intflag +${tao} * EV_large1  //updated value for st=1
	
	replace st_optimal_small=1 if value1_small>=value0_small
	replace st_optimal_small=0 if value1_small<value0_small
	replace max_value_small=value1_small if  value1_small>=value0_small
	replace max_value_small=value0_small if  value1_small<value0_small

	replace st_optimal_large=1 if value1_large>=value0_large
	replace st_optimal_large=0 if value1_large<value0_large
	replace max_value_large=value1_large if  value1_large>=value0_large
	replace max_value_large=value0_large if  value1_large<value0_large


	//compare the updated value with initial value
	replace diff_small=abs(V_small-max_value_small)
	replace diff_large=abs(V_large-max_value_large)
	egen diff=rowmax(diff_small diff_large)
	sum diff
	loc diff = r(max)
	
	replace V_small=max_value_small
	replace V_large=max_value_large

	
}
    local iter = `iter' + 1

}

//display results
qui{
collapse (mean) st_optimal*, by(stm1 P100* P*)
reshape wide st_optimal*, i(P100_large P100_small P_large P_small P) j(stm1)

rename st_optimal_small0 pr_small0 //prob of st=1 conditional on stm1=0
rename st_optimal_small1 pr_small1 //prob of st=1 conditional on stm1=1
rename st_optimal_large0 pr_large0
rename st_optimal_large1 pr_large1

gen k_small=P_small*pr_small1+(1-P_small)*pr_small0
gen k_large=P_large*pr_large1+(1-P_large)*pr_large0
gen Pdiff_small=abs(P_small-k_small)
gen Pdiff_large=abs(P_large-k_large)
gen Pdiff=Pdiff_small+Pdiff_large
sort Pdiff

keep in 1
loc Pequi=P in 1

gen Pequi_small=P_small in 1
gen Pequi_large=P_large in 1
gen Pequi=P in 1
gen rid=1
gen east=${j}
gen interest="`int'"
gen tao=${tao}
gen u_large=${u_large}
gen u_small=${u_small}

keep rid east interest tao u_large u_small  Pequi_small Pequi_large Pequi
order rid east interest tao u_large u_small  Pequi_small Pequi_large Pequi

append using Pequi_infinite_size.dta
save Pequi_infinite_size.dta, replace


}
di "east=`j', int=`int', tao100=`tao100', u_small=`u_small', u_large=`u_large': Pequi=`Pequi'"

}
}
}
}
}


use $dt/Pequi_infinite_size.dta, clear
save $do/Pequi_infinite_size.dta, replace


*---------------------------export tables (table S17)-----------------------------------*
*****tables
cd $do

use Pequi_infinite_size.dta, clear
sort east interest  u*  Pequi*
asdoc list east   u* interest  Pequi*  , save(Pequi_size.doc) replace dec(2) tzok 


