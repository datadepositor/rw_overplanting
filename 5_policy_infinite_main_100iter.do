********************************************************************************
*				program: simulation, main model, 100 iterations    			   *
********************************************************************************
global root "..."

global data "$root/1_data"
global ds "$data/source"
global dd "$data/data"
global dp "$data/proc"
global dt "$data/temp"
global do "$data/out"



*------produce main results using value iteration (Fig.4; table S13)----------------------------*
gl R=100
gl size=100
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
save $dt/equi_infinite_main_r`rid'.dta, replace

local intlist "private group"
foreach j of numlist  0 1 {
foreach int of local intlist{
foreach tao100 of numlist  85 90 95 {
foreach u100 of numlist  5 25 50 {

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

append using equi_infinite_main_r`rid'.dta
save equi_infinite_main_r`rid'.dta, replace
}


}
}
}
}

}

loc R=${R}
cd $dt
use equi_infinite_main_r1.dta, clear
forv rid=2/`R'{
	append using equi_infinite_main_r`rid'.dta
}
save  $do/equi_infinite_main_100iter.dta, replace



//export table S13
cd $do
use equi_infinite_main_100iter.dta, clear
collapse (mean) Pequi EV0 EV1, by(east interest tao100 u100)
sort east interest tao u Pequi 
asdoc list east interest tao u Pequi , save(Pequi.doc) replace dec(2) tzok 
//the table produced from R=100 is used to make table S13, the "optimum Bt planting rate" column


*---------------------------export figures-----------------------------------*
*------Fig.3 regional heterogeneity figure----------------------------*
graph set window fontface "Helvetica"
global psize "10pt" //panel text size
global tsize "7pt" //figure test size

//CDF figure
use $dt/ecdf.dta, clear
graph twoway (line cdf0 rw0 , sort(cdf0) lc("162 92 115") lp(longdash)  lw(0.6)) ///
	(line cdf1 rw1 , sort(cdf1) lc("119 143 165") lp(shortdash)  lw(0.6)), ///
	xlabel(0 (0.2) 1.0, format(%03.1f) labsize("$tsize"))  ylabel(0 (0.2) 1.0, format(%03.1f)  labsize("$tsize")) ///
	leg(region(lwidth(none)) on order(1 2) ring(0) position(5) col(1)  size("$tsize") label(1 "West") label(2 "East" )) ///
	xtitle("{bf:Root injury}", size("$tsize") margin(t=3)) /// 
	ytitle("{bf:Cumulative distribution function}", margin(r=3) size("$tsize"))  ///
	title("{bf:A}",size("8pt") margin(b=3) ring(1) position(11)) ///
	saving($do/ecdf.gph, replace) 


//yield potential figure
use $do/yield_state.dta, clear
append using $do/yield_east.dta
save $dt/temp.dta, replace

keep if state=="West"|state=="East"

forv i=1/2{
	gen hi`i'=mean`i'+invttail(n`i'-1,0.025)*(sd`i' / sqrt(n`i'))
	gen lo`i'=mean`i'-invttail(n`i'-1,0.025)*(sd`i' / sqrt(n`i'))
}

reshape long mean sd n hi lo, i(state) j(yieldvar) //yieldvar=1 for yield, =2 for yield_potential
sort statefips

save $dt/temp.dta, replace


loc bwidth=3
loc xlab1=1+`bwidth'/2
loc xlab2=`xlab1'+`bwidth'
sort state yieldvar
loc m1=round(mean[1],0.1) //yield, east
if round(`m1', 0.1) == int(`m1') {
    local m1 `m1'.0
}
di "`m1'"
loc m2=round(mean[2],0.1) //yield potential, east
loc m3=round(mean[3],0.1) //yield, west
loc m4=round(mean[4],0.1) //yield potential, west

loc d1=round(mean[2]-mean[1],0.1)
di "`d1'"
loc d0=round(mean[4]-mean[3],0.1)
di "`d0'"

use $dt/temp.dta, clear
gen statevar=yieldvar+1 if state=="West"
replace statevar=yieldvar+1+`bwidth' if state=="East"
graph twoway ///
	(bar mean statevar if yieldvar==1 & state=="West",  barw(1) fcolor("195 149 164") lcolor("162 92 115") lwidth(0.6)) ///
	(bar mean statevar if yieldvar==1 & state=="East",  barw(1) fcolor("169 184 198") lcolor("119 143 165") lwidth(0.6) ) ///
	(bar mean statevar if yieldvar==2 & state=="West", barw(1) fcolor("250 230 240") lcolor("195 149 164") lwidth(0.6)) ///
	(bar mean statevar if yieldvar==2 & state=="East", barw(1) fcolor("231 239 250")  lcolor("169 184 198") lwidth(0.6) ) ///
	(rcap hi lo statevar if yieldvar==1 & state=="West", color("162 92 115") lwidth(0.6)) ///
	(rcap hi lo statevar if yieldvar==2 & state=="West", color("195 149 164") lwidth(0.6)) ///
	(rcap hi lo statevar if yieldvar==1 & state=="East", color("119 143 165") lwidth(0.6)) ///
	(rcap hi lo statevar if yieldvar==2 & state=="East", color("169 184 198") lwidth(0.6)) ///
	, ///
	text(285 2.5 "{bf:Loss=`d0'}",color("162 92 115") size("$tsize") box bcolor("250 230 240") margin(medsmall)) ///
	text(285 5.5 "{bf:Loss=`d1'}",color("119 143 165") size("$tsize") box bcolor("231 239 250") margin(medsmall)) ///
	text(200 1.8 "R=`m3'", color("162 92 115") size("$tsize")) ///
	text(257 3.0 "P=`m4'", color("195 149 164") size("$tsize")) ///
	text(202 4.8 "R=`m1'", color("119 143 165") size("$tsize")) ///	
	text(220 6.0 "P=`m2'", color("169 184 198") size("$tsize")) ///	
	xscale(range(1 7)) yscale(range(100 300)) ///
	legend(region(lstyle(none)) off) ///
	xlabel( `xlab1' "West"  `xlab2' "East",labsize("$tsize") ) ///
	ylabel(100 150 200 250,labsize("$tsize") ) ///
	xtitle("{bf:Region}", size("$tsize") margin(t=3)) /// 
	ytitle("{bf:Yield (bushels/acre)}", size("$tsize") margin(r=2)) ///
	title("{bf:B}",size("8pt") margin(b=3) ring(1) position(11)) ///
	saving($do/yield_compare.gph, replace)

cd $do
graph combine ecdf.gph yield_compare.gph, row(1)  iscale(1) ///
	xsize(4.75) ysize(2.5)
	
graph export $do/Fig3.tif, replace width(2000) 
graph export $do/Fig3.png, replace width(2000)



*------Fig.4 all-in-one figure----------------------------*

set scheme s1color
graph set window fontface "Helvetica"
global psize "10pt" //panel text size
global tsize "7pt" //figure test size
loc lw 0.3

import excel "$root/4_submission/R2/241023_costbenefit_analysis.xlsx", sheet("figure") firstrow clear
drop bstar bbar
replace bdiff=0 if bdiff==.

reshape wide bdiff pstar pbar, ///
	i(u100 tao interest) j(east) 
reshape wide bdiff* pstar* pbar*, ///
	i(u100 tao) j(interest) string		
gen tao_private=tao-0.008
gen tao_group=tao+0.008
gen zero=0


//west
sort u100 tao
local bar=pbar0private[1]
local up1=pstar0private[1]
local up2=pstar0private[2]
local up3=pstar0private[3]
local lp1=pstar0private[7]
local lp2=pstar0private[8]
local lp3=pstar0private[9]

local ug1=pstar0group[1]
local ug2=pstar0group[2]
local ug3=pstar0group[3]
local lg1=pstar0group[7]
local lg2=pstar0group[8]
local lg3=pstar0group[9]

twoway	(scatteri `bar' 0.80 `bar' 1.0, recast(line) lpattern(shortdash) lcolor(gray)) ///
		(scatteri `up1' 0.842 `lp1' 0.842,recast(line)  lc(black) lw(`lw')) ///
		(scatteri `ug1' 0.858 `lg1' 0.858,recast(line)  lc(red) lw(`lw')) ///
		(scatteri `up2' 0.892 `lp2' 0.892,recast(line) lc(black) lw(`lw')) ///
		(scatteri `ug2' 0.908 `lg2' 0.908,recast(line) lc(red) lw(`lw')) ///
		(scatteri `up3' 0.942 `lp3' 0.942,recast(line) lc(black) lw(`lw')) ///
		(scatteri `ug3' 0.958 `lg3' 0.958,recast(line) lc(red) lw(`lw')) ///
		(scatter pstar0private tao_private if u100==5, msymbol(circle) mfcolor(white) lcolor(black) mcolor(black) msize(medium) mlw(`lw') ) ///
		(scatter pstar0private tao_private if u100==25, msymbol(X) mfcolor(white) lcolor(black) mcolor(black) mlw(`lw')) ///
		(scatter pstar0private tao_private if u100==50, msymbol(smdiamond) mfcolor(white) lcolor(black) mcolor(black) mlw(`lw'))  ///
		(scatter pstar0group tao_group if u100==5, msymbol(circle) mfcolor(white)  lcolor(red) mcolor(red) mlw(`lw')) ///
		(scatter pstar0group tao_group if u100==25, msymbol(X) mfcolor(white) lcolor(red) mcolor(red) mlw(`lw')) ///
		(scatter pstar0group tao_group if u100==50, msymbol(smdiamond) mfcolor(white) lcolor(red) mcolor(red) mlw(`lw')) ///
		(rbar  zero bdiff0private tao_private if u100==50, barw(0.012) fcolor(black%10) lw(0) yaxis(2)) ///	
		(scatter bdiff0private tao_private if u100==50, msymbol(smtriangle) mcolor(black) msize(vsmall) yaxis(2) ) ///
		(rbar  zero bdiff0group tao_group if u100==50, barw(0.012) fcolor(red%20)  lw(0) yaxis(2)) ///
		(scatter bdiff0group tao_group if u100==50, msymbol(smtriangle) mcolor(red) msize(vsmall) yaxis(2) ) ///
		, ///
		yline(0, lcolor(black)) ///
		xlabel(0.85 0.90 0.95, format(%03.2f) labsize("$tsize")) ///
		ylabel(0.0 (0.2) 0.8, format(%03.1f) labsize("$tsize")) ///
		ylabel(none, axis(2)) ///
		yscale(range(-0.4 0.8) axis(1)) yscale(range(0 900) axis(2)) ///
		leg(off) ///
		title("{bf:A. West}",size("8pt") margin(b=3) position(11)) ///
		ytitle("{bf:                           Bt planting rate}",size("$tsize") margin(r=3)) ///
		xtitle("{bf:Discount factor}",size("$tsize") margin(t=3)) ///
		saving($dt/Pwest.gph, replace) ///
		fxsize(80)
		
//east

local bar=pbar1private[1]
local up1=pstar1private[1]
local up2=pstar1private[2]
local up3=pstar1private[3]
local lp1=pstar1private[7]
local lp2=pstar1private[8]
local lp3=pstar1private[9]

local ug1=pstar1group[1]
local ug2=pstar1group[2]
local ug3=pstar1group[3]
local lg1=pstar1group[7]
local lg2=pstar1group[8]
local lg3=pstar1group[9]

twoway	(scatteri `bar' 0.80 `bar' 1.0, recast(line) lpattern(shortdash) lcolor(gray)) ///
		(scatteri `up1' 0.842 `lp1' 0.842,recast(line)  lc(black) lw(`lw')) ///
		(scatteri `ug1' 0.858 `lg1' 0.858,recast(line)  lc(red) lw(`lw')) ///
		(scatteri `up2' 0.892 `lp2' 0.892,recast(line) lc(black) lw(`lw')) ///
		(scatteri `ug2' 0.908 `lg2' 0.908,recast(line) lc(red) lw(`lw')) ///
		(scatteri `up3' 0.942 `lp3' 0.942,recast(line) lc(black) lw(`lw')) ///
		(scatteri `ug3' 0.958 `lg3' 0.958,recast(line) lc(red) lw(`lw')) ///
		(scatter pstar1private tao_private if u100==5, msymbol(circle) lw(vvthin) mfcolor(white) lcolor(black) mcolor(black)  mlw(`lw')) ///
		(scatter pstar1private tao_private if u100==25, msymbol(X) mfcolor(white) lcolor(black) mcolor(black) mlw(`lw') ) ///
		(scatter pstar1private tao_private if u100==50, msymbol(smdiamond) mfcolor(white) lcolor(black) mcolor(black)  mlw(`lw'))  ///
		(scatter pstar1group tao_group if u100==5, msymbol(circle) mfcolor(white)  lcolor(red) mcolor(red)  mlw(`lw')) ///
		(scatter pstar1group tao_group if u100==25, msymbol(X) mfcolor(white) lcolor(red) mcolor(red)  mlw(`lw')) ///
		(scatter pstar1group tao_group if u100==50, msymbol(smdiamond) mfcolor(white) lcolor(red) mcolor(red)  mlw(`lw')) ///
		(rbar  zero bdiff1private tao_private if u100==50, barw(0.012) fcolor(black%10) lw(0)  yaxis(2)) ///	
		(scatter bdiff1private tao_private if u100==5, msymbol(plus) mcolor(black) msize(small) yaxis(2) ) ///
		(scatter bdiff1private tao_private if u100==25, msymbol(P) mcolor(black) msize(tiny) yaxis(2)) ///
		(scatter bdiff1private tao_private if u100==50, msymbol(smtriangle) mcolor(black) msize(vsmall) yaxis(2) ) ///
		(rbar  zero bdiff1group tao_group if u100==50, barw(0.012) fcolor(red%20)  lw(0)  yaxis(2)) ///
		(scatter bdiff1group tao_group if u100==5, msymbol(plus) mcolor(red) msize(small) yaxis(2) ) ///
		(scatter bdiff1group tao_group if u100==25, msymbol(P) mcolor(red) msize(tiny) yaxis(2)) ///
		(scatter bdiff1group tao_group if u100==50, msymbol(smtriangle) mcolor(red) msize(vsmall) yaxis(2) ) ///
			, ///
		yline(0, lcolor(black)) ///
		xlabel(0.85 0.90 0.95, format(%03.2f) labsize("$tsize")) ///
		ylabel(none, axis(1)) ///
		ylabel(0 (100) 270, format(%2.0f) axis(2) labsize("$tsize")) ///
		yscale(range(-0.4 0.8) axis(1)) yscale(range(0 900) axis(2)) ///
		leg(region(lwidth(none) color(none)) on ring(0) position(6) col(1)  size("$tsize") ///
		order(1 "{bf:Status quo}" ///
		2 "{bf:The optimum under}" "{bf:a dispersal proportion of:}" ///
		8 "{it:u}=5%" 9 "{it:u}=25%" 10 "{it:u}=50%" ///
		14 "{bf:Benefits of shifting to}" "{bf:the optimum under}" "{bf:a dispersal proportion of:}" ///
		15 "{it:u}=5%" 16 "{it:u}=25%" 17 "{it:u}=50%")) ///
		title("{bf:B. East}",size("8pt") margin(b=3) position(11)) ///
		ytitle("{bf:Benefits ($/acre)                                            }",size("$tsize") margin(l=3) axis(2)) ///
		xtitle("{bf:Discount factor}",size("$tsize") margin(t=3)) ///
		saving($dt/Peast.gph, replace) ///
		fxsize(80)		
		
grc1leg2  $dt/Pwest.gph  $dt/Peast.gph  ,  rows(1) legendfrom( $dt/Peast.gph) ///
	position(3) ring(1) saving($do/onefigure.gph, replace)  imargin(tiny) ///
	noauto iscale(1) 
graph dis,	xsize(5) ysize(2.5)

graph use $do/onefigure.gph
graph export $do/Fig4.tif, replace width(2000) //horizontal pixels=4.75inch*300dots per inch
graph export $do/Fig4.png, replace width(2000)

*------fig. S5-7 equilibrium visualization figures----------------------------*
gl fpanel 10pt
gl faxis 8pt
gl marg 2

use $do/equi_infinite_main_100iter.dta, clear

collapse (mean) k, by(east interest tao100 u100 P)
sort east interest tao100 u100 P
save $dt/function.dta, replace

local intlist "private group"
foreach int of local intlist{
foreach j of numlist 0 1{
foreach tao100 of numlist 85 90 95 {
foreach u100 of numlist 5 25 50{
di "working on interest==`int' "
di "& tao100==`tao100' & u100==`u100'"
	use $dt/function.dta, clear
	keep if interest=="`int'" & tao100==`tao100' & u100==`u100'
	keep east P k 
	reshape wide k, i(P) j(east)
	
	graph twoway (line P P , sort lc(black) lw(0.5)) ///
		(line k1 P , sort lc(red*1.5)  lp(shortdash)  lw(0.5)) ///
		(line k0 P , sort lc(green*1.5) lp(longdash)  lw(0.5)) , ///
		xlabel(0 (0.2) 1.0, format(%03.1f) labsize(${faxis}))  ///
		ylabel(0 (0.2) 1.0, format(%03.1f) labsize(${faxis}))  ///
		title("{bf:u=`u100'%}", size(10pt) margin(b=${marg})) /// 
		leg(region(lwidth(none) color(none)) on order(3 2 1) ring(0) position(12) row(1)  size(7pt) label(1 "g(P)") label(2 "k(P) East") label(3 "k(P) West" )) ///
		xtitle("Bt planting rate (P)", size(${fpanel}) margin(t=${marg})) /// 
		ytitle("Functions", size(10pt) margin(r=${marg}))   ///
		saving($dt/functions_tao`tao100'_u`u100'_`int'.gph, replace) nodraw
}
}
}
}

cd $dt
foreach tao100 of numlist 85 90 95 {
grc1leg ///
	functions_tao`tao100'_u5_private.gph functions_tao`tao100'_u25_private.gph functions_tao`tao100'_u50_private.gph ///
	functions_tao`tao100'_u5_group.gph functions_tao`tao100'_u25_group.gph functions_tao`tao100'_u50_group.gph, ///
	row(2) legendfrom(functions_tao`tao100'_u5_private.gph) imargin(small)
graph display, xsize(5.5) ysize(4.2) 
graph export $do/equi_tao`tao100'.png, replace width(2000) 
graph export $do/equi_tao`tao100'.tif, replace width(2000) 

}	





