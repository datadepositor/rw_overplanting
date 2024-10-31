********************************************************************************
*				program: estimation			   	   							   *
********************************************************************************

global root "..."

global data "$root/1_data"
global ds "$data/source"
global dd "$data/data"
global dp "$data/proc"
global dt "$data/temp"
global do "$data/out"

*====================2005-2016 main analysis data, plot level================================*
gl vce scrd 
gl rw rw
gl pattern pattern //crop rotation pattern

use $dd/estim.dta, clear
gen pattern=lag1_cornshare_corn/100 
gen hist_scrd=lag2_rw_rate_scrd/100
gen hist_countyfips=lag2_rw_rate_countyfips/100
gen hist2_scrd=lag3_rw_rate_scrd/100
gen hist2_countyfips=lag3_rw_rate_countyfips/100

keep if year>=2005 & year<=2016
drop if soilinsecticide==1 

egen miss=rowmiss($rw  bt hist_scrd prec4cm-prec8cm  tMin1 tMin2 tMin3 seedtrtrate  locid)
drop if miss

save $dd/estim_use.dta, replace

*--------------------regression analysis (table S3,6-12)--------------------------------------------------*
use $dd/estim_use.dta, clear

//main regression: fractional response model
//scrd level Bt rate
qui fracreg probit $rw  i.bt##c.hist_scrd prec4cm-prec8cm i.seedtrtrate  i.locid , vce(cluster $vce)
est store f010
margins, dydx(*) post
est store f01

qui fracreg probit $rw  i.bt##c.hist_scrd prec4cm-prec8cm i.seedtrtrate $pattern i.locid, vce(cluster $vce)
est store f020
margins, dydx(*) post
est store f02

qui fracreg probit $rw  i.bt##c.hist_scrd prec4cm-prec8cm i.seedtrtrate tMin1 tMin2 tMin3 i.locid, vce(cluster $vce)
est store f030
margins, dydx(*) post
est store f03

//county level Bt rate
qui fracreg probit $rw  i.bt##c.hist_countyfips prec4cm-prec8cm i.seedtrtrate  i.locid, vce(cluster $vce)
est store f110
margins, dydx(*) post
est store f11

qui fracreg probit $rw  i.bt##c.hist_countyfips prec4cm-prec8cm i.seedtrtrate  $pattern i.locid, vce(cluster $vce)
est store f120
margins, dydx(*) post
est store f12

qui fracreg probit $rw i.bt##c.hist_countyfips prec4cm-prec8cm i.seedtrtrate tMin1 tMin2 tMin3 i.locid, vce(cluster $vce)
est store f130
margins, dydx(*) post
est store f13

esttab  f010 f020 f030 using $do\results.rtf, ///
	star(* 0.10 ** 0.05 *** 0.01) ///
	keep(1.bt hist_scrd 1.bt#c.hist_scrd prec4cm prec5cm prec6cm prec7cm prec8cm 1.seedtrtrate $pattern tMin1 tMin2 tMin3) ///
	title ("Table S3. Fractional response model results (coefficient estimates) for root injury, using Crop Reporting District level Bt planting rate") ///
	se b(%9.2f) label ///
	replace

esttab  f110  f120  f130  using $do\results.rtf, ///
	star(* 0.10 ** 0.05 *** 0.01) ///
	keep(1.bt hist_countyfips 1.bt#c.hist_countyfips prec4cm prec5cm prec6cm prec7cm prec8cm 1.seedtrtrate $pattern tMin1 tMin2 tMin3) ///
	title ("Table S6. Robustness regression analysis: using county-level Bt planting rate (fractional response model coefficient estimates)") ///
	se b(%9.2f) label ///
	append

//robust: linear comparison
qui reg $rw  i.bt##c.hist_scrd prec4cm-prec8cm i.seedtrtrate  i.locid, vce(cluster $vce)
est store n01

qui reg $rw  i.bt##c.hist_scrd prec4cm-prec8cm i.seedtrtrate $pattern  i.locid, vce(cluster $vce)
est store n02

qui reg $rw  i.bt##c.hist_scrd prec4cm-prec8cm i.seedtrtrate tMin1 tMin2 tMin3 i.locid, vce(cluster $vce)
est store n03

qui reg $rw  i.bt##c.hist_countyfips prec4cm-prec8cm i.seedtrtrate  i.locid, vce(cluster $vce)
est store n11

qui reg $rw  i.bt##c.hist_countyfips prec4cm-prec8cm i.seedtrtrate $pattern i.locid, vce(cluster $vce)
est store n12

qui reg $rw  i.bt##c.hist_countyfips prec4cm-prec8cm i.seedtrtrate tMin1 tMin2 tMin3 i.locid, vce(cluster $vce)
est store n13

esttab  n01 n02  n03 using $do\results.rtf, ///
	star(* 0.10 ** 0.05 *** 0.01) ///
	keep(1.bt hist_scrd 1.bt#c.hist_scrd prec4cm prec5cm prec6cm prec7cm prec8cm 1.seedtrtrate $pattern tMin1 tMin2 tMin3) ///
	title ("Table S7. Robustness regression analysis: linear regression using Crop Reporting District level Bt planting rate") ///
	se b(%9.2f) label ///
	append

esttab  n11 n12  n13 using $do\results.rtf, ///
	star(* 0.10 ** 0.05 *** 0.01) ///
	keep(1.bt hist_countyfips 1.bt#c.hist_countyfips prec4cm prec5cm prec6cm prec7cm prec8cm 1.seedtrtrate $pattern tMin1 tMin2 tMin3) ///
	title ("Table S8. Robustness regression analysis: linear regression using county level Bt planting rate") ///
	se b(%9.2f) label ///
	append	

	
//robust: multi-year effects of Bt planting history
use $dd\estim_use.dta, clear
 fracreg probit $rw  i.bt hist_scrd hist2_scrd  prec4cm-prec8cm i.seedtrtrate  i.locid , vce(cluster $vce)
est store f310
 
 fracreg probit $rw  i.bt hist_scrd hist2_scrd i.bt#c.hist_scrd i.bt#c.hist2_scrd prec4cm-prec8cm i.seedtrtrate  i.locid , vce(cluster $vce)
est store f320

esttab  f310 f320 using $do/results.rtf, ///
	star(* 0.10 ** 0.05 *** 0.01) ///
	keep(1.bt hist_scrd hist2_scrd 1.bt#c.hist_scrd 1.bt#c.hist2_scrd prec4cm prec5cm prec6cm prec7cm prec8cm 1.seedtrtrate ) ///
	title ("Table S9. Robustness regression analysis: multi-year effects of Bt planting history (fractional response model coefficient estimates)") ///
	se b(%9.2f) label ///
	replace


//robust:  multi-year effects of precipitation
use $dd\estim_use.dta, clear

fracreg probit $rw  i.bt##c.hist_scrd i.seedtrtrate prec4cm-prec8cm lag1_prec*  i.locid , vce(cluster $vce)
est store p1
fracreg probit $rw  i.bt##c.hist_scrd i.seedtrtrate prec4cm-prec8cm lag1_prec* lag2_prec* i.locid , vce(cluster $vce)
est store p2
fracreg probit $rw  i.bt##c.hist_scrd i.seedtrtrate prec4cm-prec8cm lag1_prec* lag2_prec* lag3_prec* i.locid , vce(cluster $vce)
est store p3

esttab  p1 p2 p3 p4 using $do/prec.rtf, ///
	star(* 0.10 ** 0.05 *** 0.01) ///
	keep(1.bt hist_scrd 1.bt#c.hist_scrd 1.seedtrtrate prec4cm prec5cm prec6cm prec7cm prec8cm lag1_prec* lag2_prec* lag3_prec* ) ///
	title ("Table S10. Robustness regression analysis: multi-year effects of precipitation (fractional response model coefficient estimates)") ///
	se b(%9.2f) label ///
	replace	

*-------------predicted root injury table in the main manuscript (Table 1)--------------*
est restore f010
margins bt, at(hist_scrd=(0 1)) post
lincom _b[1._at#1.bt]-_b[1._at#0.bt] 
lincom _b[2._at#1.bt]-_b[2._at#0.bt] 
lincom _b[2._at#1.bt]-_b[1._at#1.bt] 
lincom _b[2._at#0.bt]-_b[1._at#0.bt]
lincom _b[2._at#1.bt]-_b[2._at#0.bt]-(_b[1._at#1.bt]-_b[1._at#0.bt])

*-------------analyses by protein type (table S11-12, fig. S4)--------------*

use  $dd\estim_use.dta, clear

//identify trial bt type
gen bttype=0
replace bttype=1 if protein=="Cry3Bb1"
replace bttype=2 if protein=="Cry34/35Ab1"
replace bttype=3 if protein=="mCry3A"
replace bttype=4 if protein=="Cry3Bb1+Cry34/35Ab1"
replace bttype=5 if protein=="mCry3A+Cry34/35Ab1"
replace bttype=6 if protein=="mCry3A+eCry3.1Ab"

forv i=1/6{
	gen bt`i'=(bttype==`i')
}

fracreg probit $rw  ///
	i.bt1 i.bt2 i.bt3 i.bt4 i.bt5 i.bt6 hist_scrd ///
	i.bt1#c.hist_scrd i.bt2#c.hist_scrd i.bt3#c.hist_scrd i.bt4#c.hist_scrd i.bt5#c.hist_scrd i.bt6#c.hist_scrd ///
	prec4cm-prec8cm i.seedtrtrate  i.locid, vce(cluster $vce)
est store sf010

fracreg probit $rw  ///
	i.bt1 i.bt2 i.bt3 i.bt4 i.bt5 i.bt6 hist_scrd ///
	i.bt1#c.hist_scrd i.bt2#c.hist_scrd i.bt3#c.hist_scrd i.bt4#c.hist_scrd i.bt5#c.hist_scrd i.bt6#c.hist_scrd ///
	prec4cm-prec8cm i.seedtrtrate  $pattern i.locid, vce(cluster $vce)
est store sf020

fracreg probit $rw  ///
	i.bt1 i.bt2 i.bt3 i.bt4 i.bt5 i.bt6 hist_scrd ///
	i.bt1#c.hist_scrd i.bt2#c.hist_scrd i.bt3#c.hist_scrd i.bt4#c.hist_scrd i.bt5#c.hist_scrd i.bt6#c.hist_scrd ///
	prec4cm-prec8cm i.seedtrtrate  $pattern tMin1 tMin2 tMin3 i.locid, vce(cluster $vce)
est store sf030

esttab sf010 sf020 sf030 using $do\results.rtf, ///
	star(* 0.10 ** 0.05 *** 0.01) ///
	keep(1.bt1 1.bt2 1.bt3 1.bt4 1.bt5 1.bt6 hist_scrd ///
	1.bt1#c.hist_scrd 1.bt2#c.hist_scrd 1.bt3#c.hist_scrd 1.bt4#c.hist_scrd 1.bt5#c.hist_scrd 1.bt6#c.hist_scrd ///
	prec4cm prec5cm prec6cm prec7cm prec8cm 1.seedtrtrate $pattern tMin1 tMin2 tMin3) ///
	title ("Table S11. Robustness regression analysis: analyze for each Bt protein type separately (fractional response model coefficient estimates)") ///
	se b(%9.2f) label ///
	replace

//produce fig. S4 and table 12
*obtain estimates of parameter a for each protein type
est restore f010
margins bt, at(hist_scrd=(0 1)) post
lincom _b[1._at#0.bt]-_b[1._at#1.bt]
scalar a=r(estimate)
scalar a_se=r(se)

est restore sf010
margins bt1, at(hist_scrd=(0 1) bt2=0 bt3=0 bt4=0 bt5=0 bt6=0) post
lincom _b[1._at#0.bt1]-_b[1._at#1.bt1]
scalar a1=r(estimate)
scalar a1_se=r(se)

est restore sf010
margins bt2, at(hist_scrd=(0 1) bt1=0 bt3=0 bt4=0 bt5=0 bt6=0) post
lincom _b[1._at#0.bt2]-_b[1._at#1.bt2]
scalar a2=r(estimate)
scalar a2_se=r(se)

est restore sf010
margins bt3, at(hist_scrd=(0 1) bt1=0 bt2=0 bt4=0 bt5=0 bt6=0) post
lincom _b[1._at#0.bt3]-_b[1._at#1.bt3]
scalar a3=r(estimate)
scalar a3_se=r(se)

est restore sf010
margins bt4, at(hist_scrd=(0 1) bt1=0 bt2=0 bt3=0 bt5=0 bt6=0) post
lincom _b[1._at#0.bt4]-_b[1._at#1.bt4]
scalar a4=r(estimate)
scalar a4_se=r(se)

est restore sf010
margins bt5, at(hist_scrd=(0 1) bt1=0 bt2=0 bt3=0 bt4=0 bt6=0) post
lincom _b[1._at#0.bt5]-_b[1._at#1.bt5]
scalar a5=r(estimate)
scalar a5_se=r(se)

est restore sf010
margins bt6, at(hist_scrd=(0 1) bt1=0 bt2=0 bt3=0 bt4=0 bt5=0) post
lincom _b[1._at#0.bt6]-_b[1._at#1.bt6]
scalar a6=r(estimate)
scalar a6_se=r(se)

*obtain estimates of parameter b for each protein type
est restore f010
margins bt, at(hist_scrd=(0 1)) post
lincom _b[1._at#0.bt]-_b[1._at#1.bt]-_b[2._at#0.bt]+_b[2._at#1.bt]
scalar b=r(estimate)
scalar b_se=r(se)

est restore sf010
margins bt1, at(hist_scrd=(0 1) bt2=0 bt3=0 bt4=0 bt5=0 bt6=0) post
lincom _b[1._at#0.bt1]-_b[1._at#1.bt1]-_b[2._at#0.bt1]+_b[2._at#1.bt1]
scalar b1=r(estimate)
scalar b1_se=r(se)

est restore sf010
margins bt2, at(hist_scrd=(0 1) bt1=0 bt3=0 bt4=0 bt5=0 bt6=0) post
lincom _b[1._at#0.bt2]-_b[1._at#1.bt2]-_b[2._at#0.bt2]+_b[2._at#1.bt2]
scalar b2=r(estimate)
scalar b2_se=r(se)

est restore sf010
margins bt3, at(hist_scrd=(0 1) bt1=0 bt2=0 bt4=0 bt5=0 bt6=0) post
lincom _b[1._at#0.bt3]-_b[1._at#1.bt3]-_b[2._at#0.bt3]+_b[2._at#1.bt3]
scalar b3=r(estimate)
scalar b3_se=r(se)


est restore sf010
margins bt4, at(hist_scrd=(0 1) bt1=0 bt2=0 bt3=0 bt5=0 bt6=0) post
lincom _b[1._at#0.bt4]-_b[1._at#1.bt4]-_b[2._at#0.bt4]+_b[2._at#1.bt4]
scalar b4=r(estimate)
scalar b4_se=r(se)


est restore sf010
margins bt5, at(hist_scrd=(0 1) bt1=0 bt2=0 bt3=0 bt4=0 bt6=0) post
lincom _b[1._at#0.bt5]-_b[1._at#1.bt5]-_b[2._at#0.bt5]+_b[2._at#1.bt5]
scalar b5=r(estimate)
scalar b5_se=r(se)

est restore sf010
margins bt6, at(hist_scrd=(0 1) bt1=0 bt2=0 bt3=0 bt4=0 bt5=0) post
lincom _b[1._at#0.bt6]-_b[1._at#1.bt6]-_b[2._at#0.bt6]+_b[2._at#1.bt6]
scalar b6=r(estimate)
scalar b6_se=r(se)

*coefpot input matrix
matrix par_a=(a1,a2,a3,a4,a5,a6/*,a*/ \ a1_se,a2_se,a3_se,a4_se,a5_se,a6_se/*,a_se*/)
matrix par_b=(b1,b2,b3,b4,b5,b6/*,b*/ \b1_se,b2_se,b3_se,b4_se,b5_se,b6_se/*,b_se*/)

*check input matrix
mat list par_a
mat list par_b

loc a=a
loc b=b
coefplot ///
	(mat(par_a), msymbol(T) mcolor(red) ciopt(recast(rcap) lcolor(red))) ///	
	(mat(par_b), msymbol(O) mcolor(white) mlcolor(black) ciopt(recast(rcap) lcolor(black))) ///
	,se(2) ///
	ylabel(1 "Cry3Bb1" 2 "Cry34/35Ab1" 3 "mCry3A" 4 "mCry3A+Cry34/35Ab1" ///
		5 "mCry3A+Cry34/35Ab1" 6 "mCry3A+eCry3.1Ab") ///
	legend(off) ///
	xlabel(,format(%9.1f)) ///
	xline(`a', lc(red) lp(dot)) xline(`a', lc(red%10) lw(5))  ///
	xline(`b', lc(black) lp(dot)) xline(`b', lc(black%10) lw(11)) 
	//the two lw numbers (5 and 11) are obtained by adding the pooled results in the figure and calibrate the number and then removing the pooled results; so change of graph size needs to change the lw numbers.

graph dis, xsize(5) ysize(5)
graph export $do\fig4.tif, replace width(1425) //horizontal pixels=4.75inch*300dots per inch
graph export $do\fig4.png, replace width(1425)	


*------------------------scrd level analysis (table S4-5)---------------------*
use $dd/estim_use.dta, clear
drop if seedtrtrate==1
sort  stateabbr locid id year
foreach v of varlist prec4cm-prec8cm{
	bys year scrd: egen `v'_scrd=mean(`v')
} 
collapse (mean) rw hist_scrd prec*_scrd,by(bt year stateabbr statefips scrd)
reshape wide rw, i(year hist_scrd statefips scrd prec*_scrd) j(bt)

xtset scrd year
gen rw_d=rw0-rw1
gen efficacy=(rw0-rw1)/rw0
gen l_rw0=l.rw0
gen l_rw1=l.rw1
gen l_rw_d=l.rw_d
gen l_efficacy=l.efficacy
drop if rw_d<0

egen miss=rowmiss(rw0 rw_d efficacy hist_scrd l_rw0 l_efficacy prec4cm_scrd-prec8cm_scrd)
drop if miss

reg rw_d l_rw_d prec4cm_scrd-prec8cm_scrd i.scrd , vce(robust)
gen byte used=e(sample)
keep if used==1

reg rw0 hist_scrd l_rw0 l_efficacy i.scrd , vce(robust) 
est store m1
reg rw0 hist_scrd i.scrd , vce(robust)
est store m2
reg rw0 l_rw0 l_efficacy i.scrd , vce(robust)
est store m3


reg rw0 hist_scrd l_rw0 l_efficacy i.scrd prec4cm_scrd-prec8cm_scrd, vce(robust) 
est store m4
reg rw0 hist_scrd i.scrd prec4cm_scrd-prec8cm_scrd, vce(robust)
est store m5
reg rw0 l_rw0 l_efficacy i.scrd prec4cm_scrd-prec8cm_scrd, vce(robust)
est store m6


reg rw_d hist_scrd l_rw0 l_efficacy i.scrd , vce(robust)
est store n1
reg rw_d hist_scrd i.scrd , vce(robust)
est store n2
reg rw_d l_rw0 l_efficacy i.scrd , vce(robust)
est store n3


reg rw_d hist_scrd l_rw0 l_efficacy i.scrd prec4cm_scrd-prec8cm_scrd , vce(robust)
est store n4
reg rw_d hist_scrd i.scrd prec4cm_scrd-prec8cm_scrd , vce(robust)
est store n5
reg rw_d l_rw0 l_efficacy prec4cm_scrd-prec8cm_scrd i.scrd , vce(robust)
est store n6


reg efficacy hist_scrd l_efficacy i.scrd , vce(robust)
est store e1
reg efficacy hist_scrd i.scrd , vce(robust)
est store e2
reg efficacy l_efficacy i.scrd , vce(robust)
est store e3


reg efficacy hist_scrd l_efficacy prec4cm_scrd-prec8cm_scrd i.scrd , vce(robust)
est store e4
reg efficacy hist_scrd prec4cm_scrd-prec8cm_scrd i.scrd , vce(robust)
est store e5
reg efficacy l_efficacy prec4cm_scrd-prec8cm_scrd i.scrd , vce(robust)
est store e6

esttab m1  m4 e1 e4 n1 n4 using $do/results.rtf, ///
	star(* 0.10 ** 0.05 *** 0.01) ///
	keep(hist_scrd l_rw0 l_efficacy) ///
	order(hist_scrd l_rw0 l_efficacy ) ///
	title ("Table S4. Estimation results of equations of motion in terms of root injury measures") ///
	se r2 b(%9.2f) label ///
	replace
	
esttab m2 m5 e2 e5 n2 n5  using $do/results.rtf, ///
	star(* 0.10 ** 0.05 *** 0.01) ///
	keep(hist_scrd ) ///
	order(hist_scrd  ) ///
	title ("Table S5. Estimation results of equations of motion in terms of root injury measures, without lagged variables") ///
	se r2 b(%9.2f) label ///
	append

		

