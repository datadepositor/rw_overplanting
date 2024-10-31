********************************************************************************
*				program: descriptives			   	   						   *											 
********************************************************************************
global root "..."

global data "$root/1_data"
global ds "$data/source"
global dd "$data/data"
global dp "$data/proc"
global dt "$data/temp"
global do "$data/out"

set scheme s1color, permanent
graph set window fontface "Helvetica"


*------------------Fig.2 & table S2 raw comparison and Bt efficacy----------------------------* 
global psize "16pt" //panel text size
global tsize "16pt" //figure text size

***panel A	
use $dd/estim_use.dta, clear
keep if soilinsecticide==0 &seedtrtrate==0  
drop if stateabbr=="OH" | stateabbr=="MN"

duplicates drop year scrd, force
collapse (sum) *acres_scrd, by(year)

foreach v of varlist only_* two_*{
	renvarlab `v',postdrop(11) display
}

foreach v of varlist only_* two_*{
	gen `v'_rate=`v'/acres_scrd
}

gen h1=only_Cry3Bb1_rate
loc t=1
foreach v of varlist only_Cry3435A_rate only_mCry3A_rate two_Cry3Bb1_Cry3435A_rate  two_mCry3A_Cry3435A_rate two_mCry3A_eCry31Ab_rate{
	loc s=`t'
	loc t=`t'+1
	gen h`t'=h`s'+`v'
	loc s=`s'+1
}

sort year


//blue to yellow pallette
loc r1 "000 158 115"
loc r2 "000 114 178"
loc r3 "086 180 233"
loc b3 "240 228 066"
loc b2 "230 159 000"
loc b1 "213 094 000"

loc barw=0.8
loc lw=0.2
graph twoway (bar h6 h5 h4 h3 h2 h1 year,barw(`barw' `barw' `barw' `barw' `barw' `barw') ///
	lw(`lw' `lw' `lw' `lw' `lw' `lw') lcolor(black black black black black black) ///
	color("`b1'" "`b2'" "`b3'" "`r3'"  "`r2'"  "`r1'")), ///
	ylabel(0 (0.2) 1.0 ,format(%03.1f)  labsize("$tsize") ) ///
	xlabel(2004 (4) 2016, labsize("$tsize")) ///
	title("{bf:A}",size("18pt") ring(1) position(11) margin(b=3)) /// 
	xtitle("{bf:Year}", size("$psize") margin(t=3)) ///
	ytitle("{bf:Bt planting rate}", margin(r=3) size("$psize") )  ///
	leg(region(lwidth(none) color(none)) ring(0) position(11) col(2) size("$tsize") symx(*0.2) bmargin(zero) ///
	order(4 "mCry3A" 3 "Cry3Bb1+Cry34/35A" 5 "Cry34/35A" 2 "Cry34/35A+mCry3A" 6 "Cry3Bb1" 1 "mCry3A+eCry3.1Ab")) 
	
graph save $dt/rwstacking.gph, replace
graph export $dt/rwstacking.png, replace

	
//panel B: by trait	
use $dd/estim_use.dta, clear
keep if soilinsecticide==0 &seedtrtrate==0  
drop if stateabbr=="OH" | stateabbr=="MN"

collapse (mean) mean= rw  ///
	, by(bt year locid protein)
sort year locid bt protein 

gen bttype=0
replace bttype=1 if protein=="Cry3Bb1"
replace bttype=2 if protein=="Cry34/35Ab1"
replace bttype=3 if protein=="mCry3A"
replace bttype=4 if protein=="Cry3Bb1+Cry34/35Ab1"
replace bttype=5 if protein=="mCry3A+Cry34/35Ab1"
replace bttype=6 if protein=="mCry3A+eCry3.1Ab"

keep year bttype mean locid
reshape wide mean, i(year locid) j(bttype)
forv i=1/6{
	gen mean_abott`i'=100*(mean0-mean`i')/mean0
}

for var mean_abott*: drop if X<0 

egen mean_pooled=rowmean(mean1 mean2 mean3 mean4 mean5 mean6)
gen mean_abott=100*(mean0-mean_pooled)/mean0

collapse (mean) mean_abott*, by(year)

//blue to yellow pallette
loc r1 "000 158 115"
loc r2 "000 114 178"
loc r3 "086 180 233"
loc b3 "240 228 066"
loc b2 "230 159 000"
loc b1 "213 094 000"

graph twoway ///
	(lowess mean_abott1 year,lcolor("`r1'") lw(0.3) lpattern(solid)) ///
	(lowess mean_abott2 year,lcolor("`r2'") lw(0.3) lpattern(solid)) ///
	(lowess mean_abott3 year,lcolor("`r3'") lw(0.3) lpattern(solid)) ///	
	(lowess mean_abott4 year,lcolor("`b3'") lw(0.3) lpattern(solid)) ///
	(lowess mean_abott5 year,lcolor("`b2'") lw(0.3) lpattern(solid)) ///
	(scatter mean_abott1 year,mcolor("`r1'") msize(medium) msymbol(plus)) ///
	(scatter mean_abott2 year,mcolor("`r2'") msize(small) msymbol(T)) ///
	(scatter mean_abott3 year,mcolor("`r3'") msize(large) msymbol(x)) ///
	(scatter mean_abott4 year,mcolor("`b3'") msize(small) msymbol(V)) ///
	(scatter mean_abott5 year,mcolor("`b2'") msize(small) msymbol(|)) ///
	(lowess mean_abott year,lcolor(black) lw(0.5) lpattern(solid)) ///
	(scatter mean_abott year,mcolor(black)  mfcolor(white) msize(small) msymbol(O)) ///	
	, ///
	ylabel(20 (20) 100, labsize("$tsize")) ///
	xlabel(2004 (4) 2016.5, labsize("$tsize")) ///
leg(region(lwidth(none)) on  ring(0) position(6) col(2) size("$tsize") ///
	order(6 "Cry3Bb1" 9 "Cry3Bb1+Cry34/35A"  7 "Cry34/35A" 10 "Cry34/35A+mCry3A" 8 "mCry3A"  12 "pooled") ) ///
	title("{bf:B}",size("18pt") ring(1) position(11) margin(b=3)) /// 
	xtitle("{bf:Trial year}", size("$psize") margin(t=3)) ///
	ytitle("{bf:Bt efficacy (%)}", margin(r=3) size("$psize") ) 

graph save $dt/efficacy.gph, replace

//table S2
asdoc list, dec(0) replace
	
//panel C
use $dd/estim_use.dta, clear
keep if soilinsecticide==0 &seedtrtrate==0  
drop if stateabbr=="OH" | stateabbr=="MN"

preserve
duplicates drop year scrd, force
collapse (sum) lag2_rw_acres_scrd lag2_acres_scrd, by(year)
gen lag2_rw_rate=lag2_rw_acres_scrd/lag2_acres_scrd
save $dt/bt_year.dta, replace
restore

collapse (mean) mean= rw  ///
	(sd) sd= rw ///
	(count) n=rw ///
	, by(bt year)
gen hi=mean+invttail(n-1,0.025)*(sd / sqrt(n))
gen lo=mean-invttail(n-1,0.025)*(sd / sqrt(n))

reshape wide mean sd n hi lo, i(year) j(bt)
gen mean_d=mean0-mean1
merge 1:1 year using $dt/bt_year.dta

loc lw 0.6
graph twoway ///
(rarea hi1 lo1 year , fcolor("255 243 200") lw(0)) ///
 (connected mean1 year ,msymbol(i)  lw(`lw') color("230 159 000")) ///
 (rarea hi0 lo0 year , fcolor("231 239 250") lw(0)) ///
 (connected mean0 year ,msymbol(i)  lw(`lw') color("000 114 178")) ///
 (connected mean_d year,msymbol(O) lpattern(solid) lcolor(black) mcolor(black) mfcolor(white) lw(1) yaxis(1)) ///
(bar lag2_rw_rate year,bcolor("green%15") lw(0) yaxis(1)) , ///
	xlabel(2004 (4) 2016,labsize("$tsize"))  ///
	ylabel(0 (0.2) 1, labsize("$tsize") axis(1) format(%03.1f))  ///
	legend(ring(0) position(11) region(lwidth(none)) size("$tsize") col(2) ///
	order( 2 "Bt injury" 4 "non-Bt injury" 5 "Bt protection" 6 "Bt planting" "history")) ///
	xtitle("{bf:Trial year}",size("$psize") margin(t=3)) ///
	ytitle("{bf:Bt planting history / Root injury}", axis(1) margin(r=3) size("$psize")) ///	
	title("{bf:C}",size("18pt") ring(1) position(11)  margin(b=3))  

graph save $dt/bteffect_year.gph, replace

	
//panel D
use $dd/estim_use.dta, clear
keep if soilinsecticide==0 &seedtrtrate==0  
drop if stateabbr=="OH" | stateabbr=="MN"

local interv=10/100
local begin=`interv'/2
egen rwcut = cut(hist_scrd), at(0 (`interv') 1)
collapse (mean) rw, by(bt rwcut)
reshape wide rw, i(rwcut) j(bt)
gen rw_d=rw0-rw1
gen xlab=rwcut+`interv'/2
keep if rwcut<0.75

loc lw 0.6
graph twoway ///
	(lfitci rw0 xlab,  lcolor("000 114 178") clw(`lw') lpattern(shortdash) acolor("231 239 250") alwidth(0)  ) ///
	(lfitci rw1 xlab,  lcolor("230 159 000") clw(`lw') lpattern(shortdash) acolor("255 243 200") alwidth(0)) ///
	(lfitci rw_d xlab,  lcolor(black) clw(1) lpattern(solid)  acolor(grey%10) alwidth(0)  ) ///
	(scatter rw_d xlab,msymbol(O)  mcolor(black) lpattern(solid)  mfcolor(white) lw(1) msize(medium)) ///
	(scatter rw0 xlab, msymbol(plus) mcolor("000 114 178") lpattern(solid) lw(2) msize(large) ) ///
	(scatter rw1 xlab,msymbol(X) mcolor("230 159 000") lpattern(solid)  lw(2) msize(large) ), ///
	xlabel(0 (0.2) 0.8,labsize("$tsize") format(%03.1f)) ///
	ylabel(0 (0.2) 1,labsize("$tsize")format(%03.1f) ) ///
	title("{bf:D}",size("18pt") ring(1) position(11) margin(b=3)) /// 
	leg(region(lwidth(none)) on  ring(0) position(11) col(1)  size("$tsize") ///
	order(9 "Bt injury" 8 "non-Bt injury" 7 "Bt protection")) ///
	xtitle("{bf:Bt planting history}", size("$psize") margin(t=3)) ///
	ytitle("{bf:Root injury}", margin(r=3) size("$psize") )  ///
	saving($dt/bteffect_coverage_CC.gph, replace) 

graph save $dt/bteffect_coverage_CC.gph, replace
	
graph combine $dt/rwstacking.gph $dt/efficacy.gph $dt/bteffect_year.gph $dt/bteffect_coverage_CC.gph ///
	, iscale(*0.7) imargin(tiny) ///
	col(2) rows(2) ///
	xsize(6) ysize(6)

graph export $do/Fig2.tif, replace width(2000) //horizontal pixels=4.75inch*300dots per inch
graph export $do/Fig2.png, replace width(2000)	


*-----------------Fig.1B RW trial plot map-------------------------------------------------------*
use $dd/estim_use.dta, clear
keep year location stateabbr _X _Y
bys location: gen nobs_avg=_N/12 

collapse (mean) nobs_avg, by(stateabbr location  _X _Y)

egen size=cut(nobs_avg), at(0,2,4,6,8,10,12,14,16)

geo2xy _Y _X,proj(albers) replace
save $dt/RWmap.dta, replace

*keep Midwest area only
cd "/Users/eleven/Dropbox/0_Graduate Program/计量/stata/graph/examples"
//cd "D:\Dropbox\0_Graduate Program\计量\stata\graph\examples"
//use usa_state, clear
use usa_state_shp, clear
merge m:1 _ID using usa_state  // to get the identifiers
destring STATEFP, replace 
keep if inlist(STATEFP,19,17,18,39,38,46,31,27,55,26)

//twoway (scatter _Y _X) //get the map of the mainland states

*get the mainland map
drop _CX- _merge 
sort _ID shape_order
geo2xy _Y _X, proj(albers) replace
save MidWest_shp_clean.dta, replace //the basemap dataset to be used.


use usa_state_shp, clear

merge m:1 _ID using usa_state  // to get the identifiers
destring STATEFP, replace 
keep if inlist(STATEFP,19,17,18,39,38,46,31,27,55,26)

//twoway (scatter _Y _X) //get the map of the mainland states
use usa_state, clear //the master data must contain the basemap identifier
spmap using MidWest_shp_clean, id(_ID) fcolor(white) ///basemap defined as mainland shapefile in stata format
	point(data($dt/RWmap.dta) x(_X) y(_Y) fcolor(Greens) by(size) ///
	shape(triangle triangle triangle triangle triangle triangle triangle) ///
	ocolor(green*1.7 green*1.7 green*1.7 green*1.7 green*1.7 green*1.7 green*1.7) ///
	legenda(off))
graph dis, xsize(3) ysize(2.5)
//use RWmap recording.

graph use $do/RWmap_v7.gph 
graph export $do/RWmap_v7.png, replace width(900)
graph export $do/RWmap_v7.tif, replace width(900)
graph export $do/RWmap_v7.pdf, replace 

		

*--------------------fig. S1 state maps------------------------------------------------------*
*****state level
//prepare background map data for cornbelt region
use "/Users/eleven/Dropbox/1_Research/8_climate_biotech/1_data/data/v2022/dataprep/region.dta"
//use "D:\Dropbox\1_Research\8_climate_biotech\1_data\data\v2022\dataprep\region.dta", clear
keep statefips cornbelt
duplicates drop
save $dt/cornbeltstates.dta, replace

use "/Users/eleven/Dropbox/1_Research/9_halo/1_data/source/maptile_geographies/state_database_clean.dta"
//use "D:\Dropbox\1_Research\9_halo\1_data\source\maptile_geographies\state_database_clean.dta",clear
merge m:1 statefips using $dt/cornbeltstates.dta
keep if cornbelt==1
save $dt/state_database_clean_cornbelt.dta, replace

use "/Users/eleven/Dropbox/1_Research/9_halo/1_data/source/maptile_geographies/state_coords_clean.dta"
//use "D:\Dropbox\1_Research\9_halo\1_data\source\maptile_geographies\state_coords_clean.dta",clear
gen oid=_n
rename _ID _polygonid
merge m:1 _polygonid using $ds/maptile_geographies/state_database_clean.dta, keepusing(statefips) nogen
merge m:1 statefips using $dt/cornbeltstates.dta, nogen keepusing(cornbelt)
sort oid

keep if cornbelt==1
rename _polygonid _ID
drop statefips cornbelt
save $dt/state_coords_clean_cornbelt.dta, replace

use $dd/estim_use.dta, clear
keep if bt==0
keep if soilinsecticide==0 &seedtrtrate==0 
keep if year<=2016 & year>=2014

gen rw_scrd=rw_rate_scrd/100
gen rw_statefips=rw_rate_statefips/100
collapse (mean) rw rw_scrd rw_statefips, by(statefips)

capture drop _merge
merge 1:1 statefips using $dt/cornbeltstates.dta, keepusing(cornbelt) nogen

keep if cornbelt==1

merge 1:1 statefips using $dt/state_database_clean_cornbelt.dta, nogen
drop _merge


save $dt/rwxy.dta, replace

spmap rw using $ds/maptile_geographies/state_coords_clean.dta, id(_polygonid)  ///
	clmethod(custom) clbreaks(0 0.15 0.3 0.45 1) fcolor(Blues) ///
	legend(lab(1 "No data") lab(2 "0.00 - 0.15") lab(3 "0.15 - 0.30") lab(4 "0.30 - 0.45") lab(5 "0.45 - 1.00"))  ///
	legend(ring(1) position(3) size(medium))  saving($do/threeyear_rw, replace) ///
	title("{bf:A. Root injury}",size(medium) pos(11))

spmap rw_scrd using $ds/maptile_geographies/state_coords_clean.dta, id(_polygonid)  ///
	clmethod(custom) clbreaks(0 0.20 0.4 0.6 1) fcolor(Reds) ///
	legend(lab(1 "No data") lab(2 "0.0 - 0.2") lab(3 "0.2 - 0.4") lab(4 "0.4 - 0.6") lab(5 "0.6 - 1.00"))  ///
	legend(ring(1) position(3) size(medium))  saving($do/threeyear_rw_scrd, replace) ///
	title("{bf:B. Bt planting rate}",size(medium) pos(11))


graph combine $do/threeyear_rw.gph $do/threeyear_rw_scrd.gph, col(1)
graph display, xsize(4) ysize(5.5) scale(1.4)
graph export "$do/figS1.png", wid(2000) replace  
graph export "$do/figS1.tif", wid(2000) replace  
	

*----------------fig.S2 crop rotation-----------------------------------------*
use $dd\estim_use, clear
collapse (mean) cornshare_corn  share_corncorn share_corn, by(year stateabbr)
for var cornshare_corn  share_corncorn share_corn: replace X=X/100
graph twoway (bar cornshare_corn year,bcolor(gray%10) lw(0) yaxis(1)) ///
	(scatter share_corncorn year, msymbol(O) mcolor(black) mfcolor(white) msize(1)) ///
	(scatter share_corn year,msymbol(O) mcolor(black) msize(1)) ///
	,by(stateabbr, col(5)) ///
	xlabel(2004 (6) 2016, labsize(medium)) ylabel(0 (0.25) 0.50, labsize(medium) format(%03.1f) ) ///
	legend(region(lwidth(none)) size(small) col(1) ///
	lab(1 "Ratio of corn-to-corn acres to corn acres") ///
	lab(2 "Ratio of corn-to-corn acres to total acres") ///
	lab(3 "Ratio of corn acres to total acres")) ///
	xtitle("Year", size(small) margin(t=3)) ytitle("Ratio", size(small) margin(r=3)) 
//manually edit box color and so on; the "bystate" recording
graph export "$do\figS2.tif", wid(2000) replace  
graph export "$do\figS2.png", wid(2000) replace  


*----------------fig.S3 root injury and Bt planting history by state---------------*

use $dd\estim_use.dta, clear
keep if bt==0
keep if soilinsecticide==0 &seedtrtrate==0  &year>=2005 
collapse (mean) rw hist_scrd, by(stateabbr year)

graph twoway (bar hist_scrd year, bcolor(green%30))(scatter rw year, yaxis(2)) ///
	(lfit rw year, yaxis(2)  ) ///
	, by(stateabbr, col(5)) ///
	legend(region(lwidth(none)) size(small) col(3) label(1 "Bt planting rate") label(2 "Root injury") label(3 "Fitted value") ) ///
	ylabel(0 (0.4) 0.8, labsize(medium) format(%03.1f) ) ylabel(0 (0.3) 0.6, labsize(medium) format(%03.1f) axis(2) )  ///
	xlabel(2004 (6) 2016, labsize(medium))  ///
	xtitle("Trial year", size(small)) ytitle("Bt planting rate",axis(1) size(small)) ytitle("Root injury",axis(2) size(small)) 
//manually edit box color and so on; the "bystate" recording
graph export "$do\figS3.tif", wid(2000) replace  
graph export "$do\figS3.png", wid(2000) replace  


*------------------------fig.S8 below-ground trait stacking---------------------------*
use $dt\bt_nation.dta, replace
gen h1=rwonly_rate_nation/100
gen h2=h1+rwcbonly_rate_nation/100
gen h3=h2+rwhtonly_rate_nation/100
gen h4=h3+rwcbht_rate_nation/100


sort year

//red to blue pallette
loc r1 "219 049 036"
loc r2 "252 140 090"
loc r3 "255 223 146"
loc b3 "230 241 243"
loc b2 "144 190 224"
loc b1 "075 116 178"

loc barw=0.8
loc lw=0.2
graph twoway (bar   h4 h3 h2 h1 year,barw(`barw' `barw' `barw' `barw') lw(`lw' `lw' `lw' `lw') lc(black black black black) color("`b3'" "`b2'" "`b1'" "`r1'")), ///
	yscale(range(0 0.6))  ///
	ylabel(0 (0.2) 0.6 ,format(%03.1f)  angle(horizontal)) ///
	xlabel(2000 (4) 2016) ///
	xtitle("Year",margin(t=3)) ytitle("Planting rate",margin(r=3)) ///
	saving($dt\gestacking.gph, replace)   ///
	leg(region(lwidth(none) color(none)) ring(0) position(11) col(1)  size(medsmall) ///
	lab(4 "RW") lab(3 "RW+CB") lab(2 "RW+HT") lab(1 "RW+CB+HT"))

graph dis,	xsize(6) ysize(5.5) scale(0.9)
graph export "$do\figS8.tif", wid(2000) replace  
graph export "$do\figS8.png", wid(2000) replace  	


*--------fig. S9 more recent root injury data----------------*
use $dd/estim.dta, clear
keep if soilinsecticide==0 &seedtrtrate==0  
keep if stateabbr=="IA"|stateabbr=="IL"|stateabbr=="IN"|stateabbr=="NE"|stateabbr=="SD"
drop if bt==.

collapse (mean) mean= rw  ///
	, by(bt stateabbr year locid protein)
sort year locid bt protein 

drop if protein=="mCry3A+eCry3.1Ab"
gen bttype=0
replace bttype=1 if protein=="Cry3Bb1"
replace bttype=2 if protein=="Cry34/35Ab1"|protein==" Cry34/35Ab1"
replace bttype=3 if protein=="mCry3A"
replace bttype=4 if protein=="Cry3Bb1+Cry34/35Ab1"
replace bttype=5 if protein=="mCry3A+Cry34/35Ab1"

keep stateabbr year bttype mean locid
reshape wide mean, i(year stateabbr locid) j(bttype)
forv i=1/5{
	gen mean_abott`i'=100*(mean0-mean`i')/mean0
}

for var mean_abott*: drop if X<0 

egen mean_pooled=rowmean(mean1 mean2 mean3 mean4 mean5)
gen mean_abott=100*(mean0-mean_pooled)/mean0

collapse (mean) mean_abott* mean0, by(year)

loc r1 "219 049 036"
loc r2 "252 140 090"
loc r3 "255 223 146"
loc b3 "203 229 246"
loc b2 "34 128 191"
loc b1 "24 93 139"


graph twoway ///
	(scatter mean_abott1 year,mcolor("`r1'") msize(medium) msymbol(X)) ///
	(scatter mean_abott2 year,mcolor("`r2'") msize(medium) msymbol(T)) ///
	(scatter mean_abott3 year,mcolor("`r3'") msize(medium) msymbol(D)) ///
	(scatter mean_abott4 year,mcolor("`b3'") msize(medium) msymbol(D)) ///
	(scatter mean_abott5 year,mcolor("`b2'") msize(medium) msymbol(T)) ///
	(lowess mean_abott1 year,lcolor("`r1'") lw(0.6) lpattern(solid)) ///
	(lowess mean_abott2 year,lcolor("`r2'") lw(0.6) lpattern(solid)) ///
	(lowess mean_abott3 year,lcolor("`r3'") lw(0.6) lpattern(solid)) ///	
	(lowess mean_abott4 year,lcolor("`b3'") lw(0.6) lpattern(solid)) ///
	(lowess mean_abott5 year,lcolor("`b2'") lw(0.6) lpattern(solid)) ///
	, ///
	ylabel(0 (20) 100, labsize("7pt")) ///
	xlabel(2004 (4) 2020.5, labsize("7pt")) ///
	leg(region(lwidth(none)) on  ring(0) position(6) col(2) size("7pt") ///
	order(1 "Cry3Bb1" 4 "Cry3Bb1+Cry34/35A" 2 "Cry34/35A" 5 "Cry34/35A+mCry3A" 3 "mCry3A")) ///
	title("{bf:A}",size("7pt") ring(1) position(11) margin(b=3)) /// 
	xtitle("{bf:Trial year}", size("7pt") margin(t=3)) ///
	ytitle("{bf:Bt efficacy (%)}", margin(r=3) size("7pt") ) ///
	 nodraw
graph save $dt/efficacy_full.gph, replace


graph twoway ///
	(scatter mean0 year, mcolor(black) mfcolor(white) msize(medium) msymbol(O)) ///
	(lowess mean0 year,lcolor(black) lw(0.6) lpattern(solid) ) ///
	, ///
	ylabel(0 (0.2) 1, labsize("7pt") format(%03.1f)) ///
	xlabel(2004 (4) 2020.5, labsize("7pt")) ///
	leg(region(lwidth(none)) off  ring(0) position(6) col(2) size("7pt")) ///
	title("{bf:B}",size("7pt") ring(1) position(11) margin(b=3)) /// 
	xtitle("{bf:Trial year}", size("7pt") margin(t=3)) ///
	ytitle("{bf:Non-Bt injury}", margin(r=3) size("7pt") )	///
	 nodraw
graph save $dt/injury_full.gph, replace


graph combine  $dt/efficacy_full.gph $dt/injury_full.gph, row(1)  iscale(1) ///
	xsize(5) ysize(2.6)

graph export $do/figS9.tif, replace width(1425) //horizontal pixels=4.75inch*300dots per inch
graph export $do/figS9.png, replace width(1425)	


*-----------------table. S1 cross validation between TraitTrak and USDA NASS---------------------*
//obtain TraitTrak rates
import delimited "$ds\ISUTraitTrakData_NEW.csv", clear
keep if crop=="Corn"

gen bt=(strpos(seedtraitcollapsed, "RW")>0|strpos(seedtraitcollapsed, "CB")>0)
gen ht=(strpos(seedtraitcollapsed, "Gly Tol")>0|strpos(seedtraitcollapsed, "LL")>0|strpos(seedtraitcollapsed, "HT other")>0)
gen btonly=(bt==1&ht==0)
gen htonly=(bt==0&ht==1)
gen stacked=(bt==1&ht==1)

gl vlist "bt ht btonly htonly stacked "
for var $vlist :gen X_acres=acres*X

rename yearcode year
rename county countyfips
merge m:1 countyfips using $ds\state_crd_county.dta, keep(match master) nogen

collapse (sum) *acres, by(year statefips state)

foreach v of global vlist{
	gen `v'_rate_TT=`v'_acres/acres //Bt rate in Kynetec TraitTrak data.
}

save $dt/bt_TraitTrak.dta, replace

//obtain NASS rates
import delimited "$ds/BiotechCropsAllTables2022.csv", clear
gen variety="stacked" if attribute=="Stacked gene varieties (percent of all corn planted)"
replace variety="btonly" if attribute=="Insect-resistant (Bt) only (percent of all corn planted) "
replace variety="htonly" if attribute=="Herbicide-tolerant (HT) only (percent of all corn planted) "
keep if variety!=""
destring value, force replace
gen rate=value/100
keep year state rate variety
reshape wide rate, i(year state) j(variety) string
gen ratebt=ratebtonly+ratestacked
gen rateht=ratehtonly+ratestacked

//merge
merge 1:1 year state using $dt/bt_TraitTrak.dta, keep(matched)

//correlation
corr ratebt bt_rate_TT
corr rateht ht_rate_TT



