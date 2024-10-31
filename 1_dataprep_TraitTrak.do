********************************************************************************
*				program: TraitTrak data prep  	   							   *
********************************************************************************
clear all

global root "..."

global data "$root/1_data"
global ds "$data/source"
global dd "$data/data"
global dp "$data/proc"
global dt "$data/temp"
global do "$data/out"



*--------------------------below-ground (RW) trait data-----------------------------*
import delimited "$ds\ISUTraitTrakData_NEW.csv", clear
keep if crop=="Corn"


gen rw=(strpos(seedtraitcollapsed, "RW")>0)
gen cb=(strpos(seedtraitcollapsed, "CB")>0)
gen ht=(strpos(seedtraitcollapsed, "Gly Tol")>0|strpos(seedtraitcollapsed, "LL")>0|strpos(seedtraitcollapsed, "HT other")>0)
gen rwonly=(seedtraitcollapsed=="RW")

gen rwcbonly=(rw==1&cb==1&ht==0)
gen rwhtonly=(rw==1&cb==0&ht==1)
gen rwcbht=(rw==1&cb==1&ht==1) //rw cb and ht
gen cbht=(rw==0&cb==1&ht==1) // cb and ht, no rw

global bt "rw cb ht rwonly rwcbonly rwhtonly rwcbht cbht"
for var $bt :gen X_acres=acres*X

rename yearcode year
rename county countyfips
merge m:1 countyfips using $ds\state_crd_county.dta, keep(match master) nogen
gen nation=1 //national level flag

local admin "countyfips scrd statefips nation"

//coverage 
foreach a of local admin{
	preserve
	collapse (sum) *acres, by(year `a')
	for var *acres: rename X X_`a'
	
	foreach v of global bt{
		gen `v'_rate_`a'=100*`v'_acres_`a'/acres_`a'
	}

	label var rw_rate_`a' "planting rate of seed containing RW trait (%), `a' level"

	save $dt\bt_`a'.dta, replace
	restore
}

//seed price: expenditures per acre
gen east=1 if stateabbr=="IN"|stateabbr=="MI"|stateabbr=="OH"
replace east=0 if stateabbr=="ND"|stateabbr=="SD"|stateabbr=="NE"|stateabbr=="MN"|stateabbr=="IA"|stateabbr=="IL"

collapse (sum) expenditures acres units, by(rw year statefips)
gen ps_acre_statefips=expenditures/acres
reshape wide ps_acre_statefips, i(year statefips) j(rw)
save $dt\ps_statefips.dta, replace


*--------------------------lag for rw trait-------------------------------------*
local admin "countyfips scrd statefips"
foreach a of local admin{
use $dt/bt_`a'.dta, clear
keep year `a' rw_rate_`a' rw_acres_`a' acres_`a'

reshape wide rw_rate_`a' rw_acres_`a' acres_`a', i(`a') j(year)
forv i=2000/2002{
	replace rw_rate_`a'`i'=0
	replace rw_acres_`a'`i'=0

}


gen rw_rate_`a'2017=. //placeholder for 2-year lags, because the TraitTraik data end up in 2016
gen rw_rate_`a'2018=.
gen rw_rate_`a'2019=.
reshape long rw_rate_`a' rw_acres_`a' acres_`a', i(`a') j(year)

xtset `a' year
forv i=1/3{
	gen lag`i'_rw_rate_`a'=l`i'.rw_rate_`a'
	gen lag`i'_rw_acres_`a'=l`i'.rw_acres_`a'
	gen lag`i'_acres_`a'=l`i'.acres_`a'
}

keep if year>=2001
save $dt/lag_rw_rate_`a'.dta, replace
}


*--------------------------decompose RW trait data-----------------------------*
import delimited "$ds\ISUTraitTrakData_NEW.csv", clear
keep if crop=="Corn"
tab seedtrait if strpos(seedtraitcollapsed, "RW")>0
asdoc tab seedtrait if strpos(seedtraitcollapsed, "RW")>0

global traits "Cry3Bb1 mCry3A Cry3435A eCry31Ab"
foreach t of global traits {
	gen tr_`t'=0
}

replace tr_Cry3Bb1=1 if strpos(seedtraitcollapsed, "RW")>0 & ///
	(strpos(seedtrait, "SmartStax")>0 |strpos(seedtrait, "VT3")>0|strpos(seedtrait, "YG")>0) 

replace tr_mCry3A=1 if strpos(seedtraitcollapsed, "RW")>0 & ///
	(strpos(seedtrait, "Agrisure")>0 |seedtrait=="Optimum AcreMax TRIsect"| ///
	seedtrait=="Optimum AcreMax Xtreme"|seedtrait=="Optimum TRIsect") 
	
replace tr_Cry3435A=1 if strpos(seedtraitcollapsed, "RW")>0 & ///
	(strpos(seedtrait, "HX")>0 ///
	|seedtrait=="Agrisure 3122 E-Z Refuge"|seedtrait=="Optimum AcreMax Xtreme" ///
	|seedtrait=="DAS SmartStax"|seedtrait=="Genuity SmartStax"|seedtrait=="Optimum AcreMax 1" ///
	|seedtrait=="Optimum AcreMax RW"|seedtrait=="Optimum AcreMax Xtra" ///
	|seedtrait=="Optimum Intrasect Xtra"|seedtrait=="SmartStax RA"|seedtrait=="SmartStax RIB" ///
	) 

replace tr_eCry31Ab=1 if strpos(seedtraitcollapsed, "RW")>0 & ///
	(strpos(seedtrait, "Agrisure Duracade")>0) 

foreach t of global traits {
	tab year if tr_`t'==1
}

//single and stacking
gen only_Cry3Bb1=(tr_Cry3Bb1==1&tr_mCry3A==0&tr_Cry3435A==0&tr_eCry31Ab==0)
gen only_mCry3A=(tr_Cry3Bb1==0&tr_mCry3A==1&tr_Cry3435A==0&tr_eCry31Ab==0)
gen only_Cry3435A=(tr_Cry3Bb1==0&tr_mCry3A==0&tr_Cry3435A==1&tr_eCry31Ab==0)
gen only_eCry31Ab=(tr_Cry3Bb1==0&tr_mCry3A==0&tr_Cry3435A==0&tr_eCry31Ab==1)

gen two_Cry3Bb1_mCry3A=(tr_Cry3Bb1==1&tr_mCry3A==1&tr_Cry3435A==0&tr_eCry31Ab==0)
gen two_Cry3Bb1_Cry3435A=(tr_Cry3Bb1==1&tr_mCry3A==0&tr_Cry3435A==1&tr_eCry31Ab==0)
gen two_Cry3Bb1_eCry31Ab=(tr_Cry3Bb1==1&tr_mCry3A==0&tr_Cry3435A==0&tr_eCry31Ab==1)

gen two_mCry3A_Cry3435A=(tr_Cry3Bb1==0&tr_mCry3A==1&tr_Cry3435A==1&tr_eCry31Ab==0)
gen two_mCry3A_eCry31Ab=(tr_Cry3Bb1==0&tr_mCry3A==1&tr_Cry3435A==0&tr_eCry31Ab==1)
gen two_Cry3435A_eCry31Ab=(tr_Cry3Bb1==0&tr_mCry3A==0&tr_Cry3435A==1&tr_eCry31Ab==1)

gen tr_single=(only_Cry3Bb1==1|only_mCry3A==1|only_Cry3435A==1|only_eCry31Ab==1)
gen tr_double=0
replace tr_double=1 if strpos(seedtraitcollapsed, "RW")>0 & ///
	(seedtrait=="Optimum AcreMax Xtreme" |seedtrait=="Agrisure 3122 E-Z Refuge"| ///
	seedtrait=="Agrisure Duracade 5122"|seedtrait=="Agrisure Duracade 5222"| ///
	seedtrait=="Agrisure Duracade 5122A"|seedtrait=="Agrisure Duracade 5222A"| ///
	seedtrait=="SmartStax RIB" |seedtrait=="Genuity SmartStax"| ///
	seedtrait=="SmartStax RA" |seedtrait=="DAS SmartStax") //dataprep_v2.do included "HX" which is incorrect.


 
//calculate year average 
rename yearcode year
rename county countyfips
merge m:1 countyfips using $ds\state_crd_county.dta, keep(match master) nogen

for var only_* two_* tr_* :gen X_acres=acres*X

collapse (sum) *acres, by(year scrd)

for var acres only_*_acres two_*_acres tr_*_acres : rename X X_scrd


save $dt\rw_traits.dta, replace






