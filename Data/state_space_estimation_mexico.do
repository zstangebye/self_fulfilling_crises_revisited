
clear all
set more off

import excel using mexico_data, sh("Sheet2") cellrange(A2:B164) 

rename B gdp_us

gen year=substr(A,1,4)
destring year, replace
gen quarter=substr(A,-1,1)
destring quarter, replace

drop A

gen time=yq(year,quarter)

format time %tq

tsset time

gen ln_y_us=ln(gdp_us)

keep time ln_y_us


*Benchmark model:  g is AR(1), z is iid

constraint 1 [D.ln_y_us]g=1
constraint 2 [D.ln_y_us]z=1
constraint 3 [D.ln_y_us]lagz=-1
constraint 4 [lagz]L.z=1

 
 sspace (g L.g e.g, state ) (z e.z, state noconstant) (lagz L.z, state noconstant) ///
 (D.ln_y_us g z lagz, noconstant ), constraints(1/4) covstate(identity)
