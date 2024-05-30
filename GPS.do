**************************************************
* Chinese Equity Markets
* GPS model
* Rafael Perez Ribas
* Version: 01/27/2012
**************************************************

*************************************
* Packages
*************************************

net install st0120.pkg

*************************************
* System setting
*************************************

clear all
set matsize 5000
set mem 400m
set more off

cap log close
log using log\GPS.log, replace

log off

*************************************
* Database
*************************************

u data\china if (year==2004 | year==2003) & batch==0, clear 

xtset listcode year

drop if year==2003

*****************************************************
* DURATION SETUP
*****************************************************

cap drop dur
g dur = ym(year_reform, month_reform)
format dur %tmMonYY

g fail = 1

stset dur, nos o(time 543)

log on

*****************************************************
* Generalized Propensity Score
*****************************************************

mvrs stcox l_cape l_employ l_salesK w_roe w_mb synch nrpt l_orec l_sharehold ///
    manager_share nontradable state_owner l_share state_share inst_share age ///
	l_asset w_cfta l_Klabor w_debt bankloan w_pe w_tobinq repr_ind l_pgdp ///
	l_ind_sale ind_herf cash, ///
	nohr degree(3) dfdefault(5) alpha(.1) cycles(10)
est store pscox

log off

* PROBABILITY OF BEING IN THE CONTROL GROUP (TREATED AFTER JUNE 2006)

cap drop xb
predict xb, xb
cap predict S0, bases
cap predict h0, basehc


* PROBABILITY OF JOING THE REFORM IN ITS OWN PERIOD

g GPS = exp(xb)*h0*(S0^exp(xb))

g t2 = 1/(GPS^(1/2))
winsor t2, p(0.01) g(w2)


* GPS Index

g sGPS = exp(xb)

mata

    A = st_data(., "dur")

    J = J(st_nobs(), 1, 1)

    B = abs((A*J') - (A*J')')

    I = B:>6
    G = st_data(., "sGPS")
    Gt = (G*J')'
    D = I:*Gt
    T = G:/I
    C = G:>=rowmin(T) :& G:<=rowmax(D)
    st_matrix("C6", C)

    I = B:>12
    G = st_data(., "sGPS")
    Gt = (G*J')'
    D = I:*Gt
    T = G:/I
    C = G:>=rowmin(T) :& G:<=rowmax(D)
    st_matrix("C12", C)

end

svmat C6
svmat C12

qui su sGPS if C61==1
replace C61 = sGPS<r(max)

qui su sGPS if C121==1
replace C121 = sGPS<r(max)

loc m = ym(2007,2)
loc i = ym(2005,9)

tw sc sGPS dur if C61==1 & dur<=ym(2007,2), mcol(red) mstyle(p11) || ///
   sc sGPS dur if C61==0 & dur<=ym(2007,2), mcol(blue) mstyle(p4) || ///
   , scheme(s1mono) ylab(,nogrid angle(0)) ///
   leg(lab(1 "On the common support") lab(2 "Out of common support") ///
   pos(1) ring(0) col(1)) ///
   yti("GPS Index") xti("Treatment Assignment") ///
   xlab(`i'(2)`m', angle(90) labs(small)) ///
   saving(log\GPS, replace) name(GPS, replace)
graph export GPS.eps, as(eps) preview(off) replace name(GPS)


* SAVING TEMPORARY FILE

sort listcode

tempfile ps
sa `ps', replace

u data\china if batch==0, clear

sort listcode
merge listcode using `ps', keep(GPS sGPS C61 C121 w2 t2 t3)
ta _merge
drop _merge

g dose = max(0, ym(year, 12) - ym(year_reform, month_reform))
g sel = max(0, ym(year_reform, month_reform) - ym(2004, 12))


xtset listcode year

foreach v in l_fcapital l_employ l_netinc l_salesK w_mb w_roe w_tobinq {
    g d`v' = `v' - l.`v'
}

log on

*****************************************************
* Balance Property
*****************************************************

* GPS - Bandwidth = 6 months

mat B = .,.,.,.,.,.,.,.
mat coln B = unmatched p param p IPW p both p

mkspline G = sGPS, cubic nk(4)

qui foreach var in l_fcapital dl_fcapital l_employ dl_employ l_salesK dl_salesK ///
    l_netinc dl_netinc w_roe dw_roe w_mb dw_mb lr l_turnover synch w_div2 am_eqiss ///
    ma_num nrpt l_orec herfin_5 l_sharehold manager_share l_ooe ceoturn ///
	nontradable state_owner l_share state_share inst_share age l_asset l_sales ///
	w_cfta l_Klabor w_debt bankloan cash w_pe w_tobinq repr_ind repr_prov ///
	l_pgdp l_ind_sale ind_herf {

    reg `var' sel if year==2004, r
    loc b1 = _b[sel]
    loc p1 = ttail(e(df_r), abs(_b[sel]/_se[sel]))*2

    reg `var' sel G1-G3 if year==2004 & C61==1, r
    loc b2 = _b[sel]
    loc p2 = ttail(e(df_r), abs(_b[sel]/_se[sel]))*2

    reg `var' sel [aw=w2] if year==2004 & C61==1, r
    loc b3 = _b[sel]
    loc p3 = ttail(e(df_r), abs(_b[sel]/_se[sel]))*2

    reg `var' sel G1-G3 [aw=w2] if year==2004 & C61==1, r
    loc b4 = _b[sel]
    loc p4 = ttail(e(df_r), abs(_b[sel]/_se[sel]))*2

    mat A = `b1', `p1', `b2', `p2', `b3', `p3', `b4', `p4'
    mat rown A = `var'
    mat B = B \ A
}

mat B = B[2...,1...]

drop G1-G3

ta C61 if year==2004

mat list B, f(%9.3f)


*****************************************************
* Estimating Continuous Treatment Effect Before Treatment
*****************************************************

g cluster = yh(year_reform,halfyear(mdy(month_reform,1,year_reform)))
	
preserve

	* Cubic function

	cap drop Q*
	g Q0 = dose>0


	* Mean-centered spline GPS

	mkspline G = sGPS, cubic nk(4)

	qui forvalues i = 1/3 {
		su G`i' if year==2004 & C61==1
		replace G`i' = G`i' - r(mean)
	}

	qui forvalues i = 1/3 {
		g H0_`i' = Q0*G`i'
	}

	ta year, g(y)
	drop y1

	qui foreach var of varlist y2-y8 {
		su `var'
		replace `var' = `var' - r(mean)
	}

	* Panel setting

	xtset listcode year

	cap g dl_sales = l_sales - l.l_sales
	cap g dl_fcapital = l_fcapital - l.l_fcapital
	cap g dl_employ = l_employ - l.l_employ

	replace ma_num = ma_num>0 if ma_num!=.
	replace am_eqiss = am_eqiss>0 & am_eqiss!=.

	mat R = J(1,6,.z)
	mat coln R = FE se p IPW2 se p

	local titles "DeltaK DeltaL Sales/K NetIncome ROE M/B LiqRatio ShareTurnover PriceSync Dividend Leverage Issuance M\&A RPT InterLoans OwnerConcent ShareHolders ManagerShares OtherExp CEOTurnover"

	local t = 0

	qui foreach var in l_fcapital l_employ l_salesK l_netinc w_roe ///
		w_mb lr l_turnover synch w_div2 w_debt am_eqiss ma_num nrpt ///
		l_orec herfin_5 l_sharehold manager_share l_ooe ceoturn {

		local t = `t' + 1

		local title = word("`titles'",`t')


		* Fixed-Effect

		xtreg `var' f.Q0 y2-y8 if dose==0, vce(cl listcode) fe

		test f.Q0
		mat A = _b[f.Q0], _se[f.Q0], r(p)
		
		mat rown A = "`title'"

		
		* IPW

		reg d.`var' f.Q0 y2-y8 G1-G3 if C61==1 & dose==0 & year>=2004 ///
			[aw=w2], vce(cl listcode)

		test f.Q0
		mat B = _b[f.Q0], _se[f.Q0], r(p)

		mat A = A, B
		
		
		* Result Matrix

		mat R = R \ A
	}

	mat R = R[2...,1...]

	log on

	mat list R, format(%9.4f)

restore

*****************************************************
* Estimating Continuous Treatment Effect
*****************************************************

log off

preserve

	tempfile gps

	sort listcode year

	sa `gps', replace

	u data\quarter3, clear

	sort listcode year quarter

	drop lr l_turnover

	merge listcode year using `gps', keep(year_reform month_reform ///
		sGPS GPS C61 w2 l_employ w_div2 lr l_turnover ///
		herfin_5 synch l_sharehold manager_share l_ooe cluster)

	keep if _merge==3
	drop _merge

	qui for @ in var l_employ w_div2 herfin_5 lr l_turnover ///
		synch l_sharehold manager_share l_ooe: replace @ = . if quarter!=4


	g dose = max(0, ym(year, quarter*3) - ym(year_reform, month_reform))

	ta year, g(y)
	drop y1

	ta quarter, g(q)
	drop q1

	qui foreach var of varlist y2-y8 q2-q4 {
		su `var' if C61==1
		replace `var' = `var' - r(mean)
	}


	* Cubic function

	cap drop Q*
	mkspline Q = dose, cubic knots(6 12 24 30 36)


	* Reference values

	qui forvalues i = 12(12)36 {
		mean Q* if dose==`i'
		forvalues j = 1/4 {
			loc Q`j'_`i' = _b[Q`j']
		}
	}

	mat qb = J(1,4,.)
	qui forvalues i = 1/39 {
		mean Q* if dose==`i'
		mat A = e(b)
		mat rown A = r`i'
		mat qb = qb \ A
	}
	mat qb = qb[2...,1...]


	* Mean-centered spline GPS

	mkspline G = sGPS, cubic nk(4)

	qui forvalues i = 1/3 {
		su G`i' if year==2004 & quarter==4 & C61==1
		replace G`i' = G`i' - r(mean)
	}

	qui forvalues j = 1/4 {
		forvalues i = 1/3 {
			g H`j'_`i' = Q`j'*G`i'
		}
	}

	* Panel setting

	g period = yq(year, quarter)

	xtset listcode period

	cap drop dl_*
	
	g dl_fcapital = f4.l_fcapital - l_fcapital
	g dl_employ = f4.l_employ - l_employ

	mat R = J(1,2,.)
	mat coln R = IPW se

	local titles "{&Delta}K {&Delta}L Sales/K NetIncome ROE M/B LiqRatio ShareTurnover PriceInfo Dividend InterLoans OwnerConcent ShareHolders ManagerShares OtherExp Leverage"

	local t = 0

	qui foreach var in dl_fcapital dl_employ l_salesK l_netinc w_roe ///
	    w_mb lr l_turnover synch w_div2 l_orec herfin_5 l_sharehold ///
		manager_share l_ooe w_debt {

		local t = `t' + 1

		local title = word("`titles'",`t')

		local ylab "#6"


		* IPW + Parametric Adjustment

		xtreg `var' Q1-Q4 G1-G3 H* y2-y8 q2-q4 if C61==1 [aw=w2], ///
			vce(cl listcode) fe

		mat qv = e(V)
		mat qv = qv[1..4,1..4]

		mat b = e(b)
		mat b = b[1,1..4]

		forvalues i = 12(12)36 {
			lincom Q1*`Q1_`i'' + Q2*`Q2_`i'' + Q3*`Q3_`i'' + Q4*`Q4_`i''
			loc E`i' = r(estimate)
			loc S`i' = r(se)
		}
		mat B = `E12', `S12' \ `E24', `S24' \ `E36', `S36'
		
		mat rown B = `var'12 `var'24 `var'36

		mat R = R \ B

		mat GV = vecdiag(qb*qv*qb')'
		mat GB = qb*b'

		svmat GV
		svmat GB

		g LC = GB1 - 1.645*sqrt(GV1)
		g UC = GB1 + 1.645*sqrt(GV1)
		g TQ = _n if LC!=.

		tw line GB1 TQ if TQ<=30, lc(blue) || ///
		   line LC TQ if TQ<=30, lp(shortdash) lc(blue) || ///
		   line UC TQ if TQ<=30, lp(shortdash) lc(blue) || ///
           , scheme(s1mono) ylab(,nogrid) ///
           leg(off) xti("Z") ysize(4) xsize(5.5) ///
           xlab(1 6(6)30) ylab(`ylab') yline(0) ti("`title'") ///
           saving(log\g`var', replace) name(`var', replace) nodraw
        graph export `var'.eps, as(eps) preview(off) replace name(`var')

        drop GV1-TQ
	}

	mat R = R[2...,1...]

	log on

	mat list R, format(%9.4f)

restore


*****************************************************
* Other Financial Outcomes
*****************************************************

preserve

	ta year, g(y)
	drop y1

	qui foreach var of varlist y2-y8 {
		su `var'
		replace `var' = `var' - r(mean)
	}

	* Cubic function

	mkspline Q = dose, cubic knots(6 12 24 30 36)


	* Reference values

	qui forvalues i = 12(12)36 {
		mean Q* if dose==`i'
		forvalues j = 1/4 {
			loc Q`j'_`i' = _b[Q`j']
		}
	}

	mat qb = .,.,.,.
	qui forvalues i = 1/39 {
		mean Q* if dose==`i'
		mat A = e(b)
		mat rown A = r`i'
		mat qb = qb \ A
	}
	mat qb = qb[2...,1...]


	* Mean-centered spline GPS

	mkspline G = sGPS, cubic nk(4)

	qui forvalues i = 1/3 {
		su G`i' if year==2004 & C61==1
		replace G`i' = G`i' - r(mean)
	}

	qui forvalues j = 1/4 {
		forvalues i = 1/3 {
			g H`j'_`i' = Q`j'*G`i'
		}
	}

	* Panel setting

	xtset listcode year
	
	mat qbb = qb, J(39,22,0), J(39,1,1)

	local titles "M&A RPTs Issuance CEOTurnover"

	loc t = 0

	qui foreach var of varlist ma_num nrpt am_eqiss ceoturn {

		local t = `t' + 1

		local title = word("`titles'",`t')

		if "`var'" == "am_eqiss" | "`var'" == "ceoturn" | ///
			"`var'" == "ma_num" {
			logit `var' Q* G1-G3 H* y2-y8 if C61==1 [pw = w2] ///
			    , vce(cl listcode)

			mat qv = e(V)

			mat b = e(b)

			mat GB = qbb*b'

			svmat GB

			g fGB = invlogit(GB1)*(1-invlogit(GB1))
			g FGB = invlogit(GB1) - invlogit(_b[_cons])

			mkmat fGB, mat(fGB) nomis

			mat T = qbb#fGB

			mat A = T["r1:r1",1...]

			forvalues i = 2/39 {
				mat A = A \ T["r`i':r`i'",1...]
			}

			mat C = J(39,26,0), J(39,1,invlogit(_b[_cons])*(1-invlogit(_b[_cons])))

			mat A = A - C

			mat GV = vecdiag(A*qv*A')'

			svmat GV

			g LC = FGB - 1.645*sqrt(GV1)
			g UC = FGB + 1.645*sqrt(GV1)
			g TQ = _n if LC!=.

			drop GB1 fGB

		}
		else if "`var'" == "nrpt" {
			poisson `var' Q* G1-G3 H* y2-y8 if C61==1 [pw = w2] ///
			    , vce(cl listcode)

			mat qv = e(V)

			mat b = e(b)

			mat GB = qbb*b'

			svmat GB

			g fGB = exp(GB1)
			g FGB = exp(GB1) - exp(_b[_cons])

			mkmat fGB, mat(fGB) nomis

			mat T = qbb#fGB

			mat A = T["r1:r1",1...]

			forvalues i = 2/39 {
				mat A = A \ T["r`i':r`i'",1...]
			}

			mat C = J(39,26,0), J(39,1,exp(_b[_cons]))

			mat A = A - C

			mat GV = vecdiag(A*qv*A')'

			svmat GV

			g LC = FGB - 1.645*sqrt(GV1)
			g UC = FGB + 1.645*sqrt(GV1)
			g TQ = _n if LC!=.

			drop GB1 fGB

	   }

		tw line FGB TQ if TQ<=30, lc(blue) || ///
		   line LC TQ if TQ<=30, lp(shortdash) lc(blue) || ///
		   line UC TQ if TQ<=30, lp(shortdash) lc(blue) || ///
		   , scheme(s1mono) ylab(,nogrid) ///
		   leg(off) xti("Z") ysize(4) xsize(5.5) ///
		   xlab(1 6(6)30) ylab(#6 0) yline(0) ti("`title'") ///
		   saving(log\g`var', replace) name(`var', replace) nodraw
		graph export `var'.eps, as(eps) preview(off) replace name(`var')

		drop FGB-TQ
	}

restore

gr combine lr l_turnover, iscale(1) ///
    scheme(s1mono) c(2) ysize(1) xsize(2) ///
    saving(log\liquidity, replace) name(liquidity, replace)
graph export liquidity.eps, as(eps) preview(off) replace name(liquidity)

gr combine dl_fcapital dl_employ l_salesK l_netinc ///
    w_roe w_mb, scheme(s1mono) c(2) ysize(7.5) xsize(5) ///
    saving(log\effects, replace) name(effects, replace)
graph export effects.eps, as(eps) preview(off) replace name(effects)

gr combine am_eqiss w_debt w_div2, iscale(1.1) ///
    scheme(s1mono) c(3) ysize(1) xsize(3) ///
    saving(log\financial, replace) name(financial, replace)
graph export financial.eps, as(eps) preview(off) replace name(financial)

gr combine manager_share ceoturn ///
    , iscale(1) scheme(s1mono) c(2) ysize(1) xsize(2) ///
    saving(log\manager, replace) name(manager, replace)
graph export manager.eps, as(eps) preview(off) replace name(manager)

gr combine herfin_5 l_sharehold nrpt l_orec ///
    , scheme(s1mono) c(2) ysize(2) xsize(2) iscale(.5) ///
    saving(log\others, replace) name(others, replace)
graph export others.eps, as(eps) preview(off) replace name(others)

graph export synch.eps, as(eps) preview(off) replace name(synch)
graph export ma_num.eps, as(eps) preview(off) replace name(ma_num)

*****************************************************
* End of Do file
*****************************************************

log close
