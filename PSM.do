**************************************************
* Chinese Equity Markets
* Propensity Score Matching
* Rafael Perez Ribas
* Version: 01/27/2012
**************************************************

*************************************
* Packages
*************************************

net install psmatch2.pkg
net install st0120.pkg
net install st0072.pkg

*************************************
* System setting
*************************************

clear all
set matsize 5000
set mem 400m
set more off

cap log close
log using log\PSM.log, replace

log off

*************************************
* Database
*************************************

u data\china if (year==2003 | year==2004), clear


* TREATMENT VARIABLE (2005)

g treat = batch==0 & ym(year_reform, month_reform)<=ym(2006,6)


* TAKING DIFFERENCES

xtset listcode year

foreach v in l_fcapital l_employ l_netinc l_salesK w_mb w_roe w_tobinq {
    g d`v' = `v' - l.`v'
}

drop if year==2003


* DURATION OUT OF TREATMENT

cap drop dur
g dur = ym(year_reform, month_reform)
format dur %tmMonYY

g fail = 1


log on


*************************************
* Summary Statistics
*************************************

* TOTAL

mat T = .,.,.
mat coln T = total se n

qui foreach var in l_fcapital dl_fcapital l_employ dl_employ l_salesK dl_salesK l_netinc dl_netinc ///
    w_roe dw_roe w_mb dw_mb lr l_turnover am_eqiss w_div2 synch l_sharehold ma_num manager_share ceoturn ///
	herfin_5 nrpt l_orec l_ooe nontradable state_owner l_share state_share inst_share age l_asset l_sales ///
	w_cfta l_Klabor w_debt bankloan cash w_pe w_tobinq repr_ind repr_prov l_pgdp l_ind_sale ind_herf {

    mean `var'
    loc b = _b[`var']
    loc s = _se[`var']
    loc n = e(N)

    mat A = `b', `s', `n'
    mat rown A = `var'
    mat T = T \ A
}

mat T = T[2...,1...]

mat list T


* DIFFERENCES IN 2004: Pilot vs. others

mat Sb = .,.,.,.,.,.,.,.
mat coln Sb = batch se n control se n diff p-value

qui foreach var in l_fcapital dl_fcapital l_employ dl_employ l_salesK dl_salesK l_netinc dl_netinc ///
    w_roe dw_roe w_mb dw_mb lr l_turnover am_eqiss w_div2 synch l_sharehold ma_num manager_share ceoturn ///
	herfin_5 nrpt l_orec l_ooe nontradable state_owner l_share state_share inst_share age l_asset l_sales ///
	w_cfta l_Klabor w_debt bankloan cash w_pe w_tobinq repr_ind repr_prov l_pgdp l_ind_sale ind_herf {

    ttest `var' if batch==1 | ym(year_reform, month_reform)>ym(2006,6), by(batch)
    loc b1 = r(mu_2)
    loc s1 = r(sd_2)/sqrt(r(N_2))
    loc n1 = r(N_2)
    loc b0 = r(mu_1)
    loc s0 = r(sd_1)/sqrt(r(N_1))
    loc n0 = r(N_1)
    loc bd = r(mu_2) - r(mu_1)
    loc sd = r(p)

    mat A = `b1', `s1', `n1', `b0', `s0', `n0', `bd', `sd'
    mat rown A = `var'
    mat Sb = Sb \ A
}

mat Sb = Sb[2...,1...]

mat list Sb


* DIFFERENCES IN 2004: treated in 2005 vs. others

mat St = .,.,.,.,.,.,.,.
mat coln St = treated se n control se n diff p-value

qui foreach var in l_fcapital dl_fcapital l_employ dl_employ l_salesK dl_salesK l_netinc dl_netinc ///
    w_roe dw_roe w_mb dw_mb lr l_turnover am_eqiss w_div2 synch l_sharehold ma_num manager_share ceoturn ///
	herfin_5 nrpt l_orec l_ooe nontradable state_owner l_share state_share inst_share age l_asset l_sales ///
	w_cfta l_Klabor w_debt bankloan cash w_pe w_tobinq repr_ind repr_prov l_pgdp l_ind_sale ind_herf {

    ttest `var' if batch==0, by(treat)
    loc b1 = r(mu_2)
    loc s1 = r(sd_2)/sqrt(r(N_2))
    loc n1 = r(N_2)
    loc b0 = r(mu_1)
    loc s0 = r(sd_1)/sqrt(r(N_1))
    loc n0 = r(N_1)
    loc bd = r(mu_2) - r(mu_1)
    loc sd = r(p)

    mat A = `b1', `s1', `n1', `b0', `s0', `n0', `bd', `sd'
    mat rown A = `var'
    mat St = St \ A
}

mat St = St[2...,1...]

mat list St


log off

cap drop D*


*****************************************************
* DURATION SETUP
*****************************************************

stset dur, nos o(time 543)

replace _t = _t + 543
format _t %tmMonYY

sts, fail ti("") xti("") scheme(s1mono) ylab(,nogrid angle(0)) ///
    xlab(544(4)600, angle(90) labs(small)) noorig xline(548) ///
    saving(log\duration, replace) name(duration, replace)
graph export duration.eps, as(eps) preview(off) replace name(duration)

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
cap predict xb, xb
cap predict S0, bases

qui su S0 if _t==15

g PrC = r(mean)^exp(xb)


* PROBABILITY OF BEING IN THE BATCH

qui su S0 if _t==5

g Pr5 = 1 - (r(mean)^exp(xb))

g PS5 = Pr5/(Pr5 + PrC)

qui su PS5 if batch==1
loc min = r(min)
qui su PS5 if (year_reform==2006 & month_reform>6) | year_reform>2006
loc max = r(max)

cap g cPS5 = batch==0 | PS5<`max'

tw hist PS5 if ((year_reform==2006 & month_reform>6) | year_reform>2006), lc(green) freq fc(green) || ///
   hist PS5 if batch==1, lc(blue) fc(blue) freq fc(none) || ///
    , scheme(s1mono) ylab(,nogrid angle(0)) leg(lab(2 "Pilot Firms") lab(1 "Control Firms") ///
    order(2 1) pos(1) ring(0) col(1)) xline(`min' `max', lc(black) lp(dot)) yti("Frequency") ///
    xti("Propensity Score") ///
    saving(log\PS5, replace) name(PS5, replace)
graph export PSover.eps, as(eps) preview(off) replace name(PS5)

* SAVING TEMPORARY FILE

sort listcode

tempfile ps
sa `ps', replace


* ORIGINAL DATA

u data\china, clear

sort listcode
merge listcode using `ps', keep(treat PrC Pr5 PS5 cPS5)
ta _merge
drop _merge

g dose = max(0, ym(year, 12) - ym(year_reform, month_reform))
g sel = max(0, ym(year_reform, month_reform) - ym(2004, 12))


xtset listcode year

foreach v in l_fcapital l_employ l_netinc l_salesK ///
    w_mb w_roe w_tobinq {
    g d`v' = `v' - l.`v'
}


log on


*****************************************************
* PS Matching Estimation
*****************************************************

*****************************************************
* Balance Property
*****************************************************

qui {

	mat BP1 = J(1,5,.)
	mat BP2 = J(1,5,.)

	foreach v in l_fcapital dl_fcapital l_employ dl_employ l_salesK dl_salesK ///
	    l_netinc dl_netinc w_roe dw_roe w_mb dw_mb lr l_turnover synch w_div2 ///
		am_eqiss ma_num nrpt l_orec herfin_5 l_sharehold manager_share l_ooe ///
		ceoturn nontradable state_owner l_share state_share inst_share age l_asset ///
		l_sales w_cfta l_Klabor w_debt bankloan cash w_pe w_tobinq repr_ind ///
		repr_prov l_pgdp l_ind_sale ind_herf {

		* Pilot

		psmatch2 batch if cPS5==1 & year==2004 & ///
			(batch==1 | (year_reform==2006 & month_reform>6) | year_reform>2006) ///
			, out(`v') p(PS5) n(1)

		loc att = r(att)
		loc se = r(seatt)

		su `v' if batch==1 & cPS5==1 & year==2004

		loc m1 = r(mean)
		loc m0 = `m1' - `att'

		loc p = (1-normal(abs(`att'/`se')))*2

		mat A = `m1', `m0', `att', `se', `p'
		mat rown A = `v'

		mat BP1 = BP1 \ A

	}
	mat BP1 = BP1[2...,1...]
}

* Pilot

mat list BP1, f(%9.3f)


*****************************************************
* Estimating effect on Pilot firms
*****************************************************

* BETWEEN 2004/2005

preserve

keep if (year==2004 | year==2005) & ///
    (batch==1 | (year_reform==2006 & month_reform>6) | year_reform>2006)

qui {
    g t = year==2004

    mat R = J(1,6,.)
    mat coln R = OLS1 se OLS2 se NNM se
}

qui foreach var in l_fcapital l_employ l_salesK l_netinc w_roe w_mb ///
    lr l_turnover synch w_div2 w_debt am_eqiss dose {

    xtset listcode t
    g D`var' = l.`var' - `var' if year==2004

    * OLS
    reg D`var' batch, cl(listcode)
    loc b1 = _b[batch]
    loc s1 = _se[batch]

    * OLS with covariates
    reg D`var' batch nontradable state_owner l_share state_share inst_share age ///
        l_asset w_cfta l_Klabor w_debt bankloan cash w_pe w_tobinq repr_ind ///
		repr_prov l_pgdp ind_herf, cl(listcode)
    loc b2 = _b[batch]
    loc s2 = _se[batch]

    * NNM
    nnmatch D`var' batch PS5 if cPS5==1 & year==2004, tc(att) m(1) bias(bias) r(1)
    loc b3 = _b[SATT]
    loc s3 = _se[SATT]

    mat A = `b1', `s1', `b2', `s2', `b3', `s3'
    mat rown A = `var'
    mat R = R \ A
}

mat R = R[2...,1...]

ta batch cPS5 if year==2004 & w_roe!=. & PS5!=.

mat list R, f(%9.3f)

restore


* BETWEEN 2004/2006

preserve


keep if (year==2004 | year==2006) & ///
    (batch==1 | (year_reform==2006 & month_reform>6) | year_reform>2006)

qui {
    g t = year==2004

    mat R = J(1,6,.)
    mat coln R = OLS1 se OLS2 se NNM se
}

qui foreach var in l_fcapital l_employ l_salesK l_netinc w_roe w_mb ///
    lr l_turnover synch w_div2 w_debt am_eqiss dose {

    xtset listcode t
    g D`var' = l.`var' - `var' if year==2004

    * OLS
    reg D`var' batch, cl(listcode)
    loc b1 = _b[batch]
    loc s1 = _se[batch]

    * OLS with covariates
    reg D`var' batch nontradable state_owner l_share state_share inst_share age ///
        l_asset w_cfta l_Klabor w_debt bankloan cash w_pe w_tobinq repr_ind ///
		repr_prov l_pgdp ind_herf, cl(listcode)
    loc b2 = _b[batch]
    loc s2 = _se[batch]

    * NNM
    nnmatch D`var' batch PS5 if cPS5==1 & year==2004, tc(att) m(1) bias(bias) r(1)
    loc b3 = _b[SATT]
    loc s3 = _se[SATT]

    mat A = `b1', `s1', `b2', `s2', `b3', `s3'
    mat rown A = `var'
    mat R = R \ A
}

mat R = R[2...,1...]

ta batch cPS5 if year==2004 & w_roe!=. & PS5!=.

mat list R, f(%9.3f)

restore


* BETWEEN 2004/2007

preserve

keep if (year==2004 | year==2007) & ///
    (batch==1 | (year_reform==2006 & month_reform>6) | year_reform>2006)

qui {
    g t = year==2004

    mat R = J(1,6,.)
    mat coln R = OLS1 se OLS2 se NNM se
}

qui foreach var in l_fcapital l_employ l_salesK l_netinc w_roe w_mb ///
    lr l_turnover synch w_div2 w_debt am_eqiss dose {

    xtset listcode t
    g D`var' = l.`var' - `var' if year==2004

    * OLS
    reg D`var' batch, cl(listcode)
    loc b1 = _b[batch]
    loc s1 = _se[batch]

    * OLS with covariates
    reg D`var' batch nontradable state_owner l_share state_share inst_share age ///
        l_asset w_cfta l_Klabor w_debt bankloan cash w_pe w_tobinq repr_ind ///
		repr_prov l_pgdp ind_herf, cl(listcode)
    loc b2 = _b[batch]
    loc s2 = _se[batch]

    * NNM
    nnmatch D`var' batch PS5 if cPS5==1 & year==2004, tc(att) m(1) bias(bias) r(1)
    loc b3 = _b[SATT]
    loc s3 = _se[SATT]

    mat A = `b1', `s1', `b2', `s2', `b3', `s3'
    mat rown A = `var'
    mat R = R \ A
}

mat R = R[2...,1...]

ta batch cPS5 if year==2004 & w_roe!=. & PS5!=.

mat list R, f(%9.3f)

restore


* BETWEEN 2003/2004 (PLACEBO)

preserve

keep if (year==2003 | year==2004) & ///
    (batch==1 | (year_reform==2006 & month_reform>6) | year_reform>2006)

qui {
    g t = year==2003

    mat R = J(1,6,.)
    mat coln R = OLS1 se OLS2 se NNM se
}

qui foreach var in l_fcapital l_employ l_salesK l_netinc w_roe w_mb ///
    lr l_turnover synch w_div2 w_debt am_eqiss dose {

    xtset listcode t
    g D`var' = l.`var' - `var' if year==2003

    * OLS
    reg D`var' batch, cl(listcode)
    loc b1 = _b[batch]
    loc s1 = _se[batch]

    * OLS with covariates
    reg D`var' batch nontradable state_owner l_share state_share inst_share age ///
        l_asset w_cfta l_Klabor w_debt bankloan cash w_pe w_tobinq repr_ind ///
		repr_prov l_pgdp ind_herf, cl(listcode)
    loc b2 = _b[batch]
    loc s2 = _se[batch]

    * NNM
    nnmatch D`var' batch PS5 if cPS5==1 & year==2003, tc(att) m(1) bias(bias) r(1)
    loc b3 = _b[SATT]
    loc s3 = _se[SATT]

    mat A = `b1', `s1', `b2', `s2', `b3', `s3'
    mat rown A = `var'
    mat R = R \ A
}

mat R = R[2...,1...]

ta batch cPS5 if year==2003 & w_roe!=. & PS5!=.

mat list R, f(%9.3f)

restore


*****************************************************
* End of Do file
*****************************************************

log close
