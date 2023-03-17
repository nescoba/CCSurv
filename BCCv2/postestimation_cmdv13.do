 *3-3-2023: Transforming pprograms into a custom postestimation command
/*Remember to run this before running the program, after all it is a custom postestimation command for canon
 cd "C:\Users\raulcruz\OneDrive - Indiana University\Documents\students\carlos\stata"
  insheet using fit.csv, clear
  canon (weight waist pulse) (chins situps jumps) [weight=weight]
*/
*need to re-run mata drop to make stata forget the definition of the function, i.e. every time that the function is modified this command should be run
*  mata: mata drop calcpval() 
capture program drop csdcanon  
program csdcanon
   *ereturn list
   * preserving user's data
   preserve
   version 17
   args svysetcommand

   *saving original data set
    mkmat _all, matrix(All)
	
   tokenize `e(cmdline)', parse("()=]")
   /*display "1=|`1'|, 2=|`2'|, 3=|`3'|, 4=|`4'|, 5=|`5'|, 6=|`6'|, 7=|`7'|, 8=|`8'|, 9=|`9'|,10=|`10'|"*/
   /*tokens 3 and 6 also work with NYTS data*/
   /*
   di "`3'"
   di "`6'"
   di "`10'"
   */
   display("here")
   scalar varcount1 = `:word count `3''
   scalar varcount2 = `:word count `6''
   display("here5")
   /*First set of variables is the longest*/
   if varcount1 > varcount2 {
   	   display("here2")
       local X `3'
       local Y `6'
   }
   /*Second set of variables is the longest*/
  else {
  	   display("here3")
       local X `6'
       local Y `3'
   }
   display("here4")
   local finwgt `10'

   mkmat `X', matrix(OgX) 
   mkmat `Y', matrix(OgY) 
   mkmat `finwgt', matrix(W) 
   /*Creating diagonal matrix with the survey Weights*/
   matrix diag_W = I(rowsof(OgX))
   local n_f =rowsof(OgX)
   forvalues i = 1/`n_f' {
		 matrix diag_W [`i',`i']=  W[`i',1]
   }
 
   *standarixing X variables
   clear
   capture svmat double OgX
   foreach var of varlist _all {
   	egen aux_`var'= std(`var')
    drop `var'
   }
	*Not standarizing we don't get the correct canonical correlations
   mkmat  _all, matrix(X)

   *standarixing Y variables
   clear
   capture svmat double OgY
   foreach var of varlist _all {
   	egen aux_`var'= std(`var')
    drop `var'
   }
 	*Not standarizing we don't get the correct canonical correlations
   mkmat  _all, matrix(Y)  

  * mkmat sends the matrices to mata so there is no need to send them as arguments, e(rawcoef_var1) and e(rawcoef_var1) are sent to mata automatically when they are created by canon
  clear
  *capturing number of canonical correlations
  local n_cc  = e(n_cc)
  *creating matrix where all new p-values will be stored
   matrix Results = J(`n_cc',10,.) 
   mata: calcpval() 
   matrix colnames Results = CCIndex CC WilksLambda WilksLambdaFREQ RoysGreatestRoot PillaisTrace HotellingLawleyTrace WeightedReg CSDReg WeightedReg2
  
  capture svmat OGStataUVW
  capture svmat All, names(col) 
  *ereturn list
  *e(n_cc) =  3 has the number of canonical correlations; it is a scalar
 local n_cc  = e(n_cc)
 local weightindex = (2*e(n_cc)) + 1

 *=====Weighthed Simple Regression============
 
 quietly svyset  [pweight= OGStataUVW`weightindex'],

 *Notice coefficients of the simple linear regression are equal to the canonical correlations because the variances of the canonical variates are equal to 1
forval i = 1/`n_cc' { 
	local     secondindex = `i' + `n_cc'
	*summarize OGStataUVW`i'
	*summarize OGStataUVW`secondindex'
	quietly svy: regress OGStataUVW`i'  OGStataUVW`secondindex'
	matrix Results[`i',8]= `e(p)'
	quietly svy: regress OGStataUVW`secondindex' OGStataUVW`i'  
	matrix Results[`i',10]= `e(p)'
} 
 *=====End Weighthed Simple Regression============
 
 *=====Complex Survey Design Simple Regression============

* Simple linear regression this time using all the svyset factors
 *How do I make `1' a command that I can run?
 *display "The 1st argument you typed is: `svysetcommand'"
quietly `svysetcommand'
local rowfmt ="&-"
forval i = 1/`n_cc' { 
	local     secondindex = `i' + `n_cc'
	*summarize OGStataUVW`i'
	*summarize OGStataUVW`secondindex'
	quietly svy: regress OGStataUVW`i'  OGStataUVW`secondindex'
	matrix Results[`i',9]= `e(p)'
	/*
	quietly svy: regress OGStataUVW`secondindex' OGStataUVW`i'  
	matrix Results[`i',?]= `e(p)'
	*/
	if(`i'==`n_cc') {
		            local rowfmt="`rowfmt'" +"-"
                }
   else{
	     local rowfmt="`rowfmt'" +"&"
        }

} 
*display "`rowfmt'"
 *=====Complex Survey Design Simple Regression============
 *Present results from a matrix
 *Selecting only results that we want to present
 mat Results =Results[1...,1..9]
 *Using complex column names https://www.statalist.org/forums/forum/general-stata-discussion/general/1559138-%F0%9F%98%95-matrix-colnames-doesn-t-respect-the-space-character-in-the-names-of-the-columns
  matrix colnames Results = "CC Index" "CC Magnitude" "Wilk's Lambda" "Wilk's Lambda (n)" "Roy's Greatest Root" "Pillai's Trace" "Hotelling-Lawley Trace" "Weighted Regression" "CSD Regression"
  matlist Results,title("Complete List of Canonical Correlation p-values") /*border(rows)*//*rowtitle(rows)*/ /*format(%5.4f)*/ names(columns) cspec(o4& %12.0f | %9.4f & %12.4f & %12.4f & %12.4f & %12.4f & %12.4f & %12.4f & %12.4f  o2&) rspec("`rowfmt'")
 *Present results from data
 /*
   clear
   capture svmat double Results
   table
 */
 * restore user's data
 restore 
end 
/*====================end program==========================*/

/*============================================================================*/
/*=======================Start MATA Function===================================*/
version 17.0
mata:
void calcpval()
{
    myX=st_matrix("X")
	myY=st_matrix("Y")
	mydiag_W=st_matrix("diag_W")
    myW=st_matrix("W")
    myccorr=st_matrix("e(ccorr)")
	myResults=st_matrix("Results")
	myrawcoef_var1=st_matrix("e(rawcoef_var1)")
	myrawcoef_var2=st_matrix("e(rawcoef_var2)")

	myOgX=st_matrix("OgX")
	myOgY=st_matrix("OgY")
 	
	U=myOgX*myrawcoef_var1
	V=myOgY*myrawcoef_var2
	/*Now U and V are not complex matrices*/
	/*Sending U and V from mata to stata, in stata they will be called OGStataU and OGStataV*/
	st_matrix("OGStataU", U)
	st_matrix("OGStataV", V)
	UV = U,V
	UVW=UV,myW
	st_matrix("OGStataUVW", UVW)
	/*Using Eq. (8) and (9) to find the correlations and variances of the canonical matrices*/
	/*This calculations do not lead to correct canonical correlations or variances of the canonical variates as it did when using all my process above*/
	weighted_varU= 97* J(1, cols(U), 1)
	weighted_varV= 98*J(1, cols(V), 1)
	weighted_rho=  99*J(1, cols(U), 1)
	for (i=1; i<=cols(U); i++) {
		meanU=mean(U[,i])
		meanV=mean(V[,i])
		for (j=1; j<=rows(U); j++) {
			U[j,i]=U[j,i]-meanU
			V[j,i]=V[j,i]-meanV
		}/*end for j*/
        /*Now the means should be zero*/
		meanU=mean(U[,i])
		meanV=mean(V[,i])
		weighted_varU[1,i]=(U'[i,])*mydiag_W*(U[,i])/trace(mydiag_W)
		weighted_varV[1,i]=(V'[i,])*mydiag_W*(V[,i])/trace(mydiag_W)
		weighted_rho[1,i]=((U'[i,])*mydiag_W*(V[,i])/trace(mydiag_W))/sqrt(weighted_varU[1,i]*weighted_varV[1,i])
	    myResults[i,1]=i
	    myResults[i,2]=weighted_rho[1,i]
	}/*end for loop*/
	/*Notice that in this case the variances are (approximately) equal to 1 because with out normalization of the eigen values the constrain used during the calculation of the caonical correlations holds*/
    /*Now the means are zero and variances are 1 because they have been standarized*/
	/*
	mystring = "=======Weighted Mean of U========="
    display(mystring)
	meanU
	mystring = "=======Weighted Mean of U======="
    display(mystring)
	meanU
	
	mystring = "=======Weighted Variances of U======="
    display(mystring)
	weighted_varU
	mystring = "=======Weighted Variances of V======="
	display(mystring)
	weighted_varV
    */
/*The final correlation is equal to what we got form canon*/
	mystring = "=======Weighted Canonical Correlations======="
	display(mystring)
	weighted_rho
/*Finding sum of requencies*/	
 sumFREQ=sum(trunc(myW))
 /*sumFREQ*//*It works*/
/*=================== Start Wilks' Lambda================*/	
/*Gittins calls it Likelihood ratio test, pages 60-61*/
	/*trace(mydiag_W)*/
    Lambda= J(1,cols(weighted_rho),1)
	df= J(1,cols(weighted_rho),1)
	ChiSq_Wilks_Lambda=J(1,cols(weighted_rho),1)
	ChiSq_Wilks_LambdaFREQ=J(1,cols(weighted_rho),1)
	p_val_Wilks_Lambda=J(1,cols(weighted_rho),1)
	p_val_Wilks_LambdaFREQ=J(1,cols(weighted_rho),1)
	
	for (i=1; i<=cols(weighted_rho); i++) {
	   Lambda[1,i]=(1-(weighted_rho[1,i]^2))
       for (j=(i+1); j<=cols(weighted_rho); j++) {
			Lambda[1,i]=Lambda[1,i]*(1-(weighted_rho[1,j]^2)) /*Eq. (3.19) of Gittins*/
		}/*end for j*/
	   df[1,i]=(cols(myX)+1-i)*(cols(myY)+1-i)
	   ChiSq_Wilks_Lambda[1,i]=(((rows(myX)-1)-.5*(cols(myX)+cols(myY)+1))*ln(Lambda[1,i]))*-1 /*Eq. 3.20 of Gittins*/
	   ChiSq_Wilks_LambdaFREQ[1,i]=(((sumFREQ-1)-.5*(cols(myX)+cols(myY)+1))*ln(Lambda[1,i]))*-1 /*Eq. 3.20 of Gittins*/
	   p_val_Wilks_Lambda[1,i]=chi2tail(df[1,i], ChiSq_Wilks_Lambda[1,i])
	   p_val_Wilks_LambdaFREQ[1,i]=chi2tail(df[1,i], ChiSq_Wilks_LambdaFREQ[1,i])
	   myResults[i,3]=p_val_Wilks_Lambda[1,i]
	   myResults[i,4]=p_val_Wilks_LambdaFREQ[1,i]
	}/*end for loop*/
	mystring = "=======Wilks' Lambda Statistic======="
	display(mystring)
	Lambda
	/*
	mystring = "=======Wilks' Lambda Degrees of Freedom======="
	display(mystring)
	df
	*/
	mystring = "=======Wilks' Lambda Chi-Square======="
	display(mystring)
    ChiSq_Wilks_Lambda
	mystring = "=======Wilks' Lambda Chi-Square(FREQ)======="
	display(mystring)
    ChiSq_Wilks_LambdaFREQ
	
	mystring = "=======Wilks' Lambda p-value (From Chi-Square instead of F distribution)======="
	display(mystring)
	p_val_Wilks_Lambda
	mystring = "=======Wilks' Lambda p-value (FREQ) (From Chi-Square instead of F distribution)======="
	display(mystring)
	p_val_Wilks_LambdaFREQ
	/*.0415381409 vs 0.044,  .9582477164 vs 0.9592,  .8388719508  vs 0.839*/
/*These p-values are correct!!!!!*/		
     /*Working for the FREQ Option??? who knows*/
	/*p_val_Wilks_LambdaFREQ*/
	/*  0 vs <.0001,   0  vs <.0001,  0  vs 0.0021  |*/
/*=================== End Wilks' Lambda================*/		
	
/*===================Start Roy's Greatest Root================*/	
/*Gittins calls it Union-intersection test, page 61*/	
	/*Approximating Roy's Greatest Root using Eq. (12) of APPROXIMATE NULL DISTRIBUTION OF THE LARGEST ROOT IN MULTIVARIATE ANALYSIS1 BY IAIN M. JOHNSTONE*/
	/*Pθ(s,m,n)>θα)>PF(θα)=1−Fν1,ν2(ν2θα/ν1(1−θα) where ν1 =s+2m+1andν2 =s+2n+1 denote the hypothesis and error degrees of freedom respectively.
	The definitions of s=min(p,q), m=(|p-q|-1)/2 and n=(N-p-q-2)/2 is given in Gittins right above Eq. (3.17)
	*/
	p=cols(myX)
	q=cols(myY)
	pq=J(1,2,1)
	pq[1,1]=cols(myX)
	pq[1,2]=cols(myY)
	s=min(pq)
	m=(abs(p-q)-1)/2 
	n=(rows(myX)-p-q-2)/2
	v1 =s + (2*m) + 1
	v2 =s + (2*n) + 1
	largest_root=weighted_rho[1,1]^2
	Roys_Greatest_Root_stat = J(1,cols(weighted_rho),1)
	Roys_Greatest_Root_stat[1,1]=(v2*largest_root)/(v1*(1.0 - largest_root))
	Roys_Greatest_Root_p_value = J(1,cols(weighted_rho),1)
	Roys_Greatest_Root_p_value[1,1]=Ftail(v1, v2, Roys_Greatest_Root_stat[1,1])
	myResults[1,5]=Roys_Greatest_Root_p_value[1,1]
	/*10.48376021 vs 10.42 in SAS
	.0004667007 vs 0.0005 from SAS*/
	/*For the roots k where k>=2 we use the formulas at the bottom of page 61 of Gittins*/
	for (j=2; j<=cols(weighted_rho); j++) {
		s=min(pq)-j-1
	    m=(abs(p-q-j)-1)/2 
	    n=(rows(myX)-p-q-j)/2
	    v1 =s + (2*m) + 1
	    v2 =s + (2*n) + 1
	    largest_root=weighted_rho[1,j]^2
	    Roys_Greatest_Root_stat[1,j]=(v2*largest_root)/(v1*(1.0 - largest_root))
	    Roys_Greatest_Root_p_value[1,j]=Ftail(v1, v2, Roys_Greatest_Root_stat[1,j])
		myResults[j,5]=Roys_Greatest_Root_p_value[1,j]
	}/*end for j*/
	
	mystring = "=======IAIN M. JOHNSTONE Approximation of Roy's Greatest Root Statistic/F-value======="
	display(mystring)
	Roys_Greatest_Root_stat
	
	mystring = "=======IAIN M. JOHNSTONE Approximation of Roy's Greatest p-value======="
	display(mystring)
	Roys_Greatest_Root_p_value
/*Very bottom of page 61 "Like the corresponding likelihood ratio tests, these tests pf the r_k^2 for K>=2 are very conservative (Mardia, Kent & Bibby, 1979, p. 147)."*/		
/*These p-values are correct!!!!! At least the first one which is the only one provided by SAS*/	
/*===================End Roy's Greatest Root================*/	

/*===================Start Pillai's Trace=====================*/
	/*Pillai's Trace*//*This is T fin Eq(7) of All formulas p-values.pdf Comparison of Test Statistics of Nonnormal and Unbalanced Samples for Multivariate Analysis of Variance in terms of Type-I Error Rates Can Ates¸
	But for some reason T  in the Can Ates paper is kind of the inverse of V/Lambda in Gittins and SAS Cancorr manual the lambda in one is lambda/(1+lambda) or lambda/(1-lambda) in the other
	hence the formulas for Pillai's Trace and Roy's Largest Root here in this program are simple but in equations (7) and (8) of Can Ate are "complicated" (divisions) but here the
	formula for Hotellings Trace is complicated (a division) but eq. (9) in Can Ates is simple.
	All comes from the different definitions in Eq. (3.19) in Gitting vs. Eq. (5) in Can Ates  which includes and ^-1
	Althought Can Ates deal with MANOVA, not CCA, so he has "g is the number of groups, p is the number of variables in each group, N is the number of observations... and s=min(g − 1, p)"
	So for Pillai's his degrees of freedom are s(2m + s + 1) and (2n + s + 1), while we have s(2m + s + 1) and s*(2n + s + 1) from right above Roy's Maximum Root in the SAS webstie for multivariae tests
	*/
	V=0
	for (i=1; i<=cols(weighted_rho); i++) {
	    V=V+(weighted_rho[1,i]^2)/*Notice that in SAS rho_i =sqrt(Xi_i) "is the ith canonical correlation." hence Xi_i is the square of the ith canonical corrrelation*/
	}/*end for loop*/
	/*From https://documentation.sas.com/doc/en/pgmsascdc/9.4_3.3/statug/statug_introreg_sect038.htm#statug_introreg001883*/
	s=min(pq)
	m=(abs(p-q)-1)/2 
	n=(rows(myX)-p-q-2)/2 /*In the SAS help n=(v-p-1)/2 where v= error degrees of freedom (this formula for v it is right above Wilks' Lambda )*/
	df1=s*((2*m)+s+1)
	df2=s*((2*n)+s+1)
	/*The degrees of freedom are the same as in SAS 9 & 48*/
   	Pillais_Trace_stat = J(1,cols(weighted_rho),1)
	Pillais_Trace_stat[1,1]=(((2*n)+s+1)/((2*m)+s+1))*V/(s-V)
	Pillais_Trace_p_value = J(1,cols(weighted_rho),1)
	Pillais_Trace_p_value[1,1]=Ftail(df1, df2, Pillais_Trace_stat[1,1])
	myResults[1,6]=Pillais_Trace_p_value[1,1]
	mystring = "=======Pillais_Trace  Degrees of Freedom 1======="
	display(mystring)
	df1
	mystring = "=======Pillais_Trace Degrees of Freedom 2======="
	display(mystring)
	df2
	
	for (j=2; j<=cols(weighted_rho); j++) {
		V=0
		for (i=j; i<=cols(weighted_rho); i++) {
	        V=V+(weighted_rho[1,i]^2) /*Notice that in SAS rho_i =sqrt(Xi_i) "is the ith canonical correlation." hence Xi_i is the square of the ith canonical corrrelation*/
	     }/*end for i loop*/
		 /*The degrees of freedom are the same as in SAS 9 & 48,  in Eq (8) of of All formulas p-values.pdf */
		Pillais_Trace_stat[1,j]=(((2*n)+s+1)/((2*m)+s+1))*V/(s-V)
	    Pillais_Trace_p_value[1,j]=Ftail(df1, df2, Pillais_Trace_stat[1,j])
	    myResults[j,6]=Pillais_Trace_p_value[1,j]
	}/*end for j*/
		
	mystring = "=======SAS' Pillai's Trace Statistic/F-value======="
	display(mystring)
    Pillais_Trace_stat
	
	mystring = "=======SAS' Pillai's Trace p-value======="
	display(mystring)
    Pillais_Trace_p_value
	/* 1.633808758 vs 1.63 in SAS,  .1324569442 vs. 0.1341 in SAS */
    /*These p-values are correct!!!!! At least the first one which is the only one provided by SAS*/
/*===================End Pillai's Trace=====================*/

/*===================Start Hotelling-Lawley Trace=====================*/	
	/*Notice that we recycle the values of s, m and n from Pillais Trace*/
	/*These formulas come from https://documentation.sas.com/doc/en/pgmsascdc/9.4_3.3/statug/statug_introreg_sect038.htm#statug_introreg001883 which I also have in PDF*/
	/*Hence we get the SAS results which are different from the Stata Results Den DF in SAS =/= df2 in Stata, e.g. for the Middle-Aged Men in a Health Fitness Club with weights
	I get 19.05263158 which is similar to the 19.053 from SAS by very different from the 38 of Stata, the F-Values are a lot closer: SAS=2.97, Stata= 2.8078, me =2.991206843*/
	
	Hotelling_Lawley_Trace_stat = J(1,cols(weighted_rho),1)
	Hotelling_Lawley_Trace_p_value = J(1,cols(weighted_rho),1)
	
	U=0
	for (i=1; i<=cols(weighted_rho); i++) {
	    U=U+(weighted_rho[1,i]^2)/(1-(weighted_rho[1,i]^2)) /*Notice that in SAS rho_i =sqrt(Xi_i) "is the ith canonical correlation." hence Xi_i is the square of the ith canonical corrrelation*/
	}/*end for loop*/
	if (n>0){
		display("n>0")
		b=(p+2*n)*(q+2*n)/(2*(2*n+1)*(n-1))
		c=(2+(p*q+2)/(b-1))/(2*n)
		mydf1=p*q
		mydf1
		/*mydf1 is correct and so is the stat, just mydf is incorrect should much larger*/
		mydf2=4+((p*q+2)/(b-1)) 
		mydf2
		Hotelling_Lawley_Trace_stat[1,1]=(U/c)*((4+(p*q+2)/(b-1))/(p*q))
		Hotelling_Lawley_Trace_p_value[1,1]= Ftail(mydf1, mydf2, Hotelling_Lawley_Trace_stat[1,1])	
	    myResults[1,7]=Hotelling_Lawley_Trace_p_value[1,1]
		/*myResults[1,7]=1-chi2tail(mydf1, Hotelling_Lawley_Trace_stat[1,1])*/
	}else{
		display("n<=0")
		mydf1=s*(2*m+s+1)
		mydf2=2*(s*n+1)
		Hotelling_Lawley_Trace_stat[1,1]= (2*(s*n + 1)*U)/((s^2)*(2*m+s+1))
		Hotelling_Lawley_Trace_p_value[1,1]=	Ftail(mydf1, mydf2, Hotelling_Lawley_Trace_stat[1,1])
		myResults[1,7]=Hotelling_Lawley_Trace_p_value[1,1]
	}
	mystring = "=======Hotelling-Lawley Trace  Degrees of Freedom 1======="
	display(mystring)
	mydf1
	mystring = "=======Hotelling-Lawley Trace  Degrees of Freedom 2 Same as SAS different from Stata!!!======="
	display(mystring)
	mydf2
	
	for (j=2; j<=cols(weighted_rho); j++) {
		U=0
		for (i=j; i<=cols(weighted_rho); i++) {
	           U=U+(weighted_rho[1,i]^2)/(1-(weighted_rho[1,i]^2)) /*Notice that in SAS rho_i =sqrt(Xi_i) "is the ith canonical correlation." hence Xi_i is the square of the ith canonical corrrelation*/
        }/*end for i loop*/
		if (n>0){
		b=(p+2*n)*(q+2*n)/(2*(2*n+1)*(n-1))
		c=(2+(p*q+2)/(b-1))/(2*n)
		mydf1=p*q
		mydf2=4+(p*q+2)/(b-1)
		Hotelling_Lawley_Trace_stat[1,j]=(U/c)*((4+(p*q+2)/(b-1))/(p*q))
		Hotelling_Lawley_Trace_p_value[1,j]=	Ftail(mydf1, mydf2, Hotelling_Lawley_Trace_stat[1,j])
		myResults[j,7]=Hotelling_Lawley_Trace_p_value[1,j]
	    }else{
		mydf1=s*(2*m+s+1)
		mydf2=2*(s*n+1)
		Hotelling_Lawley_Trace_stat[1,j]= (2*(s*n + 1)*U)/((s^2)*(2*m+s+1))
		Hotelling_Lawley_Trace_p_value[1,j]=	Ftail(mydf1, mydf2, Hotelling_Lawley_Trace_stat[1,j])
		myResults[j,7]=Hotelling_Lawley_Trace_p_value[1,j]
	   }
	}/*end for j*/
	
	/*
	mystring = "=======Hotelling-Lawley Trace  Degrees of Freedom 1======="
	display(mystring)
	mydf1
		mystring = "=======Hotelling-Lawley Trace  Degrees of Freedom 2======="
	display(mystring)
	mydf2
	*/
	mystring = "=======SAS' Hotelling-Lawley Trace  Statistic/F-value======="
	display(mystring)
	Hotelling_Lawley_Trace_stat
	mystring = "=======SAS' Hotelling-Lawley Trace  p-value======="
	display(mystring)
	Hotelling_Lawley_Trace_p_value
	/*These p-values are correct!!!!! At least the first one which is the only one provided by canon*/	
/*==================End Hotelling-Lawley Trace=====================*/
/*Returning the variables that are needed to calculate the simple linear regressions with complex survey design: not needed because we sent them from mata to stata*/
st_matrix("Results", myResults)
}	
end
/*========== End MATA function ==========*/