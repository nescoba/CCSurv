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
   version 17
   tokenize `e(cmdline)', parse("()=]")
   /*display "1=|`1'|, 2=|`2'|, 3=|`3'|, 4=|`4'|, 5=|`5'|, 6=|`6'|, 7=|`7'|, 8=|`8'|, 9=|`9'|,10=|`10'|"*/
   /*tokens 3 and 6 also work with NYTS data*/
   /*
   di "`3'"
   di "`6'"
   di "`10'"
   */
   local X `3'
   local Y `6'
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

   clear
   capture svmat double OgX
   foreach var of varlist _all {
   	egen aux_`var'= std(`var')
    drop `var'
   }
	*Not standarizing we don't get the correct canonical correlations
   mkmat  _all, matrix(X)

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
  mata: calcpval(5)
  capture svmat OGStataUVW
  *ereturn list
  *e(n_cc) =  3 has the number of canonical correlations; it is a scalar
 
 local weightindex = (2*e(n_cc)) + 1
 local n_cc  = e(n_cc)
 
 svyset  [pweight= OGStataUVW`weightindex'],

 *Notice coefficients of the simple linear regression are equal to the canonical correlations because the variances of the canonical variates are equal to 1
forval i = 1/`=e(n_cc)' { 
	local     secondindex = `i' + `n_cc'
	*summarize OGStataUVW`i'
	*summarize OGStataUVW`secondindex'
	svy: regress OGStataUVW`i'  OGStataUVW`secondindex'
	svy: regress OGStataUVW`secondindex' OGStataUVW`i'  
} 

end 
/*====================end program==========================*/

version 17.0
mata:
void calcpval(real scalar N)
{
    myX=st_matrix("X")
	myY=st_matrix("Y")
	mydiag_W=st_matrix("diag_W")
    myW=st_matrix("W")
    	
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
/*The final correlation is equal to what we got form svy: regress*/
	mystring = "=======Weighted Canonical Correlations======="
	display(mystring)
	weighted_rho
	
/*=================== Start Wilks' Lambda================*/	
	/*trace(mydiag_W)*/
    Lambda= J(1,cols(weighted_rho),1)
	df= J(1,cols(weighted_rho),1)
	ChiSq_Wilks_Lambda=J(1,cols(weighted_rho),1)
	p_val_Wilks_Lambda=J(1,cols(weighted_rho),1)
	/*p_val_Wilks_LambdaFREQ=J(1,cols(weighted_rho),1)*//*Remeber: I haven't been able to replicate the results of SAS using FREQ instead of WEIGHT*/
	
	for (i=1; i<=cols(weighted_rho); i++) {
	   Lambda[1,i]=(1-(weighted_rho[1,i]^2))
       for (j=(i+1); j<=cols(weighted_rho); j++) {
			Lambda[1,i]=Lambda[1,i]*(1-(weighted_rho[1,j]^2))
		}/*end for j*/
	   df[1,i]=(cols(myX)+1-i)*(cols(myY)+1-i)
	   ChiSq_Wilks_Lambda[1,i]=(((rows(myX)-1)-.5*(cols(myX)+cols(myY)+1))*ln(Lambda[1,i]))*-1
	   p_val_Wilks_Lambda[1,i]=chi2tail(df[1,i], ChiSq_Wilks_Lambda[1,i])
/*===== NOT WORKING the DF in SAS with FREQ statement, it seems that trace(diagW) does not do the work=========================================*/
	  /* p_val_Wilks_LambdaFREQ[1,i]=chi2tail(trace(mydiag_W), ChiSq_Wilks_Lambda[1,i])  */
	}/*end for loop*/
	mystring = "=======Wilks' Lambda Statistic======="
	display(mystring)
	Lambda
	/*
	mystring = "=======Wilks' Lambda Degrees of Freedom======="
	display(mystring)
	df
	mystring = "=======Wilks' Lambda Chi-Square======="
	display(mystring)
    ChiSq_Wilks_Lambda
	*/
	mystring = "=======Wilks' Lambda p-value (From Chi-Square instead of F distribution)======="
	display(mystring)
	p_val_Wilks_Lambda
	/*.0415381409 vs 0.044,  .9582477164 vs 0.9592,  .8388719508  vs 0.839*/
/*These p-values are correct!!!!!*/		
     /*Working for the FREQ Option??? who knows*/
	/*p_val_Wilks_LambdaFREQ*/
	/*  0 vs <.0001,   0  vs <.0001,  0  vs 0.0021  |*/
/*=================== End Wilks' Lambda================*/		
	
/*===================Start Roy's Greatest Root================*/		
	/*Approximating Roy's Greatest Root using Eq. (12) of APPROXIMATE NULL DISTRIBUTION OF THE LARGEST ROOT IN MULTIVARIATE ANALYSIS1 BY IAIN M. JOHNSTONE*/
	/*Pθ(s,m,n)>θα)>PF(θα)=1−Fν1,ν2(ν2θα/ν1(1−θα) where ν1 =s+2m+1andν2 =s+2n+1denotethehypothesisanderrordegrees of freedom respectively.
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
	/*Pillai's Trace*/
	V=0
	for (i=1; i<=cols(weighted_rho); i++) {
	    V=V+(weighted_rho[1,i]^2) /*Notice that in SAS rho_i =sqrt(Xi_i) "is the ith canonical correlation." hence Xi_i is the square of the ith canonical corrrelation*/
	}/*end for loop*/
	/*From https://documentation.sas.com/doc/en/pgmsascdc/9.4_3.3/statug/statug_introreg_sect038.htm#statug_introreg001883*/
	s=min(pq)
	m=(abs(p-q)-1)/2 
	n=(rows(myX)-p-q-2)/2 /*In the SAS help n=(v-p-1)/2 where v= error degrees of freedom*/
	df1=s*((2*m)+s+1)
	df2=s*((2*n)+s+1)
	/*The degrees of freedom are the same as in SAS 9 & 48*/
	Pillais_Trace_stat=(((2*n)+s+1)/((2*m)+s+1))*V/(s-V)
	Pillais_Trace_p_value=Ftail(df1, df2, Pillais_Trace_stat)
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
	U=0
	for (i=1; i<=cols(weighted_rho); i++) {
	    U=U+(weighted_rho[1,i]^2)/(1-(weighted_rho[1,i]^2)) /*Notice that in SAS rho_i =sqrt(Xi_i) "is the ith canonical correlation." hence Xi_i is the square of the ith canonical corrrelation*/
	}/*end for loop*/
	if (n>0){
		b=(p+2*n)*(q+2*n)/(2*(2*n+1)*(n-1))
		c=(2+(p*q+2)/(b-1))/(2*n)
		mydf1=p*q
		mydf2=4+(p*q+2)/(b-1)
		Hotelling_Lawley_Trace_stat=(U/c)*((4+(p*q+2)/(b-1))/(p*q))
		Hotelling_Lawley_Trace_p_value=	Ftail(mydf1, mydf2, Hotelling_Lawley_Trace_stat)
		
	}else{
		mydf1=s*(2*m+s+1)
		mydf2=2*(s*n+1)
		Hotelling_Lawley_Trace_stat= (2*(s*n + 1)*U)/((s^2)*(2*m+s+1))
		Hotelling_Lawley_Trace_p_value=	Ftail(mydf1, mydf2, Hotelling_Lawley_Trace_stat)
	}
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
	/*2.991206843 vs. 2.97  in SAS, 9 vs. 9 in SAS,  19.05263158 vs 19.053 in SAS , .0211499602 vs 0.0218 in SAS */
	/*These p-values are correct!!!!! At least the first one which is the only one provided by SAS*/	
/*==================End Hotelling-Lawley Trace=====================*/
/*Returning the variables that are needed to calculate the simple linear regressions with complex survey design: not needed because we sent them from mata to stata*/
}	
end
/*========== End MATA function ==========*/