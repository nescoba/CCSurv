
matrix diag_W = I(rowsof(X))
local n_f =rowsof(X)
forvalues i = 1/`n_f' {
		 matrix diag_W [`i',`i']=  W[`i',1]
}

matrix factorYX = Y'*diag_W*X
matrix factorInvXX = inv(X'*diag_W*X)
matrix factorXY = X'*diag_W*Y
matrix factorInvYY = inv(Y'*diag_W*Y)

matrix A=factorYX*factorInvXX*factorXY*factorInvYY
mata
    myX=st_matrix("X")
	myX
	myY=st_matrix("Y")
	myY
	mydiag_W=st_matrix("diag_W")
	mydiag_W
	myfactorYX=st_matrix("factorYX")
	myfactorYX
	myfactorInvXX=st_matrix("factorInvXX")
	myfactorInvXX
	myA=st_matrix("A")
	myA
	bT = .
	lambdaSq = .
	
	lefteigensystem(myA, bT, lambdaSq)
	bT
	
	lambdaSq
	lambda=lambdaSq
	for (i=1; i<=length(lambda); i++) {
      lambda[i,1]=lambdaSq[i,1]^.5 /*https://www.ssc.wisc.edu/sscc/pubs/4-26.htm*/
    }/*end for loop*/
	lambda
	aT=bT*myfactorYX*(myfactorInvXX)
	for (i=1; i<=length(lambda); i++) {
      aT=(1/lambda[i,1])*aT 
	}/*end for loop*/
	aT	
	
	UT=aT*myX'
	UT
	U=UT'
	VT=bT*myY'
	VT
	V=VT'
	
	weighted_varU=lambda
	weighted_varV=lambda
	weighted_rho=lambda
	for (i=1; i<=length(lambda); i++) {
      weighted_varU[i,1]=UT[i,]*mydiag_W*U[,i]/trace(mydiag_W)
	  weighted_varV[i,1]=VT[i,]*mydiag_W*V[,i]/trace(mydiag_W)
	  weighted_rho[i,1]=(UT[i,]*mydiag_W*V[,i]/trace(mydiag_W))/sqrt(weighted_varU[i,1]*weighted_varV[i,1])
	}/*end for loop*/
	weighted_varU
	weighted_varV
/*The final correlation should be equalt o what we got form the eigenvalues*/
	weighted_rho
	/*
	Canonical Correlation
0.813267
0.193757
0.051574

 weighted_rho
                1
    +--------------+
  1 |  .814133001  |
  2 |  .195005619  |
  3 |  .051612239  |
    +--------------+

	 lambda
                1
    +--------------+
  1 |  .814133001  |
  2 |  .195005619  |
  3 |  .051612239  |
    +--------------+

	Canonical correlations from SAS:  0.813267  0.193757 0.051574
	Canonical correlations:
  0.8133  0.1938  0.0516
*/
/* The matrice have Imaginary parts (probably ffrom the calculation eigenvalues) but are all zeros
mata describe
myIm=Im(U)
myIm
myRe=Re(U)
myRe
*/
	st_matrix("StataU", Re(U))
	st_matrix("StataV", Re(V))
end







/*matrix list StataU*/
svmat StataU
svmat StataV
/*For correlate "aweights and fweights are allowed; see [U] 11.1.6 weight" see rcorrelate.pdf*/
correlate StataU1 StataU2 StataU3 StataV1 StataV2 StataV3 [w=weight] 
svyset  [pweight=weight],
svy: regress StataU1  StataV1   /* r= (80.74207  * SQRT(.563960908) /SQRT(5567.36064)=0.81264282  vs  .814133001 from my work & 0.8133 from canon  Prob > F  = 0.0013*/
svy: regress StataU2  StataV2   /* r= (4.587849  * SQRT(0.696245958)/SQRT(394.336132)=0.192777868 vs  .195005619 from my work & 0.1938 from canon  Prob > F  = 0.2587*/
svy: regress StataU3  StataV3   /* r= (0.3248543 * SQRT(0.267148451)/SQRT(10.5990683)=0.051574034 vs  .051612239 from my work & .0516 from canon   Prob > F  = 0.8096*/

/*We don't get the same results if we just shift the variables!!!*/
svy: regress StataV1  StataU1   /* r= (.0081916  * SQRT(5567.36064)/ SQRT(.563960908)=0.813896    vs  .814133001 from my work & 0.8133 from canon  Prob > F  =  0.0000*/
svy: regress StataV2  StataU2   /* r= (.0081835  * SQRT(394.336132)/SQRT(0.696245958)=0.194756082 vs  .195005619 from my work & 0.1938 from canon  Prob > F  = 0.4117*/
svy: regress StataV3  StataU3   /* r= (.0081882 * SQRT(10.5990683)/SQRT(0.267148451)/=0.051575798 vs  .051612239 from my work & .0516 from canon   Prob > F  = 0.8070*/


canon (weight waist pulse) (chins situps jumps) [weight=weight]


 matrix list e(rawcoef_var1)
 matrix list e(rawcoef_var2)
 mata
    myrawcoef_var1=st_matrix("e(rawcoef_var1)")
	myrawcoef_var1
	myrawcoef_var2=st_matrix("e(rawcoef_var2)")
	myrawcoef_var2
/*The matrices X and Y are calculated in Stata so they should still be available*/	
	myOgX=st_matrix("OgX")
	myOgX
	myOgY=st_matrix("OgY")
	myOgY
	/*
	U=myX*myrawcoef_var1
	V=myY*myrawcoef_var2
	*/
	U=myOgX*myrawcoef_var1
	U'
	V=myOgY*myrawcoef_var2
	V'
	/*Now U and V are not complex matrices*/
	st_matrix("OGStataU", U)
	st_matrix("OGStataV", V)
	
	/*Using Eq. (8) and (9) to find the correlations and variances of the canonical matrices*/
	/*This calculations do not lead to correct canonical correlations or variances of the canonical variates as it did when using all my process above*/
	weighted_varU= J(1, cols(U), 1)
	weighted_varV= J(1, cols(V), 1)
	weighted_rho=J(1, cols(U), 1)
	for (i=1; i<=cols(U); i++) {
		meanU=mean(U[,i])
		meanU
		meanV=mean(V[,i])
		meanV
		for (j=1; j<=rows(U); j++) {
			U[j,i]=U[j,i]-meanU
			V[j,i]=V[j,i]-meanV
		}/*end for j*/
		meanU=mean(U[,i])
		meanU
		weighted_varU[1,i]=(U'[i,])*mydiag_W*(U[,i])/trace(mydiag_W)
	    weighted_varV[1,i]=(V'[i,])*mydiag_W*(V[,i])/trace(mydiag_W)
	    weighted_rho[1,i]=((U'[i,])*mydiag_W*(V[,i])/trace(mydiag_W))/sqrt(weighted_varU[1,i]*weighted_varV[1,i])
	   
	}/*end for loop*/
	/*Notice that in this case the variances are (approximately) equal to 1 because with out normalization of the eigen values the constrain used during the calculation of the caonical correlations holds*/
	weighted_varU
	weighted_varV
/*The final correlation is equal to what we got form svy: regress*/
	weighted_rho
end 

/*matrix list StataU*/
svmat OGStataU
svmat OGStataV
svyset  [pweight=weight],
summarize OGStataU1  OGStataV1 [aweight=weight]











/*=====================GETTING SAME P-VALUES AS WITH MY CANONICAL VARIATES & same correlations but here the variance of all canonical variates is =1 due to the constrained in the calculations of the canonical coefficients====================================*/

svy: regress OGStataU1  OGStataV1   /* r= 0.8132676  vs  .814133001 from my work & 0.8133 from canon  Prob > F  = 0.0013 vs  0.001 from svy: regress */
svy: regress OGStataU2  OGStataV2   /* r= .1937568   vs  .195005619 from my work & 0.1938 from canon  Prob > F  = 0.2587 vs  0.258 from svy: regress*/
svy: regress OGStataU3  OGStataV3   /* r=  .051574  vs  .051612239 from my work & .0516 from canon   Prob > F  = 0.8096 vs   0.810 from svy: regress*/

/*We don't get the same results if we just shift the variables!!!*/
svy: regress OGStataV1  OGStataU1   /* r=  .8132674  vs  .814133001 from my work & 0.8133 from canon  Prob > F  =  0.0000 vs.  0.000 from svy: regress */
svy: regress OGStataV2  OGStataU2   /* r=  .1937568  vs  .195005619 from my work & 0.1938 from canon  Prob > F  = 0.4117 vs  0.412  from svy: regress*/
svy: regress OGStataV3  OGStataU3   /* r=   .0515741 vs  .051612239 from my work & .0516 from canon   Prob > F  = 0.8070 vs. 0.807 from svy: regress*/


 mata
    trace(mydiag_W)
    Lambda= J(1,cols(weighted_rho),1)
	df= J(1,cols(weighted_rho),1)
	ChiSq_Wilks_Lambda=J(1,cols(weighted_rho),1)
	p_val_Wilks_Lambda=J(1,cols(weighted_rho),1)
	/*p_val_Wilks_LambdaFREQ=J(1,cols(weighted_rho),1)*/
	
	for (i=1; i<=cols(weighted_rho); i++) {
	   Lambda[1,i]=(1-(weighted_rho[1,i]^2))
       for (j=(i+1); j<=cols(weighted_rho); j++) {
			Lambda[1,i]=Lambda[1,i]*(1-(weighted_rho[1,j]^2))
		}/*end for j*/
	   df[1,i]=(cols(myX)+1-i)*(cols(myY)+1-i)
	   ChiSq_Wilks_Lambda[1,i]=(((rows(myX)-1)-.5*(cols(myX)+cols(myY)+1))*ln(Lambda[1,i]))*-1
	   p_val_Wilks_Lambda[1,i]=chi2tail(df[1,i], ChiSq_Wilks_Lambda[1,i])
/*===== NOT WORKING the DF in SAS variafrom CC to CC while I always use trace(diagW)=========================================*/
	  /* p_val_Wilks_LambdaFREQ[1,i]=chi2tail(trace(mydiag_W), ChiSq_Wilks_Lambda[1,i])  */
	}/*end for loop*/
	/*Notice that in this case the variances are (approximately) equal to 1 because with out normalization of the eigen values the constrain used during the calculation of the caonical correlations holds*/
/*The final correlation is equal to what we got form svy: regress*/
	Lambda
	df
    ChiSq_Wilks_Lambda
	p_val_Wilks_Lambda
	

	p=cols(myX)
	q=cols(myY)
	pq=J(1,2,1)
	pq[1,1]=cols(myX)
	pq[1,2]=cols(myY)
	s=min(pq)
	m=(abs(p-q)-1)/2 
	n=(rows(myX)-p-q-2)/2
	v1 =s + (2*m) + 1
	v2 =s+ (2*n) + 1
	largest_root=weighted_rho[1,1]^2
	Roys_Greatest_Root_stat = J(1,cols(weighted_rho),1)
	Roys_Greatest_Root_stat[1,1]=(ν2*largest_root)/(ν1*(1.0 - largest_root))
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
	    v2 =s+ (2*n) + 1
	    largest_root=weighted_rho[1,j]^2
	    Roys_Greatest_Root_stat[1,j]=(ν2*largest_root)/(ν1*(1.0 - largest_root))
	    Roys_Greatest_Root_p_value[1,j]=Ftail(v1, v2, Roys_Greatest_Root_stat[1,j])
	}/*end for j*/
	Roys_Greatest_Root_stat
	Roys_Greatest_Root_p_value
/*Very bottom of page 61 "Like the corresponding likelihood ratio tests, these tests pf the r_k^2 for K>=2 are very conservative (Mardia, Kent & Bibby, 1979, p. 147)."*/		
	
	/*Pillai's Trace*/
	V=0
	for (i=1; i<=cols(weighted_rho); i++) {
	    V=V+(weighted_rho[1,i]^2) /*Notice that in SAS rho_i =sqrt(Xi_i) "is the ith canonical correlation." hence Xi_i is the square of the ith canonical corrrelation*/
	}/*end for loop*/
	/*From https://documentation.sas.com/doc/en/pgmsascdc/9.4_3.3/statug/statug_introreg_sect038.htm#statug_introreg001883*/
	V
	s=min(pq)
	m=(abs(p-q)-1)/2 
	n=(rows(myX)-p-q-2)/2 /*In the SAS help n=(v-p-1)/2 where v= error degrees of freedom*/
	df1=s*((2*m)+s+1)
	df1
	df2=s*((2*n)+s+1)
	df2
	/*The degrees of freedom are the same as in SAS 9 & 48*/
	Pillais_Trace_stat=(((2*n)+s+1)/((2*m)+s+1))*V/(s-V)
	Pillais_Trace_p_value=Ftail(df1, df2, Pillais_Trace_stat)
    Pillais_Trace_stat
    Pillais_Trace_p_value
	/* 1.633808758 vs 1.63 in SAS,  .1324569442 vs. 0.1341 in SAS */
	/*Hotelling-Lawley Trace*/
	/*Notice that we recycle the values of s, m and n from Pillais Trace*/
	U=0
	for (i=1; i<=cols(weighted_rho); i++) {
	    U=U+(weighted_rho[1,i]^2)/(1-(weighted_rho[1,i]^2)) /*Notice that in SAS rho_i =sqrt(Xi_i) "is the ith canonical correlation." hence Xi_i is the square of the ith canonical corrrelation*/
	}/*end for loop*/
	U
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
	mydf1
	mydf2
	Hotelling_Lawley_Trace_stat
	Hotelling_Lawley_Trace_p_value
	/*2.991206843 vs. 2.97  in SAS, 9 vs. 9 in SAS,  19.05263158 vs 19.053 in SAS , .0211499602 vs 0.0218 in SAS */
end 

