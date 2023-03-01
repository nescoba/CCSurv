*1-30-2023: Using subset of variable sent by Erin on 1-30-3023
*1-31-2023: Using subset of variable sent by Erin on 1-31-3023
*2-1-2023: Keeping only the persons who answers more than 0 for "QN89: During the past 30 days, on how many days did you use any tobacco product(s)?
*2-2-2023: Sending Ering coefficients for interpretation of selected results
 *Home
 * cd "C:\Users\Raul Cruz-Cano\OneDrive - Indiana University\Documents\students\carlos\stata"
  *Office
  cd "C:\Users\raulcruz\OneDrive - Indiana University\Documents\students\carlos\R21\data\2021-nyts-dataset-and-codebook-microsoft-excel\"
  insheet using nyts2021.csv, clear
  cls 
  replace qn9=0 if qn6 ==2  
  replace qn38=0 if qn35==2  
  replace qn40=0 if qn38==0  
  replace qn53=0 if qn51==2  
  replace qn54=0 if qn53==0  
  replace qn64=0 if qn62==2  
  replace qn69=0 if qn67==2  
  replace qn74=0 if qn73==2  
  replace qn76=0 if qn75==2  
  replace qn78=0 if qn77==2  
  replace qn80=0 if qn79==2  
  replace qn82=0 if qn81==2  
  replace qn85=0 if qn84==2  
  replace qn88=0 if qn87==2  
  replace qn89=0 if qn9==0 & qn38==0 & qn53==0 & qn64==0 & qn69==0 & qn74==0 & qn76==0 & qn78==0 & qn80==0 & qn82==0 & qn85==0 & qn88==0 
  
  
*E-cigarette TOBACCO MARKETING (excluding social media) vs. PERCEIVED ADDICTION RELATIVE TO CIGARETTES
 keep qn4a qn4b qn5b qn9 qn38 qn40 qn53 qn54 qn64 qn69 qn74 qn76 qn78 qn80 qn82 qn85 qn88 qn89 finwgt psu2 v_stratum2 qn128 qn129 qn130 qn131 qn132 qn134 
 *Getting rid of missing values in the dependent variables
  drop if qn9==. |  qn38==. |  qn40==. |  qn53==. |  qn54==. |  qn64==. |  qn69==. |  qn74==. |  qn76==. |  qn78==. |  qn80==. |  qn82==. |  qn85==. |  qn88==. |  qn89==.  
 *Getting rid of missing values in the independent variables
  drop if  qn128==. |  qn129==. |  qn130==. |  qn131 ==. |   qn132==. |  qn134==.
 
 /*Keeping only the persons who answers more than 0 for "QN89: During the past 30 days, on how many days did you use any tobacco product(s)?*/
*keep if qn89 > 0
/*Slecting people that uses at least 2 products*/
/*Notice qn40 is not in the list in email Wednesday, February 1, 2023 9:10 AM*/
gen      qn9bin=0 if qn9==0 
replace  qn9bin=1 if qn9>0
gen      qn38bin=0 if qn38==0 
replace  qn38bin=1 if qn38>0
gen      qn53bin=0 if qn53==0 
replace  qn53bin=1 if qn53>0
gen      qn64bin=0 if qn64==0 
replace  qn64bin=1 if qn64>0
gen      qn69bin=0 if qn69==0 
replace  qn69bin=1 if qn69>0
gen      qn74bin=0 if qn74==0 
replace  qn74bin=1 if qn74>0
gen      qn76bin=0 if qn76==0 
replace  qn76bin=1 if qn76>0
gen      qn78bin=0 if qn78==0 
replace  qn78bin=1 if qn78>0
gen      qn80bin=0 if qn80==0 
replace  qn80bin=1 if qn80>0
gen      qn82bin=0 if qn82==0 
replace  qn82bin=1 if qn82>0
gen      qn85bin=0 if qn85==0 
replace  qn85bin=1 if qn85>0
gen      qn88bin=0 if qn88==0 
replace  qn88bin=1 if qn88>0

gen number_of_products= qn9bin + qn38bin + qn53bin + qn64bin + qn69bin + qn74bin + qn76bin + qn78bin + qn80bin + qn82bin + qn85bin + qn88bin
*keep if number_of_products >= 2
 keep if qn5b==1
 
canon (qn9 qn38 qn40 qn53 qn54 qn64 qn69 qn74 qn76 qn78 qn80 qn82 qn85 qn88 qn89) (qn128 qn129 qn130 qn131 qn132 qn134) [weight=finwgt]
 
 mkmat qn9 qn38 qn40 qn53 qn54 qn64 qn69 qn74 qn76 qn78 qn80 qn82 qn85 qn88 qn89, matrix(OgX) 
 mkmat qn128 qn129 qn130 qn131 qn132 qn134, matrix(OgY) 
 mkmat finwgt, matrix(W) 
/*Creating diagonal matrix with the survey Weights*/
matrix diag_W = I(rowsof(OgX))
local n_f =rowsof(OgX)
forvalues i = 1/`n_f' {
		 matrix diag_W [`i',`i']=  W[`i',1]
}
/*
 From "stata canon.pdf" in the "Stored results" section:
e(rawcoef_var1): raw coefficients for varlist1
e(rawcoef_var2): raw coefficients for varlist2
 */
 matrix list e(rawcoef_var1)
 matrix list e(rawcoef_var2)
 *Let's use the raw coeff, after all, all variables are in the same scale (?)
 matrix list e(stdcoef_var1)  
 matrix list e(stdcoef_var2)   

 mata
    myrawcoef_var1=st_matrix("e(rawcoef_var1)")
	myrawcoef_var2=st_matrix("e(rawcoef_var2)")
/*The matrices X and Y are calculated in Stata so they should still be available*/	
	myOgX=st_matrix("OgX")
	myOgY=st_matrix("OgY")
    mydiag_W=st_matrix("diag_W")
	U=myOgX*myrawcoef_var1
	V=myOgY*myrawcoef_var2
	
	/*Now U and V are not complex matrices*/
	st_matrix("OGStataU", U)
	st_matrix("OGStataV", V)
	
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
		meanU
		meanV=mean(V[,i])
		meanV
		weighted_varU[1,i]=(U'[i,])*mydiag_W*(U[,i])/trace(mydiag_W)
		weighted_varV[1,i]=(V'[i,])*mydiag_W*(V[,i])/trace(mydiag_W)
		weighted_rho[1,i]=((U'[i,])*mydiag_W*(V[,i])/trace(mydiag_W))/sqrt(weighted_varU[1,i]*weighted_varV[1,i])
	   
	}/*end for loop*/
	/*Notice that in this case the variances are (approximately) equal to 1 because with out normalization of the eigen values the constrain used during the calculation of the caonical correlations holds*/
    /*Now the means are zero and variances are 1 because they have been standarized*/
	weighted_varU
	weighted_varV
/*The final correlation is equal to what we got form svy: regress*/
	weighted_rho
end 

/*matrix list StataU*/
svmat OGStataU
svmat OGStataV
/*Notice that the mean zero and std. dev. 1 of the canonical variates means that the regression coefficients are eqaul to them!!!:-)*/
svyset  [pweight=finwgt],
svy: regress OGStataU1  OGStataV1   
svy: regress OGStataU2  OGStataV2   
svy: regress OGStataU3  OGStataV3
svy: regress OGStataU4  OGStataV4  
svy: regress OGStataU5  OGStataV5   
svy: regress OGStataU6  OGStataV6  

/*Remember We don't get the same results if we just shift the variables!!!*/
/*Complex survey elements taken form  2021-NYTS-Codebook.pdf page 30 (33 of 41)*/
/*Notice that the mean zero and std. dev. 1 of the canonical variates means that the regression coefficients are eqaul to them!!!:-)*/
svyset  [pweight=finwgt], psu(psu2) strata(v_stratum2)
svy: regress OGStataU1  OGStataV1  
svy: regress OGStataU2  OGStataV2  
svy: regress OGStataU3  OGStataV3
svy: regress OGStataU4  OGStataV4  
svy: regress OGStataU5  OGStataV5   
svy: regress OGStataU6  OGStataV6  

/*-----------------Calculating other Statistics (more than anything to make sure that we get the same results as in SAS)---------------*/
egen z2qn9 = std(qn9)
egen z2qn38 = std(qn38) 
egen z2qn40 = std(qn40)
egen z2qn53 = std(qn53) 
egen z2qn54 = std(qn54)
egen z2qn64 = std(qn64) 
egen z2qn69 = std(qn69)
egen z2qn74 = std(qn74) 
egen z2qn76 = std(qn76)
egen z2qn78 = std(qn78) 
egen z2qn80 = std(qn80)
egen z2qn82 = std(qn80) 
egen z2qn85 = std(qn85)
egen z2qn88 = std(qn88) 
egen z2qn89 = std(qn89) 

egen z2qn128 = std(qn128)
egen z2qn129 = std(qn129) 
egen z2qn130 = std(qn130)
egen z2qn131 = std(qn131)
egen z2qn132 = std(qn132)
egen z2qn134 = std(qn134)

*Not standarizing we don't get the correct canonical correlations
mkmat z2qn9 z2qn38 z2qn40 z2qn53 z2qn54 z2qn64 z2qn69 z2qn74 z2qn76 z2qn78 z2qn80 z2qn82 z2qn85 z2qn88 z2qn89, matrix(X)
mkmat z2qn128 z2qn129 z2qn130 z2qn131 z2qn132 z2qn134, matrix(Y) 

mata
    trace(mydiag_W)
	myX=st_matrix("X")
	myY=st_matrix("Y")

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
	/*.0415381409 vs 0.044,  .9582477164 vs 0.9592,  .8388719508  vs 0.839*/
	/*p_val_Wilks_LambdaFREQ*/
	/*  0 vs <.0001,   0  vs <.0001,  0  vs 0.0021  |*/
end 