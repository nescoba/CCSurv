 *11-14-2022: Passing my canonical vaiates from mata to stata and using svy: regress to calculate the canonical correlations and their p-values
 *11-16-2022: Get the canonical raw coefficients from canon so that the we can calculate the canonical variates and the use svy: regress to calculate the canonical correlations and their p-values
 *11-17-2022: Replicating SAS p-values for canonical correlation with weights using  Eq.(3.19) lines 436-460 (approx) More specifically we replicate SAS' Wilks' lambda distribution using the Chi-square approximation in Eq. (3.19)-(3.20) for the overall hypothesis of independence, i.e. the Likeihood ratio test of dimmensionality of association which is what SAS appears to be testing
 *11-28-2022:Replicating SAS p-values for canonical correlation with weights using FREQ instead of WEIGTH Eq.(3.19) lines 436-460 (just text). NOT WORKING the DF in SAS variafrom CC to CC while I always use trace(diagW)
 *11-28-2022: Also adding thw work needed to calculate Roy's Greatest Root and Pillai's Trace
 *Home
 * cd "C:\Users\Raul Cruz-Cano\OneDrive - Indiana University\Documents\students\carlos\stata"
  *Office
  cd "C:\Users\raulcruz\OneDrive - Indiana University\Documents\students\carlos\stata"
  insheet using fit.csv, clear
  canon (weight waist pulse) (chins situps jumps) 
/*Canonical correlation analysis                      Number of obs =         20
Raw coefficients for the first variable set
                 |        1         2         3 
    -------------+------------------------------
          weight |   0.0314   -0.0763    0.0077 
           waist |  -0.4932    0.3687   -0.1580 
           pulse |   0.0082   -0.0321   -0.1457 
    --------------------------------------------
Raw coefficients for the second variable set
                 |        1         2         3 
    -------------+------------------------------
           chins |   0.0661   -0.0710    0.2453 
          situps |   0.0168    0.0020   -0.0198 
           jumps |  -0.0140    0.0207    0.0082 
    --------------------------------------------
----------------------------------------------------------------------------
Canonical correlations:
  0.7956  0.2006  0.0726
----------------SAS Results from CANCORR.PDF---------------------------
  Output 30.1.3 Raw and Standardized Canonical Coefficients
Raw Canonical Coefficients for the Physiological Measurements
           Physiological1     Physiological2      Physiological3
Weight      -0.031404688      -0.076319506       -0.007735047
Waist        .4932416756      0.3687229894        0.1580336471
Pulse       -0.008199315      -0.032051994        0.1457322421
Raw Canonical Coefficients for the Exercises
         Exercises1     Exercises2       Exercises3
Chins  -0.066113986    -0.071041211     -0.245275347
Situps -0.016846231    0.0019737454     0.0197676373
Jumps  0.0139715689    0.0207141063     -0.008167472
*/
egen z2weight = std(weight)
egen z2waist = std(waist) 
egen z2pulse = std(pulse)

egen z2chins = std(chins)
egen z2situps = std(situps) 
egen z2jumps = std(jumps)

*Not standarizing we don't get the correct canonical correlations
mkmat z2weight z2waist z2pulse, matrix(X)
*mkmat weight waist  pulse, matrix(X) 
matrix list X 
mkmat z2chins z2situps z2jumps, matrix(Y) 
*mkmat chins  situps jumps, matrix(Y) 
matrix list Y 

mkmat weight waist  pulse, matrix(OgX) 
matrix list OgX 
mkmat chins  situps jumps, matrix(OgY) 
matrix list OgY 

mkmat weight, matrix(W) 
/*Creating diagonal matrix with the survey Weights*/
matrix diag_W = I(rowsof(X))
local n_f =rowsof(X)
forvalues i = 1/`n_f' {
		 matrix diag_W [`i',`i']=  W[`i',1]
}
*For now let's use the Identity, i.e. all wieghts equal to 1 see if we can replicate the SAS and Stata results
matrix diag_W = I(rowsof(X))
/*My Eq. (28): A= (ydiag(W)X)(x diag(W)X)^(-1)∙(xdiag(W)Y )∙ (ydiag(W)Y)^(-1)*/
*matrix list Y
*matrix list diag_W
*matrix list X
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
	/*Notice bT is on the left in Eq. (26)*/
	lefteigensystem(myA, bT, lambdaSq)
	bT
	/*Instead of haveing [1 x l] for bT we have [l x l] because we have several eigen solutions
 Notice that now the sum of squares of the elements in the row add up to 1 for example see the first row
	0.264470452	0.789758859	0.542132664	
square =0.06994462	0.623719055	0.293907825	0.9875715
This proves that the eigenvectors are the rows of the matrix
	*/
	lambdaSq
	lambda=lambdaSq
	for (i=1; i<=length(lambda); i++) {
      lambda[i,1]=lambdaSq[i,1]^.5 /*https://www.ssc.wisc.edu/sscc/pubs/4-26.htm*/
    }/*end for loop*/
	lambda
	/* My Eq. (29) to calcultate aT = λ^(-1) (b^T ydiag(W)X)(x diag(W)X)^(-1) 
	better do it here because with the normalization of the b's done by the function eigensystem it is not clear that if done separately both 
	a and b will be normalized using the same number to make their length add up to 1 hence I can't trust using Eq. (30)*/
	aT=bT*myfactorYX*(myfactorInvXX)
	for (i=1; i<=length(lambda); i++) {
      aT=(1/lambda[i,1])*aT 
	}/*end for loop*/
	aT	
	/*"The eigenvectors returned by the above routines are scaled to have length (norm) 1." from page 5 of stata eigensystem.pdf.
	That is not he case for the raw coefficients in Stata, R or SAS so our coefficits a and b won't be the same as those provided by 
	the standard function of these system but will be internally consistent. Just make sure that we get the correct canonical correlations
	*/
	/*Calculating Canonical Variate U and V: Eq (2)     U=a^T * x and V=b^T * y   */
	UT=aT*myX'
	UT
	U=UT'
	VT=bT*myY'
	VT
	V=VT'
	/*Variance might not be 1 anymore due to the normalization performed by the eigensystem function on the eigenvectors*/
	/*Using Eq. (8) and (9) to find the correlations and variances of the canonical matrices*/
	/*just to crete a vector with the correct size*/
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
	  weighted_rho
                1
    +--------------+
  1 |  .795608148  |
  2 |  .200556046  |
  3 |  .072570299  |
    +--------------+
  lambda
                1
    +--------------+
  1 |  .795608148  |
  2 |  .200556046  |
  3 |  .072570299  |
    +--------------+
	
	canon funtion from Stata
	Canonical correlations:
  0.7956  0.2006  0.0726
	*/
	
	/*Calculation of standarized coefficients: https://rdrr.io/cran/candisc/src/R/cancor.R: stdb <- diag(sqrt(diag(cov(object$X)))) %*% coef$X*/
	/*In Stata: DiagSqrtDiagCovX = diag(sqrt(diagonal(myCovX)))*/
	/*Do not confuse diagonal() with its functional inverse, diag(); see [M-5] diag( ). diagonal()
    extracts the diagonal of a matrix into a vector; diag() creates a diagonal matrix from a vector.*/
	/*This diagonal matrix is correct (same values as in R) and the calculations are correctly implemented we are just starting with the "wrong" raw coefficients b
	which are different from to those in R, Stata and SAS I prove this in cc_example-11-9-2022.R in which the standarized coefficiets for cancor from CCA and 
	cancor from candisc are different because their raw coefficients are different.
	*/
end

/*-----------Now using survey weights----------*/
/*Creating diagonal matrix with the survey Weights*/
matrix diag_W = I(rowsof(X))
local n_f =rowsof(X)
forvalues i = 1/`n_f' {
		 matrix diag_W [`i',`i']=  W[`i',1]
}
/*My Eq. (28): A= (ydiag(W)X)(x diag(W)X)^(-1)∙(xdiag(W)Y )∙ (ydiag(W)Y)^(-1)*/
*matrix list Y
*matrix list diag_W
*matrix list X
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
	/*Notice bT is on the left in Eq. (26)*/
	lefteigensystem(myA, bT, lambdaSq)
	bT
	/*Instead of haveing [1 x l] for bT we have [l x l] because we have several eigen solutions*/
	lambdaSq
	lambda=lambdaSq
	for (i=1; i<=length(lambda); i++) {
      lambda[i,1]=lambdaSq[i,1]^.5 /*https://www.ssc.wisc.edu/sscc/pubs/4-26.htm*/
    }/*end for loop*/
	lambda
	/* My Eq. (29) to calcultate aT = λ^(-1) (b^T ydiag(W)X)(x diag(W)X)^(-1) 
	better do it here because with the normalization of the b's done by the function eigensystem it is not clear that if done separately both 
	a and b will be normalized using the same number to make their length add up to 1 hence I can't trust using Eq. (30)*/
	aT=bT*myfactorYX*(myfactorInvXX)
	for (i=1; i<=length(lambda); i++) {
      aT=(1/lambda[i,1])*aT 
	}/*end for loop*/
	aT	
	/*"The eigenvectors returned by the above routines are scaled to have length (norm) 1." from page 5 of stata eigensystem.pdf.
	That is not he case for the raw coefficients in Stata, R or SAS so our coefficits a and b won't be the same as those provided by 
	the standard function of these system but will be internally consistent. Just make sure that we get the correct canonical correlations
	*/
	/*Calculating Canonical Variate U and V: Eq (2)     U=a^T * x and V=b^T * y   */
	UT=aT*myX'
	UT
	U=UT'
	VT=bT*myY'
	VT
	V=VT'
	/*Variance might not be 1 anymore due to the normalization performed by the eigensystem function on the eigenvectors*/
	/*Using Eq. (8) and (9) to find the correlations and variances of the canonical matrices*/
	/*just to crete a vector with the correct size*/
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
correlate StataU1 StataU2 StataU3 StataV1 StataV2 StataV3 [w=weight] /*Leads to the correct correlations, i.e. those from MATA, Stata canon & SAS PROC CANCORR

             |  StataU1  StataU2  StataU3  StataV1  StataV2  StataV3
-------------+------------------------------------------------------
     StataU1 |   1.0000
     StataU2 |   0.0079   1.0000
     StataU3 |   0.0006  -0.0008   1.0000
     StataV1 |   0.8133   0.0067   0.0005   1.0000
     StataV2 |   0.0017   0.1938  -0.0002   0.0015   1.0000
     StataV3 |   0.0004  -0.0005   0.0516   0.0003  -0.0001   1.0000
*/
/*
correlate StataU1 StataU2 StataU3 StataV1 StataV2 StataV3 [w=weight],sig
correlate StataU1 StataU2 StataU3 StataV1 StataV2 StataV3 [aw=weight], sig
correlate StataU1 StataU2 StataU3 StataV1 StataV2 StataV3 [fw=weight], sig
option sig not allowed
https://www.stata.com/support/faqs/statistics/estimate-correlations-with-survey-data/
But what about a p-value?
use  svy: regress x y 
r=Coefficient *sqrt(weighted_varU)/sqrt(weighted_varV)
*/
svyset  [pweight=weight],
svy: regress StataU1  StataV1   /* r= (80.74207  * SQRT(.563960908) /SQRT(5567.36064)=0.81264282  vs  .814133001 from my work & 0.8133 from canon  Prob > F  = 0.0013*/
svy: regress StataU2  StataV2   /* r= (4.587849  * SQRT(0.696245958)/SQRT(394.336132)=0.192777868 vs  .195005619 from my work & 0.1938 from canon  Prob > F  = 0.2587*/
svy: regress StataU3  StataV3   /* r= (0.3248543 * SQRT(0.267148451)/SQRT(10.5990683)=0.051574034 vs  .051612239 from my work & .0516 from canon   Prob > F  = 0.8096*/

/*We don't get the same results if we just shift the variables!!!*/
svy: regress StataV1  StataU1   /* r= (.0081916  * SQRT(5567.36064)/ SQRT(.563960908)=0.813896    vs  .814133001 from my work & 0.8133 from canon  Prob > F  =  0.0000*/
svy: regress StataV2  StataU2   /* r= (.0081835  * SQRT(394.336132)/SQRT(0.696245958)=0.194756082 vs  .195005619 from my work & 0.1938 from canon  Prob > F  = 0.4117*/
svy: regress StataV3  StataU3   /* r= (.0081882 * SQRT(10.5990683)/SQRT(0.267148451)/=0.051575798 vs  .051612239 from my work & .0516 from canon   Prob > F  = 0.8070*/
/*p-value for canonical correlations from SAS with WEIGTH statement
	   CC       Pr > F		
1	0.813267	0.044
2	0.193757	0.9592
3	0.051574	0.839
Test of H0: The canonical correlations in the current row and all that follow are zero
Statistic	Pr > F
Wilks' Lambda	0.044
Pillai's Trace	0.1341
Hotelling-Lawley Trace	0.0218
Roy's Greatest Root	0.0005
These statistics, as described in the section "Multivariate Tests" on page 94 in Chapter 4, "Introduction to Regression Procedures

Not the same as those form Stata using svy: regress as recommended

The difference is even bigger using FREQ instead of WEIGHT:
Correlation	Pr > F
0.813267	<.0001
0.193757	<.0001
0.051574	0.0021

*/

/*In Stata's canon aweights and fweights are allowed; see [U] 11.1.6 weight. (stata canon.pdf)
u11.pdf
weightword Meaning
weight default treatment of weights
fweight or frequency frequency weights
pweight sampling weights
aweight or cellsize analytic weights
iweight importance weights*/
 canon (weight waist pulse) (chins situps jumps) [weight=weight]
 /*Now we can take the canonical variates and apply what syas here:https://www.stata.com/support/faqs/statistics/estimate-correlations-with-survey-data/
 Clearly, this population value can be estimated by replacing the various sums with weighted sums over the sample. Effectively, this is what the correlate command will compute for you when you specify aweights.

But what about a p-value?

The above comments about the equivalence of the hypotheses rho = 0 and beta = 0 might make one think,
 "Great, I can just use svy: regress to get a p-value for the correlation coefficient." Well, this is indeed the case, and it is indeed what I recommend, 
 */
 /*
 From "stata canon.pdf" in the "Stored results" section:
e(rawcoef_var1): raw coefficients for varlist1
e(rawcoef_var2): raw coefficients for varlist2
 */
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
/*p-value for canonical correlations from SAS's PROC CANCORR with WEIGTH statement are different from Stata canon with [weight=weight], they are:
	   CC       Pr > F		 Wilks' Lambda Likelihood Ratio
1	0.813267	0.044                        0.32501779
2	0.193757	0.9592
3	0.051574	0.839
but svy: regress is also prepared to accept cluster and strata information as opposed to these procedure/command which can't
Notice that the hypothesis are difference for SAS Test of H0: The canonical correlations in the current row and all that follow are zero
For my work using svy: regress that null hypothesis is about each individual correlation
Is SAS  Wilks' Lambda the same as Wilks' lambda statistic in Gittins Eq. (3.12)?
Can I approximate the Wilks' lambda distribution using the Chi-square in Eq. (3.13) for the overall hypothesis of independence and 
Eq. (3.19)-(3.20) Table 3.3 Likeihood ratio test of dimmensionality of association which is what SAS appears to be testing?

lambdaj Eq.(3.19)=(1-0.813267^2)*(1-0.193757^2)*(1-0.051574^2)=0.325018447 for j=0 is the Likelihood Ratio Likelihood Ratio in the SAS table with weights =0.32501779
using df=pq=3*3=9 then Chi-Square=-m*ln(lmbadaj)=-15.5*ln(0.325018447)=-15.5*-1.123873337=17.42003672 hence a p-value of 0.957469567 or 0.042530433

lambdaj Eq.(3.19)=(1-0.193757^2)*(1-0.051574^2)=0.959898204 for j=1 is the Likelihood Ratio Likelihood Ratio in the SAS table with weights= 0.95989828
using df=(p-1)(q-1)=(3-1)*(3-1)=2*2=4 then Chi-Square=-m*ln(lmbadaj)=-15.5*ln(0.959898204)=-15.5*-0.040928038=0.634384584 hence a p-value of 0.040832801 or 0.959167199

lambdaj Eq.(3.19)=(1-0.051574^2)=0.997340123 for j=2 is the Likelihood Ratio Likelihood Ratio in the SAS table with weights =0.99734011
using df=(p-2)(q-2)=(3-2)*(3-2)=1*1=1 then Chi-Square=-m*ln(lmbadaj)=-15.5*ln(0.997340123)=-15.5*-0.002663421=0.041283029 hence a p-value of 0.161007419 or 0.838992581

These p-values are equal/similar to what I get from SAS PROC CANCORR with weights
so I have the correct likelihood ratio and theChi-Square distribution of Gittins can approximate the F-Dist of SAS quite well
*/
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
	/*.0415381409 vs 0.044,  .9582477164 vs 0.9592,  .8388719508  vs 0.839*/
	/*p_val_Wilks_LambdaFREQ*/
	/*  0 vs <.0001,   0  vs <.0001,  0  vs 0.0021  |*/
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

