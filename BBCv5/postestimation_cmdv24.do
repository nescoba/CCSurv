 *3-3-2023: Transforming pprograms into a custom postestimation command
 *v15: Appplying conclusions from Fizing indexes-4-3-2023.docx and allowing user to provide smaller data set first
 *v16: Capturing all partial results to display 
 *v17: Adding argument so that we don't display all CC's
 *v18: Implementing the formulas in Calinski A Comparison of Some Tests for Determining the Number of Nonzero Canonical Correlations 
 *v19: Swithching to Lawley's and Fujikoshi's approximation
 *v20: Changing order of results to match SAS and Stata output
 *v21: Adding "classic" capability 
 *v22: Add UV/VU results https://www.stata.com/support/faqs/statistics/estimate-correlations-with-survey-data/
 /*For the linearization variance estimator used by svy: regress and regress, robust, the p-value for the slope of the regression of Y on X is not the same as the p-value for the regression of X on Y, unlike the case for the OLS regression estimator. (Try it!) It's really difficult to put into words why this is the case.This is my advice: use correlate with aweights for point estimates of the correlation (as you are, no doubt, doing) and use svy: regress for p-values. Do svy: regress y x and svy: regress x y and take the biggest p-value—that is the conservative thing to do. This is a fine (and perhaps superior) p-value for the test of rho = 0.*/
 *v23: Add graphs
 *v24: Maek graphs really optional "ry using args or syntax but not both; these are alternative methods of parsing command arguments." https://www.statalist.org/forums/forum/general-stata-discussion/general/1418935-passing-a-local-or-scalar-to-a-program-as-additional-argument
/*Remember to run this before running the program, after all it is a custom postestimation command for canon*/
*need to re-run mata drop to make stata forget the definition of the function, i.e. every time that the function is modified this command should be run
*  mata: mata drop calcpval() 
capture program drop csdcanon  
program csdcanon
   *ereturn list
   * preserving user's data
   preserve
   version 17
   
   args svysetcommand howmany dim1 dim2
  * syntax  varlist [, dim1(numeric int min=1) dim2(numeric int min=1) ]
   *saving original data set
    mkmat _all, matrix(All)
	
   tokenize `e(cmdline)', parse("()=]")
  * display "1=|`1'|, 2=|`2'|, 3=|`3'|, 4=|`4'|, 5=|`5'|, 6=|`6'|, 7=|`7'|, 8=|`8'|, 9=|`9'|,10=|`10'|"
   /*tokens 3 and 6 also work with NYTS data*/
   /*
   di "`3'"
   di "`6'"
   di "`10'"
   display("here")
   */
   scalar varcount1 = `:word count `3''
   scalar varcount2 = `:word count `6''
   /*Reading the variables in their original position in canon becaus ethe coefficients expected in this order*/
   mkmat `3', matrix(auxOgX) 
  * matlist auxOgX
   mkmat `6', matrix(auxOgY) 
  * matlist auxOgY
   
   /*Captuing all the variables names*/
   local names "`3' `6'"
  * display("`names'")
   /*First set of variables is the longest*/
   if varcount1 >= varcount2 {
   	  * display("here2")
       local X `3'
       local Y `6'
	   local totalvar= varcount2 +  1
	   
   }
   /*Second set of variables is the longest*/
  else {
  	 *  display("here3")
       local X `6'
       local Y `3'
	   local totalvar= varcount1 + 1
   }
  /*
   di "`X'"
   di "`Y'"
   display("here4")
   */
   /*Reading the weights*/
   local finwgt `10'
   mkmat `X', matrix(OgX) 
   *matlist OgX
   mkmat `Y', matrix(OgY) 
   *matlist OgY
   /*The variable varnames has the values of the variables in the order that canon received them*/	
   quietly gen varname = ""
   tokenize "`names'" 
   
  * local n: word count `names'
   local n=varcount1 + varcount2
   
   gen counter = _n
   
   forvalues i = 1/`n' {
   	/*getting name i*/
    local v1 : word `i' of `names'
	* display("`i'")
	* display("`v1'")
	/*assign name i to value i of the variable varname*/
    quietly replace varname="`v1'" if counter == `i'
   }
   
/*===============Graphs===============*/
 if "`dim1'"!=""{ 
/*Notice that the order of the std coeff is still the same as in canon hence the variables names are in the correct order*/
   matrix mycoef1 =e(stdcoef_var1)
   matrix stdU=  auxOgX * mycoef1
   svmat double stdU
   matrix mycoef1=mycoef1,J(rowsof(mycoef1),1,0) /*Adding a columns of zeros to aid with the painting of the markers*/
   matrix mycoef2 =e(stdcoef_var2)
   matrix stdV=  auxOgY * mycoef2
   svmat double stdV
   matrix mycoef2=mycoef2,J(rowsof(mycoef2),1,1) /*Adding a columns of ones to aid with the painting of the markers*/
 
/*=========Variables=================*/
   mat graphcoef=mycoef1\mycoef2
   svmat double graphcoef
  * matlist graphcoef
  * list in 1/5
   *Needs to add different color and shape to the markers depending on which set of variables they came from
   quietly twoway scatter graphcoef`dim2' graphcoef`dim1' if graphcoef`totalvar' == 0, mc(red)  ms(Oh) mlabel(varname) yline(0) xline(0)  ysize(8) xsize(10) || ///
          scatter graphcoef`dim2' graphcoef`dim1'  if graphcoef`totalvar' == 1, mc(blue) ms(+)   mlabel(varname) || ///
		  function y = sqrt(1 - (x)^2), range(-1 1) lwidth(thin) lcolor(black) || function y = -sqrt(1 - (x)^2), range(-1 1) lwidth(thin) lcolor(black) aspect(1) legend(off) ytitle("CC`dim2'") xtitle("CC`dim1'") saving(variables,replace)  nodraw

/*=========Units=================*/
	 *list in 1/7
	quietly scatter  stdV`dim1' stdU`dim1' , ytitle("V`dim1'") xtitle("U`dim1'") mlabel(counter) ysize(8) xsize(10) saving(units,replace) nodraw
	gr combine variables.gph units.gph
	
 }/*end if graphs*/	
/*===============End Graphs===============*/  

  if ("`svysetcommand'"!="classic"){
  	 mkmat `finwgt', matrix(W) 
     /*Creating diagonal matrix with the survey Weights*/
     matrix diag_W = I(rowsof(OgX))
     local n_f =rowsof(OgX)
     forvalues i = 1/`n_f' {
   	 		 matrix diag_W [`i',`i']=  W[`i',1]
	 } /*end for i*/ 
	 
   }  /*end if*/
   else{ /*Creating a weight matrix with all weights equal to 1 for the classic case*/
   	       matrix W = J(rowsof(OgX),1,1) 
			matrix diag_W = I(rowsof(OgX))
            local n_f =rowsof(OgX)
            forvalues i = 1/`n_f' {
   	 		       matrix diag_W [`i',`i']=  1
	        } /*end for i*/
	}/*end else*/
	  
   *===standarizing X variables===
   clear
   capture svmat double OgX
   foreach var of varlist _all {
   	egen aux_`var'= std(`var')
    drop `var'
   }/*end for var*/
  *if we do not standarize we don't get the correct canonical correlations
   mkmat  _all, matrix(X)

   *standarixing Y variables
   clear
   capture svmat double OgY
   foreach var of varlist _all {
   	egen aux_`var'= std(`var')
    drop `var'
   }
 	*if we do not standarize we don't get the correct canonical correlations
   mkmat  _all, matrix(Y)  

  * mkmat sends the matrices to mata so there is no need to send them as arguments, e(rawcoef_var1) and e(rawcoef_var1) are sent to mata automatically when they are created by canon
  clear
  *capturing number of canonical correlations, e.g. e(n_cc) =  3 has the number of canonical correlations; it is a scalar
  local n_cc  = e(n_cc)
  *creating matrix where all new p-values will be stored
   matrix Results = J(7*(`n_cc'),6,.) 
  *Runing functions that calculates existing methods p-values
   mata: calcpval() 

   matrix colnames Results = "Statistic" "df1" "df2" "Chi-Sq/F" "p-val" "Index"
   capture svmat OGStataUVW
   capture svmat All, names(col) 
  *ereturn list
*=====Weighthed Simple Regression============ 
 local weightindex = (2*e(n_cc)) + 1 /*Reading the weiths as they where stored in this column in the mata function*/
  if ("`svysetcommand'"!="classic"){
   quietly svyset  [pweight= OGStataUVW`weightindex'],
   *Notice coefficients of the simple linear regression are equal to the canonical correlations because the variances of the canonical variates are equal to 1
   forval i = 1/`n_cc' { 
	local     secondindex = `i' + `n_cc'
    quietly svy: regress OGStataUVW`i'  OGStataUVW`secondindex'
	matrix p1 = e(p)
    quietly svy: regress OGStataUVW`secondindex' OGStataUVW`i' 
	matrix p2 = e(p)
	if (p1[1,1]>=p2[1,1]){
		quietly svy: regress OGStataUVW`i'  OGStataUVW`secondindex'
	    matrix Results[(7*(`i'-1))+6,1]= e(b)
        matrix Results[(7*(`i'-1))+6,2]= `e(df_m)'
	    matrix Results[(7*(`i'-1))+6,3]= `e(df_r)'
	    matrix Results[(7*(`i'-1))+6,4]= `e(F)'
	    matrix Results[(7*(`i'-1))+6,5]= `e(p)'
	    matrix Results[(7*(`i'-1))+6,6]= `i'
	}
	else{
		quietly svy: regress OGStataUVW`secondindex' OGStataUVW`i' 
	    matrix Results[(7*(`i'-1))+6,1]=  e(b)
	    matrix Results[(7*(`i'-1))+6,2]= `e(df_m)'
	    matrix Results[(7*(`i'-1))+6,3]= `e(df_r)'
	    matrix Results[(7*(`i'-1))+6,4]= `e(F)'
	    matrix Results[(7*(`i'-1))+6,5]= `e(p)'
	    matrix Results[(7*(`i'-1))+6,6]= `i'
	}/*end else p-values comparision*/
   } /*end for i*/
  } /*end if*/
  else{
  	/*For classic there is no need to check the p-values, they are the same*/
  	forval i = 1/`n_cc' { 
	local     secondindex = `i' + `n_cc'
	*summarize OGStataUVW`i'
    quietly regress OGStataUVW`i'  OGStataUVW`secondindex'
	matrix aux = e(b)
	matrix Results[(7*(`i'-1))+6,1]= aux[1,1]
	matrix Results[(7*(`i'-1))+6,2]= e(df_m)
	matrix Results[(7*(`i'-1))+6,3]= e(df_r)
	matrix Results[(7*(`i'-1))+6,4]= e(F)
	matrix Results[(7*(`i'-1))+6,5]= Ftail(e(df_m),e(df_r),e(F))
   }/*end for i*/
  }/*end else "classic"*/
 *=====End Weighthed Simple Regression============
 
 *=====Complex Survey Design Simple Regression============
* Simple linear regression this time using all the svyset factors
 *How do I make `1' a command that I can run?
 *display "The 1st argument you typed is: `svysetcommand'"
 if ("`svysetcommand'"!="classic"){
 quietly `svysetcommand'
 forval i = 1/`n_cc' { 
	local     secondindex = `i' + `n_cc'
	quietly svy: regress OGStataUVW`i'  OGStataUVW`secondindex'
	matrix p1 = e(p)
    quietly svy: regress OGStataUVW`secondindex' OGStataUVW`i' 
	matrix p2 = e(p)
	if (p1[1,1]>=p2[1,1]){
		quietly svy: regress OGStataUVW`i'  OGStataUVW`secondindex'
	    matrix Results[(7*(`i'-1))+7,1]=  e(b)
	    matrix Results[(7*(`i'-1))+7,2]= `e(df_m)'
	    matrix Results[(7*(`i'-1))+7,3]= `e(df_r)'
	    matrix Results[(7*(`i'-1))+7,4]= `e(F)'
	    matrix Results[(7*(`i'-1))+7,5]= `e(p)'
	    matrix Results[(7*(`i'-1))+7,6]= `i'
	}
	else{
		quietly svy: regress OGStataUVW`secondindex' OGStataUVW`i' 
	    matrix Results[(7*(`i'-1))+7,1]= e(b)
	    matrix Results[(7*(`i'-1))+7,2]= `e(df_m)'
	    matrix Results[(7*(`i'-1))+7,3]= `e(df_r)'
	    matrix Results[(7*(`i'-1))+7,4]= `e(F)'
	    matrix Results[(7*(`i'-1))+7,5]= `e(p)'
	    matrix Results[(7*(`i'-1))+7,6]= `i'
	}/*end else p-values comparision*/
  }/*end for i*/ 
 }/*end if svyset*/
 *=====End Complex Survey Design Simple Regression============
 *Present results from a matrix
/*In case that this argument was not provided we display all canonical correlations*/
    if "`howmany'"==""{
	local	howmany=`n_cc'
	}
forval i = 1/`howmany' { 
	 /*Selecting rows in information for Canonical correlation i*/
	 matrix aux= Results[((7*(`i' -1))+1),1]
	 local mag  = round(aux[1,1], .0001)
	 mat ResultsAux =Results[((7*(`i' -1))+2)..((7*(`i' -1))+7),1..5]
	 matrix rownames ResultsAux = "Wilks' Lambda"  "Pillai's Trace"  "Hotelling-Lawley Trace" "Roy's Greatest Root" "Weighted Reg" "Complex Survey Design Reg"
     matlist ResultsAux,title("Statistics for Canonical Correlation: "`i')  twidth(30)  format(%10.4f)   rowtitle( "Canonical Correlation="`mag')  border(top bottom)
} 

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
    myOgX=st_matrix("OgX")
	myOgY=st_matrix("OgY")
	myrawcoef_var1=st_matrix("e(rawcoef_var1)")
	myrawcoef_var2=st_matrix("e(rawcoef_var2)")
	/*Inverting position of coefficients because the smallest data set was provided first*/
	if (rows(myrawcoef_var2)>rows(myrawcoef_var1)){
       aux=myrawcoef_var2
	   myrawcoef_var2=myrawcoef_var1
       myrawcoef_var1=aux
 	}
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
		myResults[(7*(i-1))+1,1]=weighted_rho[1,i]    /*row 1, 7,... */
		myResults[(7*(i-1))+1,2]=.                    /*row 1, 7,... */
		myResults[(7*(i-1))+1,3]=.                    /*row 1, 7,... */
		myResults[(7*(i-1))+1,4]=.                  /*row 1, 7,... */
		myResults[(7*(i-1))+1,5]=.                   /*row 1, 7,... */
		
	}/*end for loop*/

/*Finding sum of frequencies*/	
 sumFREQ=sum(trunc(myW))
 /*sumFREQ*//*It works*/
 /*=======Lawley's approximation term========*/
 Lawley= J(1,cols(weighted_rho),0)
 for (j=1; j<=cols(weighted_rho); j++) {
 	for (i=1; i<=j; i++) { /*Notice Eq. (8) starts in 1, i.e. when d=0 then i=0+1=1*/
	        Lawley[1,j]=Lawley[1,j]+(1/(weighted_rho[1,i]^2)) /*Maybe Gittins formula below Table 3.3 is wrong and it should be summation as in Calinski's formulas*/
	}/*end for i loop*/
 }/*end j for loop*/
 /*=======End Lawley's approximation term========*/

 /*=================== Start Wilks' Lambda (Barlett approximation)================*/	
/*Gittins calls it Likelihood ratio test*/
    Lambda= J(1,cols(weighted_rho),1)
	df= J(1,cols(weighted_rho),1)
	ChiSq_Wilks_Lambda=J(1,cols(weighted_rho),1)
	ChiSq_Wilks_LambdaFREQ=J(1,cols(weighted_rho),1)
	p_val_Wilks_Lambda=J(1,cols(weighted_rho),1)
	p_val_Wilks_LambdaFREQ=J(1,cols(weighted_rho),1)
	/*Independence Likelihood ratio test Eq. (3.12) and (3.13) in page 57 when i=1*/
	/*Dimensionality Likelihood ratio test Eq. (3.19) and (3.20) in page 60 when i=2,...,s*/
 	for (i=1; i<=cols(weighted_rho); i++) {
	   Lambda[1,i]=(1-(weighted_rho[1,i]^2))
       for (j=(i+1); j<=cols(weighted_rho); j++) {
			Lambda[1,i]=Lambda[1,i]*(1-(weighted_rho[1,j]^2)) /*Eq. (3.19) of Gittins and Eq. (7) of Calinski*/
		}/*end for j*/
	   df[1,i]=(cols(myX)+1-i)*(cols(myY)+1-i) /*Same degrees of freedom for Batlett/Gittins and Lawley*/
	   /*ChiSq_Wilks_Lambda[1,i]=(((rows(myX)-1)-.5*(cols(myX)+cols(myY)+1))*ln(Lambda[1,i]))*-1 *//*Eq. 3.20 of Gittins and Eq. (10) of Calinski*/
	  /* ChiSq_Wilks_LambdaFREQ[1,i]=(((sumFREQ-1)-.5*(cols(myX)+cols(myY)+1))*ln(Lambda[1,i]))*-1 *//*Eq. 3.20 of Gittins and Eq. (10) of Calinski*/
	   ChiSq_Wilks_Lambda[1,i]=(((rows(myX)-1)  -.5*(cols(myX)+cols(myY)+1) + Lawley[1,i] -i)*ln(Lambda[1,i]))*-1 /*Eq. below Table 3.3 of Gittins and Eq. (12) of Calinski*/
	   ChiSq_Wilks_LambdaFREQ[1,i]=(((sumFREQ-1)-.5*(cols(myX)+cols(myY)+1) + Lawley[1,i] -i)*ln(Lambda[1,i]))*-1 /*Eq. below Table 3.3 of Gittins and Eq. (12)  and Eq. (10) of Calinski*/
	   p_val_Wilks_Lambda[1,i]=chi2tail(df[1,i], ChiSq_Wilks_Lambda[1,i])
	   p_val_Wilks_LambdaFREQ[1,i]=chi2tail(df[1,i], ChiSq_Wilks_LambdaFREQ[1,i])
	   myResults[(7*(i-1))+2,1]= Lambda[1,i]	         /*row 2, 10,... */
	   myResults[(7*(i-1))+2,2]= df[1,i]	         /*row 2, 10,... */
	   myResults[(7*(i-1))+2,3]=.	                /*row 2, 10,... */
	   myResults[(7*(i-1))+2,4]= ChiSq_Wilks_Lambda[1,i]           /*row 2, 10,... */
	   myResults[(7*(i-1))+2,5]= p_val_Wilks_Lambda[1,i]             /*row 2, 10,... */
	   /*------FREQ------*/
	   /*myResults[(7*(i-1))+6,1]= Lambda[1,i]*/	         /*row 6, 14,... */
	  /* myResults[(7*(i-1))+6,2]= df[1,i]*/	         /*row 6, 14,... */
	  /* myResults[(7*(i-1))+6,3]=.*/	                 /*row 6, 14,... */
	  /* myResults[(7*(i-1))+6,4]= ChiSq_Wilks_LambdaFREQ[1,i] */      /*row 6, 14,... */
	   /*myResults[(7*(i-1))+6,5]= p_val_Wilks_LambdaFREQ[1,i] */        /*row 6, 14,... */
	}/*end for loop*/
/*=================== End Wilks' Lambda================*/	
	
/*===================Start Pillai's Trace =====================*/
/*Implementing the formulas in Calinski A Comparison of Some Tests for Determining the Number of Nonzero Canonical Correlations Eq(8) & (11) -(14)*/
    V=J(1,cols(weighted_rho),0)
	Pillais_Trace_stat = J(1,cols(weighted_rho),1)
	Pillais_Trace_p_value = J(1,cols(weighted_rho),1)
	for (j=1; j<=cols(weighted_rho); j++) { 
		for (i=j; i<=cols(weighted_rho); i++) { /*Notice Eq. (8) starts in 1, i.e. when d=0 then i=0+1=1*/
	        V[1,j]=V[1,j]+(weighted_rho[1,i]^2) /*Eq(8)of Calinski*/
	     }/*end for i loop*/
		 /*The degrees of freedom are the same as in SAS 9 & 48,  in Eq (8) of of All formulas p-values.pdf */
		/*Pillais_Trace_stat[1,j]=(rows(myX)-1*V[1,j]*/ /*Eq. (11) of Calinski*/
		Pillais_Trace_stat[1,j]=(rows(myX)-1-(2*j) + Lawley[1,j])*V[1,j] /*Eq. (14) of Calinski*/
	    Pillais_Trace_p_value[1,j]=chi2tail((cols(myX)+1-j)*(cols(myY)+1-j), Pillais_Trace_stat[1,j])
		myResults[(7*(j-1))+3,1]= V[1,j]	             /*row 3, 11,... */
		myResults[(7*(j-1))+3,2]= (cols(myX)+1-j)*(cols(myY)+1-j)	/*Notice the -1 to transform the i into d*/ /*row 3, 11,... */
	    myResults[(7*(j-1))+3,3]= .	             /*row 3, 11,... */
		myResults[(7*(j-1))+3,4]= Pillais_Trace_stat[1,j] 	       /*row 3, 11,... */
		myResults[(7*(j-1))+3,5]= Pillais_Trace_p_value[1,j]	         /*row 3, 11,... */
	}/*end for j*/
/*===================End Pillai's Trace=====================*/	

/*===================Start Hotelling-Lawley Trace =====================*/
/*Implementing the formulas in Calinski A Comparison of Some Tests for Determining the Number of Nonzero Canonical Correlations Eq(6) & (9) -(13)*/	
	U=J(1,cols(weighted_rho),0)
	Hotelling_Lawley_Trace_stat = J(1,cols(weighted_rho),1)
	Hotelling_Lawley_Trace_p_value = J(1,cols(weighted_rho),1)
    for (j=1; j<=cols(weighted_rho); j++) {
		for (i=j; i<=cols(weighted_rho); i++) {
	           U[1,j]=U[1,j]+(weighted_rho[1,i]^2)/(1-(weighted_rho[1,i]^2)) /*Eq(6) of Calinski*/
        }/*end for i loop*/
		/*Hotelling_Lawley_Trace_stat[1,j]= (rows(myX)-cols(myX)-cols(myY)-2)*U[1,j] *//*Eq. (9) of Calinski Bartlett's approximatio*/
		Hotelling_Lawley_Trace_stat[1,j]= (rows(myX)-cols(myX)-cols(myY)-2 + Lawley[1,j])*U[1,j] /*Eq. (13) of Calinski Fujikoshi approximation*/
		Hotelling_Lawley_Trace_p_value[1,j]=chi2tail((cols(myX)+1-j)*(cols(myY)+1-j), Hotelling_Lawley_Trace_stat[1,j])
	    myResults[(7*(j-1))+4,1]= U[1,j]	             /*row 4, 12,... */
	    myResults[(7*(j-1))+4,2]= (cols(myX)+1-j)*(cols(myY)+1-j)             /*row 4, 12,... */
	    myResults[(7*(j-1))+4,3]= .	             /*row 4, 12,... */
		myResults[(7*(j-1))+4,4]= Hotelling_Lawley_Trace_stat[1,j]  /*row 4, 12,... */
		myResults[(7*(j-1))+4,5]= Hotelling_Lawley_Trace_p_value[1,j] /*row 4, 12,... */
	}/*end for j*/
/*==================End Hotelling-Lawley Trace=====================*/

/*===================Start Roy's Greatest Root================*/	
/*Gittins calls it Union-intersection test, page 61*/	
	/*Approximating Roy's Greatest Root using Eq. (12) of APPROXIMATE NULL DISTRIBUTION OF THE LARGEST ROOT IN MULTIVARIATE ANALYSIS1 BY IAIN M. JOHNSTONE*/
	/*Pθ(s,m,n)>θα)>PF(θα)=1−Fν1,ν2(ν2θα/ν1(1−θα) where ν1 =s+2m+1 and ν2 =s+2n+1 denote the hypothesis and error degrees of freedom respectively.
	The definitions of s=min(p,q), m=(|p-q|-1)/2 and n=(N-p-q-2)/2 is given in Gittins right above Eq. (3.17)
	*/
	p=cols(myX)
	q=cols(myY)
	pq=J(1,2,1)
	pq[1,1]=cols(myX)
	pq[1,2]=cols(myY)
	largest_root=J(1,cols(weighted_rho),1)
	v1=J(1,cols(weighted_rho),1)
	v2=J(1,cols(weighted_rho),1)
	Roys_Greatest_Root_stat = J(1,cols(weighted_rho),1)
	Roys_Greatest_Root_p_value = J(1,cols(weighted_rho),1)
	/*First we have the so called "Overall test" which is described in Eq. (3.17) (Independence Union-Intersection test) 
	and is testing the hypothesis in Eq. (3.15) which is equivalent to Eq. (3.11), i.e., Eq. (3.9). */
	s=min(pq)
	m=(abs(p-q)-1)/2 
	n=(rows(myX)-p-q-2)/2
	v1[1,1] =s + (2*m) + 1
	v2[1,1] =s + (2*n) + 1
	largest_root[1,1]=weighted_rho[1,1]^2
	Roys_Greatest_Root_stat[1,1]=(v2[1,1]*largest_root[1,1])/(v1[1,1]*(1.0 - largest_root[1,1]))
    Roys_Greatest_Root_p_value[1,1]=Ftail(v1[1,1], v2[1,1], Roys_Greatest_Root_stat[1,1])
	
	myResults[5,1]= largest_root[1,1]
	myResults[5,2]= v1[1,1]	                 /*row 5, 13,... */
	myResults[5,3]= v2[1,1]	                 /*row 5, 13,... */
	myResults[5,4]= Roys_Greatest_Root_stat[1,1]        /*row 5, 13,... */
	myResults[5,5]= Roys_Greatest_Root_p_value[1,1]        /*row 5, 13,... */
	/*10.48376021 vs 10.42 in SAS (F Value)
	.0004667007 vs 0.0005 from SAS*/
	/*For the roots k where k>=2 we use the formulas at the bottom of page 61 of Gittins, 
	According to Fixing indexes-3-31-2023v2.docx the formulas in Eq. (3.17) are not reduced to the formulas at the bottom of 
	page 61 because we are testing different hypotheses
	*/
	for (k=1; k<=(cols(weighted_rho)-1); k++) {
		/*formula at the bottom of page 61 of Gittins*/
		s=min(pq)-(k)-1
		/*formulas to match Eq(3.17) when k=1*/
		/*s=min(pq)-(k)+1*/
		/*formula at the bottom of page 61 of Gittins*/
		m=(abs(p-q-(k))-1)/2
		/*formulas to match Eq(3.17) when k=1*/
		/*m=(abs(p-q-(k)+1)-1)/2*/
		/*formula at the bottom of page 61 of Gittins*/
		n=(rows(myX)-p-q-(k))/2
		/*formulas to match Eq(3.17) when k=1*/
		/*n=(rows(myX)-p-q-(k+1))/2*/
		v1[1,k+1] =s + (2*m) + 1
		v2[1,k+1] =s + (2*n) + 1
		largest_root[1,k+1]=weighted_rho[1,k]^2 /*Remember that the the cycle starts with r_1 and ends at r_s-1, see Fizing indexes-4-3-2023.docx*/
		/*
		if(v1==0){
			v1=.000001
		}
		*/
		Roys_Greatest_Root_stat[1,k+1]=(v2[1,k+1]*largest_root[1,k+1])/(v1[1,k+1]*(1.0 - largest_root[1,k+1])) /*We are storing the results in k+1 or we will delete the independence union-intersection test*/
	    Roys_Greatest_Root_p_value[1,k+1]=Ftail(v1[1,k+1], v2[1,k+1], Roys_Greatest_Root_stat[1,k+1])
		myResults[(7*(k+1-1))+5,1]= largest_root[1,k+1]	 /*row 5, 13,... */
		myResults[(7*(k+1-1))+5,2]= v1[1,k+1]	                 /*row 5, 13,... */
		myResults[(7*(k+1-1))+5,3]= v2[1,k+1]                 /*row 5, 13,... */
		myResults[(7*(k+1-1))+5,4]= Roys_Greatest_Root_stat[1,k+1]        /*row 5, 13,... */
		myResults[(7*(k+1-1))+5,5]= Roys_Greatest_Root_p_value[1,k+1]        /*row 5, 13,... */
	}/*end for k*/
/*Very bottom of page 61 "Like the corresponding likelihood ratio tests, these tests pf the r_k^2 for K>=2 are very conservative (Mardia, Kent & Bibby, 1979, p. 147)."*/		
/*These p-values are correct!!!!! At least the first one which is the only one provided by SAS*/	
/*===================End Roy's Greatest Root================*/	
for (i=1;i<=cols(weighted_rho); i++ ){/*i will control the rows*/
     /*Row 1: CC, Row 2: Wilk's Lambda, Row 3: Pillai's Trace , Row 4: Hotelling Trace , Row 5:Roy's Greatest Root, Row 6: Wilk's Lambda FREQ , Row 7 Weighted Regression, Row 8: CSD Reg*/
	 /*Column 6: Just an index that keeps track of which correlation number each row belongs to*/
	 myResults[(7*(i-1))+1,6]= i        /*row 1, 9,... */
	 myResults[(7*(i-1))+2,6]= i        /*row 2, 10,... */
	 myResults[(7*(i-1))+3,6]= i        /*row 3, 11,... */
	 myResults[(7*(i-1))+4,6]= i        /*row 4, 12,... */
	 myResults[(7*(i-1))+5,6]= i        /*row 5, 13,... */
	 myResults[(7*(i-1))+6,6]= i        /*row 6, 14,... */
}/*end for i*/
st_matrix("Results", myResults)
}	
end
/*========== End MATA function ==========*/