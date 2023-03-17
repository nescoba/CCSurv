
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
/*Selecting people that uses at least 2 products*/
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
* keep if number_of_products >= 1
 keep if qn5b==1
 *outsheet using selectednyts2021.csv, comma 
 /*
QN128	Are using the Internet, how often do you see ads or promotions for e-cigarettes?
QN129	Read newspapers or magazines, how often do you see ads or promotions for e-cigarettes?
QN130	Go to a convenience store, supermarket, gas station, kiosk/storefront, or shopping center, how often do you see ads or promotions for e-cigarettes?
QN131	Watch T.V. or streaming services (such as Netflix, Hulu, or Amazon Prime), or go to the movies, how often do you see ads or promotions for e-cigarettes?
QN132	Watch T.V. or streaming services (such as Netflix, Hulu, or Amazon Prime), or go to the movies, how often do you see people or characters using e-cigarettes?
QN134	Use social media, how often do you see posts or content (pictures, videos, or text) related to e-cigarettes?
*/
 canon (qn9 qn38 qn40 qn53 qn54 qn64 qn69 qn74 qn76 qn78 qn80 qn82 qn85 qn88 qn89) (qn128 qn129 qn130 qn131 qn132 qn134) [weight=finwgt]
 *canon (qn9 qn38 qn40 qn53 qn54 qn64 qn69 qn74 qn76 qn78 qn80 qn82 qn85 qn88 qn89) (qn128 qn131 qn132 qn134) [weight=finwgt]
  /*
  return list
  ereturn list
  sreturn list
  */
  *Double quotations must be used or the spaces will cut each pat of the svyset command into a different argument
 csdcanon "svyset [pweight=finwgt], psu(psu2) strata( v_stratum2)"