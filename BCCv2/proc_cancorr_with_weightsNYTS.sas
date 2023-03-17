dm 'log;clear;output;clear;';
ods html close; /* close previous */
ods html; /* open new */
proc import datafile="C:\Users\raulcruz\OneDrive - Indiana University\Documents\students\carlos\stata\cmd\selectednyts2021.csv"
        out=nyts
        dbms=csv
        replace;
       getnames=YES;
run;
title "PROC CANCORR with WEIGHT";
proc cancorr data=nyts all
vprefix=MTEPUse vname='MTEP Use'
wprefix=ECigMarketing wname='ECigMarketing';
var qn9 qn38 qn40 qn53 qn54 qn64 qn69 qn74 qn76 qn78 qn80 qn82 qn85 qn88 qn89;
with qn128 qn129 qn130 qn131 qn132 qn134;
weight finwgt;
run;
title "PROC CANCORR with FREQ";
proc cancorr data=nyts all
vprefix=MTEPUse vname='MTEP Use'
wprefix=ECigMarketing wname='ECigMarketing';
var qn9 qn38 qn40 qn53 qn54 qn64 qn69 qn74 qn76 qn78 qn80 qn82 qn85 qn88 qn89;
with qn128 qn129 qn130 qn131 qn132 qn134;
freq finwgt;
run;

