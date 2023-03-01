data Fit;
input Weight Waist Pulse Chins Situps Jumps;
datalines;
191 36 50 5 162 60
189 37 52 2 110 60
193 38 58 12 101 101
162 35 62 12 105 37
189 35 46 13 155 58
182 36 56 4 101 42
211 38 56 8 101 38
167 34 60 6 125 40
176 31 74 15 200 40
154 33 56 17 251 250
169 34 50 17 120 38
166 33 52 13 210 115
154 34 64 14 215 105
247 46 50 1 50 50
193 36 46 6 70 31
202 37 62 12 210 120
176 37 54 4 60 25
157 32 52 11 230 80
156 33 54 15 225 73
138 33 68 2 110 43
;
proc cancorr data=Fit all
vprefix=Physiological vname='Physiological Measurements'
wprefix=Exercises wname='Exercises';
var Weight Waist Pulse;
with Chins Situps Jumps;
weight Weight ;
title 'Middle-Aged Men in a Health Fitness Club with weights';
run;
title 'Middle-Aged Men in a Health Fitness Club with NO weights';
proc cancorr data=Fit all
vprefix=Physiological vname='Physiological Measurements'
wprefix=Exercises wname='Exercises';
var Weight Waist Pulse;
with Chins Situps Jumps;
run;
proc cancorr data=Fit all
vprefix=Physiological vname='Physiological Measurements'
wprefix=Exercises wname='Exercises';
var Weight Waist Pulse;
with Chins Situps Jumps;
freq Weight ;
title 'Middle-Aged Men in a Health Fitness Club with frequencies';
run;
