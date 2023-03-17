  *Office
  cd "C:\Users\raulcruz\OneDrive - Indiana University\Documents\students\carlos\stata"
  insheet using fit.csv, clear
  cls
  canon (weight waist pulse) (chins situps jumps) [weight=weight]
  * run postestimation_cmdv6.do
  csdcanon "svyset  [pweight= weight],"