  *Office
  cd "C:\Users\raulcruz\OneDrive - Indiana University\Documents\students\carlos\stata"
  insheet using fit.csv, clear
  cls
  canon (weight waist pulse) (chins situps jumps) [weight=weight]
  csdcanon "svyset  [pweight= weight],"
  /*
  canon (chins situps) (weight waist pulse)  [weight=weight]
  csdcanon "svyset  [pweight= weight],"
  */