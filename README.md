Muon Signal Window
==================

It is Muon Signal window with ME0.
First, using this command, all SW and efficiency codes are clone.
```
	git clone git@github.com:ParticleChef/MuonSW.git
```

+ Region of Interest (ROI) is MakeSWMuon.C file. ME0-pixel and pixel-pixel matching SW code is in ME0\_medians and SA\_medians directory.
+ The address of sample can be managed in MakeSWMuon.h file.
+ In 'scatter' directory, you can make SW/ROI plots, scatter plots and fitting.

'offi' branch -> Study using official samples.

__run.sh and run.jds are running code for KISTI server.__

