Scatter Plots for ROI
======================

This code makes plot that has scatter, SW and fitting.
Region, nth\_SW and should be changed when you draw SW scatter plots.
For Signal window (ME0-pix, pix-pix), eta\_or\_phi should be changed.

```cpp
int nth_SW = 1;
int Region = 1; // 1 or 2
int eta_or_phi = 1; // eta = 1, phi = 2

```


Run using this.
```
root
.L plot.C+
plot a
a.Loop()
```
A png file created.
