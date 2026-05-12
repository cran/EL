# EL 1.4

* Converted all legacy function documentation to `roxygen2`.
* `EL.statistic()` now returns an object of class `"htest"`.
* Added `FDEL.acf()` for frequency-domain EL test of autocorrelation differences.
* Fixed axis label issues in `EL.plot()` and updated examples.
* Renamed `d` to `Delta` in `EL.statistic()` to align with theoretical notation; `d` is retained as a deprecated alias for backward compatibility.
