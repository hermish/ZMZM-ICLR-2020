# Changelog

GRANSO: GRadient-based Algorithm for Non-Smooth Optimization


## Version 1.6 --- 2018-09-06

Description: minor fixes and improvements

### Added

- GRANSO now prints more informative summary on termination
- GRANSO can now be set to print only every kth iteration
- GRANSO now checks that all function values and gradients at x0 are
  correctly dimensioned, finite, and real valued
- GRANSO now asserts all input arguments are the right type (n a positive
  integer, options is a struct)
- better error message when obj_fn is not a function handle
- soln now contains unscaled versions for each of soln.final, soln.best,
  and soln.most_feasible when prescaling is applied

### Fixed

- constrained optimization: soln.best_to_tol is now soln.best, as it should
  be per the documentation
- newline is now printed at start, as it should have been

### Changed

- default value of .opt_tol has been changed from 1e-6 to 1e-8
- ascii printing is forced on Windows (due to MATLAB bug with extended
  ASCII chars)

### Maintenance

- renamed nameVersionCopyrightMsg.m to copyrightNotice.m
- updated URTM routines to latest versions (tablePrinter.m,
  isFiniteValued.m, optionValidator.m)

### Documentation

- changed all instances of "Matlab" to "MATLAB" (per guidelines from The
  MathWorks)
- updated URLs in APGL copyright notices to link directly to AGPL webpage
- corrected inconsistent version numbers in the main non-private routines
- fixed typo in help doc gransoOptionsAdvanced
- moved to markdown format for CHANGELOG


## Version 1.5.2 --- 2017-07-11

Description: minor bug fix

### Fixed

- fixed opts.halt_log_fn custom halt functionality


## Version 1.5.1 --- 2017-04-17

Description: minor bug fix

### Fixed

- fixed the unable to enable regularization bug


## Version: 1.5 --- 2017-02-09

Description: new/enhanced functionality added

### Added

- new limited-memory mode (LBFGS updating) for large-scale usage
- significant reduction of memory footprint in full-memory mode
- new adaptive damping of (L)BFGS updates
- added tutorial demos
- the failure rate of quadprog is now calculated and returned to the user
- GRANSO will print a warning message if quadprog appears unreliable
- user can now disable positive definite check on H0
- new "debug mode" user option
- soln.dnorm has been corrected to be soln.stat_value
- added code.ini providing machine-readable meta data about GRANSO

### Fixed

- fixed minor bug when using opts.rel_tol
- improved/corrected some printing, documentation, and code comments


## Version: 1.0 --- 2017-01-10
Description: initial public release
