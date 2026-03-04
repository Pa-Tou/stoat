Scripts that might be useful to prepare tests.

## [unit_test_expected_values.R](unit_test_expected_values.R)

We used this test to set the expected values in the unit tests of the statistical tests. 
Basically, we used R functions and packages to run the Chi2 or regression tests (linear, logistic, and Firth logistic) and extracted the p-value that we should expect.

## [check_test_files_changes.R](check_test_files_changes.R)

We have some "*system*" tests where we run STOAT on larger files and compare the results.
However, when we change something major in the code, it might change the pvalues a bit, for example. 
This script helped double-check that the changes in pvalues were minor, or zoom in to the cases where it changed enough to warrant some investigations.
