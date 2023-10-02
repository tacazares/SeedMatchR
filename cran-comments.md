## Resubmission

### 10/02/2023 - Addressing Build Issues from CRAN Checks

> We were notified of build errors from the CRAN site, but were not able to 
address them in time. https://cran.r-project.org/web/packages/SeedMatchR/index.html

These tests are passing through GitHub actions now: https://github.com/tacazares/SeedMatchR/actions/runs/6381425186

### 06/19/2023 - Addressing Initial CRAN Recommendations

> If there are references describing the methods in your package, please
> add these in the description field of your DESCRIPTION file in the form
> authors (year) <doi:...>
> authors (year) <arXiv:...>
> authors (year, ISBN:...)
> or if those are not available: <https:...>
> with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for
> auto-linking. (If you want to add a title as well please put it in
> quotes: "Title")

* There are no references to add.

> You write information messages to the console that cannot be easily
> suppressed. -> R/deseq_fc_ecdf.R; R/get_feature_seqs.R; R/get_seed.R;
> R/load_species_anno_db.R
> It is more R like to generate objects that can be used to extract the
> information a user is interested in, and then print() that object.
> Instead of print() rather use message()/warning() or
> if(verbose)print(..) (or maybe stop()) if you really have to write text
> to the console. (except for print, summary, interactive functions)

* The `message()` function has replaced `print()` or removed when not necessary. 

> Please ensure that your functions do not write by default or in your
> examples/vignettes/tests in the user's home filespace (including the
> package directory and getwd()). This is not allowed by CRAN policies.
> Please omit any default path in writing functions. In your
> examples/vignettes/tests you can write to tempdir(). -> R/get_example_data.R

* The functions now write to a directory that is identified using `tempdir()`. 

> Please do not modifiy the .GlobalEnv. This is not allowed by the CRAN
> policies. -> R/load_example_data.R

* Variables are now returned instead of adding them to `.GlobalEnv`.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

