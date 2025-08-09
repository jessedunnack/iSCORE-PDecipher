# Task Completion Checklist for iSCORE-PDecipher

## Before Starting Any Task
- [ ] Activate correct conda environment: `conda activate seuratv4`
- [ ] Check git status: `git status`
- [ ] Load development package: `devtools::load_all()`
- [ ] Understand task scope and dependencies

## After Coding Changes
- [ ] **CRITICAL**: Run any existing tests: `devtools::test()` (once implemented)
- [ ] Update documentation: `devtools::document()`
- [ ] Check package: `devtools::check()` (at minimum warnings/notes)
- [ ] Test core functionality manually if no automated tests exist
- [ ] Verify NAMESPACE is properly updated

## For New Functions
- [ ] Add comprehensive roxygen2 documentation
- [ ] Include `@export` tag if function should be public
- [ ] **CRITICAL**: Write corresponding tests (currently missing infrastructure)
- [ ] Add examples in documentation
- [ ] Update NAMESPACE: `devtools::document()`

## For Analysis Functions (HIGH RISK)
- [ ] **MANDATORY**: Create golden standard test data
- [ ] **MANDATORY**: Test with known inputs/outputs
- [ ] Validate against existing results
- [ ] Test edge cases and error conditions
- [ ] Document expected data formats

## For Shiny/UI Changes  
- [ ] Test interactively in browser
- [ ] Check mobile/responsive behavior
- [ ] Validate data flow between modules
- [ ] Test error handling in UI

## Before Committing
- [ ] **CRITICAL**: All tests pass (once implemented)
- [ ] No new warnings in `devtools::check()`
- [ ] Documentation is updated
- [ ] Git commit message is descriptive
- [ ] Consider impact on package users

## Current Critical Missing Items
- [ ] **URGENT**: Set up formal testing infrastructure with `usethis::use_testthat()`
- [ ] **URGENT**: Fix NAMESPACE sync issues with @export functions  
- [ ] **HIGH**: Create tests for `run_mast_analysis` and `run_mixscale_analysis`
- [ ] **HIGH**: Establish golden standard test data
- [ ] **MEDIUM**: Add comprehensive input validation