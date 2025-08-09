# Unit Testing Strategy for iSCORE-PDecipher
**Created:** 2025-08-08
**Principle:** STABILITY WORSHIPPED THROUGH TESTING

## 🛡️ What Are Unit Tests?
Automated tests that verify individual functions produce expected outputs for given inputs. They are our safety net for the "stability above all" principle.

## 🎯 Why Unit Tests Are NON-NEGOTIABLE

1. **Guarantee Stability**: Catch ANY deviation from expected behavior
2. **Enable Safe Optimization**: Change code fearlessly with test protection  
3. **Document Behavior**: Tests show exactly what functions should do
4. **Prevent Regression**: Previously fixed bugs stay fixed
5. **Build Confidence**: Know immediately if changes break anything

## 📋 Implementation Plan

### Phase 1: Establish Testing Framework
```r
# Install testing infrastructure
install.packages("testthat")
usethis::use_testthat()

# Create test structure
tests/
├── testthat.R
├── testthat/
│   ├── test-mast-analysis.R
│   ├── test-mixscale-analysis.R
│   ├── test-data-import.R
│   ├── test-enrichment.R
│   └── fixtures/         # Golden standard test data
│       ├── mast_expected_results.rds
│       └── mixscale_expected_results.rds
```

### Phase 2: Create Golden Standard Test Data
```r
# Capture CURRENT WORKING outputs as golden standard
# These become our immutable reference points

# Example: Save current MAST results
current_mast <- run_MAST_analysis(test_data)
saveRDS(current_mast, "tests/testthat/fixtures/mast_golden_standard.rds")

# Any future changes MUST produce identical results
```

### Phase 3: Critical Tests to Implement

#### 1. Data Import Tests
```r
test_that("MAST data import handles all mutations correctly", {
  data <- import_mast_data("test_file.rds")
  expect_equal(ncol(data), expected_cols)
  expect_true("LRRK2" %in% data$gene)
  expect_false("A15_FPD-24" %in% data$experiment) # Bug we fixed!
})

test_that("MixScale import distinguishes experiments correctly", {
  data <- import_mixscale_data("test_file.rds")
  # Test for the C12_FPD-24 vs A15_FPD-24 confusion we fixed
  expect_true(all(data$experiment %in% c("C12_FPD-23", "C12_FPD-24", "C18_FPD-23")))
})
```

#### 2. Statistical Calculation Tests
```r
test_that("MAST p-values match golden standard", {
  result <- run_MAST(test_gene, test_cluster)
  golden <- readRDS("fixtures/mast_golden_standard.rds")
  expect_equal(result$p_val_adj, golden$p_val_adj, tolerance = 1e-10)
})

test_that("MixScale weighted p-values calculated correctly", {
  result <- calculate_weighted_pval(test_data)
  expect_equal(result, expected_weighted_pval, tolerance = 1e-10)
})
```

#### 3. Enrichment Analysis Tests
```r
test_that("Enrichment returns correct pathway counts", {
  result <- run_enrichment(test_genes)
  expect_equal(nrow(result), expected_pathways)
  expect_true(all(result$p.adjust <= 0.05))
})
```

#### 4. Memory and Performance Tests
```r
test_that("Analysis completes within memory limits", {
  mem_before <- pryr::mem_used()
  result <- run_analysis(large_test_data)
  mem_after <- pryr::mem_used()
  expect_true((mem_after - mem_before) < 5e9) # Less than 5GB increase
})
```

#### 5. Edge Case Tests
```r
test_that("Handles empty clusters gracefully", {
  expect_warning(result <- analyze_cluster(empty_cluster))
  expect_true(is.na(result) || nrow(result) == 0)
})

test_that("Handles missing experiments without crashing", {
  data_missing <- test_data[test_data$experiment != "C12_FPD-24", ]
  expect_silent(result <- process_data(data_missing))
})
```

### Phase 4: Integration with Development Workflow

#### Before ANY Code Change:
```r
# 1. Run all tests - MUST PASS
devtools::test()

# 2. Make your change

# 3. Run tests again - MUST STILL PASS
devtools::test()

# 4. If any test fails - REVERT IMMEDIATELY
```

#### Continuous Testing:
```r
# Watch for changes and auto-run tests
testthat::auto_test("R/", "tests/testthat/")
```

## 🚨 Test Coverage Requirements

### Minimum Coverage Before Optimization:
- **Critical Functions**: 100% coverage
  - All MAST analysis functions
  - All MixScale analysis functions
  - Data import functions
  - Enrichment calculations
  
- **Supporting Functions**: 80% coverage
  - Visualization functions
  - Helper utilities
  - Data transformations

- **Overall Package**: 85% coverage minimum

### Measuring Coverage:
```r
library(covr)
cov <- package_coverage()
report(cov)
# Should show >= 85% coverage
```

## 📜 Testing Commandments

1. **NEVER skip failing tests** - Fix the code or fix the test
2. **NEVER change golden standards** without team consensus
3. **ALWAYS run full test suite** before commits
4. **ALWAYS add tests** for bug fixes (regression tests)
5. **ALWAYS test edge cases** - empty data, single values, NAs
6. **NEVER optimize** without test coverage
7. **ALWAYS compare** optimized vs original outputs

## 🎯 Success Criteria

- [ ] All existing functionality has tests
- [ ] Tests run in < 5 minutes
- [ ] Zero test failures on main branch
- [ ] Golden standards established and versioned
- [ ] Coverage report shows >= 85%
- [ ] CI/CD pipeline runs tests automatically

## 🔥 Example: Testing Before Optimization

```r
# BEFORE optimization
test_that("Original function works", {
  result_original <- slow_function(test_data)
  saveRDS(result_original, "golden_result.rds")
  expect_equal(result_original$value, 42)
})

# AFTER optimization
test_that("Optimized function matches original EXACTLY", {
  result_optimized <- fast_function(test_data)
  golden <- readRDS("golden_result.rds")
  expect_identical(result_optimized, golden) # MUST BE IDENTICAL!
})
```

## 📚 Resources

- [testthat documentation](https://testthat.r-lib.org/)
- [R Packages book - Testing chapter](https://r-pkgs.org/tests.html)
- [covr for coverage reports](https://github.com/r-lib/covr)

---

**Remember:** Every test we write is a guardian of stability. No optimization without verification!