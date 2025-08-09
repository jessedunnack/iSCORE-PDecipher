# iSCORE-PDecipher Code Style and Conventions

## Documentation Standards
- **Roxygen2**: All functions documented with `#'` comments
- **Export Tags**: `@export` for public functions (currently not synced with NAMESPACE)
- **Parameters**: `@param name Description` format
- **Return Values**: `@return Description` format
- **Examples**: Limited examples in documentation

## Naming Conventions
- **Functions**: `snake_case` (e.g., `run_mast_analysis`, `launch_iscore_app`)
- **Variables**: `snake_case` 
- **Files**: `snake_case.R` (e.g., `MAST_analysis.R`, `signature_data_helpers.R`)
- **Exports**: Descriptive function names for user-facing functions

## Code Organization
- **One main function per file**: Generally followed
- **Helper functions**: Often in same file as main function
- **Modular structure**: Separate files for different analysis types
- **Shiny modules**: Organized in `inst/shiny/modules/`

## Error Handling
- **Package checks**: `requireNamespace()` with helpful error messages
- **Input validation**: Basic validation in core functions
- **Warnings**: Used for optional package dependencies

## Dependencies Management
- **Imports**: Listed in DESCRIPTION, loaded with `::`
- **Suggests**: Optional packages with graceful degradation
- **Version constraints**: Specified for critical packages

## Current Issues Needing Attention
1. **CRITICAL**: No formal testing conventions established
2. **HIGH**: NAMESPACE out of sync with @export tags  
3. **MEDIUM**: Inconsistent error handling patterns
4. **MEDIUM**: Limited input validation in complex functions

## Best Practices to Implement
- Formal testthat 3e testing structure
- Comprehensive input validation
- Golden standard test data for scientific functions  
- Consistent error messaging
- Regular NAMESPACE synchronization checks