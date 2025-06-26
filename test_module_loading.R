# Test Module Loading
# Check if signature nomination module can be loaded without errors

cat("=== Testing Module Loading ===\n")

# Test if we can source the module
tryCatch({
  source("inst/shiny/modules/mod_signature_nomination.R")
  cat("✓ Signature nomination module loaded successfully\n")
  
  # Check if functions exist
  if (exists("mod_signature_nomination_ui")) {
    cat("✓ UI function exists\n")
  } else {
    cat("✗ UI function missing\n")
  }
  
  if (exists("mod_signature_nomination_server")) {
    cat("✓ Server function exists\n")
  } else {
    cat("✗ Server function missing\n")
  }
  
}, error = function(e) {
  cat("✗ Error loading module:", e$message, "\n")
})

# Test if we can create the UI without errors
tryCatch({
  test_ui <- mod_signature_nomination_ui("test")
  cat("✓ UI creation successful\n")
}, error = function(e) {
  cat("✗ Error creating UI:", e$message, "\n")
})

cat("=== Module Loading Test Complete ===\n")