# 🚨 CRITICAL INVESTIGATION: P-Value Discrepancy Between Individual Experiments and App Display

**Date:** January 17, 2025  
**Status:** URGENT - Requires immediate attention  
**Impact:** Major implications for manuscript conclusions

## The Problem

There is a massive discrepancy between:
1. **Individual experiment p-values**: SYNJ1 cluster_7 = **6.47 × 10⁻⁴²** (C12_FPD-24 only)
2. **App combined p-values**: SYNJ1 cluster_0 = **1.70 × 10⁻³** (all experiments combined)

This is a **39 orders of magnitude** difference!

## Critical Findings

### Strongest Individual Signals (Not Shown in App)
1. **SYNJ1 cluster_7 (C12_FPD-24)**: p = 6.47 × 10⁻⁴²  ⚠️ HIDDEN
2. **ATP13A2 cluster_6 (A15_FPD-24)**: p = 1.62 × 10⁻³⁵  ⚠️ HIDDEN
3. **FBXO7 cluster_3 (C12_FPD-24)**: p = 4.60 × 10⁻³⁵  ⚠️ HIDDEN

### What App Shows Instead (Combined)
1. **PARK7 cluster_0**: p = 6.10 × 10⁻⁴¹
2. **FBXO7 cluster_1**: p = 6.10 × 10⁻⁴⁰
3. **PARK7 cluster_1**: p = 5.62 × 10⁻²⁸

## Why This Matters

### Scientific Impact
- **False negatives**: Extremely significant experiment-specific effects are being masked
- **Wrong clusters highlighted**: App shows cluster_0/1 as strongest, but cluster_7 has the real signal
- **Biological misinterpretation**: Cell type-specific effects (e.g., synaptic cluster_7) are lost

### Statistical Issues
1. **Simpson's Paradox**: Strong effect in one experiment becomes non-significant when combined
2. **Heterogeneity ignored**: Different experiments may target different cell populations
3. **Power loss**: Combining experiments with different effect sizes reduces overall significance

## Specific Examples of Concern

### SYNJ1 - The Hidden Super-Signal
- **Individual**: Cluster_7 + C12_FPD-24 = p < 10⁻⁴⁰ (synaptic dysfunction)
- **Combined**: Not even in top clusters, diluted to p = 0.0017 in cluster_0
- **Impact**: Missing critical synaptic-specific dysfunction signal

### ATP13A2 - Wrong Cluster Focus
- **Individual**: Cluster_6 shows p = 1.62 × 10⁻³⁵ (lysosomal specialization)
- **Combined**: App emphasizes cluster_1 instead
- **Impact**: Misidentifying the functionally relevant cell population

## Urgent Actions Needed

### 1. Immediate Verification
```r
# Check if experiment-specific effects are being properly preserved
# Verify the combination methodology
# Test for heterogeneity across experiments
```

### 2. App Modifications to Consider
- Add "Individual Experiments" view option
- Show heterogeneity statistics (I² or similar)
- Flag when combined p-value masks strong individual signals
- Allow filtering by specific experiments

### 3. Statistical Framework Review
- Is combining experiments always appropriate?
- Should we use meta-analysis methods instead?
- Need formal heterogeneity testing before combining

## Potential Solutions

### Option 1: Dual Reporting
- Show both combined AND best individual p-values
- Flag discrepancies > 10 orders of magnitude

### Option 2: Experiment-Aware Analysis
- Stratify by experiment in the app
- Use meta-analysis framework with heterogeneity testing
- Weight experiments by quality/cell count

### Option 3: Highlight Maximum Signal
- Always report the strongest signal (individual or combined)
- Annotate which experiment drives significance

## Questions to Answer

1. **Why are we combining experiments?** 
   - Different experiments may have different perturbation efficiencies
   - Cell populations may vary between experiments

2. **Is the current combination method valid?**
   - Simple union of DE genes may not be appropriate
   - Should use proper meta-analysis framework

3. **What are we potentially missing?**
   - Experiment-specific biology
   - Cell type-specific effects
   - Temporal dynamics

## Recommendation

**DO NOT RELY SOLELY ON APP P-VALUES FOR MANUSCRIPT**

Until this is resolved:
1. Report individual experiment results for key findings
2. Note when combined analysis masks strong signals
3. Perform heterogeneity testing before combining
4. Consider experiment-specific effects as primary evidence

## Files for Investigation
- `R/signature_analysis.R` - How experiments are combined
- `inst/shiny/modules/mod_signature_nomination.R` - What gets displayed
- Individual experiment results need separate analysis

**This is a critical issue that could affect the main conclusions of the manuscript!**