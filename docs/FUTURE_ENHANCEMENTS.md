# Future Enhancements for iSCORE-PDecipher

## Data-Driven PD Biology Focus (Priority Enhancement)

### Current Implementation
- Uses predefined PD-relevant pathway terms
- Manual curation of biological categories
- Focuses on known PD mechanisms (mitochondrial, autophagy, dopamine, etc.)

### Desired Enhancement: Unbiased Discovery Approach
Instead of manually defining which pathways to focus on, implement data-driven discovery:

#### 1. **Term Frequency Analysis**
- **Most Abundant Terms**: Terms appearing most frequently across all conditions
- **Consistent Terms**: Terms appearing across multiple gene/cluster combinations
- **Cross-Method Terms**: Terms appearing in both MAST and CRISPRi results

#### 2. **Impact-Based Prioritization** 
- **Strongest Statistical Significance**: Terms with lowest p-values
- **Highest Effect Sizes**: Terms with strongest enrichment scores
- **Most Robust Terms**: Terms with consistent significance across replicates

#### 3. **Frequency-Based Metrics**
- **Breadth Score**: Number of conditions where term appears
- **Depth Score**: Strength of enrichment when term appears
- **Consistency Score**: Variance in significance across appearances

#### 4. **Interactive Discovery Interface**
- **Dynamic Term Rankings**: Real-time sorting by different metrics
- **Threshold Controls**: Adjustable significance and frequency cutoffs
- **Visual Indicators**: Term strength, frequency, and breadth in single view
- **Export Capabilities**: Top terms for manuscript reporting

### Implementation Approach
```r
# Example functions to implement:
calculate_term_frequency_metrics(enrichment_data)
rank_terms_by_impact(enrichment_data, metric = "combined_score")
identify_cross_method_terms(mast_terms, crispri_terms)
generate_unbiased_term_report(enrichment_data, top_n = 50)
```

### Benefits
- **Hypothesis Generation**: Discover unexpected pathway patterns
- **Reduced Bias**: Not limited to known PD mechanisms
- **Manuscript Ready**: Data-driven term prioritization for publication
- **Reproducible**: Systematic ranking approach vs manual selection

### Priority
**High** - This aligns with unbiased discovery goals and manuscript needs.

---

## Other Planned Enhancements

### Multi-Dataset Comparison
- Cross-study validation capabilities
- Batch effect analysis and correction

### Advanced Statistics
- Multiple testing correction across all comparisons
- Meta-analysis capabilities across conditions

### Enhanced Visualizations
- Network analysis of term relationships
- Temporal analysis if time series data available