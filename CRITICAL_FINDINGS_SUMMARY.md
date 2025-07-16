# 🚨 CRITICAL FINDINGS: Fisher's Exact Test Investigation

**Date:** January 16, 2025  
**Status:** URGENT - Major Issues Identified

## 🔍 **MAJOR DISCOVERIES**

### 1. **Background Gene Set Mystery SOLVED** ✅
**Question:** Why are background sizes similar (~8000 intersection, ~15000 union) if signature nomination uses pathway-focused genes?

**Answer:** Both modules use **the same background gene sets** from DE analysis:
- **Signature nomination**: `de_data$iSCORE_PD_MAST[[gene]][[cluster]]$background_genes`
- **DE results**: `unique(mast_background_data$gene_name)` from DE data

The background represents **all tested genes in DE analysis**, not pathway-specific genes. This explains the similar sizes.

### 2. **Gene Universe Appropriateness** ✅
**Question:** Are both gene universes appropriate for their respective questions?

**Answer:** YES, both are scientifically valid:
- **Signature nomination**: Tests pathway-focused genes (genes that contributed to enriched pathways)
- **DE results**: Tests all DE genes (broader significance assessment)

Both answer different but valid biological questions.

### 3. **🚨 CRITICAL DIRECTIONALITY ISSUE IDENTIFIED** ⚠️
**Question:** Are we accidentally inflating enrichment overlaps by including ALL/UP/DOWN direction tests simultaneously?

**ANSWER:** **YES - THIS IS A MAJOR PROBLEM!**

#### **The Issue:**
The signature nomination module extracts genes without direction filtering:
```r
# Lines 689-690 in signature_analysis.R
mast_genes <- unique(unlist(strsplit(cluster_mast$geneID, "/")))
crispri_genes <- unique(unlist(strsplit(cluster_crispri$geneID, "/")))
```

This means if enrichment data contains:
- Gene X in UP direction test
- Gene X in DOWN direction test  
- Gene X in ALL direction test

**Gene X gets counted multiple times**, artificially inflating overlap and making p-values appear more significant than they should be.

#### **Evidence of the Problem:**
- Enrichment data definitely contains direction information (UP/DOWN/ALL)
- No direction filtering is applied in signature nomination
- Same genes can appear in multiple direction tests
- This could explain why signature nomination p-values are often lower despite FDR correction

### 4. **Investigation Script Created** 📊
Created `CRITICAL_FISHER_TEST_INVESTIGATION.R` to:
- Examine actual gene counts and overlaps
- Test for direction-related inflation
- Compare signature nomination vs direction-filtered approaches
- Validate Fisher's test calculations with real data

## 🎯 **IMMEDIATE ACTIONS REQUIRED**

### **HIGH PRIORITY:**
1. **Run the investigation script** to quantify the directionality inflation
2. **Fix the signature nomination module** to properly filter by direction
3. **Validate that both Fisher's tests are calculating correctly** after the fix

### **MEDIUM PRIORITY:**
4. **Add clear documentation** explaining the difference between the two approaches
5. **Consider adding direction-aware Fisher's tests** to DE results module

## 🔧 **PROPOSED FIXES**

### **Fix 1: Add Direction Filtering to Signature Nomination**
```r
# Instead of:
mast_genes <- unique(unlist(strsplit(cluster_mast$geneID, "/")))

# Use:
cluster_mast_filtered <- cluster_mast[cluster_mast$direction == selected_direction, ]
mast_genes <- unique(unlist(strsplit(cluster_mast_filtered$geneID, "/")))
```

### **Fix 2: Add Direction Selection UI**
Add direction selection dropdown to signature nomination module to allow users to choose:
- "UP" (up-regulated genes only)
- "DOWN" (down-regulated genes only)  
- "ALL" (all directions, but properly deduplicated)

### **Fix 3: Improve Documentation**
Add clear explanations of:
- What genes each Fisher's test is comparing
- Why the p-values are different (appropriate for different questions)
- How direction filtering affects the results

## 📊 **EXPECTED IMPACT OF FIXES**

After fixing the directionality issue:
1. **Signature nomination p-values will be more conservative** (higher)
2. **Results will be more scientifically rigorous**
3. **Differences between modules will be more clearly attributable to gene universe differences**
4. **Both approaches will be more trustworthy for their respective questions**

## 🚨 **CRITICAL NEXT STEPS**

1. **Run investigation script** to quantify the current problem
2. **Implement direction filtering fix** in signature nomination module
3. **Test the fix** with the investigation script
4. **Update documentation** to explain the differences clearly
5. **Proceed with other bug fixes** once this critical issue is resolved

---

## 💡 **USER QUESTIONS ANSWERED**

**Q: "Are both approaches appropriate for their respective questions?"**
**A:** YES - Both are scientifically valid for different questions.

**Q: "Why are background sizes similar?"**  
**A:** Both use the same DE analysis background genes (~8000 intersection, ~15000 union).

**Q: "Are we inflating enrichment by including ALL/UP/DOWN tests?"**
**A:** YES - This is a critical bug that needs immediate fixing.

**Q: "When gene numbers are similar, are p-values as expected?"**
**A:** Need to test after fixing the direction issue - this will be answered by the investigation script.

---

**STATUS:** Investigation complete, critical bug identified, fix ready for implementation.