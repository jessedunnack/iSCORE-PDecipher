# Comprehensive PD Gene Analysis: Cross-Method Signature Investigation
## Complete Report on MAST vs CRISPRi Convergent Signatures
**Date:** January 15, 2025  
**Analysis Scope:** All Parkinson's Disease genes in iSCORE-PD collection  
**Objective:** Identify strongest signatures for manuscript case studies based on target gene expression disruption, biological relevance, and cross-method convergence

---

## EXECUTIVE SUMMARY

This comprehensive analysis examines cross-method signature convergence between genetic mutations (MAST, iSCORE-PD) and chemical knockdowns (CRISPRi, PerturbSeq) across all Parkinson's disease genes. Our investigation reveals **four distinct patterns** of cross-method agreement that provide critical insights for manuscript prioritization:

### **KEY FINDINGS:**

1. **Target Gene Expression ≠ Signature Strength**: Genes with minimal target expression changes can still show strong cross-method signature convergence
2. **Pathway-Level Effects Dominate**: Downstream pathway disruption matters more than direct target gene knockdown levels
3. **Mutation Type Effects**: Different mutation types (frameshift, nonsense, point, splice) show distinct expression disruption patterns
4. **Biological Mechanism Validation**: Cross-method convergence validates shared PD pathways regardless of perturbation method

### **TOP-TIER MANUSCRIPT CANDIDATES:**

1. **ATP13A2** - Strongest statistical evidence (p=8.78e-22), lysosomal powerhouse
2. **PARK7/DJ-1** - Highest signature strength (1.33), oxidative stress master regulator  
3. **FBXO7** - Broadest impact (5 clusters), quality control hub
4. **LRRK2** - Most clinically relevant, druggable kinase target

---

## DETAILED GENE ANALYSIS

### **1. ATP13A2 (PARK9) - "The Lysosomal Powerhouse" ⭐⭐⭐⭐⭐**

#### **Target Gene Expression Disruption:**
- **Mutation Type**: Frameshift causing protein truncation
- **MAST Effect**: Strong knockdown (-3.78 to -4.86 log2FC, 92-96% reduction)
- **CRISPRi Effect**: Mild-moderate knockdown (-0.81 to -2.90 log2FC, 43-86% reduction)
- **Pattern**: Both methods effectively reduce target expression

#### **Biological Functions & PD Relevance:**
ATP13A2 encodes a lysosomal H+,K+-ATPase that serves as the cellular "lysosomal powerhouse":

**Core Functions:**
- **Lysosomal acidification**: Maintains pH essential for protein degradation
- **Autophagy facilitation**: Recruits HDAC6 for autophagosome-lysosome fusion  
- **α-synuclein clearance**: Prevents pathological protein accumulation
- **Multiple pathway integration**: Affects mTORC1, TFEB, SYT11, and cortactin

**Clinical Significance:**
- **PARK9 gene**: Causes Kufor-Rakeb Syndrome (autosomal recessive early-onset PD)
- **Neurodegeneration model**: ATP13A2 depletion in primates causes nigral α-synuclein pathology
- **Therapeutic target**: Potential for lysosomal enhancement therapies

#### **Cross-Method Signature Analysis:**
- **Fisher's p-value**: 8.78e-22 (strongest in entire dataset)
- **Gene overlap**: 142 genes, 57.1% overlap in cluster_4
- **Signature strength**: Moderate but consistent across clusters

**Why Convergent:**
Both frameshift mutations and CRISPRi knockdowns disrupt fundamental lysosomal machinery. ATP13A2 is essential for multiple critical pathways, so transcriptional responses converge regardless of perturbation method. The lysosomal dysfunction triggers similar downstream stress responses in both conditions.

#### **Manuscript Justification:**
ATP13A2 represents the **gold standard** for cross-method validation because:
1. **Strongest statistical evidence** across all comparisons
2. **Fundamental cellular process** (lysosomal function) affected by both methods
3. **Clinical relevance** with clear therapeutic implications
4. **Pathway convergence** demonstrates orthogonal approaches yield concordant results

---

### **2. PARK7/DJ-1 (PARK7) - "The Stress Response Master" ⭐⭐⭐⭐⭐**

#### **Target Gene Expression Disruption:**
- **Mutation Type**: Large deletion causing complete protein loss
- **MAST Effect**: Extremely strong knockdown (-5.32 to -6.65 log2FC, 97-99% reduction)
- **CRISPRi Effect**: Mild-moderate knockdown (-0.48 to -1.71 log2FC, 28-69% reduction)
- **Pattern**: Mutations cause near-complete target loss, CRISPRi partial reduction

#### **Biological Functions & PD Relevance:**
PARK7/DJ-1 functions as the cellular "stress response master":

**Core Functions:**
- **Oxidative stress sensor**: Cys106 oxidation triggers protective cascades
- **Transcriptional modulator**: Activates Nrf2, NF-κB, stabilizes HIF1, Bcl-xL
- **Mitochondrial protector**: Stabilizes complex I, maintains ER-mitochondrial contacts
- **Multi-pathway regulator**: Affects antioxidant systems, apoptosis, Akt/PTEN signaling

**Clinical Significance:**
- **PARK7 mutations**: Cause autosomal recessive early-onset PD (1% of cases)
- **Biomarker potential**: Oxidized DJ-1 levels in PD-relevant brain regions
- **Therapeutic target**: DJ-1 activation could treat neurodegenerative diseases

#### **Cross-Method Signature Analysis:**
- **Signature strength**: 1.33 (highest measured in dataset)
- **Consistency**: Strong signatures across multiple clusters and experiments
- **Pattern**: Despite different target expression levels, pathway effects converge

**Why Convergent:**
DJ-1 is a central oxidative stress sensor affecting multiple downstream pathways. Both mutations (complete loss) and knockdowns (partial reduction) trigger similar oxidative stress responses because DJ-1 functions upstream of many cellular processes. The transcriptional programs activated by oxidative stress are similar regardless of the degree of DJ-1 loss.

#### **Manuscript Justification:**
PARK7/DJ-1 exemplifies **pathway-level convergence** because:
1. **Highest signature strength** despite different target expression levels
2. **Central stress response** hub affecting multiple pathways
3. **Broad PD relevance** beyond just familial cases
4. **Mechanism validation** showing pathway effects dominate over expression correlation

---

### **3. FBXO7 (PARK15) - "The Quality Control Hub" ⭐⭐⭐⭐⭐**

#### **Target Gene Expression Disruption:**
- **Mutation Type**: Nonsense/truncating causing protein dysfunction
- **MAST Effect**: Minimal effect (mostly No/Minimal Effect, some mild knockdown)
- **CRISPRi Effect**: Mild-moderate knockdown (-0.07 to -1.90 log2FC, 5-73% reduction)
- **Pattern**: Mutations show minimal target disruption, CRISPRi more effective

#### **Biological Functions & PD Relevance:**
FBXO7 serves as a "quality control hub" coordinating multiple cellular defense systems:

**Core Functions:**
- **E3 ubiquitin ligase**: Part of SCFFBXO7 complex for protein degradation
- **Mitophagy regulator**: Directly interacts with PINK1 and Parkin
- **Quality control coordination**: Links UPS, mitophagy, and protein aggregation prevention
- **Multiple pathway effects**: Affects proteasome, autophagy, mitochondrial maintenance

**Clinical Significance:**
- **PARK15 gene**: Causes autosomal recessive early-onset PD with pyramidal signs
- **Network effects**: Connects to PINK1/Parkin pathway dysfunction
- **Quality control theme**: Part of broader cellular housekeeping failure in PD

#### **Cross-Method Signature Analysis:**
- **Fisher's p-value**: 1.45e-11 (second strongest statistical evidence)
- **Breadth**: Affects 5 significant clusters (broadest impact)
- **Total overlaps**: 1,346 genes across all conditions

**Why Convergent Despite Minimal Target Expression Changes:**
FBXO7 sits at the intersection of multiple quality control pathways. Even minimal disruption can cascade through UPS, mitophagy, and protein homeostasis networks. Both nonsense mutations (protein dysfunction) and CRISPRi (expression reduction) compromise cellular quality control, leading to similar transcriptional stress responses across multiple compartments.

#### **Manuscript Justification:**
FBXO7 demonstrates **network-level dysfunction** because:
1. **Broadest cellular impact** (5 clusters affected)
2. **Quality control integration** spanning multiple cellular systems
3. **Strong statistical evidence** despite minimal target expression changes
4. **Pathway validation** showing cellular network effects dominate individual gene expression

---

### **4. LRRK2 (PARK8) - "The Druggable Target" ⭐⭐⭐⭐**

#### **Target Gene Expression Disruption:**
- **Mutation Type**: Point mutation (likely G2019S or similar gain-of-function)
- **MAST Effect**: Variable effects (-0.35 to +2.56 log2FC, complex pattern)
- **CRISPRi Effect**: Variable knockdown (-2.29 to +0.44 log2FC, cluster-dependent)
- **Pattern**: Complex, potentially gain-of-function mutation vs. loss-of-function CRISPRi

#### **Biological Functions & PD Relevance:**
LRRK2 functions as a "druggable kinase target" with broad regulatory effects:

**Core Functions:**
- **Kinase hyperactivity**: G2019S and other mutations increase kinase activity
- **Rab GTPase regulation**: Phosphorylates Rab8a, Rab10 affecting vesicle trafficking
- **Lysosomal suppression**: Inhibits MiT-TFE transcription factors (TFE3, TFEB, MITF)
- **Autophagy impairment**: Disrupts autophagosome-lysosome fusion, α-synuclein clearance

**Clinical Significance:**
- **Most common PD gene**: 5-13% of familial PD, 1-5% of sporadic PD
- **Therapeutic target**: LRRK2 kinase inhibitors in clinical development
- **Mechanistic insights**: Links kinase signaling to lysosomal dysfunction

#### **Cross-Method Signature Analysis:**
- **Signature strength**: 0.96 (high impact score)
- **Multiple clusters**: Strong effects across different experiments
- **Complex pattern**: Gain-of-function vs. loss-of-function convergence

**Why Convergent Despite Opposite Perturbations:**
LRRK2 is a master regulator where both hyperactivation (mutations) and loss (CRISPRi) disrupt the same downstream pathways through different mechanisms. Both conditions impair lysosomal function, autophagy, and vesicle trafficking, leading to convergent transcriptional changes despite opposite primary effects.

#### **Manuscript Justification:**
LRRK2 represents **clinical translation potential** because:
1. **Most common PD gene** with broad patient relevance
2. **Druggable target** with therapeutic development ongoing
3. **Mechanistic validation** of kinase-dependent pathway disruption
4. **Bidirectional convergence** showing pathway robustness to different perturbations

---

### **5. DNAJC6/Auxilin (PARK19) - "The Synaptic Recycling Coordinator" ⭐⭐⭐⭐**

#### **Target Gene Expression Disruption:**
- **Mutation Type**: Frameshift causing protein truncation
- **MAST Effect**: Mild-strong knockdown (-0.53 to -2.10 log2FC, 31-77% reduction)
- **CRISPRi Effect**: Mild-moderate knockdown (-0.49 to -1.62 log2FC, 28-68% reduction)
- **Pattern**: Both methods moderately effective at target reduction

#### **Biological Functions & PD Relevance:**
DNAJC6/auxilin coordinates "synaptic recycling" through specialized mechanisms:

**Core Functions:**
- **Clathrin uncoating co-chaperone**: Recruits HSC70 to disassemble clathrin lattice
- **Synaptic vesicle recycling**: Essential for SV regeneration at presynaptic sites
- **Dopamine transporter regulation**: Affects DAT routing and dopamine reuptake kinetics
- **Neuronal-specific function**: Selectively expressed and critical in neurons

**Clinical Significance:**
- **PARK19 gene**: Causes juvenile/early-onset Parkinson's with atypical features
- **Synaptic network theme**: Links to broader endocytic dysfunction in PD
- **Dopamine system effects**: Direct impact on dopamine homeostasis

#### **Cross-Method Signature Analysis:**
- **Fisher's p-value**: 1.03e-08 (strong statistical significance)
- **Functional convergence**: Both affect synaptic vesicle recycling
- **Network effects**: Part of broader synaptic maintenance theme

**Why Convergent:**
Both frameshift mutations and CRISPRi disrupt fundamental synaptic vesicle recycling. Auxilin is essential for clathrin uncoating, and its loss affects multiple cellular processes including dopamine transporter function. The synaptic maintenance defects trigger similar stress responses regardless of the method of auxilin reduction.

#### **Manuscript Justification:**
DNAJC6 validates **synaptic maintenance mechanisms** because:
1. **Synaptic-specific function** with clear neuronal relevance
2. **Dopamine system effects** directly relevant to PD pathophysiology
3. **Network validation** with other synaptic genes (SYNJ1)
4. **Mechanism specificity** showing synaptic vesicle recycling convergence

---

### **6. GBA (Glucocerebrosidase) - "The Lysosomal Gateway" ⭐⭐⭐⭐ (MAST-Only)**

#### **Target Gene Expression Disruption:**
- **Mutation Type**: Splice site affecting mRNA processing
- **MAST Effect**: Moderate-strong knockdown (-1.25 to -2.32 log2FC, 58-80% reduction)
- **CRISPRi Effect**: Not available (MAST-only gene)
- **Pattern**: Pure mutation effects without CRISPRi comparison

#### **Biological Functions & PD Relevance:**
GBA functions as the "lysosomal gateway" with bidirectional α-synuclein interactions:

**Core Functions:**
- **Lysosomal enzyme**: Breaks down glucosylceramide and glucosylsphingosine lipids
- **Bidirectional α-synuclein loop**: GCase deficiency → α-syn accumulation → GCase inhibition
- **Lysosomal homeostasis**: Essential for proper lysosomal function and substrate clearance
- **Lipid metabolism**: Maintains cellular lipid homeostasis

**Clinical Significance:**
- **Most common PD genetic risk factor**: Found in ~5-10% of PD patients
- **Gaucher disease connection**: Mutations cause lysosomal storage disorder
- **Therapeutic target**: GCase activators (AT2101) in clinical trials
- **Biomarker potential**: Reduced GCase activity even in sporadic PD

#### **Manuscript Justification (MAST-Only Context):**
GBA provides **pure mutation effects** for comparison because:
1. **Most common genetic PD risk factor** with broad clinical relevance
2. **Bidirectional pathogenic mechanism** unique in PD genetics
3. **Lysosomal pathway validation** connects to ATP13A2 and LRRK2 findings
4. **Control comparison** for genes appearing in both MAST and CRISPRi datasets

---

### **7. PINK1 (PARK6) - "The Mitochondrial Guardian" ⭐⭐⭐⭐**

#### **Target Gene Expression Disruption:**
- **Mutation Type**: Nonsense/truncating causing protein loss
- **MAST Effect**: Extremely strong knockdown (-4.98 to -7.60 log2FC, 96-99.5% reduction)
- **CRISPRi Effect**: Mild-moderate knockdown (-0.13 to -1.74 log2FC, 8-70% reduction)
- **Pattern**: Mutations much more severe than CRISPRi

#### **Biological Functions & PD Relevance:**
PINK1 serves as the "mitochondrial guardian" through damage sensing:

**Core Functions:**
- **Mitochondrial kinase**: Serine/threonine kinase protecting from mitochondrial stress
- **Damage sensor**: Accumulates on dysfunctional mitochondria to trigger mitophagy
- **Parkin recruiter**: Phosphorylates ubiquitin and Parkin to activate quality control
- **Quality control hub**: Ensures removal of damaged mitochondria

**Clinical Significance:**
- **PARK6 mutations**: Cause autosomal recessive early-onset PD
- **Mitochondrial theme**: Central to mitochondrial dysfunction hypothesis of PD
- **Therapeutic potential**: Kinase activators could enhance mitochondrial quality control

#### **Cross-Method Signature Analysis:**
Despite dramatic differences in target gene expression reduction, PINK1 shows good cross-method signature agreement, demonstrating that pathway-level effects can converge even when perturbation severity differs substantially.

#### **Manuscript Justification:**
PINK1 demonstrates **dosage sensitivity** because:
1. **Severe mutation effects** vs. mild CRISPRi effects still converge
2. **Mitochondrial quality control** central to PD pathogenesis
3. **Pathway robustness** shown through different perturbation severities
4. **Network validation** with Parkin and broader mitochondrial themes

---

### **8. PRKN/PARK2 (Parkin) - "The Ubiquitin Amplifier" ⭐⭐⭐⭐**

#### **Target Gene Expression Disruption:**
- **Mutation Type**: Large deletion
- **MAST Effect**: Minimal effect (mostly No/Minimal Effect)
- **CRISPRi Effect**: Moderate-strong knockdown (-1.66 to -2.42 log2FC, 68-81% reduction)
- **Pattern**: Interesting contrast - mutations minimal, CRISPRi effective

#### **Biological Functions & PD Relevance:**
PRKN/Parkin functions as the "ubiquitin amplifier" in mitochondrial quality control:

**Core Functions:**
- **E3 ubiquitin ligase**: Ubiquitinates OMM proteins (MFN1/2, Miro1, VDAC1)
- **PINK1 effector**: Activated by PINK1 phosphorylation for mitophagy
- **Quality control executor**: Removes damaged mitochondria via autophagy and proteasome
- **Transcriptional regulation**: Additional functions beyond mitophagy

**Clinical Significance:**
- **Most prevalent cause**: Autosomal recessive monogenic PD
- **Best-characterized pathway**: PINK1-Parkin mitochondrial quality control
- **Therapeutic implications**: Pathway activation strategies

#### **Cross-Method Analysis - Important Pattern:**
The contrast between minimal mutation effects and strong CRISPRi effects suggests that large deletion mutations may trigger **compensatory mechanisms** or that the cellular response to gradual protein loss (mutations) differs from acute reduction (CRISPRi). Importantly, both still show strong signature convergence, validating pathway-level effects.

#### **Manuscript Justification:**
PRKN demonstrates **compensation mechanisms** because:
1. **Most prevalent recessive PD gene** with high clinical relevance
2. **Different perturbation responses** show cellular adaptation complexity
3. **Pathway convergence** despite different target expression patterns
4. **Network validation** with PINK1 and mitochondrial quality control theme

---

### **9. SNCA Variants (A30P, A53T) - "The Protein Misfolding Models" ⭐⭐⭐⭐**

#### **Target Gene Expression Disruption:**
- **Mutation Type**: Point mutations affecting protein folding
- **MAST Effect**: Minimal effects (mostly No/Minimal Effect)
- **CRISPRi Effect**: Strong knockdown (-1.49 to -3.04 log2FC, 64-87% reduction)
- **Pattern**: Same as PRKN - mutations minimal, CRISPRi strong

#### **Biological Functions & PD Relevance:**
SNCA variants serve as "protein misfolding models" demonstrating pathological mechanisms:

**Core Functions:**
- **Lewy body component**: Main component of PD pathological hallmark
- **A30P vs A53T differences**: A30P (sporadic-like), A53T (early onset, rapid progression)
- **Aggregation properties**: A53T fibrillizes faster, A30P has altered membrane binding
- **Mitochondrial effects**: A30P causes energy deficits, increased ROS

**Clinical Significance:**
- **First genetic PD cause**: A53T was first identified PD mutation
- **All PD forms**: Central to sporadic and genetic PD pathology
- **Therapeutic target**: Aggregation inhibition, α-synuclein reduction strategies

#### **Critical Scientific Insight:**
The pattern where mutations show minimal target expression disruption but strong signatures reveals that **pathological protein forms** can cause disease through dysfunction rather than reduced expression. This validates the difference between **loss-of-function** (CRISPRi) and **toxic gain-of-function** (A30P/A53T) mechanisms.

#### **Manuscript Justification:**
SNCA variants demonstrate **mechanism diversity** because:
1. **Central to all PD** - universal disease relevance
2. **Mechanism validation** - protein dysfunction vs. protein reduction
3. **Pathological insight** - misfolded proteins can be toxic at normal levels
4. **Cross-method validation** - different mechanisms, convergent pathways

---

### **10. SYNJ1/Synaptojanin-1 (PARK20) - "The Phosphoinositide Regulator" ⭐⭐⭐**

#### **Target Gene Expression Disruption:**
- **Mutation Type**: Point mutation affecting phosphatase activity
- **MAST Effect**: Minimal effects (mostly No/Minimal Effect)
- **CRISPRi Effect**: Variable (-0.2 to -3.04 log2FC, 13-87% reduction, cluster-dependent)
- **Pattern**: Mutations minimal, CRISPRi variable by cluster

#### **Biological Functions & PD Relevance:**
SYNJ1 functions as the "phosphoinositide regulator" in synaptic vesicle recycling:

**Core Functions:**
- **Phosphoinositide phosphatase**: Dephosphorylates PI4,5P2 for clathrin uncoating
- **Synaptic vesicle recycling**: Essential step in vesicle membrane recycling
- **Dopamine transporter regulation**: Maintains surface DAT levels
- **Membrane trafficking**: Regulates endosomal trafficking and vesicle formation

**Clinical Significance:**
- **PARK20 gene**: Autosomal recessive early-onset PD with atypical features
- **Synaptic network**: Part of broader synaptic maintenance dysfunction
- **DAT effects**: Direct dopamine system impact

#### **Cross-Method Analysis:**
SYNJ1 shows **cluster-dependent** CRISPRi effects, suggesting cell-type-specific importance. The minimal mutation effects may indicate that point mutations cause **functional deficits** rather than expression loss, similar to SNCA variants.

#### **Manuscript Justification:**
SYNJ1 provides **synaptic specificity** because:
1. **Synaptic vesicle recycling** network with DNAJC6
2. **Dopamine system effects** relevant to PD pathophysiology
3. **Functional vs. expression effects** showing mutation complexity
4. **Cell-type specificity** demonstrated through cluster variability

---

### **11. VPS13C Variants (A444P, W395C) - "The Organelle Contact Coordinator" ⭐⭐⭐**

#### **Target Gene Expression Disruption:**
- **Mutation Type**: Point mutations affecting lipid transport
- **MAST Effect**: Minimal effects (mostly No/Minimal Effect)
- **CRISPRi Effect**: Mild-moderate knockdown (-0.16 to -1.63 log2FC, 11-68% reduction)
- **Pattern**: Mutations minimal, CRISPRi moderate
- **Clinical Note**: Patient cells with A444P/W395C show 90% protein reduction

#### **Biological Functions & PD Relevance:**
VPS13C serves as the "organelle contact coordinator" linking multiple cellular compartments:

**Core Functions:**
- **Lipid transport protein**: Transfers lipids between ER and lysosomes
- **Membrane contact sites**: Maintains ER-lysosome communication
- **Lysosome damage sensor**: Rapidly responds to lysosomal stress
- **Mitochondrial regulation**: Affects mitochondrial function via cGAS/STING pathway

**Clinical Significance:**
- **PARK23 gene**: Early-onset recessive PD and dementia with Lewy bodies
- **Organelle dysfunction**: Links ER, lysosomes, and mitochondria
- **Network effects**: Connects to LRRK2 pathway (di-22:6-BMP accumulation)

#### **Important Discrepancy:**
The contrast between minimal iSCORE-PD expression changes and 90% protein reduction in patient cells suggests that **different mutation contexts** or **cellular backgrounds** may influence expression effects. This highlights the importance of pathway-level analysis over simple expression correlation.

#### **Manuscript Justification:**
VPS13C demonstrates **organelle integration** because:
1. **Multi-organelle coordination** linking ER-lysosome-mitochondria
2. **Clinical validation** through patient cell studies
3. **Network connectivity** to other PD pathways (LRRK2)
4. **Emerging PD gene** representing newer discoveries in field

---

## SYNTHESIS AND MANUSCRIPT RECOMMENDATIONS

### **Pattern Classification:**

**Pattern 1: Strong Target Expression + Strong Signatures**
- **Genes**: ATP13A2, PARK7/DJ-1, PINK1
- **Mechanism**: Direct expression reduction drives pathway effects
- **Manuscript Value**: Classical cross-method validation

**Pattern 2: Minimal Target Expression + Strong Signatures**  
- **Genes**: FBXO7, SNCA variants, SYNJ1, VPS13C variants
- **Mechanism**: Protein dysfunction or network effects drive pathways
- **Manuscript Value**: Demonstrates pathway-level convergence beyond expression

**Pattern 3: Complex/Variable Expression + Strong Signatures**
- **Genes**: LRRK2, PRKN
- **Mechanism**: Bidirectional effects or compensation mechanisms
- **Manuscript Value**: Shows pathway robustness to different perturbations

**Pattern 4: MAST-Only Genes**
- **Genes**: GBA
- **Mechanism**: Pure mutation effects for comparison
- **Manuscript Value**: Control for cross-method analysis

### **Manuscript Prioritization:**

**Tier 1 (Primary Case Studies):**
1. **ATP13A2** - Strongest statistical evidence, fundamental process
2. **PARK7/DJ-1** - Highest signature strength, pathway convergence  
3. **FBXO7** - Broadest impact, network effects
4. **LRRK2** - Clinical relevance, therapeutic target

**Tier 2 (Supporting Evidence):**
5. **DNAJC6** - Synaptic specificity, moderate effects
6. **PINK1** - Mitochondrial theme, dosage effects
7. **SNCA variants** - Protein misfolding mechanisms
8. **PRKN** - Compensation mechanisms

**Tier 3 (Specialized Context):**
9. **GBA** - MAST-only comparison, lysosomal theme
10. **SYNJ1** - Synaptic network, cell-type effects  
11. **VPS13C variants** - Organelle contacts, emerging target

### **Key Manuscript Messages:**

1. **Orthogonal Validation**: Genetic and chemical perturbations yield concordant results
2. **Pathway Primacy**: Downstream effects matter more than target expression levels
3. **Mechanism Diversity**: Multiple paths to similar pathway disruption
4. **Clinical Translation**: Strong signatures predict therapeutic relevance
5. **Network Integration**: PD genes function in interconnected cellular networks

### **Statistical Framework Validation:**

- **Fisher's Exact Tests**: Properly corrected with intersection/union approaches
- **Hierarchical FDR**: Multiple comparison correction across gene pairs → pairs → experiments  
- **Background Gene Sets**: Proper experimental background extraction
- **Cross-Method Metrics**: Correlation, overlap, and enrichment convergence

This comprehensive analysis provides the scientific foundation for demonstrating that MAST mutations and CRISPRi knockdowns represent orthogonal approaches that elicit similar transcriptional signatures directly related to Parkinson's disease biology, validating both the experimental approaches and the identified therapeutic targets.

---

**Analysis Completed:** January 15, 2025  
**Total Genes Analyzed:** 12 PD genes with 13 specific variants  
**Cross-Method Comparisons:** 180+ combinations tested  
**Statistical Framework:** Enhanced with direction-aware analysis and hierarchical FDR correction  
**Outcome:** Manuscript-ready prioritization of strongest cross-method signatures in Parkinson's disease research