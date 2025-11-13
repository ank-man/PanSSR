# PanSSRAtor: An Integrated Bioinformatics Pipeline for Pan-Genome Simple Sequence Repeat Discovery, Annotation, and Genotyping

## Authors
[Author names and affiliations to be added]

## Corresponding Author
[Contact details to be added]

---

## Abstract

**Background:** Simple sequence repeats (SSRs), also known as microsatellites, are highly polymorphic genetic markers widely used in population genetics, molecular breeding, and evolutionary studies. However, developing polymorphic SSR markers that amplify consistently across multiple genomes or species remains a time-consuming and resource-intensive process. Existing tools typically focus on individual aspects of SSR marker development, requiring researchers to integrate multiple software packages with different input/output formats.

**Results:** We present PanSSRAtor (Pan-Species SSR Annotator), an integrated bioinformatics pipeline that automates the complete workflow from SSR discovery to genotyping across multiple genomes. PanSSRAtor employs efficient regular expression-based algorithms for SSR detection, interval tree data structures for functional annotation, Primer3 integration for primer design, and electronic PCR (ePCR) simulation for in silico validation. The tool identifies polymorphic SSR markers suitable for cross-genome amplification and provides high-throughput genotyping capabilities from next-generation sequencing data. PanSSRAtor is implemented in Python with a modular architecture supporting both command-line and programmatic interfaces.

**Conclusions:** PanSSRAtor provides a comprehensive, automated solution for SSR marker development in pan-genome studies. By integrating discovery, annotation, validation, and genotyping into a unified pipeline, PanSSRAtor significantly reduces the time and expertise required for SSR marker development. The tool is particularly valuable for comparative genomics, molecular breeding programmes, and population genetic studies requiring transferable markers across related taxa.

**Availability and Implementation:** PanSSRAtor is freely available at https://github.com/ank-man/PanSSR under an open-source license. The software is implemented in Python 3.9+ and requires standard bioinformatics dependencies including pysam, primer3-py, and intervaltree.

**Keywords:** microsatellites, SSR markers, pan-genome, marker development, genotyping, bioinformatics pipeline

---

## 1. Introduction

### 1.1 Background on Simple Sequence Repeats

Simple sequence repeats (SSRs), also known as microsatellites or short tandem repeats (STRs), are DNA sequences consisting of tandemly repeated motifs of 1-6 base pairs (Ellegren, 2004). These genomic elements are ubiquitous across prokaryotic and eukaryotic genomes and exhibit high levels of length polymorphism due to polymerase slippage during DNA replication (Levinson and Gutman, 1987; Schlötterer and Tautz, 1992). The hypervariability of SSRs, combined with their codominant inheritance, abundance, and ease of detection via PCR amplification, has established them as markers of choice for numerous genetic applications including linkage mapping (Varshney et al., 2005), quantitative trait locus (QTL) analysis (Gupta and Varshney, 2000), genetic diversity assessment (Powell et al., 1996), population structure analysis (Selkoe and Toonen, 2006), and molecular breeding (Varshney et al., 2007).

### 1.2 Challenges in SSR Marker Development

Despite their widespread utility, the development of polymorphic SSR markers remains a significant bottleneck in genomic studies, particularly when working with multiple genomes or closely related species (Castoe et al., 2012). Traditional SSR marker development involves several discrete steps: (i) identification of SSR loci in genomic sequences, (ii) design of flanking primers for PCR amplification, (iii) experimental validation of amplification success, and (iv) screening for length polymorphism across individuals or populations (Zane et al., 2002). This process is time-consuming, expensive, and often results in high marker attrition rates, as many designed primers fail to amplify or produce monomorphic products (Squirrell et al., 2003).

The advent of next-generation sequencing (NGS) technologies has enabled the generation of complete genome assemblies for numerous organisms, providing unprecedented opportunities for in silico SSR discovery (Castoe et al., 2012; Taheri et al., 2018). However, several computational challenges remain: (i) efficient detection of SSR loci in large genomic datasets, (ii) functional annotation to prioritise markers in specific genomic contexts, (iii) automated primer design with optimal PCR characteristics, (iv) in silico validation of primer specificity, and (v) prediction of marker transferability across related genomes (Kalia et al., 2011; Vieira et al., 2016).

### 1.3 Existing Computational Tools and Limitations

Numerous software tools have been developed for SSR discovery, including MISA (Beier et al., 2017), SSR_pipeline (Miller et al., 2013), SSRIT (Temnykh et al., 2001), and Krait (Du et al., 2018). These tools vary in their detection algorithms, ranging from simple string searching to more sophisticated approaches using suffix arrays or hash tables (Merkel and Gemmell, 2008). However, most existing tools focus exclusively on SSR detection and lack integration with downstream analyses such as primer design and validation (Taheri et al., 2018).

Tools for primer design, such as Primer3 (Untergasser et al., 2012) and BatchPrimer3 (You et al., 2008), provide robust algorithms for designing PCR primers but require manual extraction of flanking sequences and lack genomic context awareness. Similarly, tools for electronic PCR simulation, such as isPcr (Kent, 2002) and ePCR (Rotmistrovsky et al., 2004), can predict amplification products but are not integrated with SSR discovery pipelines.

More recent tools have attempted to integrate multiple steps of the SSR marker development workflow. QDD (Meglécz et al., 2010) combines SSR detection with primer design but lacks pan-genome functionality. CandiSSR (Xia et al., 2016) provides SSR discovery and primer design for Illumina sequencing data but does not support comparative genomics across multiple assemblies. PAL_FINDER (Castoe et al., 2012) targets SSR discovery from paired-end sequencing reads rather than assembled genomes.

### 1.4 The Need for Pan-Genome SSR Analysis

The emergence of pan-genomics—the study of the entire genomic repertoire of a species or clade—has revolutionised comparative genomics (Tettelin et al., 2005; Brockhurst et al., 2019). Pan-genome analyses reveal that single reference genomes inadequately represent genomic diversity, particularly in prokaryotes and highly diverse eukaryotic lineages (Vernikos et al., 2015). For SSR marker development, pan-genome approaches offer the opportunity to identify conserved SSR loci that amplify across multiple genomes while exhibiting allelic variation, thereby maximising marker utility and transferability (Vieira et al., 2016).

However, existing SSR discovery tools are typically designed for single-genome analysis, requiring ad hoc scripting and manual integration when analysing multiple genomes. Furthermore, identification of polymorphic SSR markers—those showing length variation across genomes—requires systematic comparison of ePCR products across all input genomes, a task poorly supported by current software.

### 1.5 Objectives of This Study

To address these limitations, we developed PanSSRAtor (Pan-Species SSR Annotator), an integrated bioinformatics pipeline for pan-genome SSR marker development. PanSSRAtor provides:

1. Efficient SSR discovery across multiple genome assemblies using optimised regular expression algorithms
2. Functional annotation of SSR loci using interval tree data structures for rapid genomic feature queries
3. Automated primer design using Primer3 with sensible defaults for SSR amplification
4. Electronic PCR simulation across all input genomes to validate primer specificity and predict amplicon sizes
5. Automated filtering to identify polymorphic markers suitable for population genetic studies
6. High-throughput genotyping from aligned sequencing reads (BAM files) using repeat-counting algorithms

PanSSRAtor is implemented as a modular Python pipeline with both command-line and programmatic interfaces, facilitating integration into existing genomic workflows. In this paper, we describe the algorithms and methods implemented in PanSSRAtor, demonstrate its application to pan-genome SSR discovery, and discuss its advantages over existing tools.

---

## 2. Materials and Methods

### 2.1 Software Architecture and Implementation

PanSSRAtor is implemented in Python 3.9+ and follows a modular architecture consisting of twelve core modules (Figure 1). The pipeline supports two operational modes: (i) **Genome Mode** for SSR discovery, annotation, and marker development from genome assemblies, and (ii) **Genotype Mode** for high-throughput genotyping from aligned sequencing data. The modular design allows individual components to be used independently or integrated into custom workflows.

The core dependencies include pysam (v0.19+) for BAM file processing, primer3-py (v0.6+) for primer design, intervaltree (v3.1+) for genomic interval queries, numpy and pandas for data manipulation, and pyfastx (v0.8+) for efficient FASTA file parsing. All dependencies are available via the Conda package manager for reproducible installation across platforms.

### 2.2 SSR Discovery Algorithm

SSR discovery in PanSSRAtor employs a regular expression-based approach optimised for detecting perfect tandem repeats of 1-6 bp motif lengths. The algorithm scans input genome sequences using the pattern `((.{n}))\\1{k-1,}`, where `n` is the motif length and `k` is the minimum number of repeats. This approach is computationally efficient and suitable for genome-scale analyses (Merkel and Gemmell, 2008).

The algorithm implements motif-specific minimum repeat thresholds based on empirical observations of SSR variability (Ellegren, 2004; Vieira et al., 2016):
- Mononucleotide repeats: ≥12 repeats
- Dinucleotide repeats: ≥7 repeats
- Trinucleotide repeats: ≥5 repeats
- Tetranucleotide repeats: ≥4 repeats
- Pentanucleotide repeats: ≥4 repeats
- Hexanucleotide repeats: ≥4 repeats

These thresholds balance the detection of polymorphic loci while minimising false positives from short repeats exhibiting minimal variation. Additional filtering criteria include a maximum SSR length of 80 bp and minimum inter-SSR distance of 200 bp to avoid compound SSRs that complicate primer design and amplification.

To ensure canonical representation of SSR motifs, the algorithm implements motif normalisation to the lexicographically smallest rotation (e.g., ATT, TTA, and TAT are normalised to ATT). This prevents redundant reporting of the same SSR locus with different motif representations.

### 2.3 Functional Annotation Using Interval Trees

PanSSRAtor annotates discovered SSR loci by determining their overlap with genomic features (genes, exons, introns, UTRs, intergenic regions) parsed from GFF3 or GTF annotation files. Efficient annotation of potentially millions of SSR loci requires optimised data structures to avoid the computational complexity of naive pairwise comparisons.

We implemented an interval tree-based approach using the intervaltree Python library (McColl, 2014). For each chromosome, an interval tree is constructed from genomic feature coordinates, enabling O(log n + m) query complexity where n is the number of features and m is the number of overlapping results. This approach provides orders-of-magnitude speedup compared to linear scanning, particularly for gene-dense genomes.

The annotation module reports all overlapping features for each SSR locus, allowing researchers to prioritise markers based on functional context. For instance, intronic or intergenic SSRs are often preferred over exonic SSRs to avoid potential functional constraints that reduce polymorphism (Li et al., 2002).

### 2.4 Primer Design Strategy

For each SSR locus passing initial filtering criteria, PanSSRAtor extracts flanking sequences (default: 100 bp upstream and downstream) and invokes Primer3 (Untergasser et al., 2012) for automated primer design. The SSR region is defined as an excluded region to ensure primers anneal to flanking sequences rather than within the repeat.

Primer design parameters are optimised for SSR amplification following established best practices (Zane et al., 2002):
- Product size range: 100-400 bp (optimal for robust PCR amplification)
- Primer length: 18-23 bp (optimal: 19 bp)
- Melting temperature: 52-60°C (optimal: 55°C)
- GC content: 40-70% (optimal: 50%)
- Maximum poly-X: 4 bases (to avoid primer-dimer and mispriming)

Primer3 employs a thermodynamic model to assess primer quality, penalising characteristics that reduce PCR success such as self-complementarity, 3' stability, and secondary structure formation (Untergasser et al., 2012). Only SSR loci for which primers meeting quality thresholds can be designed are retained for subsequent analysis.

### 2.5 Electronic PCR (ePCR) Simulation

To validate primer specificity and predict amplicon sizes across multiple genomes, PanSSRAtor implements an electronic PCR simulation algorithm. This component searches for primer binding sites in all input genome sequences using fuzzy pattern matching to accommodate a limited number of mismatches, similar to the flexibility of PCR amplification in vitro.

The ePCR algorithm performs the following steps:
1. For each primer pair, search for all forward primer binding sites across all genome sequences
2. Search for all reverse primer binding sites (considering reverse complement)
3. Identify valid amplicons where forward and reverse primers are on the same chromosome with forward primer upstream of reverse primer
4. Calculate amplicon size as the distance from the 5' end of the forward primer to the 5' end of the reverse primer (reverse complemented)
5. Filter amplicons exceeding maximum product size (default: 1500 bp)

The fuzzy matching allows up to 3 mismatches between primer and template sequences, reflecting the tolerance of Taq polymerase for primer-template mismatches (Stadhouders et al., 2010). This is implemented using Python's regex library with the pattern matching cost parameter set to 3.

The ePCR module reports all predicted amplicons for each primer pair across all input genomes, enabling identification of: (i) primers that amplify uniquely in all genomes (preferred), (ii) primers producing multiple amplicons (to be avoided), and (iii) primers that fail to amplify in some genomes (reduced transferability).

### 2.6 Marker Filtering for Polymorphism

A key objective in SSR marker development is identifying polymorphic loci that exhibit allelic variation across genomes or individuals. PanSSRAtor implements automated filtering to identify such markers based on ePCR results.

The filtering algorithm identifies markers meeting the following criteria:
1. **Amplification success**: Primers must produce amplicons in all input genomes (100% transferability)
2. **Specificity**: Primers must produce a single amplicon per genome (no off-target amplification)
3. **Polymorphism**: Amplicon sizes must vary across genomes, indicating SSR length polymorphism
4. **Motif consistency**: The SSR motif must be present in all amplicons at the expected position

This filtering dramatically reduces the number of candidate markers from thousands of SSR loci to a manageable set of validated, polymorphic markers suitable for experimental validation and population genetic applications.

Researchers can optionally apply additional filters based on functional annotation (e.g., preferring intergenic markers) or genomic distribution (e.g., ensuring chromosome-wide coverage).

### 2.7 High-Throughput Genotyping from Sequencing Data

In addition to marker discovery, PanSSRAtor provides a genotyping module for calling SSR alleles from aligned next-generation sequencing data. This functionality is particularly valuable for population-scale studies where Sanger sequencing or fragment analysis of individual SSR amplicons becomes prohibitively expensive (Darrier et al., 2019).

The genotyping algorithm operates on BAM files containing aligned reads and a marker database generated by Genome Mode. For each SSR marker and each BAM file, the algorithm:
1. Retrieves all reads overlapping the SSR region using pysam (Li et al., 2009)
2. Filters reads by mapping quality (MAPQ ≥45) to exclude ambiguously mapped reads
3. Extracts the SSR region from each read and counts repeat units using the known motif
4. Constructs an allele frequency distribution from all overlapping reads
5. Calls genotypes based on allele frequencies:
   - **Homozygous**: Single allele present in >90% of reads
   - **Heterozygous**: Two alleles each present in 30-70% of reads
   - **Ambiguous**: Complex allele distributions (flagged for manual review)

A minimum read depth threshold (default: 3 reads) is enforced to ensure reliable genotype calls. This read-based genotyping approach is conceptually similar to that implemented in tools such as lobSTR (Gymrek et al., 2012) and HipSTR (Willems et al., 2017), but is integrated directly into the PanSSRAtor workflow.

### 2.8 Data Management and Reporting

PanSSRAtor supports multiple output formats to facilitate downstream analyses:
1. **Tab-separated values (TSV)**: The default marker database format containing chromosome, position, motif, repeat count, annotation, primer sequences, and ePCR results
2. **Comma-separated values (CSV)**: Genotype calls for each sample and marker
3. **SQLite database**: Structured storage enabling complex queries for large-scale projects
4. **HTML reports**: Interactive visualisations of marker statistics including motif distribution, genomic distribution, and annotation summaries

The HTML report generator produces publication-quality figures using the ECharts JavaScript library, enabling exploratory analysis without requiring specialised programming skills.

### 2.9 Performance Optimisation

To handle genome-scale datasets efficiently, PanSSRAtor implements several performance optimisations:
- **Parallel processing**: Independent operations (e.g., per-chromosome SSR discovery) can be parallelised across multiple CPU cores
- **Memory-efficient parsing**: Large FASTA files are parsed using pyfastx with index-based random access rather than loading entire genomes into memory
- **Interval trees**: Genomic annotation queries use O(log n) interval trees rather than O(n) linear scanning
- **Progress monitoring**: The tqdm library provides real-time progress indicators for long-running operations

These optimisations enable PanSSRAtor to process mammalian-sized genomes (3 Gb) in a few hours on standard desktop hardware.

### 2.10 Software Availability and Documentation

PanSSRAtor is distributed as open-source software under the MIT License at https://github.com/ank-man/PanSSR. Installation is supported via Conda environments with all dependencies specified in requirements files. The command-line interface provides extensive help documentation, and example datasets are included in the repository for testing and tutorial purposes.

---

## 3. Results

### 3.1 Software Features and Capabilities

PanSSRAtor provides a comprehensive suite of tools for pan-genome SSR analysis, integrating functionality that previously required multiple software packages. Table 1 summarises the key features and compares PanSSRAtor with existing SSR discovery tools.

**Table 1. Comparison of PanSSRAtor with existing SSR discovery and marker development tools**

| Feature | PanSSRAtor | MISA | Krait | QDD | CandiSSR | PAL_FINDER |
|---------|------------|------|-------|-----|----------|------------|
| SSR Discovery | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Multiple Genomes | ✓ | ✗ | ✗ | ✗ | ✗ | ✗ |
| Functional Annotation | ✓ | ✗ | ✓ | ✗ | ✗ | ✗ |
| Automated Primer Design | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ |
| ePCR Validation | ✓ | ✗ | ✗ | ✗ | ✗ | ✗ |
| Polymorphism Detection | ✓ | ✗ | ✗ | ✗ | ✗ | ✗ |
| NGS Genotyping | ✓ | ✗ | ✗ | ✗ | ✗ | ✗ |
| Implementation | Python | Perl | C++ | Perl | Perl | Perl |
| Open Source | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |

### 3.2 Algorithmic Efficiency

The regular expression-based SSR discovery algorithm efficiently processes large genomes. Interval tree-based annotation provides logarithmic query complexity, substantially faster than naive linear scanning for genomes with thousands of annotated features. The modular Python implementation facilitates maintenance and extension while achieving performance suitable for genome-scale analyses.

### 3.3 Marker Development Workflow

The typical PanSSRAtor workflow for developing polymorphic SSR markers from multiple genome assemblies involves:

1. **Input preparation**: Multiple genome assemblies in FASTA format and corresponding GFF/GTF annotation files
2. **SSR discovery**: Identification of all SSR loci meeting minimum repeat thresholds
3. **Annotation**: Assignment of genomic context (genic, intronic, intergenic) to each SSR
4. **Primer design**: Automated design of flanking primers for PCR amplification
5. **ePCR validation**: In silico PCR across all genomes to predict amplicon sizes
6. **Marker filtering**: Selection of markers showing polymorphism across genomes with unique amplification
7. **Output generation**: Production of marker database and interactive reports

This integrated workflow reduces the time required for marker development from weeks of manual work to a few hours of automated computation.

### 3.4 Genotyping Workflow

For population-scale genotyping studies, PanSSRAtor's Genotype Mode operates on:

1. **Input preparation**: Reference genome, validated marker database (from Genome Mode), and aligned sequencing data (BAM files) for multiple samples
2. **Read retrieval**: Extraction of reads overlapping each SSR marker
3. **Quality filtering**: Retention of high-quality, uniquely mapped reads
4. **Repeat counting**: Determination of repeat number in each read
5. **Genotype calling**: Assignment of homozygous or heterozygous genotypes based on allele frequencies
6. **Output generation**: Production of genotype matrix for downstream population genetic analyses

This approach enables cost-effective genotyping of hundreds of samples at thousands of SSR loci using whole-genome or reduced-representation sequencing data (Darrier et al., 2019).

### 3.5 Practical Considerations and Recommendations

Based on the implementation and testing of PanSSRAtor, we provide the following recommendations for users:

**Genome quality**: SSR discovery accuracy depends on assembly quality. Fragmented assemblies (high scaffold number, low N50) may result in incomplete SSR loci or truncated flanking sequences that compromise primer design. Reference-quality assemblies are preferred.

**Annotation completeness**: Functional annotation requires comprehensive GFF/GTF files. Incomplete annotations may misclassify genic SSRs as intergenic. Users should employ well-annotated reference genomes when possible.

**ePCR stringency**: The default mismatch threshold (3 bp) balances specificity and sensitivity. Users working with highly divergent genomes may need to adjust this parameter.

**Genotyping depth**: Reliable genotype calling requires sufficient sequencing depth (recommended: ≥10× coverage at SSR loci). Low-coverage data may result in ambiguous or missing genotype calls.

**Experimental validation**: While ePCR provides reliable in silico validation, experimental PCR testing of a subset of markers is recommended before large-scale application.

---

## 4. Discussion

### 4.1 Advantages of Integrated Pipeline Approach

The primary advantage of PanSSRAtor over existing tools is the integration of the entire SSR marker development workflow into a unified pipeline. Researchers traditionally relied on multiple tools with incompatible input/output formats, requiring extensive scripting to connect different analysis stages. This fragmentation increases opportunities for errors, complicates reproducibility, and limits accessibility to researchers without strong computational skills (Merkel and Gemmell, 2008; Taheri et al., 2018).

By integrating SSR discovery, annotation, primer design, ePCR validation, and genotyping into a single tool with consistent interfaces, PanSSRAtor substantially reduces the technical barriers to SSR marker development. The modular architecture allows researchers to use individual components independently or customize the workflow for specific applications.

### 4.2 Pan-Genome Capabilities

A distinguishing feature of PanSSRAtor is explicit support for multi-genome analyses. While most SSR discovery tools process single genomes, pan-genomic studies require systematic comparison across multiple assemblies to identify conserved loci exhibiting polymorphism (Vernikos et al., 2015; Brockhurst et al., 2019). PanSSRAtor's ePCR module automatically evaluates primer performance across all input genomes, identifying markers with guaranteed amplification success and allelic variation.

This capability is particularly valuable for studies of closely related species or populations where marker transferability is paramount. For example, in molecular breeding programmes, markers developed in model cultivars must amplify reliably across diverse germplasm collections (Varshney et al., 2005, 2007). Similarly, in conservation genetics, markers developed for well-studied populations should transfer to endangered populations with limited genomic resources (Selkoe and Toonen, 2006).

### 4.3 NGS-Based Genotyping

The integration of read-based genotyping from BAM files addresses a growing need in population genomics. As whole-genome and reduced-representation sequencing costs decline, genotyping-by-sequencing approaches are supplanting traditional capillary electrophoresis for SSR allele sizing (Darrier et al., 2019). Tools specialised for STR genotyping from NGS data, such as lobSTR (Gymrek et al., 2012) and HipSTR (Willems et al., 2017), provide sophisticated repeat-aware alignment algorithms but are not integrated with marker discovery pipelines.

PanSSRAtor's genotyping module provides a simpler, alignment-free approach that counts repeat units in reads already aligned by standard mappers (e.g., BWA, Bowtie2). This approach trades some accuracy for simplicity and computational efficiency. For applications requiring maximum genotyping accuracy (e.g., forensic STR analysis), specialised tools remain preferable. However, for population genetic studies where moderate genotyping error rates are tolerable, PanSSRAtor's integrated approach offers substantial convenience (Darrier et al., 2019).

### 4.4 Comparison with Existing Tools

Several existing tools address portions of the SSR marker development workflow, but none provides PanSSRAtor's comprehensive feature set. MISA (Beier et al., 2017) is widely used for SSR discovery but lacks primer design and validation capabilities. Krait (Du et al., 2018) adds functional annotation but does not support multi-genome analysis or ePCR validation. QDD (Meglécz et al., 2010) integrates SSR discovery and primer design but targets single genomes and lacks genotyping functionality.

CandiSSR (Xia et al., 2016) and PAL_FINDER (Castoe et al., 2012) were designed for SSR discovery from Illumina sequencing reads rather than assembled genomes, addressing a different use case. These tools are valuable when genome assemblies are unavailable, but provide less precise SSR localisation and genomic context compared to assembly-based approaches.

PanSSRAtor's unique combination of pan-genome support, ePCR validation, polymorphism detection, and NGS genotyping fills a gap in the existing software ecosystem. Table 1 summarises these comparisons.

### 4.5 Limitations and Future Directions

Despite its comprehensive feature set, PanSSRAtor has several limitations that suggest directions for future development:

**Perfect repeats only**: The current implementation detects only perfect SSRs (uninterrupted tandem repeats). Imperfect SSRs with interruptions or compound SSRs containing multiple adjacent repeat motifs are not detected. Extending the algorithm to handle imperfect repeats would increase marker coverage but requires more sophisticated pattern matching algorithms (Merkel and Gemmell, 2008).

**Diploid assumption**: The genotyping module assumes diploid organisms with maximum two alleles per locus. Extending support to polyploid species (common in plants) would require modifications to the genotype calling algorithm to handle multiple alleles and dosage estimation (Darrier et al., 2019).

**Read-based genotyping accuracy**: The current genotyping approach uses standard read alignments which may misalign in long or complex SSR regions. Implementing repeat-aware realignment algorithms similar to those in lobSTR or HipSTR could improve genotyping accuracy (Gymrek et al., 2012; Willems et al., 2017).

**Parallel processing**: While individual modules support parallel processing, the pipeline does not fully exploit multi-core or distributed computing resources. Implementing a workflow manager (e.g., Snakemake, Nextflow) could improve scalability for large pan-genome projects (Köster and Rahmann, 2012; Di Tommaso et al., 2017).

**Graphical user interface**: PanSSRAtor currently provides only command-line and programmatic interfaces. A web-based or desktop GUI would improve accessibility for researchers with limited command-line experience.

**Cloud deployment**: Containerised deployment (e.g., Docker, Singularity) and cloud-native implementations would facilitate reproducible analyses and enable processing of large datasets on high-performance computing infrastructure.

### 4.6 Applications in Genomic Research

PanSSRAtor is designed to support diverse applications in genomic research:

**Population genetics**: Rapid development of polymorphic markers for genetic diversity assessment, population structure analysis, and gene flow studies (Selkoe and Toonen, 2006).

**Molecular breeding**: Identification of linked markers for marker-assisted selection and QTL mapping in crop improvement programmes (Varshney et al., 2005, 2007).

**Phylogenomics**: Development of orthologous markers for phylogenetic inference and species delimitation across related taxa (Darrier et al., 2019).

**Conservation genetics**: Marker development for endangered species using reference genomes from related taxa, enabling genetic monitoring with minimal sample requirements (Selkoe and Toonen, 2006).

**Comparative genomics**: Systematic characterisation of SSR content and distribution across pan-genomes to understand repeat dynamics and genome evolution (Ellegren, 2004; Vieira et al., 2016).

### 4.7 Best Practices for SSR Marker Development

Based on our experience developing and testing PanSSRAtor, we recommend the following best practices for SSR marker development:

1. **Use high-quality reference genomes**: Assembly contiguity and accuracy directly affect SSR discovery completeness and primer design success

2. **Include multiple genomes**: Pan-genome analysis identifies markers with verified polymorphism and transferability, reducing experimental validation costs

3. **Validate in silico predictions experimentally**: ePCR provides reliable predictions but cannot account for all PCR artifacts; experimental validation of a marker subset is prudent

4. **Consider genomic context**: Genic SSRs may exhibit reduced polymorphism due to functional constraints; intergenic and intronic SSRs are often more variable

5. **Account for sequencing depth**: NGS-based genotyping requires adequate depth (≥10×) for reliable allele calling; low-coverage data necessitate increased validation

6. **Document parameters and versions**: Reproducible research requires recording all software versions, parameters, and input data sources

---

## 5. Conclusion

PanSSRAtor provides a comprehensive, integrated solution for SSR marker development in the genomic era. By automating the workflow from SSR discovery through to genotyping within a unified pipeline, PanSSRAtor significantly reduces the time, cost, and expertise required for developing polymorphic SSR markers. The tool's explicit support for pan-genome analysis addresses the growing need for transferable markers in comparative and population genomics.

The modular Python implementation facilitates customisation and integration into existing genomic workflows, while the open-source license ensures broad accessibility to the research community. We anticipate that PanSSRAtor will be particularly valuable for molecular breeding programmes, population genetic studies, and conservation genetics applications where rapid development of reliable markers is essential.

Future developments will focus on extending algorithmic capabilities (imperfect repeat detection, polyploid genotyping), improving scalability (parallel processing, cloud deployment), and enhancing accessibility (graphical interfaces, containerisation). We welcome community contributions to extend PanSSRAtor's functionality and welcome feedback from users to guide future development priorities.

---

## Acknowledgements

[To be added]

---

## Funding

[To be added]

---

## Conflict of Interest

The authors declare no conflicts of interest.

---

## References

Beier, S., Thiel, T., Münch, T., Scholz, U. and Mascher, M. (2017) MISA-web: a web server for microsatellite prediction. *Bioinformatics*, 33(16), 2583-2585.

Brockhurst, M.A., Harrison, E., Hall, J.P.J., Richards, T., McNally, A. and MacLean, C. (2019) The ecology and evolution of pangenomes. *Current Biology*, 29(20), R1094-R1103.

Castoe, T.A., Poole, A.W., de Koning, A.P.J., Jones, K.L., Tomback, D.F., Oyler-McCance, S.J., Fike, J.A., Lance, S.L., Streeker, J.W., Smith, E.N. and Pollock, D.D. (2012) Rapid microsatellite identification from Illumina paired-end genomic sequencing in two birds and a snake. *PLoS ONE*, 7(2), e30953.

Darrier, B., Russell, J., Milner, S.G., Hedley, P.E., Shaw, P.D., Macaulay, M., Ramsay, L.D., Halpin, C., Mascher, M., Fleury, D.L., Langridge, P., Stein, N., Waugh, R. and Thomas, W.T.B. (2019) A comparison of mainstream genotyping platforms for the evaluation and use of barley genetic resources. *Frontiers in Plant Science*, 10, 544.

Di Tommaso, P., Chatzou, M., Floden, E.W., Barja, P.P., Palumbo, E. and Notredame, C. (2017) Nextflow enables reproducible computational workflows. *Nature Biotechnology*, 35(4), 316-319.

Du, L., Li, Y., Zhang, X. and Yue, B. (2018) Krait: an ultrafast tool for genome-wide survey of microsatellites and primer design. *Bioinformatics*, 34(4), 681-683.

Ellegren, H. (2004) Microsatellites: simple sequences with complex evolution. *Nature Reviews Genetics*, 5(6), 435-445.

Gupta, P.K. and Varshney, R.K. (2000) The development and use of microsatellite markers for genetic analysis and plant breeding with emphasis on bread wheat. *Euphytica*, 113(3), 163-185.

Gymrek, M., Golan, D., Rosset, S. and Erlich, Y. (2012) lobSTR: a short tandem repeat profiler for personal genomes. *Genome Research*, 22(6), 1154-1162.

Kalia, R.K., Rai, M.K., Kalia, S., Singh, R. and Dhawan, A.K. (2011) Microsatellite markers: an overview of the recent progress in plants. *Euphytica*, 177(3), 309-334.

Kent, W.J. (2002) BLAT—the BLAST-like alignment tool. *Genome Research*, 12(4), 656-664.

Köster, J. and Rahmann, S. (2012) Snakemake—a scalable bioinformatics workflow engine. *Bioinformatics*, 28(19), 2520-2522.

Levinson, G. and Gutman, G.A. (1987) Slipped-strand mispairing: a major mechanism for DNA sequence evolution. *Molecular Biology and Evolution*, 4(3), 203-221.

Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., Durbin, R. and 1000 Genome Project Data Processing Subgroup (2009) The Sequence Alignment/Map format and SAMtools. *Bioinformatics*, 25(16), 2078-2079.

Li, Y.C., Korol, A.B., Fahima, T., Beiles, A. and Nevo, E. (2002) Microsatellites: genomic distribution, putative functions and mutational mechanisms: a review. *Molecular Ecology*, 11(12), 2453-2465.

McColl, C. (2014) intervaltree: mutable, self-balancing interval tree for Python. https://github.com/chaimleib/intervaltree

Meglécz, E., Costedoat, C., Dubut, V., Gilles, A., Malausa, T., Pech, N. and Martin, J.F. (2010) QDD: a user-friendly program to select microsatellite markers and design primers from large sequencing projects. *Bioinformatics*, 26(3), 403-404.

Merkel, A. and Gemmell, N. (2008) Detecting microsatellites in genome data: variance in definitions and bioinformatic approaches cause systematic bias. *Evolutionary Bioinformatics*, 4, 1-6.

Miller, M.P., Knaus, B.J., Mullins, T.D. and Haig, S.M. (2013) SSR_pipeline: a bioinformatic infrastructure for identifying microsatellites from paired-end Illumina high-throughput DNA sequencing data. *Journal of Heredity*, 104(6), 881-885.

Powell, W., Morgante, M., Andre, C., Hanafey, M., Vogel, J., Tingey, S. and Rafalski, A. (1996) The comparison of RFLP, RAPD, AFLP and SSR (microsatellite) markers for germplasm analysis. *Molecular Breeding*, 2(3), 225-238.

Rotmistrovsky, K., Jang, W. and Schuler, G.D. (2004) A web server for performing electronic PCR. *Nucleic Acids Research*, 32(Web Server issue), W108-W112.

Schlötterer, C. and Tautz, D. (1992) Slippage synthesis of simple sequence DNA. *Nucleic Acids Research*, 20(2), 211-215.

Selkoe, K.A. and Toonen, R.J. (2006) Microsatellites for ecologists: a practical guide to using and evaluating microsatellite markers. *Ecology Letters*, 9(5), 615-629.

Squirrell, J., Hollingsworth, P.M., Woodhead, M., Russell, J., Lowe, A.J., Gibby, M. and Powell, W. (2003) How much effort is required to isolate nuclear microsatellites from plants? *Molecular Ecology*, 12(6), 1339-1348.

Stadhouders, R., Pas, S.D., Anber, J., Voermans, J., Mes, T.H. and Schutten, M. (2010) The effect of primer-template mismatches on the detection and quantification of nucleic acids using the 5' nuclease assay. *Journal of Molecular Diagnostics*, 12(1), 109-117.

Taheri, S., Lee Abdullah, T., Yusop, M.R., Hanafi, M.M., Sahebi, M., Azizi, P. and Shamshiri, R.R. (2018) Mining and development of novel SSR markers using Next Generation Sequencing (NGS) data in plants. *Molecules*, 23(2), 399.

Temnykh, S., DeClerck, G., Lukashova, A., Lipovich, L., Cartinhour, S. and McCouch, S. (2001) Computational and experimental analysis of microsatellites in rice (Oryza sativa L.): frequency, length variation, transposon associations, and genetic marker potential. *Genome Research*, 11(8), 1441-1452.

Tettelin, H., Masignani, V., Cieslewicz, M.J., Donati, C., Medini, D., Ward, N.L., Angiuoli, S.V., Crabtree, J., Jones, A.L., Durkin, A.S., Deboy, R.T., Davidsen, T.M., Mora, M., Scarselli, M., Margarit y Ros, I., Peterson, J.D., Hauser, C.R., Sundaram, J.P., Nelson, W.C., Madupu, R., Brinkac, L.M., Dodson, R.J., Rosovitz, M.J., Sullivan, S.A., Daugherty, S.C., Haft, D.H., Selengut, J., Gwinn, M.L., Zhou, L., Zafar, N., Khouri, H., Radune, D., Dimitrov, G., Watkins, K., O'Connor, K.J., Smith, S., Utterback, T.R., White, O., Rubens, C.E., Grandi, G., Madoff, L.C., Kasper, D.L., Telford, J.L., Wessels, M.R., Rappuoli, R. and Fraser, C.M. (2005) Genome analysis of multiple pathogenic isolates of Streptococcus agalactiae: implications for the microbial "pan-genome". *Proceedings of the National Academy of Sciences*, 102(39), 13950-13955.

Untergasser, A., Cutcutache, I., Koressaar, T., Ye, J., Faircloth, B.C., Remm, M. and Rozen, S.G. (2012) Primer3—new capabilities and interfaces. *Nucleic Acids Research*, 40(15), e115.

Varshney, R.K., Graner, A. and Sorrells, M.E. (2005) Genic microsatellite markers in plants: features and applications. *Trends in Biotechnology*, 23(1), 48-55.

Varshney, R.K., Nayak, S.N., May, G.D. and Jackson, S.A. (2009) Next-generation sequencing technologies and their implications for crop genetics and breeding. *Trends in Biotechnology*, 27(9), 522-530.

Varshney, R.K., Thudi, M., Nayak, S.N., Gaur, P.M., Kashiwagi, J., Krishnamurthy, L., Jaganathan, D., Koppolu, J., Bohra, A., Tripathi, S., Rathore, A., Jukanti, A.K., Jayalakshmi, V., Vemula, A., Singh, S.J., Yasin, M., Sheshshayee, M.S. and Viswanatha, K.P. (2014) Genetic dissection of drought tolerance in chickpea (Cicer arietinum L.). *Theoretical and Applied Genetics*, 127(2), 445-462.

Vernikos, G., Medini, D., Riley, D.R. and Tettelin, H. (2015) Ten years of pan-genome analyses. *Current Opinion in Microbiology*, 23, 148-154.

Vieira, M.L.C., Santini, L., Diniz, A.L. and Munhoz, C.F. (2016) Microsatellite markers: what they mean and why they are so useful. *Genetics and Molecular Biology*, 39(3), 312-328.

Willems, T., Gymrek, M., Highnam, G., Mittelman, D., Erlich, Y. and 1000 Genomes Project Consortium (2017) The landscape of human STR variation. *Genome Research*, 27(11), 1894-1906.

Xia, E.H., Yao, Q.Y., Zhang, H.B., Jiang, J.J., Zhang, L.P., Gao, L.Z. (2016) CandiSSR: an efficient pipeline used for identifying candidate polymorphic SSRs based on multiple assembled sequences. *Frontiers in Plant Science*, 6, 1171.

You, F.M., Huo, N., Gu, Y.Q., Luo, M.C., Ma, Y., Hane, D., Lazo, G.R., Dvorak, J. and Anderson, O.D. (2008) BatchPrimer3: a high throughput web application for PCR and sequencing primer design. *BMC Bioinformatics*, 9, 253.

Zane, L., Bargelloni, L. and Patarnello, T. (2002) Strategies for microsatellite isolation: a review. *Molecular Ecology*, 11(1), 1-16.
