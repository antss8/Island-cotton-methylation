This directory contains all scripts for the three omics analyses. All of them are job scripts run on an HPC cluster (LSF scheduler), i.e. scripts with `#BSUB` headers, or helper scripts.

## Directory structure

| Directory | Contents | Scripts |
| --- | --- | --- |
| `Methy/` | WGBS (whole-genome bisulfite sequencing): alignment, methylation levels, DMR / DMG | 15 |
| `RNA-seq/` | Transcriptome: QC, alignment, quantification, differential expression | 10 |
| `smallrna/` | Small RNA: Rfam depletion, alignment, collapsing, hc-siRNA loci and abundance | 8 |

## Environment and dependencies

- Scheduler: LSF, submitted with `bsub < script.lsf`; queues are `normal`, `q2680v2` and `smp`
- Reference genome: `Gbarbadense_genome_HAU_v2.0` (sequence `Gbarbadense_genome_HAU_v2.0.fasta`), annotation `Gbarbadense_gene_model_ensembl.gtf`
- Main software: `fastp`, `fastqc`, `bismark` (+ `bowtie2`), `deduplicate_bismark`, `bismark_methylation_extractor`, `coverage2cytosine`, `ViewBS`, `bedtools`, `samtools`, `hisat2`, `stringtie`, `trimmomatic`, `bowtie`, `sRNAminer`, `R`, `Python`
- `RNA-seq/6.1.getFPKM.py` and `RNA-seq/6.2.getTPM.py` are Python 2 scripts (they use `print` statements) and must be run with `python2`; the remaining Python scripts are compatible with Python 3
- `Methy/13.get_DMRpvalue_and_FDR.R` requires R; usage: `Rscript get_DMRpvalue_and_FDR.R <input_file> <output_file>`

## 1. Methy: WGBS methylation pipeline (`Methy/`)

Run the scripts in numerical order.

| No. | Script | Function | Main input → output |
| --- | --- | --- | --- |
| 1 | `1.fastp.lsf` | QC and adapter trimming of raw sequencing data (`--length_required 50 --qualified_quality_phred 30`) | raw `fq.gz` → paired-end clean data in `clean1/` |
| 2 | `2.combine_reads.lsf` | Merge the clean data of multiple lanes from the same library, separately for read 1 and read 2 | multi-lane `fastq.gz` → `*_1_clean1.fastq.gz`, `*_2_clean2.fastq.gz` |
| 3 | `3.bismark_genome.lsf` | Build the Bismark genome index with bowtie2 | genome directory → `mapping2/` |
| 4 | `4.compare.lsf` | Bismark alignment (`--non_directional -N 1 -L 30`, PE) | clean data → BAM in `comparison3/` |
| 5 | `5.redup.lsf` | Deduplication + removal of Scaffold sequences + sorting (`deduplicate_bismark`, `samtools`) | BAM → `rededuplicate2/` (including `sorted/`) |
| 6 | `6.extractor.lsf` | Methylation extraction (`--CX_context --cytosine_report`, with `coverage2cytosine` generating the CX report); also batch-gzips the results | deduplicated BAM → `Extractor2/*/*.CX.txt.CX_report.txt` |
| 7 | `7.dbinom.py` | Binomial test to select high-confidence mC sites: `sc >= 3` and `binom.pmf(x, n, 0.005) <= 1e-4` | `ViewBS1/Vpre/*_View1.tab` → `*_dbionom2.CX.txt.CX_report.txt` |
| 8 | `8.filter.py` | Simple coverage filtering (keep sites with `sc >= 3`) | `*.CX.txt.CX_report.txt` → `*_filiter.CX.txt.CX_report.txt` |
| 9 | `9.getmC.lsf` | Sort the binomial test results of the 16 samples by chromosome and position (`DbinomSort/`) | `Intersection2/dbinom/` → `*_Sorted2.txt` |
| 10 | `10.getCG-CHG-CHHmC.lsf` | Divide the all-site counts by the binomial test site counts to compute relative methylation levels in the three contexts CG/CHG/CHH | `All/` + `dbinom/` → `Intersection2/{CG,CHG,CHH}/*_calculate*.txt` |
| 11 | `11.get_CG-CHG-CHHtc.lsf` | Intersect sites between the two replicates; split by context and count sites within 1 Mb windows (`*_SCount_all_*.txt`), also outputting `bed3` | `Intersection3/All/` → `Intersection3/All/<sample>/FSallsort{CG,CHG,CHH}/` |
| 12 | `12.getCG-CHG-CHHDMRpvalue_pre.lsf` | Use `bedtools coverage` to count mC / tC within 200 bp sliding windows (requiring coverage ≥ 10) and `paste` them into the Fisher test input table; the windows are those shared by both samples | `bed3` + `All_win_200.txt` → `3.DMR/PerTXT/.../*_mc_*.txt`, `*_tc_*.txt` → `*_Fisherper.txt`, `*_Fisherlast.txt` |
| 13 | `13.get_DMRpvalue_and_FDR.R` | Chi-square test on each 2×2 contingency table to obtain a p value, then `p.adjust(method = "fdr")` correction; outputs the windows with FDR < 0.05 | `*_Fisherlast.txt` → `*_FDR_result.txt` |
| 14 | `14.get_DMR_result.lsf` | Batch-calls the R script of step 13, filters windows with p < 0.05, sorts by p value and outputs the final DMR results plus `bed3` (three branches: CG / CHG / CHH) | `*_Fisherlast.txt` → `ReallyDMR/`, `FDR/`, `sorted_last/bed3/` |
| 15 | `15.get_DMG.lsf` | Use `bedtools intersect -wa -wb` to intersect the DMRs with the gene annotation and deduplicate by gene to obtain the differentially methylated genes (DMGs) | `*_DMR_all.txt` + `sorted_genome_feature.bed6` → `4.DMG/DMRAnno/{CG,CHG,CHH}/*_GeneID_unique.txt` |

Main line of analysis: `raw data → fastp → Bismark alignment → deduplication/cleaning → methylation extraction` → (binomial test filtering by `7.dbinom.py`) → single-base methylation levels → intersection of the two replicates → 1 Mb window counting → 200 bp sliding-window counting → chi-square test + FDR → DMRs → gene annotation to obtain DMGs.

---

## 2. RNA-seq: transcriptome pipeline (`RNA-seq/`)

| No. | Script | Function |
| --- | --- | --- |
| 1 | `1.get_fastqc.lsf` | FastQC quality reports on the clean paired-end data |
| 2 | `2.get_cleandata.lsf` | `fastp` QC and adapter trimming, outputting clean data plus html/json reports |
| 3 | `3.get_index.lsf` | `hisat2-build` to build the genome index |
| 4 | `4.get_mapping.lsf` | `hisat2` alignment + `samtools` conversion to BAM, removal of Scaffold sequences and sorting, producing `*_clean.sorted.bam` |
| 5 | `5.getStrongtie_Result.lsf` | `stringtie -e -B -G` quantification, outputting per-sample transcript GTF and gene expression tables |
| 6.1 | `6.1.getFPKM.py` | Compute the **FPKM** matrix from the StringTie output (genes / transcripts, written as CSV) |
| 6.2 | `6.2.getTPM.py` | The **TPM** version of the above |
| 6.3 | `6.3get_txt.py` | Convert the expression matrix CSV into a tab-separated TXT |
| 8 | `8.getDEG_HN1.lsf` | `Trinity run_DE_analysis.pl --method DESeq2` for H vs N differential expression analysis; requires `samplegroup.txt` and `contrast.txt` |
| 9 | `9.getDEG_NN.lsf` | Same as above, for the other comparison (N vs N) |
| 10 | `10.get_cleanDEGs.sh` | Remove genes shared between the two differential expression results to obtain the specific DEG lists (`5HN*DEG.txt`, `9HN*DEG.txt`) |

Note: `6.1` / `6.2` are adapted from the official StringTie `prepDE.py`, with the FPKM-related logic replaced by TPM (or vice versa); the arguments and usage are identical (`-i` input directory or GTF list, `-g` / `-t` the output names of the gene / transcript matrices).

---

## 3. small RNA pipeline (`smallrna/`)

| No. | Script | Function |
| --- | --- | --- |
| 1 | `1.get_cleandata.lsf` | `trimmomatic` in single-end mode to remove adapters and trim low-quality bases (`MINLEN:18`) |
| 2 | `2.getrfamDATAbase.lsf` | Download all Rfam 14.2 `RF*.fa.gz` files and decompress them (the curl command is commented out by default) |
| 3 | `3.getRfam.lsf` | Merge all Rfam fasta files into `database/Rfamdatabase.fa` |
| 4 | `4.get_noRfam_reads.lsf` | Align against Rfam with `bowtie -v 0 --un` and keep the reads that do **not** map to Rfam for downstream analysis |
| 5 | `5.get_sRNAminer_mapping.lsf` | `sRNAminer bowtie_to_sRNAminer` + `Reads_mapping`: convert the bowtie alignment results into sRNAminer format and align them to the genome |
| 6 | `6.get_collapsing.lsf` | `sRNAminer Seq_collapsing` to collapse (deduplicate) the clean reads |
| 7 | `7.get_24nt_hcsiRNA_locu.lsf` | `phasiRNA_identification_no_align` to identify 24 nt phasiRNAs, and `hc-siRNA_locus_identification` to identify hc-siRNA loci |
| 8 | `8.getAbundance.lsf` | `sRNAlocus_abundance` to quantify the hc-siRNA locus abundance of each sample (`*_siRNA_abundance.xls`) |

---

## Usage notes

1. Run the steps in numerical order; submit to LSF with `bsub < script_name.lsf`.
2. On the first run, replace all absolute paths inside the scripts (data directory, software paths, genome directory).
3. Sample lists are mostly hard-coded; add or remove samples according to the actual sequencing batches, and note the commented-out lines for historical batches.
4. Steps `12`, `13` and `14` are chained through the fixed `3.DMR/PerTXT` directory, and step `14` depends on the intermediate results produced by step `13`; they must be run in order and the directories must stay consistent.
5. Resource requirements are declared in the `#BSUB` header of each script (cores `-n`, queue `-q`, memory `-M`); adjust them according to your cluster.
