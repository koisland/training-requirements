# RNA-Seq
* Quantify expression of coding sequences.
* Exons usually. Intron/promoter regions ignored.

## Data Structure 
* One observation
	* Gene ID
	* Sample Name
	* Seq Read Counts 

## Workflow
1. Data transformation
2. Data normalization
3. Differential expression
4. Gene selection
5. Clustering

## Diff. Expr. Gene Selection Thresholds
* `pThr`
	* Adjusted p-value thresholds
* `logFCThr`
	* Log2 fold-change
* `baseMeanThr` 
	* Minimum average expression level in counts
* `cpmThr`
	* Copy-per-million threshold

## DESeq Dataset (`DESeqDataSet()`)
* `se`
	* Structured dataset of assay (values), cols (samples), and rows (genes)
* `design`
	* DESeq Design Matrix
		* Provides args for DESeq algo.
		* Examples:
			* `~ cell + dex`
				* Is there a difference between a pair of cell types (`cell`) and treatment (`dex`)
			* `~ dex`
				* Is there a difference between treatments

## DESEq (`DESeq`) 
* Output DESeq analysis
* Uses [Benjamini-Hochberg](https://www.statisticshowto.com/benjamini-hochberg-procedure/) to reduce FP rate to get `pAdj`

## Some issues with the final report
* Column header is not aligned and is missing geneID
* Gene symbols and descriptions aren't there.
* Normalized data is not shown.
* Variance stabilization normalization?

To fix this, we use the libraries, "org.Hs.eg.db" and "AnnotationDbi", for human gene annotation.

Pick "ENSEMBL", "ENTREZID", "SYMBOL" and "GENENAME" columns.
* DESeq data uses Ensembl ID
* Might be interesting to select GO column and perform gene ontology.