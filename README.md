Testing Shambhala harmonizer.

The Shambhala method (Borisov, 2019) depends upon two datasets, the auxiliary one (P0), and the definitive, (Q).

The code for Shambhala is downloaded for the repository mentioned in (Borisov, 2019).

The mix of 1) SEQC gene expression profiles marked as A, B, C, D. A: Stratagene all but brain. B: Ambion brain. C: 75% A + 25% B. D: 75%B + 25% A.

Profiled at GPL11154 (Illumina HiSeq 2000, marked as Illumina_NGS), GPL10558 (Illumina_chip), GPL17930 (Affymetrix_HUG) and GPL16043 (Affymetrix_PRV).

and 2) the content of the Legacy GDC archive (Glioblastoma multiforme, Ovarian Serous Cystadenocarcinoma, Uterine Corpus Endometrial Carcinoma, Rectum Adenocarcinoma, Colon Adenocarcinoma, Lung Squamous Cell Carcinoma, Acute Myeloid Leukemia) profiled at Affymetrix U133 Plus 2 (Affymetrix_571), Illumina HiSeq 2000 (GPL11154, Illumina_NGS) и AgilentG4502A. 

were harmonized by the following methods:

1. QN (Q.csv)
2. DESeq2 (D.csv).
3. Shambhala, P0 = hybrid of Illumina chip GPL10558 (50%) + Agilent (GPL4133+GPL1708) (50%) and Q = GTEx Illumina GPL11154 (00.csv).
4. Shambhala, P0 = GPL10558 (50%) + (GPL4133+GPL1708) (50%), Q = GTEx Affymetrix GPL16977 (01.csv).
5. Shambhala, P0 = CustomArray, Q = GTEx Illumina GPL11154  (10.csv).
6. Shambhala, P0 = CustomArray, Q = GTEx Affymetrix GPL16977 (11.csv).
7. Shambhala, P0 = Affymetrix GPL571, Q = GTEx Illumina GPL11154 (20.csv).
8. Shambhala, P0 = Affymetrix GPL571, Q = GTEx Affymetrix GPL16977 (21.csv).
9. Shambhala, P0 = Illumina chip GP10558, Q = GTEx Illumina GPL11154 (30.csv).
10. Shambhala, P0 = Illumina chip GP10558, Q = GTEx Affymetrix GPL16977 (31.csv).
11. Shambhala, P0 = Illumina NGS GP11154, Q = GTEx Illumina GPL11154 (40.csv).
12. Shambhala, P0 = Illumina NGS GP11154, Q = GTEx Affymetrix GPL16977 (41.csv).
13. Shambhala, P0 = Agilent (GPL4133+GPL1708), Q = GTEx Ilumina GPL11154 (50.csv).
14. Shambhala, P0 = Agilent (GPL4133+GPL1708), Q = GTEx Affymetrix GPL16977 (51.csv).
15. Shambhala, P0 = hybrid of Illumina chip GPL10558 (75%) + Agilent (GPL4133+GPL1708) (25%), Q = GTEx Illumina GPL11154 (60.csv).
16. Shambhala, P0 = hybrid Illumina chip GPL10558 (75%) + Agilent (GPL4133+GPL1708) (25%), Q = GTEx Affymetrix GPL16977 (61.csv).

These data are deposited at https://drive.google.com/drive/u/0/folders/1mPaFjFGu99FBKqFZFxsIgkaEgDpMFDGS (please feel free to ask permission for downloading; totally > 8 Gb).

The code for Watermelon multisection (WM) metric (Zolotovskaya, 2020) of the dendrogram and multiscale bootstrap for the dendrogram, are also provided.

The META.csv file contains the type of metadata for sample type and experimental platform.

A good harmonization technique should group profiles according to their sample type rather than experimental platform.

Consequently, the WM metric should be higher for sample-wize rather than platform-wise division.

Note that the WM metric is entropy-based, and entropy is additive, so that the WM-metric should be normalized by the number of categories (classes) in the dendrogram (sample types vs experimental platforms, respectively).

Note also that perhaps we should group all solid cancer except leukemias as "solid tumors" (contrary to leukemias) and (A+C) as "mostly A" vs (B+D) as "mostly B".

In the resulting *.png files, we can see the examples of dendrograms, as well as the distribution (violin plots) of the ratio of sample-wise vs platform-wise WM metrics.

Note also that bot WM-metric calculation and dendrogram bootstrapping is time-consuming. Therefore, it is preferred to analyze *sampled* dendrograms, which have only about ten profiles of selected combination (sample type + platform type).

Additionally, we should try another metrics for sample uniformity suggested in (Aliev, 2018), which are based on Student t-test. For these metrics, we should also compare sample-wise vs. platform-wise values.




