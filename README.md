# Nano Sniffer
An ONT Amplicon-Based CNV detection tool using a negative control and intra-sample coverage normalization

1. Purpose

This document describes a  Python-based workflow designed to detect copy number variations (CNVs) in Oxford Nanopore Technologies (ONT) amplicon sequencing data. The method is optimized for small targeted regions, such as a single TP53 or BRCA amplicon split into multiple sub-targets, where standard CNV callers are often unsuitable due to the limited number of genomic intervals and the strong influence of PCR/amplicon coverage bias.

The workflow compares each sample against a negative-control and applies an additional intra-sample normalization step across all BED-defined targets. 

The final output is an Excel report containing both high-confidence CNV calls and target-level coverage metrics.
