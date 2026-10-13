# UCSC_Toil_TCGA_GTEx

- Canonical DATA path: `DATA/bulkRNAseq/UCSC_Toil_TCGA_GTEx/raw/` (local SSD; Argos same path for expected_count, TPM/phenotype at `DATA/bulkRNAseq/public_bulk_DB/UCSC_Xena/`)
- Source: UCSC Xena toil hub, https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/ (Vivian et al., Nat Biotechnol 2017); public
- Files: `TcgaTargetGtex_gene_expected_count.gz` (1,289,124,490 B, md5 87b9b2c409846b412396836faedf4814, acquired 2026-10-13 for this project), `TcgaTargetGtex_rsem_gene_tpm.gz` (md5 8006d92eeb199bf4701102750f45fe65, pre-existing), `TcgaTargetGTEX_phenotype.txt.gz` (md5 25cc6f2bdbc56da2989df79a2ad86308, pre-existing), `gencode.v23.annotation.gene.probemap` (md5 c9f7e859d17d4c9102c07b98ffbd753d)
- Why: uniform TCGA + GTEx processing for Step 1 DE (C1/C2/C3). Not downloaded: GDC STAR-Counts (Argos per-project Rdata lost sample-type barcodes; see hot data list §2.2)
- Status: verified (byte size = S3 Content-Length, gzip -t OK, md5 identical local/Argos)
