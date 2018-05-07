This folder contains tools for predicting open reading frames and non-sense mediated decay transcripts. 


Scripts:
tama_orf_seeker.py
tama_orf_blastp_parser.py
tama_cds_regions_bed_add.py
tama_bed_extract_cds.py


Pipeline:

1. (transcript bed file) -> tama_orf_seeker.py -> (orf fasta)
2. (orf fasta) -> blastp uniref/uniprot -> (blastp bls file)
3. (blastp bls file) -> tama_orf_blastp_parser.py -> (top hits file)
4. (top hits file) + (transcript bed file) -> tama_cds_regions_bed_add.py -> (bed file with cds)

5. (bed file with cds) -> tama_bed_extract_cds.py -> (bed file representing only cds regions)


