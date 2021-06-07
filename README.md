# DrosCellID

Masking of reference genome was performed with ncbi blast-2.2.26:

blastall -i te_consensus.fasta -d dmel-all-chromosome-r6.30.fasta -p blastn -a 10 -e 1e-100 -F "m L" -U T -K 20000 -b 20000 -m 8 | perl maskRef.pl dmel-all-chromosome-r6.30.fasta - > dmel-all-chromosome-r6.30.masked.fasta
