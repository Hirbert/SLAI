docker run --rm -v `pwd`/test_data:/data slai  python3 SLAI.py \
	--bamfile /data/sample_data/sample20220912-A-A00168-R-20220908-C93_I4.sorted.rmdup.realign.bam  \
	--fasta /data/database/ucsc.hg19.fasta \
	--outdir /data/output \
	--snv /data/sample_data/sample20220912-A-A00168-R-20220908-C93_I4.final.fp_filter.txt \
	--flank 300 \
	--cpu 8 \
	--outfile sample20220912-A-A00168-R-20220908-C93_I4
