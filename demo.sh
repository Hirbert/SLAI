docker run --rm -v `pwd`/test_data:/data slai  python3 SLAI.py \
	--bamfile /data/sample_data/P265700-WatchMaker.demo.bam  \
	--fasta /data/database/ucsc.hg19.fasta \
	--outdir /data/output \
	--snv /data/sample_data/P265700-WatchMaker.snv.chr1.txt \
	--flank 300 \
	--cpu 8 \
	--outfile P265700-WatchMaker.chr1
