docker run --rm \
  -v "/structure-dot":/data \
  quay.io/biocontainers/gepard:2.1.0--hdfd78af_0 \
  gepardcmd \
    -seq1 /data/data/sample_1.fasta \
    -seq2 /data/data/sample_2.fasta \
    -outfile /data/output/dotplot.png \
    -word 10 \
    -lower 50

#  sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
gepardcmd \
  -seq "/structure-dot/data/sample_1.fasta //structure-dot/data/sample_2.fasta" \
  -matrix /structure-dot/data/edna.mat \
  -outfile /structure-dot/output/dotplot.png \
  -word 10 \
  -lower 50
