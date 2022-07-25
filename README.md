# Viral Primal Design
Niema's scripts/data for designing viral primers (e.g. SARS-CoV-2 and Monkeypox).

# Step 1: Data Acquisition
For SARS-CoV-2, I will be using the [consensus sequences](https://github.com/andersen-lab/HCoV-19-Genomics/tree/master/consensus_sequences) from the [Andersen lab](https://andersen-lab.com/)'s [HCoV-19-Genomics GitHub repo](https://github.com/andersen-lab/HCoV-19-Genomics). There are a *lot* of consensus sequences, each in a separate FASTA file, so `git clone` would be slow. Instead, we can download them and merge them into a single FASTA as follows:

```bash
wget -O data/sarscov2/consensus_sequences_YYYY-MM-DD.fasta.gz "https://github.com/andersen-lab/HCoV-19-Genomics/releases/latest/download/consensus_sequences.fasta.gz"
```

# Step 2: Multiple Sequence Alignment
I will be using [ViralMSA](https://github.com/niemasd/ViralMSA) ([Moshiri *et al*., 2021](https://doi.org/10.1093/bioinformatics/btaa743)) wrapping around [Minimap2](https://github.com/lh3/minimap2) ([Li *et al*., 2018](https://doi.org/10.1093/bioinformatics/bty191)) to perform reference-guided Multiple Sequence Alignment. ViralMSA doesn't (currently) directly support gzipped FASTA files, so I'll temporary decompress the `fasta.gz` file, and the output files are large, so I'll compress them before uploading to GitHub:

```bash
gunzip -k data/consensus_sequences_*.fasta.gz
ViralMSA.py -s data/sarscov2/consensus_sequences_*.fasta -r SARS-CoV-2 -e niemamoshiri@gmail.com -o data/sarscov2/viralmsa_out --omit_ref
pigz -9 -p 8 data/sarscov2/viralmsa_out/*.aln data/sarscov2/viralmsa_out/*.sam
rm data/consensus_sequences_*.fasta.gz
```

# Step 3: Conserved Regions
I wrote a [script (`conserved.py`)](scripts/conserved.py) to compute the maximum base frequency and Shannon Entropy at every position of the reference genome. It takes a Multiple Sequence Alignment (one-line FASTA) via standard input and prints a TSV to standard output:

```bash
zcat data/sarscov2/viralmsa_out/*.aln.gz | ./scripts/conserved.py > data/sarscov2/conserved_2022-07-22.tsv
```
