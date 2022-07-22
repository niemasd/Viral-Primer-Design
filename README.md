# SARS-CoV-2-Primer-Design
Niema's scripts/data for designing SARS-CoV-2 primers

# Step 1: Data Acquisition
I will be using the [consensus sequences](https://github.com/andersen-lab/HCoV-19-Genomics/tree/master/consensus_sequences) from the [Andersen lab](https://andersen-lab.com/)'s [HCoV-19-Genomics GitHub repo](https://github.com/andersen-lab/HCoV-19-Genomics). There are a *lot* of consensus sequences, each in a separate FASTA file, so `git clone` would be slow. Instead, we can download them and merge them into a single FASTA as follows:

```bash
wget -O data/consensus_sequences_YYYY-MM-DD.fasta.gz "https://github.com/andersen-lab/HCoV-19-Genomics/releases/latest/download/consensus_sequences.fasta.gz"
```

# Step 2: Multiple Sequence Alignment
I will be using [ViralMSA](https://github.com/niemasd/ViralMSA) ([Moshiri *et al*., 2021](https://doi.org/10.1093/bioinformatics/btaa743)) to perform reference-guided Multiple Sequence Alignment.
