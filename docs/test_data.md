# Test data

The `test` profile (`-profile test`) expects two things to exist under
`assets/test/`:

```
assets/test/
├── samplesheet.test.csv          # one or two tiny samples
└── refs/
    └── phix.fa.gz                # any small reference (e.g. PhiX174, ~5.4 kb)
```

These are intentionally **not** committed — fetch what you need before running:

```bash
mkdir -p assets/test/refs

# PhiX reference (5,386 bp - perfect smoke-test size)
curl -L https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001422.1&rettype=fasta \
    | gzip > assets/test/refs/phix.fa.gz

# A few hundred HiFi reads of your choosing, then:
cat > assets/test/samplesheet.test.csv <<'EOF'
sample,fastq
test01,assets/test/reads/test01.hifi.fastq.gz
EOF
```

Then run:

```bash
nextflow run main.nf -profile test,docker
```

The `test` profile pins all processes to 2 CPU / 4 GB and sets
`errorStrategy = 'terminate'` so problems surface immediately.
