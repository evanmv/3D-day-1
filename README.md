# 3D-day-1
**1. Login to the Saga server:**
- Show here

**2. Setting up today's working directory**
```bash
git clone https://github.com/MBV-INF4410/3D-day-1.git
cd 3D-day-1
```

**3. Loading and setting up HiC-Pro**
```bash
module purge
module load HiC-Pro/2.11.4-foss-2019a-Python-2.7.15
```
**4. Loading and setting up HiC-Pro** 
HiC-Pro will be used to process the Hi-C data, including mapping the reads and aggregation of the contact frequencies (takes a few minutes).
```bash
curl -O https://raw.githubusercontent.com/nservant/HiC-Pro/master/config-hicpro.txt
curl -O https://raw.githubusercontent.com/nservant/HiC-Pro/master/annotation/chrom_hg19.sizes
curl -O https://raw.githubusercontent.com/nservant/HiC-Pro/master/annotation/HindIII_resfrag_hg19.bed
tar -zxvf hg19_chr18/* -C hg19_chr18
```
**5. Adapting the `config-hicpro.txt` file**

To prepare for running HiC-Pro, we will need to change two lines in the `config-hicpro.txt` file. Use a text-editor (like emacs, vim, nano, (or TextEdit [MacOS]) to:
- Add the bowtie path to line nr. 39: `BOWTIE2_IDX_PATH =` -> `BOWTIE2_IDX_PATH = [fullpath]/hg19_chr18/` where `[fullpath]` is the full path to your current working directory, i.e. the `INC-tutorial` directory (use the `pwd` command if you are uncertain about the working dir.) 
- Change the reference genome on 47: `REFERENCE_GENOME =` -> `REFERENCE_GENOME = hg19_chr18`
- Change line nr. 89: `BIN_SIZE = 20000 40000 150000 500000 1000000` -> `BIN_SIZE = 50000 1000000`

**6. Run HiC-Pro (Takes ~10 minutes)**
```bash
HiC-Pro --input fastq --output hicpro_results --conf config-hicpro.txt
```

**7. Setup the folder structure for the HiC contacts**
```bash
mkdir -p hic/bedpe/intra
mkdir -p hic/matrix
```

**8. Convert Hi-C to BEDPE and matrix format***

The output from HiC-Pro needs to be converte to [BEDPE](https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format) in order to be processed further by Chrom3D, and to a matrix format in order to be compatible with the Armatus TAD caller.
```bash
awk 'NR==FNR { map[$4] = $1"\t"$2"\t"$3; next } { print $0,map[$1],map[$2] }' hicpro_results/hic_results/matrix/chr18/raw/50000/chr18_50000_abs.bed hicpro_results/hic_results/matrix/chr18/raw/50000/chr18_50000.matrix  | awk '$4==$7' | awk '{print $4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$3}' > hic/bedpe/intra/chr18

for chr in hic/bedpe/intra/*
do
chrname=$(basename $chr)
cut -f 2,5,7 $chr > hic/matrix/$chrname
done
```
```diff
!If you are stuck at this point, you can copy nesessary files to proceed with the remaining steps by:
mkdir -p hic/bedpe/intra/
mkdir -p hic/matrix/
cp backup/9/bedpe/chr18 hic/bedpe/intra/
cp backup/9/matrix/chr18 hic/matrix/
```

**9. Running Armatus to call TADs**
```bash
mkdir hic/tads
armatus-linux-x64 -r 50000 -c chr18 -S -i hic/matrix/chr18 -g .6 -o hic/tads/chr18
```
Note: if you get error messages here, try to see if you have installed armatus such that you can execute it as `armatus` instead of `armatus-linux-x64`.

The `-r 50000` sets the bin-size to 50000 bp,  `-c chr18` specifies that only chromosome 18 should be considdered, `-S` specifies that sparse matrix format (3 column text file) is used, `-i hic/matrix/chr18` provides the input data (in matrix format), `-g .6` is the gamma-max parameter indicating the highest resolution to generate domains (often is set based on trial and error), `-o hic/tads/chr18` gives the output for the TADs.

**10. Since Armatus output is end-inclusive, convert to BED by adding 1 to end position**
```bash
awk '{printf("%s\t%i\t%i\n",$1,$2,$3+1)}' hic/tads/chr18.consensus.txt > hic/tads/chr18.consensus.bed
```
