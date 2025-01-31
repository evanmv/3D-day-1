# Mapping and processing Hi-C data and calling Topologically Associated Domains (TADs)
In this exercise, you will again log in to the Fox HPC system to perform analyses. You will be:
- Mapping and processing Hi-C sequencing reads (fastq format) from Adipose Stem Cells (using HiC-Pro)
- Using containers (with apptainer)
- Generating (raw and normalized) Hi-C contact maps from the mapped reads (using HiC-Pro)
- Generating Topologically Associating Domains (TADs) from the Hi-C contact map
- Visualizing and analyzing the TADs in the UCSC genome browser

In some places, there are orange lines
```diff
! Like this
```
Here, you should stop and note down your answers to the questions asked. We will go through these in plenary later, so you can see if you understood correctly. Also note down any questions you might have.

Good luck!

**1. Login to the Fox HPC cluster:**

```bash
ssh [ec-username]@fox.educloud.no
```

**2. Setting up today's working directory in your home directory**
```bash
cd
git clone https://github.com/MBV-INF4410/3D-day-1.git
cd 3D-day-1
```

```diff
! Do you recognize the file type in the fastq/chr18 folder?
! Are these paired-end or single-end reads?
! How long are the reads (hint: zcat fastq/chr18/SRR6657510_chr18_R1.fastq.gz | head)
! Another hint: To count the number of nucleotides in the read, you could use R-command `nchar`
```

**3. Setting up nescessary files for running HiC-Pro** 
HiC-Pro will be used to process the Hi-C data, including mapping the reads and aggregation of the contact frequencies (takes a few minutes).
```bash
curl -O https://raw.githubusercontent.com/nservant/HiC-Pro/master/config-hicpro.txt
curl -O https://raw.githubusercontent.com/nservant/HiC-Pro/master/annotation/chrom_hg19.sizes
curl -O https://raw.githubusercontent.com/nservant/HiC-Pro/master/annotation/HindIII_resfrag_hg19.bed
tar -zxvf hg19_chr18/* -C hg19_chr18
```
```diff
! Look at the file chrom_hg19.sizes (hint: use `head`) and HindIII_resfrag_hg19.bed. 
! Can you guess what these two files contain?
```

**4. Adapting the `config-hicpro.txt` file**
To prepare for running HiC-Pro, we will need to change two lines in the `config-hicpro.txt` file. 
Start by copying your current working directory to your clipboard by:
```bash
pwd
```
and copying the resulting path to your clipboard (using "Control + C" / "command + C"). We will call this '[pwd]' from now on

Use a text-editor (we have been using `nano` previously in the course) to edit the `config-hicpro.txt` file: 
```bash
nano config-hicpro.txt
```
- Add the bowtie path to line nr. 39: `BOWTIE2_IDX_PATH =` -> `BOWTIE2_IDX_PATH = [pwd]/hg19_chr18/` where `[pwd]` is the full path to your current working directory. Paste the [pwd] in by "Control + V"
- Change the reference genome on line nr. 47: `REFERENCE_GENOME =` -> `REFERENCE_GENOME = hg19_chr18`
- Change the reference genome on line nr. 48: `GENOME_SIZE = chrom_hg19.sizes` -> `GENOME_SIZE = [pwd]/chrom_hg19.sizes`
- Change restriction fragment on line nr. 67: `GENOME_FRAGMENT = HindIII_resfrag_hg19.bed` -> `GENOME_FRAGMENT = [pwd]/HindIII_resfrag_hg19.bed`
- Change line nr. 89: `BIN_SIZE = 20000 40000 150000 500000 1000000` -> `BIN_SIZE = 50000 1000000`
- Save your changes ("Control + O" / [WriteOut])
- Exit with Control + X
- Lost? To see which line you are at, you can use Control + W. If you have forgotten how to use `nano`, go back to the material from the first week

```diff
! What do you think line number 89 in the config file specifies?
```

**5. Set up your own interactive environment**
Like earlier in the week, we will use `sallic` to allocate resources in an interactive enviroment:
```bash
salloc --ntasks=1 --mem-per-cpu=4G --time=03:00:00  --account=ec34
```

**6. Run HiC-Pro using a container called hicpro_2.11.3_ubuntu.img (Takes ~15 minutes)**
```bash

apptainer exec /projects/ec34/biosin5410/HiC/hicpro_2.11.3_ubuntu.img HiC-Pro --input fastq --output hicpro_results --conf config-hicpro.txt
```
```diff
! Look at the output from HiC-Pro. 
! Compare with the slides from the presentation to see if you can understand what is happening
```

**7. Setup the folder structure for the HiC contacts**
```bash
mkdir -p hic/bedpe/intra
mkdir -p hic/matrix
```

**8. Convert Hi-C to BEDPE and matrix format**

The output from HiC-Pro needs to be converted to [BEDPE](https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format) in order to be processed further by Chrom3D, and to a matrix format in order to be compatible with the Armatus TAD caller.
```bash
awk 'NR==FNR { map[$4] = $1"\t"$2"\t"$3; next } { print $0,map[$1],map[$2] }' hicpro_results/hic_results/matrix/chr18/raw/50000/chr18_50000_abs.bed hicpro_results/hic_results/matrix/chr18/raw/50000/chr18_50000.matrix  | awk '$4==$7' | awk '{print $4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$3}' > hic/bedpe/intra/chr18

for chr in hic/bedpe/intra/*
do
chrname=$(basename $chr)
cut -f 2,5,7 $chr > hic/matrix/$chrname
done
```
```diff
! Compare the files hic/bedpe/intra/chr18 and hic/matrix/chr18 (hint: use `head`). 
! What is the difference?
```

**9. Running Armatus to call TADs**
```bash
module purge
module load Armatus/2.3-GCC-11.3.0
mkdir hic/tads
armatus -r 50000 -c chr18 -S -i hic/matrix/chr18 -g .6 -o hic/tads/chr18
```
The parameter `-r 50000` sets the bin-size to 50000 bp,  `-c chr18` specifies that only chromosome 18 should be considdered, `-S` specifies that sparse matrix format (3 column text file) is used, `-i hic/matrix/chr18` provides the input data (in matrix format), `-g .6` is the gamma-max parameter indicating the highest resolution to generate domains (often is set based on trial and error), `-o hic/tads/chr18` gives the output for the TADs.

**10. Since Armatus output is end-inclusive, convert to BED by adding 1 to end position**
```bash
awk '{printf("%s\t%i\t%i\n",$1,$2,$3+1)}' hic/tads/chr18.consensus.txt > hic/tads/chr18.consensus.bed
```

```diff
! Try to explain how the awk script above works and what it does.
! Compare the files hic/tads/chr18.consensus.txt and hic/tads/chr18.consensus.bed
! What is the difference between the two files?
```

**11. Download the `hic/tads/chr18.consensus.bed` file to your local computer**

**12. Investigating the TADs (chr18.consensus.bed) in the UCSC genome browser**
- Go to the UCSC genome browser at  http://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg19
- Upload the BED file `chr18.consensus.bed` going through step 1-8 in the following visualization:

![UCSC_upload_BED](https://user-images.githubusercontent.com/5373069/100238933-1066af00-2f31-11eb-93d1-3945f8879dd6.png)

- Turn on a "Pack" view of your uploaded "TADs" track by going through step 1-4 in the following visualization:

![UCSC_display_pack](https://user-images.githubusercontent.com/5373069/100244151-e9ab7700-2f36-11eb-9196-5f7262582f50.png)


**13. Try to answer the following questions**
```diff
! Do you see a relationship between the borders of the TADs (the parts where the line breaks up/down) and the genes? [if so, what kind(s) of relationship(s)]?
```

```diff
! Select "Regulation" -> "Rao 2014 Hi-C" -> "Full" to turn on some Hi-C data visualization. 
! Can you see a relationship between the Hi-C data and the TADs?
```
**Bonus: 14. Visualize the Hi-C data as a heatmap**

A classical scenario in bioinformatics is that a tool has been build with libraries not available on your current HPC system. Unless the developer of the tool is already providing a container for us, it may seem impossible to run the tool on the HPC system. But fear not! With apptainer, we can easily deploy environments from existing images, An Apptainer container can e.g. easily built with the 'build' command, for example from existing containers with the environment we need. For this bonus exercise, we need an older version of python (Python2.7), and a library called matplotlib. Someone has already made a container with these requirements, and this can be made into an apptainer container with: apptainer build python-healpy.sif docker://mwirtz/python2.7-healpy
Do not do this yourself, this has been done already, and a container is available at /projects/ec34/biosin5410/HiC/python-healpy.sif 

First, clone the `HiCPlotter` tool that we need for creating visualizations of Hi-C data:
```bash
module purge
git clone https://github.com/kcakdemir/HiCPlotter.git
```
The tool requires python2.7 which is not available on Fox, but we can run interactively within a container (`python-healpy.sif`) which already contains Python2.7 and matplotlib:
```bash
apptainer shell /projects/ec34/biosin5410/HiC/python-healpy.sif 
```
 Once inside the interactive apptainer shell we can use the Python2.7-program HiCPlotter:
```bash 
python HiCPlotter/HiCPlotter.py -f hicpro_results/hic_results/matrix/chr18/iced/50000/chr18_50000_iced.matrix -o example -r 50000 -tri 1 -bed hicpro_results/hic_results/matrix/chr18/raw/50000/chr18_50000_abs.bed -n chr18 -chr chr18 -ptr 1 -hmc 1 -up 1
mv "example-chr18.ofBins(0-1561).50K.png" chr18.png
```
use the command `exit` to exit the apptainer interactive shell

Download `chr18.png` to your own computer and visualize it there.
