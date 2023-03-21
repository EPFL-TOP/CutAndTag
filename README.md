# CutAndTag
using conda virtual environment
Get the latest installer from [here](https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh) and install it

```shell
bash Miniconda3-latest-MacOSX-x86_64.sh
```

We need to get `bioconda` channel on top, let's see if that is the case:

```shell
conda config --show channels
```

if not add it
```shell
conda config --add channels bioconda
```

check that the channel priority is `flexible`
```shell
conda config --describe channel_priority
```

On mac M1, pay attention to `osx-arm64` and `osx-64` order in `conda info` as `bioconda` is empty for `osx-arm64`

Let's now run the virtual env creation:

```shell
conda env create -f cutandtag_env.yaml
```

Answer yes towards the end and activate the environment each time it is needed `conda activate cutandtag`


Now let's run !
We first need to get the `fastq` files somewhere on our local computer. For example in `data`

To run the fastq quality check, just do:
```shell
mkdir -p data/fastqc/
fastqc -o  data/fastqc -f fastq data/fastq/Dr_GFPpoly_S17/Dr_GFPpoly_S17_R1_001.fastq.gz
fastqc -o  data/fastqc -f fastq data/fastq/Dr_GFPpoly_S17/Dr_GFPpoly_S17_R2_001.fastq.gz
fastqc -o  data/fastqc -f fastq data/fastq/Dr_VenusAb_S18/Dr_VenusAb_S18_R1_001.fastq.gz
fastqc -o  data/fastqc -f fastq data/fastq/Dr_VenusAb_S18/Dr_VenusAb_S18_R2_001.fastq.gz
fastqc -o  data/fastqc -f fastq data/fastq/Dr_H2AZ_S19/Dr_H2AZ_S19_R1_001.fastq.gz
fastqc -o  data/fastqc -f fastq data/fastq/Dr_H2AZ_S19/Dr_H2AZ_S19_R2_001.fastq.gz
fastqc -o  data/fastqc -f fastq data/fastq/Dr_Neg_S20/Dr_Neg_S20_R1_001.fastq.gz
fastqc -o  data/fastqc -f fastq data/fastq/Dr_Neg_S20/Dr_Neg_S20_R2_001.fastq.gz
```

We now need to build the `bowtie2` index using the `ensembl` reference genome.
First get the `top_level` reference genome:

```shell
mkdir -p data/ensembl
curl https://ftp.ensembl.org/pub/release-109/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.toplevel.fa.gz -o data/ensembl/Danio_rerio.GRCz11.dna.toplevel.fa.gz
```

and crete the bowtie2 index (this might take some time...)

```shell
mkdir -p data/bowtie2Index
bowtie2-build data/ensembl/Danio_rerio.GRCz11.dna.toplevel.fa.gz data/bowtie2Index/danio_rerio
```

Now for each of our sample, we can run the alignment  (this might take some time especially for H2AZ:

```shell
mkdir data/mapping/
bowtie2 --local --very-sensitive-local --no-unal --phred33 -I 10 -X 700 -p 10 -x data/bowtie2Index/danio_rerio -1 data/fastq/Dr_GFPpoly_S17/Dr_GFPpoly_S17_R1_001.fastq.gz -2 data/fastq/Dr_GFPpoly_S17/Dr_GFPpoly_S17_R2_001.fastq.gz -S data/mapping/Dr_GFPpoly_S17_bowtie2.sam &>> summary_alignment.txt
bowtie2 --local --very-sensitive-local --no-unal --phred33 -I 10 -X 700 -p 10 -x data/bowtie2Index/danio_rerio -1 data/fastq/Dr_VenusAb_S18/Dr_VenusAb_S18_R1_001.fastq.gz -2 data/fastq/Dr_VenusAb_S18/Dr_VenusAb_S18_R2_001.fastq.gz -S data/mapping/Dr_VenusAb_S18_bowtie2.sam &>> summary_alignment.txt
bowtie2 --local --very-sensitive-local --no-unal --phred33 -I 10 -X 700 -p 10 -x data/bowtie2Index/danio_rerio -1 data/fastq/Dr_H2AZ_S19/Dr_H2AZ_S19_R1_001.fastq.gz -2 data/fastq/Dr_H2AZ_S19/Dr_H2AZ_S19_R2_001.fastq.gz -S data/mapping/Dr_H2AZ_S19_bowtie2.sam  &>> summary_alignment.txt
bowtie2 --local --very-sensitive-local --no-unal --phred33 -I 10 -X 700 -p 10 -x data/bowtie2Index/danio_rerio -1 data/fastq/Dr_Neg_S20/Dr_Neg_S20_R1_001.fastq.gz -2 data/fastq/Dr_Neg_S20/Dr_Neg_S20_R2_001.fastq.gz -S data/mapping/Dr_Neg_S20_bowtie2.sam  &>> summary_alignment.txt
```

We now need to convert our `sam` files to `bam` files
```shell
samtools view -h -S -b -o data/mapping/Dr_GFPpoly_S17_bowtie2.bam data/mapping/Dr_GFPpoly_S17_bowtie2.sam
samtools view -h -S -b -o data/mapping/Dr_VenusAb_S18_bowtie2.bam data/mapping/Dr_VenusAb_S18_bowtie2.sam
samtools view -h -S -b -o data/mapping/Dr_H2AZ_S19_bowtie2.bam data/mapping/Dr_H2AZ_S19_bowtie2.sam
samtools view -h -S -b -o data/mapping/Dr_Neg_S20_bowtie2.bam data/mapping/Dr_Neg_S20_bowtie2.sam
```

We sort the bam files
```shell
samtools sort data/mapping/Dr_GFPpoly_S17_bowtie2.bam -o data/mapping/Dr_GFPpoly_S17_bowtie2.sorted.bam
samtools sort data/mapping/Dr_VenusAb_S18_bowtie2.bam -o data/mapping/Dr_VenusAb_S18_bowtie2.sorted.bam
samtools sort data/mapping/Dr_H2AZ_S19_bowtie2.bam -o data/mapping/Dr_H2AZ_S19_bowtie2.sorted.bam
samtools sort data/mapping/Dr_Neg_S20_bowtie2.bam -o data/mapping/Dr_Neg_S20_bowtie2.sorted.bam
```

And we index them
```shell
samtools index data/mapping/Dr_GFPpoly_S17_bowtie2.sorted.bam
samtools index data/mapping/Dr_VenusAb_S18_bowtie2.sorted.bam
samtools index data/mapping/Dr_H2AZ_S19_bowtie2.sorted.bam
samtools index data/mapping/Dr_Neg_S20_bowtie2.sorted.bam
```

We can now run the peak calling
```shell
mkdir data/peakcalling
macs2 callpeak -t data/mapping/Dr_VenusAb_S18_bowtie2.sorted.bam -c data/mapping/Dr_Neg_S20_bowtie2.sorted.bam -g hs -f BAMPE -n macs2_peak_q0.01 --outdir data/peakcalling/VenusAb_S18 -q 0.01 --keep-dup all
macs2 callpeak -t data/mapping/Dr_GFPpoly_S17_bowtie2.sorted.bam -c data/mapping/Dr_Neg_S20_bowtie2.sorted.bam -g hs -f BAMPE -n macs2_peak_q0.01 --outdir data/peakcalling/Dr_GFPpoly_S17 -q 0.01 --keep-dup all
```