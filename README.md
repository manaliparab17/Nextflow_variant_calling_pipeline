Variant Calling Pipeline (Nextflow)
Nextflow (DSL2) version of a variant calling pipeline I originally wrote in bash.
Works on any paired-end fastq sample — reference used here is mouse chr11 (mm39).

What it does

Raw fastq → aligned BAM → variants → filtered → annotated → high-impact variants.

Steps: index reference → FastQC → BWA-MEM alignment → sort/index BAM → flagstat → bcftools call → filter (QUAL/DP) → split SNPs/indels → snpEff annotation → extract HIGH impact.

Requirements

Conda env with: bwa, samtools, bcftools, fastqc, snpeff

Folder structur
Sequences/   -> fastq.gz files
Reference/   -> chr11.fa

Everything else is auto-generated.

Run
bashconda activate <env_name>
nextflow run main.nf --reads 'Sequences/*_{1,2}.fastq.gz' --reference Reference/chr11.fa

Update nextflow.config first: set process.conda to your env's full path, and workDir to where you want intermediate files stored.

Notes
On NTFS-mounted drives (/mnt/d, /mnt/e), symlink staging can fail — stageInMode = 'copy' fixes it.
process.conda needs the full env path, not just the name.
Editing the pipeline breaks -resume caching for affected steps.
