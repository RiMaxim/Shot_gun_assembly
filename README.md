# Pipeline
```
#$1 - sample name
#$2 - threads

#1. Quality contol of raw reads
fastqc -t $2 $1_1.fastq.gz $1_2.fastq.gz

#2. Filtration by trimmomatic:
java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar \
PE \
-threads $2 \
-phred33 \
$1_1.fastq.gz $1_2.fastq.gz \
$1_R1_paired.fastq.gz $1_R1_unpaired.fastq.gz \
$1_R2_paired.fastq.gz $1_R2_unpaired.fastq.gz \
ILLUMINACLIP:/opt/Trimmomatic-0.39/adapters/All_adapters.fa:2:30:10 \
LEADING:0 TRAILING:0 SLIDINGWINDOW:4:0 HEADCROP:0

#3. Quality contol of filtered reads
fastqc -t 2 $1_R1_paired.fastq.gz $1_R2_paired.fastq.gz

#4. Assembly by SPAdes
 /opt/SPAdes-3.15.4-Linux/bin/metaspades.py \
-o ./Spades_assembly/$1 \
-1 ./$1_R1_paired.fastq.gz \
-2 ./$1_R2_paired.fastq.gz \
-s ./$1_R1_unpaired.fastq.gz \
-m 1005 -t $2 -k 33,55,99

#5. Quast
/opt/quast/quast.py -t $2 -m 500 -l ./Spades_assembly/$1/contigs.fasta -o ./quast_results/$1 ./Spades_assembly/$1/contigs.fasta
/opt/quast/quast.py -t $2 -m 5000 -l ./Spades_assembly/$1/contigs.fasta_5kb -o ./quast_results/$1_5kb ./Spades_assembly/$1/contigs.fasta

#6. Filtration contigs by length >= 5kb
python3 Select_contigs_by_length.py \
--inputfile ./Spades_assembly/$1/contigs.fasta \
--outputdir ./Assembly_short_name_5000bp/ \
--outputfile $1.fasta

#7. Run prokka
prokka \
--cpus $2 \
--gcode 11 \
--metagenome \
--force \
--outdir ./Annotation_prokka/$1 \
--prefix $1 \
./Assembly_short_name_5000bp/$1.fasta

#8. Run Antismash
/opt/antismash/run_antismash \
./Annotation_prokka/$1/$1.gbk \
./Antismash/ \
-c $2 \
-v \
--clusterhmmer \
--tigrfam \
--cb-general \
--cb-subclusters \
--cb-knownclusters \
--cc-mibig \
--output-basename $1 \
--genefinding-tool prodigal \
--html-title $1

#9. Run Bigscape
/opt/BiG-SCAPE/run_bigscape \
./Antismash/$1 \
./Bigscape/$1 \
-c $2 \
-v

#10. Run metagenemark
/opt/MetaGeneMark-2/src/gmhmmp2 \
-s ./Assembly_short_name_5000bp/$1.fasta \
-f gff3 \
-M /opt/MetaGeneMark-2/src/mgm2_11.mod \
--gid_per_contig \
--verbose \
--out ./Annotation_mgm/$1.gff3 \
--AA ./Annotation_mgm/$1.fasta

#11. Run Cas_effectors and Cas_proteins
mkdir ./Cas_effectors/$1
hmmsearch --tblout ./Cas_effectors/$1/$1_Cas9.txt -E 10e-3 --cpu 20 ./Cas_profiles/TIGR01865.1.HMM ./Annotation_prokka/$1/$1.faa
hmmsearch --tblout ./Cas_effectors/$1/$1_Cas12.txt -E 10e-3 --cpu 20 ./Cas_profiles/TIGR04330.1.HMM ./Annotation_prokka/$1/$1.faa
mkdir ./Cas_proteins
python3 Read_hmmsearch_output_return_sequences.py -s ./Annotation_prokka/$1/$1.faa -m .k/Cas_effectors/$1/$1_Cas9.txt -e 0.001 -o ./Cas_proteins/$1_Cas9.fasta
python3 Read_hmmsearch_output_return_sequences.py -s ./Annotation_prokka/$1/$1.faa -m ./Cas_effectors/$1/$1_Cas12.txt -e 0.001 -o .k/Cas_proteins/$1_Cas12.fasta

#12.
#13.
#14.
#15.
#16.
#17.
#18.







