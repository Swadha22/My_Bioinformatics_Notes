# My_Bioinformatics_Notes

My Cheat Sheet

## To download data from one server to another:

    sftp chud_sfsu@sftp.genewiz.com
    sftp> lcd ~/swadha/histone/RNAseq_Data
    sftp>  cd /genewiz-us-ngs-sftp/chud_sfsu/30-386429628/00_fastq
    sftp> ls -l
    sftp> get name_of_data

## Important things about SAM format:

  https://genome.sph.umich.edu/wiki/SAM

    #To convert SAM to BAM
     samtools view -S -b sample.sam > sample.bam 
    #Count total number of reads:
      samtools view -c <filename.bam>
    # Count number of mapped reads:
      samtools view -c -F 4 <filename.bam>
    # Count number of unmapped reads:
      samtools view -c -f 4 <filename.bam>
    # To get mapped reads from a bam file:
      samtools view -F 4 <filename.bam>
    # To get unmapped reads from a bam file:
      samtools view -f 4 <filename.bam>
    # To get stats of the sam/bam file
       samtools flagstat
    # To extract unmapped reads whose mates are also unmapped:
      samtools view -f 12 in.bam > uu.sam

     # To extract mapped reads whose mates are unmapped:
      samtools view -F 4 -f 8 in.bam > mu.sam


#### To remove rows from 2 to 4 from a file
    
        sed 2,4d example.txt


### To convert fastq to fasta

         cat test.fastq | paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\t' '\n' > test..fa

        Fastq to fasta conversion
         cat SRR445174_1.fastq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > SRR445174_1.fa

### To get unmapped reads from a bam file with their headers
        samtools view -h -f 4 aln.bam > aln_only_mapped.sam

        # output back to BAM
        samtools view -h -f 4 -b aln.bam > aln_only_mapped.bam


### To count the number of reads having both itself and it’s mate mapped:
        samtools view -c -f 1 -F 12

### Truncating bam files:
        samtools view -h big.bam | head -n 1000 > little.bam

### To trim last 30 bases from a sequence
        fastx_trimmer -t 30 -i SRR445174_1.fastq -Q33 -o SRR445174_1_trimmed.fasq

### To sort a BAM file
        samtools sort sample.bam -o sample.sorted.bam
        
### BAM to FASTQ conversion (paired end reads)
        bedtools bamtofastq -i input_sorted.bam -fq output_r1.fastq -fq2 output_r2.fastq


### SAM/BAM to fasta conversion
         samtools view file.sam | awk '{print ">"$1"\n"$10}' > output.fa


### To copy PATH of each file in folder
        ls -R1  ~/swadha/phylogeny/nematodegenomes/species | while read l; do case $l in *:) d=${l%:};; "") d=;; *) echo "$d/$l";; esac; done > test.txt

### Bowtie-2 mapping command
        bowtie2-build genome.fa indexf
         bowtie2 -x index -1 file1.fastq -2 file2.fastq  -S  bowtie_output.sam

        Bowtie output:
        Field    What You Find There
        1          readname
        2          bitflag -- for SE bt2 mapping, only 0 (pos strand), 4 (unmapped), 16 (neg strand) are relevant
        3          chrName/refname
        4          1-based start on fwd strand (this is the end if on neg strand)
        5          MAPQ
        6          CIGAR string rep of alignment (most are 50M which means 50 match+mismatch)
        7          chrName/refName of mate if PE; * if SE
        8          1 based start of mate; * if SE
        9          Inferred fragment length if PE; * if SE
        10        the read's sequence
        11        the reads base call qualities (Qs)


### EXTENDED ANALYSIS
        perl longest.pl c.elegans.fa > c.elegans.fa.opg
        diamond makedb --in c.elegans.fa.opg --db c.elegans.db
        diamond blastp --query japonica.fa --db c.elegans.db -elE-20 | sort --k=12rn
        Perl get_orthologs.pl 

### Hisat2 mapping command

        hisat2-build genome.fa index

        hisat2_extract_splice_sites.py /home/roylab/swadha/test_trans_splicing/Giardia_lamblia.GL2.40.gtf > giardia_cis__splicesites.txt

        hisat2 -x  /hisat2_ref/ref --known-splicesite-infile giardia_cis__splicesites.txt  -1 1.fastq.gz -2 2.fastq.gz > hisat2_mapped_reads.sam


### BLAT mapping command
       # Index file is same as bowtie’s
        BLAT database query_in_fasta output.pslx

## BLAST
### To make the database:
        makeblastdb -in database_file -dbtype prot -out database.out
### To blast:
        blastp -query querry.fa -db database -out output -outfmt "7 qacc sacc evalue qstart qend sstart send qlen slen length" -evalue 1e-10 

        the parameter to mask low complexity reads in protein query is "-seg". When the query is nucleotide, the parameter is -dust

        tblastn -query query_His-35.fa -db C.elegans -seg no > blastout_c.ele.txt


        Mapping quality → higher = more unique the alignment is. 
        Mapping quality is related to “uniqueness.” We say an alignment is unique if it has a much higher alignment score than all the other possible alignments. The bigger the gap between the best alignment’s score and the second-best alignment’s score, the more unique the best alignment, and the higher its mapping quality should be.
        Accurate mapping qualities are useful for downstream tools like variant callers. For instance, a variant caller might choose to ignore evidence from alignments with mapping quality less than, say, 10. A mapping quality of 10 or less indicates that there is at least a 1 in 10 chance that the read truly originated elsewhere.
        1- Bowtie 2 generates MAPQ scores between 0–42
        2- BWA generates MAPQ scores between 0–37
        3- You should always take a look at your mapped sequence data to see what ranges of scores are present before doing anything else with your BAM/SAM files
        4- TopHat: the meaning of the scores are as follows:
        TopHat outputs MAPQ scores in the BAM/SAM files with possible values 0, 1, 2, or 50. The first three values indicate mappings to 5, 3–4, or 2 locations, whereas a value of 50 represents a unique match. Please note that older versions of TopHat used a value of 255 for unique matches. Further note that standalone versions of Bowtie and Bowie 2 (used by TopHat) produce a different range of MAPQ scores (0–42).

## MACS2
    #   V1: ./macs2 callpeak -t HTZ1.SSTT004_S36HTZ1EX39_TTAGGC_L005_R1_001.sam -c Control.bam  -n HTAS1 -f SAM -g ce --bdg --keep-dup=auto --broad --broad-cutoff=0.01 --nomodel --extsize=250 --SPMR --outdir macs_output

     #   V2: macs2 callpeak -t treatmentfilesam.gz -c inputfile.sam.gz -n name -g ce --bdg --keep-dup=auto --broad --broad-cutoff=0.01 --nomodel --extsize=250 –-SPMR --outdir outDir

      #   V3: Macs2  callpeak -t HTAS1.C.s39.merged.mapped.sam -c /home/roylab/swadha/histone/input_files/Input.S39B.Cov.r2_SSTT007_5_S60_L005_R1_001.sam -n HTAS1 -f SAM -g ce -B --nomodel --extsize 147 --SPMR --outdir macsOutput_V2

### NEXT Step: Fold enrichment and LR lot

        /macs2 bdgcmp -t HTAS1_treat_pileup.bdg -c HTAS1_control_lambda.bdg -o fold_enrichment.bdg -m FE
        /macs2 bdgcmp -t HTAS1_treat_pileup.bdg  -c HTAS1_control_lambda.bdg -o fold_enrichment_logLR.bdg -m logLR -p 0.00001

         * -m FE means to calculate fold enrichment. Other options can be logLR for log likelihood, subtract for subtracting noise from treatment sample.
         * -p sets pseudocount. This number will be added to 'pileup per million reads' value. You don't need it while generating fold enrichment track because control lambda will always >0. But in order to avoid log(0) while calculating log likelihood, we'd add pseudocount. Because I set precision as 5 decimals, here I use 0.00001.

## Bowtie-1
### In the 1st step you, Generate index files:
               bowtie-build c_elegans_genome.fa index_files
### In the second step, Align the ChIP reads to the reference genome: 
               bowtie -x index_files fastq_file -S output.sam

## Bowtie-2
### 1st step, generate index file:
        Bowtie2-build c_elegans_genome.fa index_files
### 2nd step, align the ChIP reads to reference genome (paired end files)
        bowtie2 -x <ref_genome_files> -1 <path_to_fastq 1> -2 <path_to_fastq 2> -S <path_to_output_sam_file>
        Align the ChIP reads to reference genome (single end)
        bowtie2 -x <ref_genome_files> -U <path_to_fastq_files> -S <path_to_output_sam_file>


## To fetch chromosome sizes
        fetchChromSizes ce11 > ce.chrom.sizes


## STAR MApping (http://chagall.med.cornell.edu/RNASEQcourse/STARmanual.pdf)

        STAR --runMode genomeGenerate --genomeDir home/swadha/star_index/  --genomeFastaFiles Giardia_lamblia.GL2.dna.nonchromosomal.fa --sjdbGTFfile Giardia_lamblia.GL2.40.gtf --genomeSAindexNbases 12

        STAR --runMode genomeGenerate --genomeDir /home/roylab/swadha/histone/RNAseq_analysis/genome_index/  --genomeFastaFiles /home/roylab/swadha/histone/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa   --sjdbGTFfile  ~/swadha/histone/Caenorhabditis_elegans.WBcel235.94.gtf  --genomeSAindexNbases 12

         STAR --genomeDir star_index/ --readFilesIn test_1.fq --outFileNamePrefix star  --outFilterMultimapNmax 1 --outReadsUnmapped unmapped_a --outSAMtype BAM SortedByCoordinate

        nohup STAR --genomeDir /home/roylab/swadha/histone/RNAseq_analysis/genome_index/   --readFilesIn /home/roylab/swadha/histone/RNAseq_analysis/RNAseq_Data/HTAS-1_XC18-12_R1_001.fastq.gz /home/roylab/swadha/histone/RNAseq_analysis/RNAseq_Data/HTAS-1_XC18-12_R2_001.fastq.gz  --readFilesCommand zcat --outFileNamePrefix HTAS-1_XC18-12  --outFilterMultimapNmax 1 --outReadsUnmapped unmapped_a --outSAMtype BAM SortedByCoordinate

#### FINAL:
        STAR --genomeDir genome_index/ --readFilesIn R1_001.fastq.gz R2_001.fastq.gz --readFilesCommand zcat --outFileNamePrefix _mapped --outFilterMultimapNmax 10 --outSAMtype BAM SortedByCoordinate

#### Star for ChIP
         nohup STAR --genomeDir /home/roylab/swadha/histone/RNAseq_analysis/genome_index  --readFilesIn /home/roylab/swadha/histone/ChIP_analysis/DATA/HTZ1.S39B.Ex49.Cov.i2.SSPG005B_S54.160510_50SR_HS4KB_L007.fastq.gz  --readFilesCommand zcat --alignIntronMax 1 --alignEndsType EndToEnd  --outSAMtype BAM SortedByCoordinate
        
        
        
## How a quality score is born
        Illumina systems create a quality score in three steps.
        1- First of all it evaluates the detected light signal for each base call for every cluster, on every tile, for every sequencing cycle simultaneous with a sequencing run. Thereby it measures various aspects correlating with the quality of the base call like single-to-noise ratio and light intensity. Based on these parameters a quality predictor value (QPV) is calculated. 
        2- In the second step the QPV has to be translated into a quality score with the help of a quality table. The Q-table based on a statistical calibration curve derived from empirical data including various well-characterized human and not-human samples mostly using a version of the so-called Phred algorithm (1). That is why the score is also called Phred quality score. In the last step the quality score (per cycle) is recorded common with the base call in a base call file (.bcl) which is later converted to FASTQ files (.fastq).


## Z-score
        Z-score= (Rawscore - mean)/standard deviation

## MULTIPLE SEQUENCE ALIGNMENT
        muscle -in sequences.fasta -out alignment
        trimal -in alignment -out alignment_trimmed
        clustalw  -type=dna -infile=gisaid_hcov-19_2020_05_02_07.fasta > testAlignment

## RAxML was called as follows:
        raxmlHPC-SSE3 -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -#10 -s alignment -n SwadhaTree.raxml


## IQ Tree
        nohup iqtree -s alignment -m VT+R9 -bb 1000 -nm 10000

## Server Jobs Information
        ps -eo pid,lstart,cmd | grep 'iqtree'


## To know about recursive filenames
        ls -R1  ~/swadha/phylogeny/nematodegenomes/species | while read l; do case $l in *:) d=${l%:};; "") d=;; *) echo "$d/$l";; esac; done > test.txt

## To download data in bulk from a website
        wget -r -A xyz.fa.gz <<URL>>


## Bhaiya’s Covid Stuff
        sed 's/>/\n>/' gisaid_hcov-19_2020_05_02_07.fasta > 1.fasta
        sed 's/>/#>/' gisaid_hcov-19_2020_05_02_07.fasta > hashInBeginning.fasta
        sed 's/^\(#.*\)$/\1%/' hashInBeginning.fasta > editEndOfHeader.fasta
        perl -pe 's/\n//'g editEndOfHeader.fasta > deleteNewline.fasta
        sed 's/#/\n/' deleteNewline.fasta > replaceHashWithNewline.fasta
        sed 's/%/\n/' replaceHashWithNewline.fasta > finalFatsaFile.fasta
        tail -n +2 file.txt > file.stdout
        zcat  genomic.fa.gz | sed 's/>/#>/' | sed 's/^\(#.*\)$/\1%/'  | perl -pe 's/\n//'g | sed 's/#/\n/' | sed 's/%/\n/' > .txt


## To delete file extension after a certain regex
        rename -v -n 's/[.PRJ].*//' *[.PRJ]*




## Bash loop to read a file line-by-line
        while IFS= read -r line
                   do
                         cp $line /home/roylab/swadha/phylogeny/nematodegenomes/GFF_files/
                         #command2 on $line

                   done < /home/roylab/swadha/phylogeny/nematodegenomes/gff_species_file_paths.txt


## Bash nested loops
        while read filename_line
        do
            while read genename_line
            do

                cat $filename_line* | grep -A 1 "$genename_line"


            done < geneNames.txt
        done < filename.txt





## Perl Script to call another script
        use strict;
        use warnings;

        my $filename = '/home/roylab/swadha/phylogeny/nematodegenomes/genome_file_names.txt';
        open(my $fh, '<:encoding(UTF-8)', $filename)
          or die "Could not open file '$filename' $!";


        while (my $row = <$fh>) {
          chomp $row;
          #print "$row\n";
          my $x = $row;
          my $y = substr($x, 0, index($x, 'WBPS11.genomic.fa'));
          #print "$y\n"


          my $gff_ext = "WBPS11.annotations.gff3";
          my $genome_ext = "WBPS11.genomic.fa";



          my $gff_file = '/home/roylab/swadha/phylogeny/nematodegenomes/GFF_files/' . $y . $gff_ext;
          my $genome_file = '/home/roylab/swadha/phylogeny/nematodegenomes/genome_files/' . $y . $genome_ext;

           print "$gff_file\n";
           #print "$genome_file\n";

           system ("exint-master.pl $gff_file $genome_file > /home/roylab/swadha/phylogeny/nematodegenomes/intron_exon_files/$y.exons-introns");
        }





## To delete file extensions
        find -type f -name '*.txt' | while read f; do mv "$f" "${f%.txt}"; done

## To delete all the text after a delimiter in a line
        cut -f1 -d":" <<filename>>
        sed 's/\.com.*/.com/' file.txt

## To include filenames as first column 
         for f in *.exons-introns; do sed "s/^/$f\t/" $f; done
        OR
        for i in *.txt; do nawk '{print FILENAME","$0}' $i ; done > test.txt

## To delete underscores from sequences without changing the headers
        awk 'BEGIN{RS=">";FS="\n"}NR>1{printf ">%s\n",$1;for (i=2;i<=NF;i++) {gsub(/.../,"",$i);printf "%s\n",$i}}' <<filename>>

## To delete lowercase characters from sequences without changing the headers
        awk 'BEGIN{RS=">";FS="\n"}NR>1{printf ">%s\n",$1;for (i=2;i<=NF;i++) {gsub(/[[:lower:]]/,"",$i);printf "%s\n",$i}}' 


## To add a character at the end of each line in fasta
        sed 's/.*/&char/' abcd.txt
## To remove spaces from the headers
        sed 's, ,_,g' file.txt > xyz.txt

## To add count numbers to the headers of the sequences
        awk '/^>/{$0=$0"_"(++i)}1' infile

## To extract a tar.gz file, use the --extract (-x) option and specify the archive file name after the f option:
        tar -xf archive.tar.gz

## To get Nth line from a text file
        cat output.txt | head -10 | tail -1

## To fetch characters from x position to y
        cat test.txt | head -c 4153809 | tail -c 744 > test2.txt

## To get a whole sequence from the genome file using the contig name
        samtools faidx c.genomic.fa Contig27 > test.txt

## outputs the content of the file starting with the 11th line
        tail -n +11 file

## To add a header in a file
        echo -e "Runs\tOpposition\tDate" | cat - table > table_with_header
        sed '1i Runs\tOpposition\tDate' table > table_with_header
        awk 'BEGIN {print "Runs\tOpposition\tDate"} {print}' table > table_with_header
        
## To remove string between two HASH characters
        sed -e 's/\(#\).*\(#\)/\1\2/' filename.txt

## To remove repeated sequences
        awk '/^>/{f=!d[$1];d[$1]=1}f' in.fa > out.fa
## To remove 1st character of lines other than header.
         sed -i 's/^[-]\{1\}//' *.txt
## To delete two consecutive dots
        perl -pe 's/\.\./\n/’g < inputfile 

## To add a character % at the end of each header in fasa
        sed 's/^\(>.*\)$/\1%/' test.txt > test2.txt

## To put a ‘>’ at the beginning of each line
        sed 's/^/>/' filename

## To remove 1st line of a file
        tail -n +2 "$FILE"

## To filter sequences having a length less than or equal to 200 aminoacid
        awk '!/^>/ { next } { getline seq } length(seq) <= 200 { print $0 "\n" seq }' <filename>

        awk '!/^>/ { next } { getline seq } length(seq) >= 90 { print $0 "\n" seq }' 14_format_change.txt



## To delete size zero files
        find . -type f -empty -delete
## To delete 4th column
        sed -i -r 's/\S+//3' file


## To get unique words in each row
        #!/usr/bin/env python
        for line in open('introns.txt', 'r'):
             seen = []
            words = line.rstrip('\n').split(',')
            for word in words:
                if not word in seen:
                    print word,
                    seen.append(word)
            print 

## To fetch only those lines which do not match the string
        while IFS= read -r line
                   do         
            perl -ni.bak -e "print unless /$line/" step_11_groupBySpecies.txt              
                   done < non_nematode_list.txt


## To make Index file from BAM file
        samtools index -b /home/roylab/swadha/histone/rep_3/bowtie_2/TMT01e_S5.bam /home/roylab/swadha/histone/rep_3/bowtie_2/TMT01e_S5.index
## To make fragment size plots
        ./bamPEFragmentSize -o /home/roylab/swadha/histone/rep_3/bowtie_2/fragment_size.png -T 75  -b /home/roylab/swadha/histone/rep_3/bowtie_2/sorted/TMT01e_S5_sorted.bam
        OR

        ./bamPEFragmentSize -o /home/roylab/swadha/histone/rep_3/bowtie_2/fragment_size.png -T 130  -b /home/roylab/swadha/histone/rep_3/bowtie_2/sorted/TMT01e_S5_sorted.bam /home/roylab/swadha/histone/rep_3/bowtie_2/sorted/TMT01a_S1.bam

## To delete the lines which have a certain string
        sed -i '/MtDNA/d' HTAS1.N_TMT01b_S2_sorted.bam.bdg
## To delete the lines which have a MtDNA and X chr string
        sed '/MtDNA/d;/X/d' $line > file.bdg



## To make bed file from bam file (index file should be in the folder ending with .bai extension)
        bamCoverage --bam a.bam -o a.SeqDepthNorm.bw --binSize 50

## To change >I to >chr1 

        cat TMT01d_S4_sorted.bam.bdg | sed 's/^/>/' |  sed 's/>I\t/chr1\t/g' | sed 's/>II\t/chr2\t/g' | sed 's/>III\t/chr3\t/g' | sed 's/>III\t/chr3\t/g' | sed 's/>IV\t/chr4\t/g' | sed 's/>V\t/chr5\t/g' | sed 's/>X\t/chrX\t/g' |  sed 's/>//' > test.txt


## To make bins in the genome file
        bedtools makewindows -g ce.chrom.sizes -w 50 > genome.1Mb.bed

## To sort all the bam files together
        for f in *.bam; do filename="${f%%.*}"; samtools sort $f -o $f.sorted.bam; done

## Sort the bed file
        bedtools sort -i TMT01c_S3_sorted.bam_V2.bdg > sorted_bed_file
       
## screen
        screen  (http://aperiodic.net/screen/quick_reference) is a command-line tool which creates a virtual session on a server, and lets you to run commands as normal while allowing you to close and reopen the connection while everything continues to run.
        
        For us, the most useful application of screen is starting and monitoring long-running processes, without having to send them to the background. This has the benefits of letting you monitor messages being send to the screen, and (if needed) killing the process via control + c

        basic usage
        Start a  new screen: screen -S [screen name]
        List all available screens: screen -ls
        Detach from the current screen: control + a, d (while holding control, hit a then d)
        Reattach to a screen: screen -r [screen name]
        Kill the current screen: control + d

        A word of caution: it is possible to get oneself into an Inception-style screen situation if one creates a new screen within an existing screen. If you find that something just isn't making sense, check to make sure you haven't inadvertently Incepted yourself.



## Bin the sorted bed file 
        bedtools intersect -a genome.1Mb.bed -b data.sorted.bed -c -sorted
##  To sort Bam file
       samtools sort my.sam > my_sorted.bam
## bamCoverage: to count ChIP reads (USED THIS IN ChIP reads)
        samtools index -b -m 50 BAMFILE output.bai (this .bai has to be in the folder to run bamCoverage)
        bamCoverage --bam a.bam -o a.SeqDepthNorm.bw  --binSize 50

## To convert GBFF to GFF3 file format
        Install bioperl
        bp_genbank2gff3 Caenorhabditis_sulstoni.gbff.gz

## To find chromosoem size
        /home/roylab/swadha/histone/bin/fetchChromSizes ce11 ce.chrom.sizes

## To sort all the bam files at once
        for f in *.bam; do  samtools sort $f  ${f/.bam/sorted.bam}; done

## To Find Z scores in bedgraph file

        import pandas as pd
        import numpy as np
        from scipy import stats
        my_data = pd.read_csv("/home/roylab/swadha/histone/bam_files/step3.1_un_normalized_bedgraph/HTAS1.C.s39.merged.mapped.sorted.rmDup.bam.sorted.bam.bedgraph", sep="\t", header=None)
        #my_data.head(2)
        my_data.columns = ["chr","start","end","ReadsMapped"]
        a = my_data["ReadsMapped"]
        Zscores = pd.Series(stats.zscore(a))
        my_data['ReadsMapped'] = Zscores
        np.savetxt(r'/home/roylab/swadha/histone/bam_files/step4.bdgNormbyZscores/test.txt', my_data.values, fmt='%s')




## To fetch and analyze contigs from blast output
        samtools faidx caenorhabditis_remanei.PRJNA248911.WBPS11.genomic.fa Contig27 | grep -v "^>" | perl -pe 's/\n//' |  cut -c1-6

## ChIP data analysis steps:
        1- Mapping

        2- samtools index -b -m 50 BAMFILE output.bai

        3- bamCoverage -of bedgraph  --normalizeUsing RPKM --bam file.bam -o $line.bedgraph --binSize 50dee

        4- Normalize using Zscores

        4-  computeMatrix scale-regions --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 -S HTAS1.bw  HTZ_1.bw input.bw  -R /home/roylab/swadha/histone/Caenorhabditis_elegans.WBcel235.94.gtf -o Replicate_S43_antibody3656

        4- computeMatrix reference-point -S <biwig file(s)> -R /home/roylab/swadha/histone/Caenorhabditis_elegans.WBcel235.94.gtf -o Replicate_S43_antibody3656 -a 5000 -b 1000 --metagene


        5-  ~/swadha/histone/bin/deepTools-master/bin/plotProfile -m s39DataMatrix -out s39DataV3.png  --perGroup --colors red yellow blue --plotTitle "S39 Replicate Data"  --plotWidth 20

## bamCompare (to normalize sample to input)
         ~/swadha/histone/bin/deepTools-master/bin/bamCompare -b1 HTAS1.C.S39B..sam -b2 Input.S39B..sam --binSize 1 -o /home/roylab/swadha/histone/All_bam_files_together/BAM_to_bigwigs_viaBamCompare_normalization/HTAS1.C.S39.bw

## computeMatrix scale-regions
        nohup ~/swadha/histone/bin/deepTools-master/bin/computeMatrix scale-regions --beforeRegionStartLength 1000 --afterRegionStartLength 5000 -R /home/roylab/swadha/histone/Caenorhabditis_elegans.WBcel235.94.gtf -S HTAS1.C.S41B..sam.bedgraph.bw HTAS1.N.S41B..sam.bedgraph.bw HTZ1.S41B..sam.bedgraph.bw Input.S41B..sam.bedgraph.bw -o /home/roylab/swadha/histone/All_bam_files_together/computeMatrix_Binsize_ONE/for_quals/scaled_regions/S_41 --binSize 1  --startLabel TSS --endLabel 5000 --samplesLabel HTAS-1.C HTZ-1 Input 

## computeMatrix reference-point
         nohup ~/swadha/histone/bin/deepTools-master/bin/computeMatrix reference-point -R /home/roylab/swadha/histone/Caenorhabditis_elegans.WBcel235.94.gtf -S HTAS1.C.S39B..sam.bedgraph.bw HTZ1.S39B..sam.bedgraph.bw Input.S39B..sam.bedgraph.bw --referencePoint TSS -b 1000 -a 5000 --binSize 1  --samplesLabel HTAS-1.C HTZ-1 Input -o /home/roylab/swadha/histone/All_bam_files_together/computeMatrix_Binsize_ONE/for_quals/reference_point/S39

## Multibamsummary:
        nohup multiBamSummary bins -bs 50 --bamfiles  HTAS1.C.S39B..bam HTAS1.C.S39B.reseq..bam HTAS1.C.S41B..bam  HTAS1.C.S41B.reseq..bam HTAS1.N.S41B..bam HTAS1.N.S41B.reseq..bam HTAS1C_Rb3656_TMT01d_S4.bam HTAS1C_Rb5657_TMT01c_S3.bam  HTAS1.C.S43_3656_TMT01d_S4_NoDups.bam HTAS1.C.S43_5657_TMT01c_S3_noDups.bam HTAS1.N.S43_TMT01b_S2.bam HTAS1.N.S43_TMT01b_S2_noDups.bam HTZ1.S39B..bam HTZ1.S39B.reseq.bam HTZ1.S41B..bam  HTZ1.S41B.reseq..bam HTZ_1_.S43.TMT01e_S5_NoDups.bam HTZ-1_TMT01e_S5sorted.bam Input.S39B..bam Input.S41B..bam InputTMT01a_S1.bam input_.S43.TMT01a_S1_noDups.bam -o all_in_one_MultiBamSummery.npz



## To remove duplicates: -s Remove duplicate for single-end reads. By default, the command works for paired-end reads only. -S Treat paired-end reads as single-end reads.
        samtools rmdup <input.bam> <output.bam>
        OR
        samtools rmdup -s <input.bam> <output.bam>

## Tell grep to use the regular expressions as defined by Perl (Perl has \t as tab):
        grep -P "\t" <file name>
## Use the literal tab character:
        grep "^V<tab>" <filename>
## Use printf to print a tab character for you:
        grep "$(printf '\t')" <filename>


## To find mean of column 4
         awk '{ total += $4 } END { print total/NR }' filename.txt

## To subtract each element of a column 4 by mean
        awk -v s=3.58702 '{print $1, $2, $3, $4-s}' filename.txt

## To divide each element of a column 4 by mean
        awk -v s=3.58702 '{print $1, $2, $3, $4/s}' filename.txt

## Observed-Expected/Expected
        awk -v s=3.58702 '{print $1,"\t", $2,"\t", $3,"\t", ($4-s)/s}' filename.txt > output.txt

## To find maximum in a column
        cut -f  4 file.bdg  | grep -Eo '[0-9]+' | sort -rn | head -n 1

## To find minimum in a column 
        cut -f  4 filename.bdg  | grep -Eo '[0-9]+' | sort -rn | tail -n 1

## To print uniq values in column 4
        cut -f  4 file.bedgraph  | grep -Eo '[0-9]+' | sort -rn | cat | uniq

## To sort by column 4 and print
        sort -k4 -n file.bdg
## To fetch all the rows whose 4th column value is greater than 10
        sort -k4 -n file.bedgraph | awk -F"\t" '$4>10' 


## Togenerate multibam comparision matrix
         nohup ~/swadha/histone/bin/deepTools-master/bin/multiBamSummary bins -bs 50 --bamfiles a,b,c,d,e -o output.npz -l sample1 sample2 sample3 sample4 sample5 --outRawCounts rawReadCounts.table

## Coorelation Plot
        ~/swadha/histone/bin/deepTools-master/bin/plotCorrelation --corData Replicate_S43_spearmanCor.npz  --corMethod spearman --colorMap RdYlBu --plotNumbers  --plotTitle "Spearman Correlation of S43 Read Counts" --whatToPlot heatmap  -o S43_heatMap.corr.png


## To remove duplicates by markdup:
        # The first sort can be omitted if the file is already name ordered
        samtools sort -n -o namesort.bam example.bam

        # Add ms (mate score) and MC (mate cigar) tags for markdup to use later
        samtools fixmate -m namesort.bam fixmate.bam

        # Markdup needs position order
        samtools sort -o positionsort.bam fixmate.bam

        # Finally remove the marked duplicates
        samtools markdup -r positionsort.bam markdup.bam
        To rename att the files at once
         rename 's/.sam.sorted.bai/noDup.bai/' *.sam.sorted.bai

## To list the filenames along with its path
        ls -d "$PWD"/*.sam > filenames.txt


## To count the number of columns in a file separated by \t
        head -2 format_change_codons.txt |tail -1 |tr '\t' '\n' |wc -l



## RNASeq
        infer_experiment.py -i HTAS-1_XC18-12Aligned.sortedByCoord.out.bam -r /home/roylab/swadha/histone/ce11_RefSeq.bed

        Reading reference gene model /home/roylab/swadha/histone/all.genes.BED ... Done
        Loading SAM/BAM file ...  Total 200000 usable reads were sampled
        This is PairEnd Data
        Fraction of reads failed to determine: 0.1060
        Fraction of reads explained by "1++,1--,2+-,2-+": 0.4041
        Fraction of reads explained by "1+-,1-+,2++,2--": 0.4900

         nohup FPKM_count.py -i HTAS-1_XC18-12Aligned.sortedByCoord.out.bam -o FPKM_values -r /home/roylab/swadha/histone/ce11_RefSeq.bed -d 1++,1--,2+-,2-+ 1+-,1-+,2++,2--

        read_distribution.py -r  /home/roylab/swadha/histone/ce11_RefSeq.bed -i HTAS-1_XC18-12Aligned.sortedByCoord.out.bam > read_distribution.txt


        Rscript to_visualizeReadDistributionasBarChart.R

        nohup geneBody_coverage.py -i /home/roylab/swadha/histone/RNAseq_analysis/mapped_reads/HIS-35_FX1328-12/HIS-35_FX1328-12Aligned.sortedByCoord.out.bam /home/roylab/swadha/histone/RNAseq_analysis/mapped_reads/HTAS-1_XC18-12/HTAS-1_XC18-12Aligned.sortedByCoord.out.bam /home/roylab/swadha/histone/RNAseq_analysis/mapped_reads/wildType_N2-12/wildType_N2-12Aligned.sortedByCoord.out.bam  -r /home/roylab/swadha/histone/ce11_RefSeq.bed -o geneBodyCoverage_Rep-1

## HTSEQ
        htseq-count /home/roylab/swadha/histone/RNAseq_analysis/1_mapped_reads/HIS-35_FX1328-12/HIS-35_FX1328-12Aligned.sortedByCoord.out.bam  /home/roylab/swadha/histone/Caenorhabditis_elegans.WBcel235.94.gtf -o /home/roylab/swadha/histone/RNAseq_analysis/2_HTSEQ/HIS-35_FX1328-12Aligned.sortedByCoord.out.bam.txt

        nohup htseq-count -m union -f bam HIS-35_FX1328-12Aligned.sortedByCoord.out.bam Caenorhabditis_elegans.WBcel235.94.gtf -r pos -c HIS-35_FX1328-12_htseq.txt -s no -o HIS-35_FX1328-12_htseq.bam

## Final HTSEQ
        htseq-count -m union -f bam mapped..bam  Caenorhabditis_elegans.WBcel235.102.gtf -r pos -c htseq.txt -s no



## To Convert bam file into bed file, using Bedtools
        bamToBed -split -i filename.bam > filename.bed
## Determine scaling factor: count mapped reads after MtDNA and rRNA removal to scale
        samtools view -c -F 4 filename.bam
        #10M divided by number of mapped reads to determine the scaling factor
        #STEP8b.Converting to bedgraph file and scaling to 10M reads
        genomeCoverageBed -scale 0.85 -i accepted_hits.bed -bg -g ~/Documents/Strome_Lab/Bioinformatics_Workshop/Generating_browser_tracks/ce10.fai.txt > pgc.10K.bg
        #STEP9b. Making bigwig files from bedgraph to display on the genome browser
        bedGraphToBigWig pgc.10K.bg ~/Documents/Strome_Lab/Bioinformatics_Workshop/Generating_browser_tracks/ce10.fai.txt pgc.10K.bw 


        Counting number of mapped reads of all the chromosomes except MtDNA and X
        samtools view -b sortedByCoord.out.bam I II III IV V | samtools view -c -F 4



        Sorting bedgraph_file
        sort -k1,1 -k2,2 unsorted.bedGraph > sorted.bedGraph

        Scaling in all the 33 files using bash loop (input file has two columns, first will be in variable 1: “filename”, and 2nd column will be used in variable 2, “scale_factor”)
        while IFS= read -r line
                   do         
            filename=$(echo "$line" | cut -f 1)
            scale_factor=$(echo "$line" | cut -f 2)

            genomeCoverageBed -scale "$scale_factor" -i "$filename" -bg -g Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.fai > "$filename"
                   done < filenames.txt


        To calculate mean fragment length
        java -jar picard.jar CollectInsertSizeMetrics I=/home/roylab/swadha/histone/RNAseq_analysis/1_all_bam_files/male_germ_WT_N2-12_mappedAligned.sortedByCoord.out.bam O= /home/roylab/swadha/histone/RNAseq_analysis/output.txt H=insert_size_histogram.pdf M=0.5


        UCSC tools
        http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/




