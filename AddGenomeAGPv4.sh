#!/bin/bash
# USER DEFINED INSTALLATION DIRECTORY


#Prompt for name of genome
echo -e "\nAGPv4"


# Check if mysql and R is installed
#my $FC_path = `which mysql
#print STDERR "\n\nMysql is not installed. \n" unless ( $FC_path );
#exit unless ( $FC_path );

#my $R_path = `which R
#print STDERR "\n\nR is not executable, please see: http://www.r-project.org/\n\n" unless ( $R_path );
#exit unless ( $R_path );

# Fetch gene annotation and run R scripts

echo -e "\n\tSTATUS:"
echo -e "\t\tMaking Temporary Folder."

export TMPDIR=/home/lschulte/iRNA-v1.1/tmp

# Make a filtered inclusion list
echo -e "\t\tFiltering transcripts.";
VARIABLE=1
Rscript $INSTALLDIR"/home/lschulte/iRNA-v1.1/bin/Filter_BED.R" $TMPDIR $VARIABLE

# Give each gene a unique identifier
awk '{ print $0"\t"NR }' $TMPDIR/Gene.Dump > $TMPDIR/Tmp
mv $TMPDIR/Tmp $TMPDIR/Gene.Dump
awk '{ print $0"\t"NR }' $TMPDIR/Gene.Dump.Inc > $TMPDIR/Tmp
mv $TMPDIR/Tmp $TMPDIR/Gene.Dump.Inc

echo -e "\t\tGenerating BED files.";
# Make gene lists
awk '{ print $3"\t"$5"\t"$6"\t"$13"\t0\t"$4 }' $TMPDIR/Gene.Dump | uniq - > $TMPDIR/Gene.All.bed
awk '{ print $3"\t"$5"\t"$6"\t"$13"\t0\t"$4 }' $TMPDIR/Gene.Dump.Inc > $TMPDIR/Gene.Inc.bed

# Make exon lists
awk '{ print $10 }' $TMPDIR/Gene.Dump | sed "s/.$//"  | sed -e 's/,/\n/g' > $TMPDIR/Exon.Start
awk '{ print $11 }' $TMPDIR/Gene.Dump | sed "s/.$//" | sed -e 's/,/\n/g' > $TMPDIR/Exon.End
awk '{i = 1; while (i <= $9) { print $4; i++ }}' $TMPDIR/Gene.Dump > $TMPDIR/Exon.Strand
awk '{i = 1; while (i <= $9) { print $3; i++ }}' $TMPDIR/Gene.Dump > $TMPDIR/Exon.Chr
awk '{i = 1; while (i <= $9) { print $13; i++ }}' $TMPDIR/Gene.Dump > $TMPDIR/Exon.Name
awk '{i = 1; while (i <= $9) { print i; i++ }}' $TMPDIR/Gene.Dump > $TMPDIR/Exon.Number
awk '{i = 1; while (i <= $9) { print $17; i++ }}' $TMPDIR/Gene.Dump > $TMPDIR/Gene.Number

paste $TMPDIR/Exon.Chr $TMPDIR/Exon.Start $TMPDIR/Exon.End $TMPDIR/Exon.Strand $TMPDIR/Exon.Name $TMPDIR/Exon.Number $TMPDIR/Gene.Number > $TMPDIR/Exon.All.bed
rm $TMPDIR/Exon.Chr
rm $TMPDIR/Exon.Start
rm $TMPDIR/Exon.End
rm $TMPDIR/Exon.Strand
rm $TMPDIR/Exon.Name
rm $TMPDIR/Exon.Number
rm $TMPDIR/Gene.Number
awk '{ print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$4"\t"$7 }' $TMPDIR/Exon.All.bed | sort -k7n,7 -k2n,2 - > $TMPDIR/Tmp
mv $TMPDIR/Tmp $TMPDIR/Exon.All.bed 

awk '{ print $10 }' $TMPDIR/Gene.Dump.Inc | sed "s/.$//"  | sed -e 's/,/\n/g' > $TMPDIR/Exon.Start
awk '{ print $11 }' $TMPDIR/Gene.Dump.Inc | sed "s/.$//" | sed -e 's/,/\n/g' > $TMPDIR/Exon.End
awk '{i = 1; while (i <= $9) { print $4; i++ }}' $TMPDIR/Gene.Dump.Inc > $TMPDIR/Exon.Strand
awk '{i = 1; while (i <= $9) { print $3; i++ }}' $TMPDIR/Gene.Dump.Inc > $TMPDIR/Exon.Chr
awk '{i = 1; while (i <= $9) { print $13; i++ }}' $TMPDIR/Gene.Dump.Inc > $TMPDIR/Exon.Name
awk '{i = 1; while (i <= $9) { print i; i++ }}' $TMPDIR/Gene.Dump.Inc > $TMPDIR/Exon.Number
awk '{i = 1; while (i <= $9) { print $19; i++ }}' $TMPDIR/Gene.Dump.Inc > $TMPDIR/Gene.Number

paste $TMPDIR/Exon.Chr $TMPDIR/Exon.Start $TMPDIR/Exon.End $TMPDIR/Exon.Strand $TMPDIR/Exon.Name $TMPDIR/Exon.Number $TMPDIR/Gene.Number > $TMPDIR/Exon.Inc.bed
rm $TMPDIR/Exon.Chr
rm $TMPDIR/Exon.Start
rm $TMPDIR/Exon.End
rm $TMPDIR/Exon.Strand
rm $TMPDIR/Exon.Name
rm $TMPDIR/Exon.Number
rm $TMPDIR/Gene.Number
awk '{ print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$4"\t"$7 }' $TMPDIR/Exon.Inc.bed | sort -k7n,7 -k2n,2 - > $TMPDIR/Tmp
mv $TMPDIR/Tmp $TMPDIR/Exon.Inc.bed 

# Make intron list
awk ' NR==1{printf("%s", $0"\t");next}{print FS $0;printf("%s", $0"\t")}END{if(!(NR%2)){print ""}}' $TMPDIR/Exon.All.bed > $TMPDIR/Intron.All.bed
awk '$7 == $14 && $5 < $12 { print $1"\t"$3"\t"$9"\t"$4"\t"$5"\t"$6 }' $TMPDIR/Intron.All.bed > $TMPDIR/Tmp
mv $TMPDIR/Tmp $TMPDIR/Intron.All.bed

awk ' NR==1{printf("%s", $0"\t");next}{print FS $0;printf("%s", $0"\t")}END{if(!(NR%2)){print ""}}' $TMPDIR/Exon.Inc.bed > $TMPDIR/Intron.Inc.bed
awk '$7 == $14 && $5 < $12 { print $1"\t"$3"\t"$9"\t"$4"\t"$5"\t"$6 }' $TMPDIR/Intron.Inc.bed > $TMPDIR/Tmp
mv $TMPDIR/Tmp $TMPDIR/Intron.Inc.bed

# Finalize exon lists
cut -f1-6 $TMPDIR/Exon.All.bed > $TMPDIR/Tmp
mv $TMPDIR/Tmp $TMPDIR/Exon.All.bed 
cut -f1-6 $TMPDIR/Exon.Inc.bed > $TMPDIR/Tmp
mv $TMPDIR/Tmp $TMPDIR/Exon.Inc.bed 

# Make mRNA list
awk '{ print $20 }' $TMPDIR/mRNA.Dump  | sed "s/.$//" | sed -e 's/,/\n/g' > $TMPDIR/Block.Length
awk '{ print $22 }' $TMPDIR/mRNA.Dump  | sed "s/.$//" | sed -e 's/,/\n/g' > $TMPDIR/Block.Start
awk '{i = 1; while (i <= $19) { print $10; i++ }}' $TMPDIR/mRNA.Dump > $TMPDIR/Block.Strand
awk '{i = 1; while (i <= $19) { print $15; i++ }}' $TMPDIR/mRNA.Dump > $TMPDIR/Block.Chr
awk '{i = 1; while (i <= $19) { print i; i++ }}' $TMPDIR/mRNA.Dump > $TMPDIR/Block.Number
awk '{i = 1; while (i <= $19) { print $11; i++ }}' $TMPDIR/mRNA.Dump > $TMPDIR/Block.Name
paste $TMPDIR/Block.Chr $TMPDIR/Block.Start $TMPDIR/Block.Length $TMPDIR/Block.Strand $TMPDIR/Block.Name $TMPDIR/Block.Number > $TMPDIR/mRNA.All.bed
rm $TMPDIR/Block.Chr
rm $TMPDIR/Block.Start
rm $TMPDIR/Block.Length
rm $TMPDIR/Block.Strand
rm $TMPDIR/Block.Name
rm $TMPDIR/Block.Number
awk '{ print $1"\t"$2"\t"$2+$3"\t"$5"\t"$6"\t"$4 }' $TMPDIR/mRNA.All.bed > $TMPDIR/Tmp
mv $TMPDIR/Tmp $TMPDIR/mRNA.All.bed

cp $TMPDIR/Gene.All.bed $TMPDIR/File1.All
cp $TMPDIR/Intron.All.bed $TMPDIR/File2.All
cp $TMPDIR/Exon.All.bed $TMPDIR/File3.All
cp $TMPDIR/Gene.Inc.bed $TMPDIR/File1.Inc
cp $TMPDIR/Intron.Inc.bed $TMPDIR/File2.Inc
cp $TMPDIR/Exon.Inc.bed $TMPDIR/File3.Inc

echo -e "\t\t\tFiltering genes."

## UNIX: Make unstranded non-overlapping gene list
# Identify, subtract and remove genes that have 100% overlap with non-self
intersectBed -f 1.0 -wa -wb -a $TMPDIR/Gene.Inc.bed -b $TMPDIR/Gene.All.bed > $TMPDIR/Overlap
awk '$4 != $10 { print $0 }' $TMPDIR/Overlap > $TMPDIR/tmp
mv $TMPDIR/tmp $TMPDIR/Overlap
awk '{ print $2"_"$3"_"$4"_"$5 }' $TMPDIR/Overlap | sort | uniq > $TMPDIR/Remove
awk 'NR==FNR{a[$0];next}!($2"_"$3"_"$4"_"$5 in a)' $TMPDIR/Remove $TMPDIR/Gene.Inc.bed  > $TMPDIR/Keep
awk 'NR==FNR{a[$0];next}$2"_"$3"_"$4"_"$5 in a' $TMPDIR/Remove $TMPDIR/Gene.All.bed  > $TMPDIR/Rem
subtractBed -a $TMPDIR/Keep -b $TMPDIR/Rem > $TMPDIR/Gene.Inc.bed
rm $TMPDIR/Keep 
rm $TMPDIR/Overlap 
rm $TMPDIR/Rem
rm $TMPDIR/Remove 

# Identify and subtract non-100% overlap with non-self
intersectBed -wa -wb -a $TMPDIR/Gene.Inc.bed -b $TMPDIR/Gene.All.bed > $TMPDIR/Overlap	
awk '$4 != $10 { print $0 }' $TMPDIR/Overlap > $TMPDIR/tmp
mv $TMPDIR/tmp $TMPDIR/Overlap
awk '{ print $8"_"$9"_"$10"_"$11 }' $TMPDIR/Overlap > $TMPDIR/Overlap1
awk '{ print $2"_"$3"_"$4"_"$5 }' $TMPDIR/Overlap | cat $TMPDIR/Overlap1 - | sort | uniq > $TMPDIR/Overlap.genes
rm $TMPDIR/Overlap1
mv $TMPDIR/Overlap.genes $TMPDIR/Overlap
awk 'NR==FNR{a[$0];next}!($2"_"$3"_"$4"_"$5 in a)' $TMPDIR/Overlap $TMPDIR/Gene.Inc.bed > $TMPDIR/NoOverlap
awk 'NR==FNR{a[$0];next}$2"_"$3"_"$4"_"$5 in a' $TMPDIR/Overlap $TMPDIR/Gene.Inc.bed > $TMPDIR/HasOverlapInc
awk 'NR==FNR{a[$0];next}$2"_"$3"_"$4"_"$5 in a' $TMPDIR/Overlap $TMPDIR/Gene.All.bed > $TMPDIR/HasOverlapAll

Count=$(wc -l $TMPDIR/HasOverlapInc | awk '{ print $1}' -)
for (( i=1; i<=$Count; i++ ))
do
awk -v VAR=$i 'NR == VAR { print $0 }' $TMPDIR/HasOverlapInc  > $TMPDIR/Tmp.bed 
ID=$(awk -v VAR=$i 'NR == VAR { print $4 }' $TMPDIR/HasOverlapInc)
awk -v VAR=$ID '$4 != VAR { print $0 }' $TMPDIR/HasOverlapAll > $TMPDIR/Tmp1.bed
subtractBed -a $TMPDIR/Tmp.bed -b $TMPDIR/Tmp1.bed >> $TMPDIR/gene.subtracted.bed			
done

cat $TMPDIR/NoOverlap $TMPDIR/gene.subtracted.bed > $TMPDIR/Gene.filt.unstranded.bed
rm $TMPDIR/HasOverlapInc 
rm $TMPDIR/HasOverlapAll
rm $TMPDIR/NoOverlap 
rm $TMPDIR/Overlap 
rm $TMPDIR/Tmp.bed 
rm $TMPDIR/Tmp1.bed 
rm $TMPDIR/gene.subtracted.bed 

## UNIX: Make stranded non-overlapping gene list
# Regen original files
cp $TMPDIR/File1.All $TMPDIR/Gene.All.bed
cp $TMPDIR/File2.All $TMPDIR/Intron.All.bed
cp $TMPDIR/File3.All $TMPDIR/Exon.All.bed
cp $TMPDIR/File1.Inc $TMPDIR/Gene.Inc.bed
cp $TMPDIR/File2.Inc $TMPDIR/Intron.Inc.bed
cp $TMPDIR/File3.Inc $TMPDIR/Exon.Inc.bed

# Identify, subtract and remove genes that have 100% overlap with non-self
intersectBed -s -f 1.0 -wa -wb -a $TMPDIR/Gene.Inc.bed -b $TMPDIR/Gene.All.bed > $TMPDIR/Overlap
awk '$4 != $10 { print $0 }' $TMPDIR/Overlap > $TMPDIR/tmp
mv $TMPDIR/tmp $TMPDIR/Overlap
awk '{ print $2"_"$3"_"$4"_"$5 }' $TMPDIR/Overlap | sort | uniq > $TMPDIR/Remove
awk 'NR==FNR{a[$0];next}!($2"_"$3"_"$4"_"$5 in a)' $TMPDIR/Remove $TMPDIR/Gene.Inc.bed  > $TMPDIR/Keep
awk 'NR==FNR{a[$0];next}$2"_"$3"_"$4"_"$5 in a' $TMPDIR/Remove $TMPDIR/Gene.All.bed  > $TMPDIR/Rem
subtractBed -s -a $TMPDIR/Keep -b $TMPDIR/Rem > $TMPDIR/Gene.Inc.bed
rm $TMPDIR/Keep 
rm $TMPDIR/Overlap 
rm $TMPDIR/Rem
rm $TMPDIR/Remove 

# Identify and subtract non-100% overlap with non-self
intersectBed -s -wa -wb -a $TMPDIR/Gene.Inc.bed -b $TMPDIR/Gene.All.bed > $TMPDIR/Overlap	
awk '$4 != $10 { print $0 }' $TMPDIR/Overlap > $TMPDIR/tmp
mv $TMPDIR/tmp $TMPDIR/Overlap
awk '{ print $8"_"$9"_"$10"_"$11 }' $TMPDIR/Overlap > $TMPDIR/Overlap1
awk '{ print $2"_"$3"_"$4"_"$5 }' $TMPDIR/Overlap | cat $TMPDIR/Overlap1 - | sort | uniq > $TMPDIR/Overlap.genes
rm $TMPDIR/Overlap1
mv $TMPDIR/Overlap.genes $TMPDIR/Overlap
awk 'NR==FNR{a[$0];next}!($2"_"$3"_"$4"_"$5 in a)' $TMPDIR/Overlap $TMPDIR/Gene.Inc.bed > $TMPDIR/NoOverlap
awk 'NR==FNR{a[$0];next}$2"_"$3"_"$4"_"$5 in a' $TMPDIR/Overlap $TMPDIR/Gene.Inc.bed > $TMPDIR/HasOverlapInc
awk 'NR==FNR{a[$0];next}$2"_"$3"_"$4"_"$5 in a' $TMPDIR/Overlap $TMPDIR/Gene.All.bed > $TMPDIR/HasOverlapAll

Count=$(wc -l $TMPDIR/HasOverlapInc | awk '{ print $1}' -)
for (( i=1; i<=$Count; i++ ))
do
awk -v VAR=$i 'NR == VAR { print $0 }' $TMPDIR/HasOverlapInc  > $TMPDIR/Tmp.bed 
ID=$(awk -v VAR=$i 'NR == VAR { print $4 }' $TMPDIR/HasOverlapInc)
awk -v VAR=$ID '$4 != VAR { print $0 }' $TMPDIR/HasOverlapAll > $TMPDIR/Tmp1.bed
subtractBed -s -a $TMPDIR/Tmp.bed -b $TMPDIR/Tmp1.bed >> $TMPDIR/gene.subtracted.bed			
done

cat $TMPDIR/NoOverlap $TMPDIR/gene.subtracted.bed > $TMPDIR/Gene.filt.stranded.bed
rm $TMPDIR/HasOverlapInc 
rm $TMPDIR/HasOverlapAll
rm $TMPDIR/NoOverlap 
rm $TMPDIR/Overlap 
rm $TMPDIR/Tmp.bed 
rm $TMPDIR/Tmp1.bed 
rm $TMPDIR/gene.subtracted.bed 

echo -e "\t\t\tFiltering exons."

## UNIX: Make unstranded non-overlapping exon list
# Regen original files
cp $TMPDIR/File1.All $TMPDIR/Gene.All.bed
cp $TMPDIR/File2.All $TMPDIR/Intron.All.bed
cp $TMPDIR/File3.All $TMPDIR/Exon.All.bed
cp $TMPDIR/File1.Inc $TMPDIR/Gene.Inc.bed
cp $TMPDIR/File2.Inc $TMPDIR/Intron.Inc.bed
cp $TMPDIR/File3.Inc $TMPDIR/Exon.Inc.bed

## UNIX: Make unstranded non-overlapping exon list
intersectBed -f 1.0 -wa -wb -a $TMPDIR/Exon.Inc.bed -b $TMPDIR/Exon.All.bed > $TMPDIR/Overlap
awk '$4 != $10 { print $0 }' $TMPDIR/Overlap > $TMPDIR/tmp
mv $TMPDIR/tmp $TMPDIR/Overlap
awk '{ print $2"_"$3"_"$4"_"$5 }' $TMPDIR/Overlap | sort | uniq > $TMPDIR/Remove
awk 'NR==FNR{a[$0];next}!($2"_"$3"_"$4"_"$5 in a)' $TMPDIR/Remove $TMPDIR/Exon.Inc.bed  > $TMPDIR/Keep
awk 'NR==FNR{a[$0];next}$2"_"$3"_"$4"_"$5 in a' $TMPDIR/Remove $TMPDIR/Exon.All.bed  > $TMPDIR/Rem
subtractBed -a $TMPDIR/Keep -b $TMPDIR/Rem > $TMPDIR/Exon.Inc.bed
rm $TMPDIR/Keep 
rm $TMPDIR/Overlap 
rm $TMPDIR/Rem
rm $TMPDIR/Remove 

# Identify and subtract non-100% overlap with non-self
intersectBed -wa -wb -a $TMPDIR/Exon.Inc.bed -b $TMPDIR/Exon.All.bed > $TMPDIR/Overlap	
awk '$4 != $10 { print $0 }' $TMPDIR/Overlap > $TMPDIR/tmp
mv $TMPDIR/tmp $TMPDIR/Overlap
awk '{ print $8"_"$9"_"$10"_"$11 }' $TMPDIR/Overlap > $TMPDIR/Overlap1
awk '{ print $2"_"$3"_"$4"_"$5 }' $TMPDIR/Overlap | cat $TMPDIR/Overlap1 - | sort | uniq > $TMPDIR/Overlap.genes
rm $TMPDIR/Overlap1
mv $TMPDIR/Overlap.genes $TMPDIR/Overlap
awk 'NR==FNR{a[$0];next}!($2"_"$3"_"$4"_"$5 in a)' $TMPDIR/Overlap $TMPDIR/Exon.Inc.bed > $TMPDIR/NoOverlap
awk 'NR==FNR{a[$0];next}$2"_"$3"_"$4"_"$5 in a' $TMPDIR/Overlap $TMPDIR/Exon.Inc.bed > $TMPDIR/HasOverlapInc
awk 'NR==FNR{a[$0];next}$2"_"$3"_"$4"_"$5 in a' $TMPDIR/Overlap $TMPDIR/Exon.All.bed > $TMPDIR/HasOverlapAll

Count=$(wc -l $TMPDIR/HasOverlapInc | awk '{ print $1}' -)
for (( i=1; i<=$Count; i++ ))
do
awk -v VAR=$i 'NR == VAR { print $0 }' $TMPDIR/HasOverlapInc  > $TMPDIR/Tmp.bed 
ID=$(awk -v VAR=$i 'NR == VAR { print $4 }' $TMPDIR/HasOverlapInc)
awk -v VAR=$ID '$4 != VAR { print $0 }' $TMPDIR/HasOverlapAll > $TMPDIR/Tmp1.bed
subtractBed -a $TMPDIR/Tmp.bed -b $TMPDIR/Tmp1.bed >> $TMPDIR/exon.subtracted.bed			
done

cat $TMPDIR/NoOverlap $TMPDIR/exon.subtracted.bed > $TMPDIR/Exon.filt.unstranded.bed
rm $TMPDIR/HasOverlapInc 
rm $TMPDIR/HasOverlapAll
rm $TMPDIR/NoOverlap 
rm $TMPDIR/Overlap 
rm $TMPDIR/Tmp.bed 
rm $TMPDIR/Tmp1.bed 
rm $TMPDIR/exon.subtracted.bed 

## UNIX: Make stranded non-overlapping exon list
# Regen original files
cp $TMPDIR/File1.All $TMPDIR/Gene.All.bed
cp $TMPDIR/File2.All $TMPDIR/Intron.All.bed
cp $TMPDIR/File3.All $TMPDIR/Exon.All.bed
cp $TMPDIR/File1.Inc $TMPDIR/Gene.Inc.bed
cp $TMPDIR/File2.Inc $TMPDIR/Intron.Inc.bed
cp $TMPDIR/File3.Inc $TMPDIR/Exon.Inc.bed

# Identify, subtract and remove genes that have 100% overlap with non-self
intersectBed -s -f 1.0 -wa -wb -a $TMPDIR/Exon.Inc.bed -b $TMPDIR/Exon.All.bed > $TMPDIR/Overlap
awk '$4 != $10 { print $0 }' $TMPDIR/Overlap > $TMPDIR/tmp
mv $TMPDIR/tmp $TMPDIR/Overlap
awk '{ print $2"_"$3"_"$4"_"$5 }' $TMPDIR/Overlap | sort | uniq > $TMPDIR/Remove
awk 'NR==FNR{a[$0];next}!($2"_"$3"_"$4"_"$5 in a)' $TMPDIR/Remove $TMPDIR/Exon.Inc.bed  > $TMPDIR/Keep
awk 'NR==FNR{a[$0];next}$2"_"$3"_"$4"_"$5 in a' $TMPDIR/Remove $TMPDIR/Exon.All.bed  > $TMPDIR/Rem
subtractBed -s -a $TMPDIR/Keep -b $TMPDIR/Rem > $TMPDIR/Exon.Inc.bed
rm $TMPDIR/Keep 
rm $TMPDIR/Overlap 
rm $TMPDIR/Rem
rm $TMPDIR/Remove 

# Identify and subtract non-100% overlap with non-self
intersectBed -s -wa -wb -a $TMPDIR/Exon.Inc.bed -b $TMPDIR/Exon.All.bed > $TMPDIR/Overlap	
awk '$4 != $10 { print $0 }' $TMPDIR/Overlap > $TMPDIR/tmp
mv $TMPDIR/tmp $TMPDIR/Overlap
awk '{ print $8"_"$9"_"$10"_"$11 }' $TMPDIR/Overlap > $TMPDIR/Overlap1
awk '{ print $2"_"$3"_"$4"_"$5 }' $TMPDIR/Overlap | cat $TMPDIR/Overlap1 - | sort | uniq > $TMPDIR/Overlap.genes
rm $TMPDIR/Overlap1
mv $TMPDIR/Overlap.genes $TMPDIR/Overlap
awk 'NR==FNR{a[$0];next}!($2"_"$3"_"$4"_"$5 in a)' $TMPDIR/Overlap $TMPDIR/Exon.Inc.bed > $TMPDIR/NoOverlap
awk 'NR==FNR{a[$0];next}$2"_"$3"_"$4"_"$5 in a' $TMPDIR/Overlap $TMPDIR/Exon.Inc.bed > $TMPDIR/HasOverlapInc
awk 'NR==FNR{a[$0];next}$2"_"$3"_"$4"_"$5 in a' $TMPDIR/Overlap $TMPDIR/Exon.All.bed > $TMPDIR/HasOverlapAll

Count=$(wc -l $TMPDIR/HasOverlapInc | awk '{ print $1}' -)
for (( i=1; i<=$Count; i++ ))
do
awk -v VAR=$i 'NR == VAR { print $0 }' $TMPDIR/HasOverlapInc  > $TMPDIR/Tmp.bed 
ID=$(awk -v VAR=$i 'NR == VAR { print $4 }' $TMPDIR/HasOverlapInc)
awk -v VAR=$ID '$4 != VAR { print $0 }' $TMPDIR/HasOverlapAll > $TMPDIR/Tmp1.bed
subtractBed -s -a $TMPDIR/Tmp.bed -b $TMPDIR/Tmp1.bed >> $TMPDIR/exon.subtracted.bed			
done

cat $TMPDIR/NoOverlap $TMPDIR/exon.subtracted.bed > $TMPDIR/Exon.filt.stranded.bed
rm $TMPDIR/HasOverlapInc 
rm $TMPDIR/HasOverlapAll
rm $TMPDIR/NoOverlap 
rm $TMPDIR/Overlap 
rm $TMPDIR/Tmp.bed 
rm $TMPDIR/Tmp1.bed 
rm $TMPDIR/exon.subtracted.bed 

echo -e "\t\t\tFiltering introns."

## UNIX: Make unstranded non-overlapping intron list
# Regen original files
cp $TMPDIR/File1.All $TMPDIR/Gene.All.bed
cp $TMPDIR/File2.All $TMPDIR/Intron.All.bed
cp $TMPDIR/File3.All $TMPDIR/Exon.All.bed
cp $TMPDIR/File1.Inc $TMPDIR/Gene.Inc.bed
cp $TMPDIR/File2.Inc $TMPDIR/Intron.Inc.bed
cp $TMPDIR/File3.Inc $TMPDIR/Exon.Inc.bed

## UNIX: Make unstranded non-overlapping intron list
intersectBed -f 1.0 -wa -wb -a $TMPDIR/Intron.Inc.bed -b $TMPDIR/Intron.All.bed > $TMPDIR/Overlap
awk '$4 != $10 { print $0 }' $TMPDIR/Overlap > $TMPDIR/tmp
mv $TMPDIR/tmp $TMPDIR/Overlap
awk '{ print $2"_"$3"_"$4"_"$5 }' $TMPDIR/Overlap | sort | uniq > $TMPDIR/Remove
awk 'NR==FNR{a[$0];next}!($2"_"$3"_"$4"_"$5 in a)' $TMPDIR/Remove $TMPDIR/Intron.Inc.bed  > $TMPDIR/Keep
awk 'NR==FNR{a[$0];next}$2"_"$3"_"$4"_"$5 in a' $TMPDIR/Remove $TMPDIR/Intron.All.bed  > $TMPDIR/Rem
subtractBed -a $TMPDIR/Keep -b $TMPDIR/Rem > $TMPDIR/Intron.Inc.bed
rm $TMPDIR/Keep 
rm $TMPDIR/Overlap 
rm $TMPDIR/Rem
rm $TMPDIR/Remove 

# Identify and subtract non-100% overlap with non-self
intersectBed -wa -wb -a $TMPDIR/Intron.Inc.bed -b $TMPDIR/Intron.All.bed > $TMPDIR/Overlap	
awk '$4 != $10 { print $0 }' $TMPDIR/Overlap > $TMPDIR/tmp
mv $TMPDIR/tmp $TMPDIR/Overlap
awk '{ print $8"_"$9"_"$10"_"$11 }' $TMPDIR/Overlap > $TMPDIR/Overlap1
awk '{ print $2"_"$3"_"$4"_"$5 }' $TMPDIR/Overlap | cat $TMPDIR/Overlap1 - | sort | uniq > $TMPDIR/Overlap.genes
rm $TMPDIR/Overlap1
mv $TMPDIR/Overlap.genes $TMPDIR/Overlap
awk 'NR==FNR{a[$0];next}!($2"_"$3"_"$4"_"$5 in a)' $TMPDIR/Overlap $TMPDIR/Intron.Inc.bed > $TMPDIR/NoOverlap
awk 'NR==FNR{a[$0];next}$2"_"$3"_"$4"_"$5 in a' $TMPDIR/Overlap $TMPDIR/Intron.Inc.bed > $TMPDIR/HasOverlapInc
awk 'NR==FNR{a[$0];next}$2"_"$3"_"$4"_"$5 in a' $TMPDIR/Overlap $TMPDIR/Intron.All.bed > $TMPDIR/HasOverlapAll

Count=$(wc -l $TMPDIR/HasOverlapInc | awk '{ print $1}' -)
for (( i=1; i<=$Count; i++ ))
do
awk -v VAR=$i 'NR == VAR { print $0 }' $TMPDIR/HasOverlapInc  > $TMPDIR/Tmp.bed 
ID=$(awk -v VAR=$i 'NR == VAR { print $4 }' $TMPDIR/HasOverlapInc)
awk -v VAR=$ID '$4 != VAR { print $0 }' $TMPDIR/HasOverlapAll > $TMPDIR/Tmp1.bed
subtractBed -a $TMPDIR/Tmp.bed -b $TMPDIR/Tmp1.bed >> $TMPDIR/intron.subtracted.bed			
done

cat $TMPDIR/NoOverlap $TMPDIR/intron.subtracted.bed > $TMPDIR/Intron.filt.unstranded.bed
rm $TMPDIR/HasOverlapInc 
rm $TMPDIR/HasOverlapAll
rm $TMPDIR/NoOverlap 
rm $TMPDIR/Overlap 
rm $TMPDIR/Tmp.bed 
rm $TMPDIR/Tmp1.bed 
rm $TMPDIR/intron.subtracted.bed 

# Subtract all exons and mRNA regions
subtractBed -a $TMPDIR/Intron.filt.unstranded.bed -b $TMPDIR/Exon.All.bed > $TMPDIR/Intron.subtracted.bed	
subtractBed -a $TMPDIR/Intron.subtracted.bed -b $TMPDIR/mRNA.All.bed > $TMPDIR/Intron.filt.unstranded.bed
rm $TMPDIR/Intron.subtracted.bed

## UNIX: Make stranded non-overlapping intron list
# Regen original files
cp $TMPDIR/File1.All $TMPDIR/Gene.All.bed
cp $TMPDIR/File2.All $TMPDIR/Intron.All.bed
cp $TMPDIR/File3.All $TMPDIR/Exon.All.bed
cp $TMPDIR/File1.Inc $TMPDIR/Gene.Inc.bed
cp $TMPDIR/File2.Inc $TMPDIR/Intron.Inc.bed
cp $TMPDIR/File3.Inc $TMPDIR/Exon.Inc.bed

## UNIX: Make stranded non-overlapping intron list
intersectBed -s -f 1.0 -wa -wb -a $TMPDIR/Intron.Inc.bed -b $TMPDIR/Intron.All.bed > $TMPDIR/Overlap
awk '$4 != $10 { print $0 }' $TMPDIR/Overlap > $TMPDIR/tmp
mv $TMPDIR/tmp $TMPDIR/Overlap
awk '{ print $2"_"$3"_"$4"_"$5 }' $TMPDIR/Overlap | sort | uniq > $TMPDIR/Remove
awk 'NR==FNR{a[$0];next}!($2"_"$3"_"$4"_"$5 in a)' $TMPDIR/Remove $TMPDIR/Intron.Inc.bed  > $TMPDIR/Keep
awk 'NR==FNR{a[$0];next}$2"_"$3"_"$4"_"$5 in a' $TMPDIR/Remove $TMPDIR/Intron.All.bed  > $TMPDIR/Rem
subtractBed -s -a $TMPDIR/Keep -b $TMPDIR/Rem > $TMPDIR/Intron.Inc.bed
rm $TMPDIR/Keep 
rm $TMPDIR/Overlap 
rm $TMPDIR/Rem
rm $TMPDIR/Remove 

# Identify and subtract non-100% overlap with non-self
intersectBed -s -wa -wb -a $TMPDIR/Intron.Inc.bed -b $TMPDIR/Intron.All.bed > $TMPDIR/Overlap	
awk '$4 != $10 { print $0 }' $TMPDIR/Overlap > $TMPDIR/tmp
mv $TMPDIR/tmp $TMPDIR/Overlap
awk '{ print $8"_"$9"_"$10"_"$11 }' $TMPDIR/Overlap > $TMPDIR/Overlap1
awk '{ print $2"_"$3"_"$4"_"$5 }' $TMPDIR/Overlap | cat $TMPDIR/Overlap1 - | sort | uniq > $TMPDIR/Overlap.genes
rm $TMPDIR/Overlap1
mv $TMPDIR/Overlap.genes $TMPDIR/Overlap
awk 'NR==FNR{a[$0];next}!($2"_"$3"_"$4"_"$5 in a)' $TMPDIR/Overlap $TMPDIR/Intron.Inc.bed > $TMPDIR/NoOverlap
awk 'NR==FNR{a[$0];next}$2"_"$3"_"$4"_"$5 in a' $TMPDIR/Overlap $TMPDIR/Intron.Inc.bed > $TMPDIR/HasOverlapInc
awk 'NR==FNR{a[$0];next}$2"_"$3"_"$4"_"$5 in a' $TMPDIR/Overlap $TMPDIR/Intron.All.bed > $TMPDIR/HasOverlapAll

Count=$(wc -l $TMPDIR/HasOverlapInc | awk '{ print $1}' -)
for (( i=1; i<=$Count; i++ ))
do
awk -v VAR=$i 'NR == VAR { print $0 }' $TMPDIR/HasOverlapInc  > $TMPDIR/Tmp.bed 
ID=$(awk -v VAR=$i 'NR == VAR { print $4 }' $TMPDIR/HasOverlapInc)
awk -v VAR=$ID '$4 != VAR { print $0 }' $TMPDIR/HasOverlapAll > $TMPDIR/Tmp1.bed
subtractBed -s -a $TMPDIR/Tmp.bed -b $TMPDIR/Tmp1.bed >> $TMPDIR/intron.subtracted.bed			
done

cat $TMPDIR/NoOverlap $TMPDIR/intron.subtracted.bed > $TMPDIR/Intron.filt.stranded.bed
rm $TMPDIR/HasOverlapInc 
rm $TMPDIR/HasOverlapAll
rm $TMPDIR/NoOverlap 
rm $TMPDIR/Overlap 
rm $TMPDIR/Tmp.bed 
rm $TMPDIR/Tmp1.bed 
rm $TMPDIR/intron.subtracted.bed 

# Subtract all exons and mRNA regions
subtractBed -s -a $TMPDIR/Intron.filt.stranded.bed -b $TMPDIR/Exon.All.bed > $TMPDIR/Intron.subtracted.bed	
subtractBed -s -a $TMPDIR/Intron.subtracted.bed -b $TMPDIR/mRNA.All.bed > $TMPDIR/Intron.filt.stranded.bed
rm $TMPDIR/Intron.subtracted.bed

echo -e "\t\tGetting information.";
VARIABLE=2
Rscript $INSTALLDIR"/home/lschulte/iRNA-v1.1/bin/Filter_BED.R" $TMPDIR $VARIABLE

sort -k1,1 -k3n,3 $TMPDIR/Intron.filt.unstranded.bed > $TMPDIR/Tmp
mv $TMPDIR/Tmp $TMPDIR/Intron.filt.unstranded.bed 
awk '{ print $1"\t"$6 }' $TMPDIR/Intron.filt.unstranded.bed | uniq | awk '{i = 1; while (i <= $2) { print $1"_intron_"i; i++ }}' - > $TMPDIR/Tmp
awk '{ print $2"\t"$3"\t"$4"\t"$5 }' $TMPDIR/Intron.filt.unstranded.bed | paste $TMPDIR/Tmp - | awk '{ print $2"\t"$3"\t"$4"\t"$1"\t0\t"$5 }' - > $TMPDIR/Tmp2
mv $TMPDIR/Tmp2 $TMPDIR/Intron.filt.unstranded.bed 

sort -k1,1 -k3n,3 $TMPDIR/Intron.filt.stranded.bed > $TMPDIR/Tmp
mv $TMPDIR/Tmp $TMPDIR/Intron.filt.stranded.bed 
awk '{ print $1"\t"$6 }' $TMPDIR/Intron.filt.stranded.bed | uniq | awk '{i = 1; while (i <= $2) { print $1"_intron_"i; i++ }}' - > $TMPDIR/Tmp
awk '{ print $2"\t"$3"\t"$4"\t"$5 }' $TMPDIR/Intron.filt.stranded.bed | paste $TMPDIR/Tmp - | awk '{ print $2"\t"$3"\t"$4"\t"$1"\t0\t"$5 }' - > $TMPDIR/Tmp2
mv $TMPDIR/Tmp2 $TMPDIR/Intron.filt.stranded.bed 

sort -k1,1 -k3n,3 $TMPDIR/Exon.filt.unstranded.bed > $TMPDIR/Tmp
mv $TMPDIR/Tmp $TMPDIR/Exon.filt.unstranded.bed 
awk '{ print $1"\t"$6 }' $TMPDIR/Exon.filt.unstranded.bed | uniq | awk '{i = 1; while (i <= $2) { print $1"_exon_"i; i++ }}' - > $TMPDIR/Tmp
awk '{ print $2"\t"$3"\t"$4"\t"$5 }' $TMPDIR/Exon.filt.unstranded.bed | paste $TMPDIR/Tmp - | awk '{ print $2"\t"$3"\t"$4"\t"$1"\t0\t"$5 }' - > $TMPDIR/Tmp2
mv $TMPDIR/Tmp2 $TMPDIR/Exon.filt.unstranded.bed 

sort -k1,1 -k3n,3 $TMPDIR/Exon.filt.stranded.bed > $TMPDIR/Tmp
mv $TMPDIR/Tmp $TMPDIR/Exon.filt.stranded.bed 
awk '{ print $1"\t"$6 }' $TMPDIR/Exon.filt.stranded.bed | uniq | awk '{i = 1; while (i <= $2) { print $1"_exon_"i; i++ }}' - > $TMPDIR/Tmp
awk '{ print $2"\t"$3"\t"$4"\t"$5 }' $TMPDIR/Exon.filt.stranded.bed | paste $TMPDIR/Tmp - | awk '{ print $2"\t"$3"\t"$4"\t"$1"\t0\t"$5 }' - > $TMPDIR/Tmp2
mv $TMPDIR/Tmp2 $TMPDIR/Exon.filt.stranded.bed 

awk '{ print $2"\t"$3"\t"$4"\t"$1"\t0\t"$5 }' $TMPDIR/Gene.filt.unstranded.bed > $TMPDIR/Tmp
mv $TMPDIR/Tmp $TMPDIR/Gene.filt.unstranded.bed
awk '{ print $2"\t"$3"\t"$4"\t"$1"\t0\t"$5 }' $TMPDIR/Gene.filt.stranded.bed > $TMPDIR/Tmp
mv $TMPDIR/Tmp $TMPDIR/Gene.filt.stranded.bed


awk '$4 == "-" { print $3"\t"$6-500"\t"$6+1000"\t"$2"\t0\t"$4 }' $TMPDIR/Gene.Dump > $TMPDIR/Minus
awk '$4 == "+" { print $3"\t"$5-1000"\t"$5+500"\t"$2"\t0\t"$4 }' $TMPDIR/Gene.Dump > $TMPDIR/Plus
cat $TMPDIR/Minus $TMPDIR/Plus > $TMPDIR/Promoter
awk '{$2=($2<0)?0:$2}1' OFS='\t' $TMPDIR/Promoter > $TMPDIR/Tmp
mv $TMPDIR/Tmp $TMPDIR/Promoter
rm $TMPDIR/Plus
rm $TMPDIR/Minus

cp $TMPDIR/Gene.filt.unstranded.bed $TMPDIR/Pol.filt.unstranded.bed
cp $TMPDIR/Gene.filt.stranded.bed $TMPDIR/Pol.filt.stranded.bed

awk '{$2=($2<0)?0:$2}1' OFS='\t' $TMPDIR/Pol.filt.stranded.bed > $TMPDIR/Tmp
mv $TMPDIR/Tmp $TMPDIR/Pol.filt.stranded.bed
awk '{$2=($2<0)?0:$2}1' OFS='\t' $TMPDIR/Pol.filt.unstranded.bed > $TMPDIR/Tmp
mv $TMPDIR/Tmp $TMPDIR/Pol.filt.unstranded.bed
awk '{$2=($2<0)?0:$2}1' OFS='\t' $TMPDIR/Promoter > $TMPDIR/Tmp
mv $TMPDIR/Tmp $TMPDIR/Promoter

subtractBed -a $TMPDIR/Pol.filt.unstranded.bed -b $TMPDIR/Promoter > $TMPDIR/Tmp
mv $TMPDIR/Tmp $TMPDIR/Pol.filt.unstranded.bed 

subtractBed -a $TMPDIR/Pol.filt.stranded.bed -b $TMPDIR/Promoter > $TMPDIR/Tmp
mv $TMPDIR/Tmp $TMPDIR/Pol.filt.stranded.bed 
rm $TMPDIR/Promoter


echo -e "\t\tConverting to GTF format and adding to script";

genome="AGPv4";

./bedToGenePred $TMPDIR/Gene.filt.stranded.bed stdout | ./genePredToGtf file stdin $TMPDIR/$genome.gene.strand.gtf
./bedToGenePred $TMPDIR/Gene.filt.unstranded.bed stdout | ./genePredToGtf file stdin $TMPDIR/$genome.gene.unstranded.gtf
./bedToGenePred $TMPDIR/Exon.filt.stranded.bed stdout | ./genePredToGtf file stdin $TMPDIR/$genome.exon.strand.gtf
./bedToGenePred $TMPDIR/Exon.filt.unstranded.bed stdout | ./genePredToGtf file stdin $TMPDIR/$genome.exon.unstranded.gtf
./bedToGenePred $TMPDIR/Intron.filt.stranded.bed stdout | ./genePredToGtf file stdin $TMPDIR/$genome.intron.strand.gtf
./bedToGenePred $TMPDIR/Intron.filt.unstranded.bed stdout | ./genePredToGtf file stdin $TMPDIR/$genome.intron.unstranded.gtf
./bedToGenePred $TMPDIR/Pol.filt.stranded.bed stdout | ./genePredToGtf file stdin $TMPDIR/$genome.pol.strand.gtf
./bedToGenePred $TMPDIR/Pol.filt.unstranded.bed stdout | ./genePredToGtf file stdin $TMPDIR/$genome.pol.unstranded.gtf

echo -e "RefSeq\tChr\tStrand\tStart\tEnd\tSymbol" > $TMPDIR/Tmp
awk '{ print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$13 }' $TMPDIR/Gene.Dump.Inc | cat $TMPDIR/Tmp - > $TMPDIR/$genome.annotation.txt 
rm $TMPDIR/Tmp

#rm $TMPDIR/*.bed
rm $TMPDIR/*.Dump
rm $TMPDIR/*.All
rm $TMPDIR/*.Inc

ANNDIR=$INSTALLDIR"/home/lschulte/iRNA-v1.1/ann/"
DATADIR=$INSTALLDIR"/home/lschulte/iRNA-v1.1/data/"
GENOMEFILE=$INSTALLDIR"/home/lschulte/iRNA-v1.1/data/Genomes"

mv $TMPDIR/$genome.pol.strand.gtf $DATADIR
mv $TMPDIR/$genome.pol.unstranded.gtf $DATADIR
mv $TMPDIR/$genome.gene.strand.gtf $DATADIR
mv $TMPDIR/$genome.gene.unstranded.gtf $DATADIR
mv $TMPDIR/$genome.exon.strand.gtf $DATADIR
mv $TMPDIR/$genome.exon.unstranded.gtf $DATADIR
mv $TMPDIR/$genome.intron.strand.gtf $DATADIR
mv $TMPDIR/$genome.intron.unstranded.gtf $DATADIR
mv $TMPDIR/$genome.annotation.txt $ANNDIR

echo $genome >> $GENOMEFILE
