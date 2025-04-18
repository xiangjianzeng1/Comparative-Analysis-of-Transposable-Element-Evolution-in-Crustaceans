#Backup
mv genome.fa.mod.EDTA.TEanno.sum genome.fa.mod.EDTA.TEanno.sum.bak
mv genome.fa.mod.EDTA.TEanno.gff3 genome.fa.mod.EDTA.TEanno.gff3.bak
mv genome.fa.mod.EDTA.final/genome.fa.mod.EDTA.TElib.fa genome.fa.mod.EDTA.final/genome.fa.mod.EDTA.TElib.fa.bak

#Extract the unknown from TElib constructed by EDTA
grep "LTR/unknown" genome.fa.mod.EDTA.TElib.fa | sed 's/>//' | seqtk subseq genome.fa.mod.EDTA.TElib.fa - >LTR_unknown.fa
grep -v "LTR/unknown" genome.fa.mod.EDTA.TElib.fa | sed 's/>//' | seqtk subseq genome.fa.mod.EDTA.TElib.fa - >LTR_known.fa

#conda activate DeepTE
#Reclassify the unknown extracted in the previous step
python /home/zxj/software/DeepTE/DeepTE.py -i LTR_unknown.fa -sp M -m_dir /home/zxj/workspace/crab/TE/TEAnalysis/TE_crab/11.software/Metazoans_model -fam LTR 1>DeepTE.log 2>&1

#Rename the transposons in the library reclassified by DeepTE to standardize the naming
sed 's/LTR\/unknown__ClassI_LTR_Copia/LTR\/Copia/' opt_DeepTE.fasta | sed 's/LTR\/unknown__ClassI_LTR_Gypsy/LTR\/Gypsy/'|sed 's/LTR\/unknown__ClassI_LTR/LTR\/unknown/' >LTR_unknown_DeepTE.fa

#Merge know and unknown and rebuild TElib
cat LTR_known.fa LTR_unknown_DeepTE.fa >genome.fa.mod.EDTA.final/genome.fa.mod.EDTA.TElib.fa

#Replace the original TElib with the built TElib
cp genome.fa.mod.EDTA.final/genome.fa.mod.EDTA.TElib.fa genome.fa.mod.EDTA.TElib.fa
# source deactivate  
