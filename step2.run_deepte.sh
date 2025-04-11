#备份
mv genome.fa.mod.EDTA.TEanno.sum genome.fa.mod.EDTA.TEanno.sum.bak
mv genome.fa.mod.EDTA.TEanno.gff3 genome.fa.mod.EDTA.TEanno.gff3.bak
mv genome.fa.mod.EDTA.final/genome.fa.mod.EDTA.TElib.fa genome.fa.mod.EDTA.final/genome.fa.mod.EDTA.TElib.fa.bak

#对EDTA构建的TElib中的unknow提取出来
grep "LTR/unknown" genome.fa.mod.EDTA.TElib.fa | sed 's/>//' | seqtk subseq genome.fa.mod.EDTA.TElib.fa - >LTR_unknown.fa
grep -v "LTR/unknown" genome.fa.mod.EDTA.TElib.fa | sed 's/>//' | seqtk subseq genome.fa.mod.EDTA.TElib.fa - >LTR_known.fa

#conda activate DeepTE
#对上一步提取出来的unknown进行重分类
python /home/zxj/software/DeepTE/DeepTE.py -i LTR_unknown.fa -sp M -m_dir /home/zxj/workspace/crab/TE/TEAnalysis/TE_crab/11.software/Metazoans_model -fam LTR 1>DeepTE.log 2>&1

#将DeepTE重新分类构建的库里面的转座子进行重命名，规范命名
sed 's/LTR\/unknown__ClassI_LTR_Copia/LTR\/Copia/' opt_DeepTE.fasta | sed 's/LTR\/unknown__ClassI_LTR_Gypsy/LTR\/Gypsy/'|sed 's/LTR\/unknown__ClassI_LTR/LTR\/unknown/' >LTR_unknown_DeepTE.fa

#合并know和unknown，重新构建TElib
cat LTR_known.fa LTR_unknown_DeepTE.fa >genome.fa.mod.EDTA.final/genome.fa.mod.EDTA.TElib.fa

#用构建好的TElib替换原本的TElib
cp genome.fa.mod.EDTA.final/genome.fa.mod.EDTA.TElib.fa genome.fa.mod.EDTA.TElib.fa
# source deactivate  
