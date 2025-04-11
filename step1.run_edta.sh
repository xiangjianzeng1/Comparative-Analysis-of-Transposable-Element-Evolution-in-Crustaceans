# source activate EDTA
ln -s ../../12.ref/genome.fa genome.fa

# 使用EDTA注释
nohup EDTA.pl --genome genome.fa --species others --sensitive 1 --anno 1 --threads 30 --curatedlib ../../11.software/SINE_LINE.fa & 

# source deactivate
