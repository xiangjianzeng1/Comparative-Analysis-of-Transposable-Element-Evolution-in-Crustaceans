# source activate EDTA
ln -s ../../12.ref/genome.fa genome.fa

# Annotation with EDTA
nohup EDTA.pl --genome genome.fa --species others --sensitive 1 --anno 1 --threads 30 --curatedlib ../../11.software/SINE_LINE.fa & 

# source deactivate
