mkdir 30x_from5G_spades_time 50x_from5G_spades_time 70x_from5G_spades_time 90x_from5G_spades_time 110x_from5G_spades_time
cp -r 30x_from5G_spades/spades 30x_from5G_spades/*.bam 30x_from5G_spades_time/
cp -r 50x_from5G_spades/spades 50x_from5G_spades/*.bam 50x_from5G_spades_time/
cp -r 70x_from5G_spades/spades 70x_from5G_spades/*.bam 70x_from5G_spades_time/
cp -r 90x_from5G_spades/spades 90x_from5G_spades/*.bam 90x_from5G_spades_time/
cp -r 110x_from5G_spades/spades 110x_from5G_spades/*.bam 110x_from5G_spades_time/
nohup bash ~/code/assembler/binning_evaluate_spades_simulate.sh `pwd`/10x_from5G_1.fq 10x_from5G_spades_time 200 400 `pwd`/10x_from5G_2.fq 10 0.95 2 > 10x_from5G_spades_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_spades_simulate.sh `pwd`/30x_from5G_1.fq 30x_from5G_spades_time 200 400 `pwd`/30x_from5G_2.fq 10 0.95 2 > 30x_from5G_spades_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_spades_simulate.sh `pwd`/50x_from5G_1.fq 50x_from5G_spades_time 200 400 `pwd`/50x_from5G_2.fq 10 0.95 2 > 50x_from5G_spades_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_spades_simulate.sh `pwd`/70x_from5G_1.fq 70x_from5G_spades_time 200 400 `pwd`/70x_from5G_2.fq 10 0.95 2 > 70x_from5G_spades_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_spades_simulate.sh `pwd`/90x_from5G_1.fq 90x_from5G_spades_time 200 400 `pwd`/90x_from5G_2.fq 10 0.95 2 > 90x_from5G_spades_time.log 2>&1 &

nohup bash ~/code/assembler/binning_evaluate_spades_simulate.sh `pwd`/110x_from5G_1.fq 110x_from5G_spades_time 200 400 `pwd`/110x_from5G_2.fq 10 0.95 2 > 110x_from5G_spades_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_spades_simulate.sh `pwd`/130x_from5G_1.fq 130x_from5G_spades_time 200 400 `pwd`/130x_from5G_2.fq 10 0.95 2 > 130x_from5G_spades_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_spades_simulate.sh `pwd`/150x_from5G_1.fq 150x_from5G_spades_time 200 400 `pwd`/150x_from5G_2.fq 10 0.95 2 > 150x_from5G_spades_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_spades_simulate.sh `pwd`/170x_from5G_1.fq 170x_from5G_spades_time 200 400 `pwd`/170x_from5G_2.fq 10 0.95 2 > 170x_from5G_spades_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_spades_simulate.sh `pwd`/190x_from5G_1.fq 190x_from5G_spades_time 200 400 `pwd`/190x_from5G_2.fq 10 0.95 2 > 190x_from5G_spades_time.log 2>&1 &

cp -r 30x_from5G_all/megahit 30x_from5G_all/*.bam 30x_from5G_megahit_time/
cp -r 50x_from5G_all/megahit 50x_from5G_all/*.bam 50x_from5G_megahit_time/
cp -r 70x_from5G_all/megahit 70x_from5G_all/*.bam 70x_from5G_megahit_time/
cp -r 90x_from5G_all/megahit 90x_from5G_all/*.bam 90x_from5G_megahit_time/
cp -r 110x_from5G_all/megahit 110x_from5G_all/*.bam 110x_from5G_megahit_time/
nohup bash ~/code/assembler/binning_evaluate_megahit_simulate.sh `pwd`/10x_from5G_1.fq 10x_from5G_megahit_time 200 400 `pwd`/10x_from5G_2.fq 10 0.95 2 > 10x_from5G_megahit_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_megahit_simulate.sh `pwd`/30x_from5G_1.fq 30x_from5G_megahit_time 200 400 `pwd`/30x_from5G_2.fq 10 0.95 2 > 30x_from5G_megahit_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_megahit_simulate.sh `pwd`/50x_from5G_1.fq 50x_from5G_megahit_time 200 400 `pwd`/50x_from5G_2.fq 10 0.95 2 > 50x_from5G_megahit_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_megahit_simulate.sh `pwd`/70x_from5G_1.fq 70x_from5G_megahit_time 200 400 `pwd`/70x_from5G_2.fq 10 0.95 2 > 70x_from5G_megahit_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_megahit_simulate.sh `pwd`/90x_from5G_1.fq 90x_from5G_megahit_time 200 400 `pwd`/90x_from5G_2.fq 10 0.95 2 > 90x_from5G_megahit_time.log 2>&1 &

nohup bash ~/code/assembler/binning_evaluate_megahit_simulate.sh `pwd`/110x_from5G_1.fq 110x_from5G_megahit_time 200 400 `pwd`/110x_from5G_2.fq 10 0.95 2 > 110x_from5G_megahit_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_megahit_simulate.sh `pwd`/130x_from5G_1.fq 130x_from5G_megahit_time 200 400 `pwd`/130x_from5G_2.fq 10 0.95 2 > 130x_from5G_megahit_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_megahit_simulate.sh `pwd`/150x_from5G_1.fq 150x_from5G_megahit_time 200 400 `pwd`/150x_from5G_2.fq 10 0.95 2 > 150x_from5G_megahit_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_megahit_simulate.sh `pwd`/170x_from5G_1.fq 170x_from5G_megahit_time 200 400 `pwd`/170x_from5G_2.fq 10 0.95 2 > 170x_from5G_megahit_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_megahit_simulate.sh `pwd`/190x_from5G_1.fq 190x_from5G_megahit_time 200 400 `pwd`/190x_from5G_2.fq 10 0.95 2 > 190x_from5G_megahit_time.log 2>&1 &