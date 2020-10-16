# real
nohup bash ~/code/assembler/binning_evaluate_megahit_real.sh `pwd`/sharon_1.fastq megahit_time 200 679 `pwd`/sharon_2.fastq 10 0.95 4 > megahit_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_spades_real.sh `pwd`/sharon_1.fastq metaspades_time 200 679 `pwd`/sharon_2.fastq 10 0.95 4 > metaspades_time.log 2>&1 &

# mock12
# spades
nohup bash ~/code/assembler/binning_evaluate_spades_mock12.sh `pwd`/SRR8073716.1_1.fastq metaspades_time 200 300 `pwd`/SRR8073716.1_2.fastq 10 0.95 4 > metaspades_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_spades_mock12.sh `pwd`/SRR8073716.1_1.fastq metaspades_time 200 300 `pwd`/SRR8073716.1_2.fastq 20 0.98 5 > metaspades_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_spades_mock12.sh `pwd`/SRR8073716.1_1.fastq metaspades_time 200 300 `pwd`/SRR8073716.1_2.fastq 10 0.95 3 > metaspades_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_spades_mock12.sh `pwd`/SRR8073716.1_1.fastq metaspades_time 200 300 `pwd`/SRR8073716.1_2.fastq 10 0.95 2 > metaspades_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_spades_mock12.sh `pwd`/SRR8073716.1_1.fastq metaspades_time 200 300 `pwd`/SRR8073716.1_2.fastq 20 0.98 2 > metaspades_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_spades_mock12.sh `pwd`/SRR8073716.1_1.fastq metaspades_time 200 300 `pwd`/SRR8073716.1_2.fastq 20 0.98 3 > metaspades_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_spades_mock12.sh `pwd`/SRR8073716.1_1.fastq metaspades_time 200 300 `pwd`/SRR8073716.1_2.fastq 20 0.98 4 > metaspades_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_spades_mock12.sh `pwd`/SRR8073716.1_1.fastq metaspades_time 200 300 `pwd`/SRR8073716.1_2.fastq 20 0.90 3 > metaspades_time.log 2>&1 &

# megahit
nohup bash ~/code/assembler/binning_evaluate_megahit_mock12.sh `pwd`/SRR8073716.1_1.fastq megahit_time 200 300 `pwd`/SRR8073716.1_2.fastq 20 0.98 2 > megahit_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_megahit_mock12.sh `pwd`/SRR8073716.1_1.fastq megahit_time 200 300 `pwd`/SRR8073716.1_2.fastq 10 0.95 2 > megahit_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_megahit_mock12.sh `pwd`/SRR8073716.1_1.fastq megahit_time 200 300 `pwd`/SRR8073716.1_2.fastq 20 0.98 3 > megahit_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_megahit_mock12.sh `pwd`/SRR8073716.1_1.fastq megahit_time 200 300 `pwd`/SRR8073716.1_2.fastq 20 0.98 4 > megahit_time.log 2>&1 &

# mock64
# spades
nohup bash ~/code/assembler/binning_evaluate_spades_mock64.sh `pwd`/SRR606249.1_1.fastq metaspades_time 200 200 `pwd`/SRR606249.1_2.fastq 10 0.95 3 > metaspades_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_spades_mock64.sh `pwd`/SRR606249.1_1.fastq metaspades_time 200 200 `pwd`/SRR606249.1_2.fastq 20 0.98 5 > metaspades_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_spades_mock64.sh `pwd`/SRR606249.1_1.fastq metaspades_time 200 200 `pwd`/SRR606249.1_2.fastq 10 0.95 2 > metaspades_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_spades_mock64.sh `pwd`/SRR606249.1_1.fastq metaspades_time 200 200 `pwd`/SRR606249.1_2.fastq 20 0.98 2 > metaspades_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_spades_mock64.sh `pwd`/SRR606249.1_1.fastq metaspades_time 200 200 `pwd`/SRR606249.1_2.fastq 20 0.98 3 > metaspades_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_spades_mock64.sh `pwd`/SRR606249.1_1.fastq metaspades_time 200 200 `pwd`/SRR606249.1_2.fastq 20 0.98 4 > metaspades_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_spades_mock64.sh `pwd`/SRR606249.1_1.fastq metaspades_time 200 200 `pwd`/SRR606249.1_2.fastq 10 0.98 4 > metaspades_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_spades_mock64.sh `pwd`/SRR606249.1_1.fastq metaspades_time 200 200 `pwd`/SRR606249.1_2.fastq 10 0.98 3 > metaspades_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_spades_mock64.sh `pwd`/SRR606249.1_1.fastq metaspades_time 200 200 `pwd`/SRR606249.1_2.fastq 10 0.90 3 > metaspades_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_spades_mock64.sh `pwd`/SRR606249.1_1.fastq metaspades_time 200 200 `pwd`/SRR606249.1_2.fastq 10 0.80 3 > metaspades_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_spades_mock64.sh `pwd`/SRR606249.1_1.fastq metaspades_time 200 200 `pwd`/SRR606249.1_2.fastq 10 0.50 3 > metaspades_time.log 2>&1 &

# megahit
nohup bash ~/code/assembler/binning_evaluate_megahit_mock64.sh `pwd`/SRR606249.1_1.fastq megahit_time 200 200 `pwd`/SRR606249.1_2.fastq 20 0.98 2 > megahit_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_megahit_mock64.sh `pwd`/SRR606249.1_1.fastq megahit_time 200 200 `pwd`/SRR606249.1_2.fastq 10 0.95 2 > megahit_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_megahit_mock64.sh `pwd`/SRR606249.1_1.fastq megahit_time 200 200 `pwd`/SRR606249.1_2.fastq 20 0.98 3 > megahit_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_megahit_mock64.sh `pwd`/SRR606249.1_1.fastq megahit_time 200 200 `pwd`/SRR606249.1_2.fastq 20 0.98 4 > megahit_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_megahit_mock64.sh `pwd`/SRR606249.1_1.fastq megahit_time 200 200 `pwd`/SRR606249.1_2.fastq 10 0.98 4 > megahit_time.log 2>&1 &
nohup bash ~/code/assembler/binning_evaluate_megahit_mock64.sh `pwd`/SRR606249.1_1.fastq megahit_time 200 200 `pwd`/SRR606249.1_2.fastq 10 0.98 3 > megahit_time.log 2>&1 &
