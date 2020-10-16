## $1 reads1 full path
## $2 out dir
## $3 thread
## $4 insert size
## $5 reads2 full path
## $6 mapping quality
## $7 identity
## $8 num of PEs

# assembly
if [ ! -d "$2" ]; then
    mkdir $2
fi
cd $2

if [ ! -d "megahit" ]; then
    megahit -1 $1 -2 $5 -t $3 -o megahit > megahit.log 2>&1
fi

if [ ! -f "final.contigs.fa" ]; then
    # seqkit seq -m 1000 megahit/final.contigs.fa > final.contigs.fa
    cp megahit/final.contigs.fa ./
fi

# maxbin
start=`date +%s`
if [ ! -d "maxbin" ]; then
    mkdir maxbin
    run_MaxBin.pl -contig final.contigs.fa -reads $1 -reads2 $5 -out maxbin/maxbin -thread $3 > maxbin.log 2>&1
fi
end=`date +%s`
dif=$[ end - start ]
echo "Time of maxbin: "$dif"s"

# ground truth alignment
if [ ! -d "align_ref" ]; then
    mkdir align_ref; cd align_ref
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/1.fna ../final.contigs.fa | samtools sort -@ 128 -o 1.bam" > align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/2.fna ../final.contigs.fa | samtools sort -@ 128 -o 2.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/3.fna ../final.contigs.fa | samtools sort -@ 128 -o 3.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/4.fna ../final.contigs.fa | samtools sort -@ 128 -o 4.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/5.fna ../final.contigs.fa | samtools sort -@ 128 -o 5.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/6.fna ../final.contigs.fa | samtools sort -@ 128 -o 6.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/7.fna ../final.contigs.fa | samtools sort -@ 128 -o 7.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/8.fna ../final.contigs.fa | samtools sort -@ 128 -o 8.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/9.fna ../final.contigs.fa | samtools sort -@ 128 -o 9.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/10.fna ../final.contigs.fa | samtools sort -@ 128 -o 10.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/11.fna ../final.contigs.fa | samtools sort -@ 128 -o 11.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/12.fna ../final.contigs.fa | samtools sort -@ 128 -o 12.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/13.fna ../final.contigs.fa | samtools sort -@ 128 -o 13.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/14.fna ../final.contigs.fa | samtools sort -@ 128 -o 14.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/15.fna ../final.contigs.fa | samtools sort -@ 128 -o 15.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/16.fna ../final.contigs.fa | samtools sort -@ 128 -o 16.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/17.fna ../final.contigs.fa | samtools sort -@ 128 -o 17.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/18.fna ../final.contigs.fa | samtools sort -@ 128 -o 18.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/19.fna ../final.contigs.fa | samtools sort -@ 128 -o 19.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/20.fna ../final.contigs.fa | samtools sort -@ 128 -o 20.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/21.fna ../final.contigs.fa | samtools sort -@ 128 -o 21.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/22.fna ../final.contigs.fa | samtools sort -@ 128 -o 22.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/23.fna ../final.contigs.fa | samtools sort -@ 128 -o 23.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/24.fna ../final.contigs.fa | samtools sort -@ 128 -o 24.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/25.fna ../final.contigs.fa | samtools sort -@ 128 -o 25.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/26.fna ../final.contigs.fa | samtools sort -@ 128 -o 26.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/27.fna ../final.contigs.fa | samtools sort -@ 128 -o 27.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/28.fna ../final.contigs.fa | samtools sort -@ 128 -o 28.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/29.fna ../final.contigs.fa | samtools sort -@ 128 -o 29.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/30.fna ../final.contigs.fa | samtools sort -@ 128 -o 30.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/31.fna ../final.contigs.fa | samtools sort -@ 128 -o 31.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/32.fna ../final.contigs.fa | samtools sort -@ 128 -o 32.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/33.fna ../final.contigs.fa | samtools sort -@ 128 -o 33.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/34.fna ../final.contigs.fa | samtools sort -@ 128 -o 34.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/35.fna ../final.contigs.fa | samtools sort -@ 128 -o 35.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/36.fna ../final.contigs.fa | samtools sort -@ 128 -o 36.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/37.fna ../final.contigs.fa | samtools sort -@ 128 -o 37.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/38.fna ../final.contigs.fa | samtools sort -@ 128 -o 38.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/39.fna ../final.contigs.fa | samtools sort -@ 128 -o 39.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/40.fna ../final.contigs.fa | samtools sort -@ 128 -o 40.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/41.fna ../final.contigs.fa | samtools sort -@ 128 -o 41.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/42.fna ../final.contigs.fa | samtools sort -@ 128 -o 42.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/43.fna ../final.contigs.fa | samtools sort -@ 128 -o 43.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/44.fna ../final.contigs.fa | samtools sort -@ 128 -o 44.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/45.fna ../final.contigs.fa | samtools sort -@ 128 -o 45.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/46.fna ../final.contigs.fa | samtools sort -@ 128 -o 46.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/47.fna ../final.contigs.fa | samtools sort -@ 128 -o 47.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/48.fna ../final.contigs.fa | samtools sort -@ 128 -o 48.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/49.fna ../final.contigs.fa | samtools sort -@ 128 -o 49.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/50.fna ../final.contigs.fa | samtools sort -@ 128 -o 50.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/51.fna ../final.contigs.fa | samtools sort -@ 128 -o 51.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/52.fna ../final.contigs.fa | samtools sort -@ 128 -o 52.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/53.fna ../final.contigs.fa | samtools sort -@ 128 -o 53.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/54.fna ../final.contigs.fa | samtools sort -@ 128 -o 54.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/55.fna ../final.contigs.fa | samtools sort -@ 128 -o 55.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/56.fna ../final.contigs.fa | samtools sort -@ 128 -o 56.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/57.fna ../final.contigs.fa | samtools sort -@ 128 -o 57.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/58.fna ../final.contigs.fa | samtools sort -@ 128 -o 58.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/59.fna ../final.contigs.fa | samtools sort -@ 128 -o 59.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/60.fna ../final.contigs.fa | samtools sort -@ 128 -o 60.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/61.fna ../final.contigs.fa | samtools sort -@ 128 -o 61.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/62.fna ../final.contigs.fa | samtools sort -@ 128 -o 62.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/63.fna ../final.contigs.fa | samtools sort -@ 128 -o 63.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/64Mock/reference/64.fna ../final.contigs.fa | samtools sort -@ 128 -o 64.bam" >> align.sh
    bash align.sh > align.log 2>&1
    ~/code/assembler/align_ref_64/parse_alignment > answer.txt
    cd ../
fi

# bwa
if [ ! -f "align.megahit.bam" ]; then
    bwa index final.contigs.fa
    bwa mem -t 100 final.contigs.fa $1 $5 | samtools sort -@ 100 -o align.megahit.bam
fi

start=`date +%s`
if [ ! -d "final.contigs.fa.metabat-bins1500" ]; then
    runMetaBat.sh -t $3 -m 1500 final.contigs.fa align.megahit.bam
fi
end=`date +%s`
dif=$[ end - start ]
echo "Time of metabat: "$dif"s"

start=`date +%s`
if [ ! -d "concoct_output" ]; then
    samtools index -@ 200 align.megahit.bam
    cut_up_fasta.py final.contigs.fa -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa
    concoct_coverage_table.py contigs_10K.bed align.megahit.bam > coverage_table.tsv
    concoct --composition_file contigs_10K.fa --coverage_file coverage_table.tsv -b concoct_output/ --threads 200
    merge_cutup_clustering.py concoct_output/clustering_gt1000.csv > concoct_output/clustering_merged.csv
    extract_fasta_bins.py final.contigs.fa concoct_output/clustering_merged.csv --output_path concoct_output/
fi
end=`date +%s`
dif=$[ end - start ]
echo "Time of concort: "$dif"s"

start=`date +%s`
if [ ! -d "solidbin" ]; then
    source ~/.bashrc
    conda activate solidbin
    ln -s /home/comp/zmzhang/software/SolidBin/auxiliary auxiliary
    mkdir input_solidbin
    mkdir input_solidbin/sr
    mkdir input_solidbin/pb
    ln -s $1 input_solidbin/sr/reads_1.fq
    ln -s $5 input_solidbin/sr/reads_2.fq
    # cp $1 $5 input_solidbin/sr
    cp final.contigs.fa input_solidbin
    python ~/software/SolidBin/scripts/gen_kmer.py final.contigs.fa 1000 4
    bash ~/software/SolidBin/scripts/gen_cov.sh
    python ~/software/SolidBin/SolidBin.py --contig_file final.contigs.fa --coverage_profiles input_solidbin/coverage_new.tsv --composition_profiles kmer_4_f1000.csv --use_sfs --output solidbin_output/solidbin.txt
    mkdir solidbin
    cd solidbin_output/good_bins/
    for file in *; do
        cp $file ../../solidbin/$file.fa
    done
    cd ../../
    conda deactivate
fi
end=`date +%s`
dif=$[ end - start ]
echo "Time of solidbin: "$dif"s"

start=`date +%s`
if [ ! -d "mycc" ]; then
    cut -f 1,3 final.contigs.fa.depth.txt > mycc.depth.txt
    sed -i '1d' mycc.depth.txt
    MyCC.py final.contigs.fa -a mycc.depth.txt -meta -keep
    mkdir mycc
    cp 20*/7_AllClusters/*.fa mycc
    cp 20*/1_Data/Alias.txt mycc
    python ~/code/assembler/parse_mycc.py
fi
end=`date +%s`
dif=$[ end - start ]
echo "Time of mycc: "$dif"s"

python ~/software/GraphBin_New/support/prepResult.py --binned maxbin --assembler megahit --output maxbin
python ~/software/GraphBin_New/support/prepResult.py --binned final.contigs.fa.metabat-bins1500 --assembler megahit --output final.contigs.fa.metabat-bins1500
python ~/software/GraphBin_New/support/prepResult.py --binned concoct_output --assembler megahit --output concoct_output
python ~/software/GraphBin_New/support/prepResult.py --binned solidbin --assembler megahit --output solidbin
python ~/software/GraphBin_New/support/prepResult.py --binned mycc --assembler megahit --output mycc

# get fastg
contig=`head -n 1 final.contigs.fa`
header=${contig%%_*}
k_mer=${header##*k}

if [ ! -f "final.contigs.fastg" ]; then
    megahit_toolkit contig2fastg $k_mer final.contigs.fa > final.contigs.fastg
    megahit_toolkit contig2fastg $k_mer final.contigs.fa > final.contigs.fastg
fi
~/software/gfa1/misc/fastg2gfa final.contigs.fastg > final.contigs.gfa

# for maxbin
if [ ! -d "assembler_mq$6_idt$7_pe$8" ]; then
    mkdir assembler_mq$6_idt$7_pe$8
fi
cd assembler_mq$6_idt$7_pe$8
start=`date +%s`
if [ ! -f "assembler.log" ]; then
    ~/code/assembler/assembler ngs_megahit ../final.contigs.fa ../final.contigs.fastg ../align_ref/answer.txt ../align.megahit.bam ../maxbin/initial_contig_bins.csv $4 $k_mer $6 $7 $8> assembler.log 2>&1
fi
end=`date +%s`
dif=$[ end - start ]
echo "Time of generating graphs: "$dif"s"

# maxbin without purify
start=`date +%s`
if [ ! -d "edge_graph_graphbin" ]; then
    ~/software/GraphBin_New/graphbin --binned ../maxbin/initial_contig_bins.csv --graph ../final.contigs.gfa --output edge_graph_graphbin --assembler megahit --contigs ../final.contigs.fa > edge_graph_graphbin.log 2>&1
fi
end=`date +%s`
dif=$[ end - start ]
echo "Time of graphbin(maxbin): "$dif"s"

start=`date +%s`
if [ ! -d "final" ]; then
    python ~/code/assembler/AMGL-SEMI/AMGL-SEMI.py --remove True --assembly_graph edge_graph.txt --PE_graph PE_graph.txt --binned ../maxbin/initial_contig_bins.csv --max_iter 100 --thresh 0.00000001 > amgl.log 2>&1
    mkdir final; mv final.cluster final
fi
end=`date +%s`
dif=$[ end - start ]
echo "Time of amgl(maxbin): "$dif"s"

~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../align_ref/answer.txt ../maxbin/initial_contig_bins.csv $k_mer > evaluate_megahit.maxbin
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../align_ref/answer.txt edge_graph_graphbin/graphbin_output.csv $k_mer > evaluate_megahit.edge_graph
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../align_ref/answer.txt final/final.cluster $k_mer > evaluate_megahit.amgl

python ~/code/assembler/parse_result.py --assembler megahit --cluster edge_graph_graphbin/graphbin_output.csv --output edge_graph_graphbin --contigs ../final.contigs.fa
python ~/code/assembler/parse_result.py --assembler megahit --cluster final/final.cluster --output final --contigs ../final.contigs.fa

# for metabat
if [ ! -d "../assembler_mq$6_idt$7_pe$8_metabat" ]; then
    mkdir ../assembler_mq$6_idt$7_pe$8_metabat
fi
cd ../assembler_mq$6_idt$7_pe$8_metabat
cp ../assembler_mq$6_idt$7_pe$8/PE_graph.txt ../assembler_mq$6_idt$7_pe$8/edge_graph.txt ./

# maxbin without purify
start=`date +%s`
if [ ! -d "edge_graph_graphbin" ]; then
    ~/software/GraphBin_New/graphbin --binned ../final.contigs.fa.metabat-bins1500/initial_contig_bins.csv --graph ../final.contigs.gfa --output edge_graph_graphbin --assembler megahit --contigs ../final.contigs.fa > edge_graph_graphbin.log 2>&1
fi
end=`date +%s`
dif=$[ end - start ]
echo "Time of graphbin(metabat): "$dif"s"

start=`date +%s`
if [ ! -d "final" ]; then
    python ~/code/assembler/AMGL-SEMI/AMGL-SEMI.py --remove True --assembly_graph edge_graph.txt --PE_graph PE_graph.txt --binned ../final.contigs.fa.metabat-bins1500/initial_contig_bins.csv --max_iter 100 --thresh 0.00000001 > amgl.log 2>&1
    mkdir final; mv final.cluster final
fi
end=`date +%s`
dif=$[ end - start ]
echo "Time of amgl(metabat): "$dif"s"

~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../align_ref/answer.txt ../final.contigs.fa.metabat-bins1500/initial_contig_bins.csv $k_mer > evaluate_megahit.maxbin
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../align_ref/answer.txt edge_graph_graphbin/graphbin_output.csv $k_mer > evaluate_megahit.edge_graph
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../align_ref/answer.txt final/final.cluster $k_mer > evaluate_megahit.amgl

python ~/code/assembler/parse_result.py --assembler megahit --cluster edge_graph_graphbin/graphbin_output.csv --output edge_graph_graphbin --contigs ../final.contigs.fa
python ~/code/assembler/parse_result.py --assembler megahit --cluster final/final.cluster --output final --contigs ../final.contigs.fa

# for concoct
if [ ! -d "../assembler_mq$6_idt$7_pe$8_concoct" ]; then
    mkdir ../assembler_mq$6_idt$7_pe$8_concoct
fi
cd ../assembler_mq$6_idt$7_pe$8_concoct
cp ../assembler_mq$6_idt$7_pe$8/PE_graph.txt ../assembler_mq$6_idt$7_pe$8/edge_graph.txt ./

# maxbin without purify
start=`date +%s`
if [ ! -d "edge_graph_graphbin" ]; then
    ~/software/GraphBin_New/graphbin --binned ../concoct_output/initial_contig_bins.csv --graph ../final.contigs.gfa --output edge_graph_graphbin --assembler megahit --contigs ../final.contigs.fa > edge_graph_graphbin.log 2>&1
fi
end=`date +%s`
dif=$[ end - start ]
echo "Time of graphbin(concort): "$dif"s"

start=`date +%s`
if [ ! -d "final" ]; then
    python ~/code/assembler/AMGL-SEMI/AMGL-SEMI.py --remove True --assembly_graph edge_graph.txt --PE_graph PE_graph.txt --binned ../concoct_output/initial_contig_bins.csv --max_iter 100 --thresh 0.00000001 > amgl.log 2>&1
    mkdir final; mv final.cluster final
fi
end=`date +%s`
dif=$[ end - start ]
echo "Time of amgl(concort): "$dif"s"

~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../align_ref/answer.txt ../concoct_output/initial_contig_bins.csv $k_mer > evaluate_megahit.maxbin
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../align_ref/answer.txt edge_graph_graphbin/graphbin_output.csv $k_mer > evaluate_megahit.edge_graph
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../align_ref/answer.txt final/final.cluster $k_mer > evaluate_megahit.amgl

python ~/code/assembler/parse_result.py --assembler megahit --cluster edge_graph_graphbin/graphbin_output.csv --output edge_graph_graphbin --contigs ../final.contigs.fa
python ~/code/assembler/parse_result.py --assembler megahit --cluster final/final.cluster --output final --contigs ../final.contigs.fa

# for solidbin
if [ ! -d "../assembler_mq$6_idt$7_pe$8_solidbin" ]; then
    mkdir ../assembler_mq$6_idt$7_pe$8_solidbin
fi
cd ../assembler_mq$6_idt$7_pe$8_solidbin
cp ../assembler_mq$6_idt$7_pe$8/PE_graph.txt ../assembler_mq$6_idt$7_pe$8/edge_graph.txt ./

# maxbin without purify
start=`date +%s`
if [ ! -d "edge_graph_graphbin" ]; then
    ~/software/GraphBin_New/graphbin --binned ../solidbin/initial_contig_bins.csv --graph ../final.contigs.gfa --output edge_graph_graphbin --assembler megahit --contigs ../final.contigs.fa > edge_graph_graphbin.log 2>&1
fi
end=`date +%s`
dif=$[ end - start ]
echo "Time of graphbin(solidbin): "$dif"s"

start=`date +%s`
if [ ! -d "final" ]; then
    python ~/code/assembler/AMGL-SEMI/AMGL-SEMI.py --remove True --assembly_graph edge_graph.txt --PE_graph PE_graph.txt --binned ../solidbin/initial_contig_bins.csv --max_iter 100 --thresh 0.00000001 > amgl.log 2>&1
    mkdir final; mv final.cluster final
fi
end=`date +%s`
dif=$[ end - start ]
echo "Time of amgl(solidbin): "$dif"s"

~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../align_ref/answer.txt ../solidbin/initial_contig_bins.csv $k_mer > evaluate_megahit.maxbin
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../align_ref/answer.txt edge_graph_graphbin/graphbin_output.csv $k_mer > evaluate_megahit.edge_graph
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../align_ref/answer.txt final/final.cluster $k_mer > evaluate_megahit.amgl

python ~/code/assembler/parse_result.py --assembler megahit --cluster edge_graph_graphbin/graphbin_output.csv --output edge_graph_graphbin --contigs ../final.contigs.fa
python ~/code/assembler/parse_result.py --assembler megahit --cluster final/final.cluster --output final --contigs ../final.contigs.fa

# for mycc
if [ ! -d "../assembler_mq$6_idt$7_pe$8_mycc" ]; then
    mkdir ../assembler_mq$6_idt$7_pe$8_mycc
fi
cd ../assembler_mq$6_idt$7_pe$8_mycc
cp ../assembler_mq$6_idt$7_pe$8/PE_graph.txt ../assembler_mq$6_idt$7_pe$8/edge_graph.txt ./

# maxbin without purify
start=`date +%s`
if [ ! -d "edge_graph_graphbin" ]; then
    ~/software/GraphBin_New/graphbin --binned ../mycc/initial_contig_bins.csv --graph ../final.contigs.gfa --output edge_graph_graphbin --assembler megahit --contigs ../final.contigs.fa > edge_graph_graphbin.log 2>&1
fi
end=`date +%s`
dif=$[ end - start ]
echo "Time of graphbin(mycc): "$dif"s"

start=`date +%s`
if [ ! -d "final" ]; then
    python ~/code/assembler/AMGL-SEMI/AMGL-SEMI.py --remove True --assembly_graph edge_graph.txt --PE_graph PE_graph.txt --binned ../mycc/initial_contig_bins.csv --max_iter 100 --thresh 0.00000001 > amgl.log 2>&1
    mkdir final; mv final.cluster final
fi
end=`date +%s`
dif=$[ end - start ]
echo "Time of amgl(mycc): "$dif"s"

~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../align_ref/answer.txt ../mycc/initial_contig_bins.csv $k_mer > evaluate_megahit.maxbin
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../align_ref/answer.txt edge_graph_graphbin/graphbin_output.csv $k_mer > evaluate_megahit.edge_graph
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../align_ref/answer.txt final/final.cluster $k_mer > evaluate_megahit.amgl

python ~/code/assembler/parse_result.py --assembler megahit --cluster edge_graph_graphbin/graphbin_output.csv --output edge_graph_graphbin --contigs ../final.contigs.fa
python ~/code/assembler/parse_result.py --assembler megahit --cluster final/final.cluster --output final --contigs ../final.contigs.fa
