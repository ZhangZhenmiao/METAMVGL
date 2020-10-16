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
    echo "minimap2 -t 128 -a /home/comp/zmzhang/12Mock/reference/Cohaesibacter.ES.047.fna ../final.contigs.fa | samtools sort -@ 128 -o Cohaesibacter.ES.047.bam" > align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/12Mock/reference/Halomonas.HL-4.fna ../final.contigs.fa | samtools sort -@ 128 -o Halomonas.HL-4.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/12Mock/reference/Halomonas.HL-93.fna ../final.contigs.fa | samtools sort -@ 128 -o Halomonas.HL-93.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/12Mock/reference/Marinobacter.LV10MA510-1.fna ../final.contigs.fa | samtools sort -@ 128 -o Marinobacter.LV10MA510-1.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/12Mock/reference/Marinobacter.LV10R510-8.fna ../final.contigs.fa | samtools sort -@ 128 -o Marinobacter.LV10R510-8.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/12Mock/reference/Micromonospora_coxensis.DSM.45161.fna ../final.contigs.fa | samtools sort -@ 128 -o Micromonospora_coxensis.DSM.45161.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/12Mock/reference/Micromonospora_echinaurantiaca.DSM.43904.fna ../final.contigs.fa | samtools sort -@ 128 -o icromonospora_echinaurantiaca.DSM.43904.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/12Mock/reference/Micromonospora_echinofusca.DSM.43913.fna ../final.contigs.fa | samtools sort -@ 128 -o Micromonospora_echinofusca.DSM.43913.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/12Mock/reference/Muricauda.ES.050.fna ../final.contigs.fa | samtools sort -@ 128 -o Muricauda.ES.050.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/12Mock/reference/Propionibacteriaceae.ES.041.fna ../final.contigs.fa | samtools sort -@ 128 -o Propionibacteriaceae.ES.041.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/12Mock/reference/Psychrobacter.LV10R520-6.fna ../final.contigs.fa | samtools sort -@ 128 -o Psychrobacter.LV10R520-6.bam" >> align.sh
    echo "minimap2 -t 128 -a /home/comp/zmzhang/12Mock/reference/Thioclava.ES.032.fna  ../final.contigs.fa | samtools sort -@ 128 -o Thioclava.ES.032.bam" >> align.sh
    bash align.sh > align.log 2>&1
    ~/code/assembler/align_ref_12/parse_alignment > answer.txt
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
