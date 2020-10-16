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
if [ ! -f "kraken2_0.9.tsv" ]; then
    kraken2 --db ~/software/kraken2_install/standard_db --confidence 0.9 --output kraken2_0.9.tsv -t 100 final.contigs.fa > kraken2_0.9.log 2>&1 &
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

# maxbin purify
if [ ! -d "maxbin_cleaned" ]; then
    cd maxbin
    source ~/.bashrc
    conda activate python3
    echo "magpurify phylo-markers \$1 purify/\$1" > magpurify.sh
    echo "magpurify clade-markers \$1 purify/\$1" >> magpurify.sh
    echo "magpurify conspecific \$1 purify/\$1 /home/comp/zmzhang/Uhgg/Uhgg.msh" >> magpurify.sh
    echo "magpurify tetra-freq \$1 purify/\$1" >> magpurify.sh
    echo "magpurify gc-content \$1 purify/\$1" >> magpurify.sh
    echo "magpurify known-contam \$1 purify/\$1" >> magpurify.sh
    echo "magpurify coverage \$1 purify/\$1 ../align.megahit.bam" >> magpurify.sh
    echo "magpurify clean-bin \$1 purify/\$1 purify/\$1/\$1" >> magpurify.sh
    for file in *.fasta; do
        nohup bash magpurify.sh $file > $file.log 2>&1 &
    done
    wait
    conda deactivate
    cd ../; mkdir maxbin_cleaned
    cp maxbin/purify/*/*.fasta maxbin_cleaned
fi

# megahit purify
if [ ! -d "final.contigs.fa.metabat-bins1500_cleaned" ]; then
    cd final.contigs.fa.metabat-bins1500
    source ~/.bashrc
    conda activate python3
    echo "magpurify phylo-markers \$1 purify/\$1" > magpurify.sh
    echo "magpurify clade-markers \$1 purify/\$1" >> magpurify.sh
    echo "magpurify conspecific \$1 purify/\$1 /home/comp/zmzhang/Uhgg/Uhgg.msh" >> magpurify.sh
    echo "magpurify tetra-freq \$1 purify/\$1" >> magpurify.sh
    echo "magpurify gc-content \$1 purify/\$1" >> magpurify.sh
    echo "magpurify known-contam \$1 purify/\$1" >> magpurify.sh
    echo "magpurify coverage \$1 purify/\$1 ../align.megahit.bam" >> magpurify.sh
    echo "magpurify clean-bin \$1 purify/\$1 purify/\$1/\$1" >> magpurify.sh
    for file in *.fa; do
        nohup bash magpurify.sh $file > $file.log 2>&1 &
    done
    wait
    conda deactivate
    cd ../; mkdir final.contigs.fa.metabat-bins1500_cleaned
    cp final.contigs.fa.metabat-bins1500/purify/*/*.fa final.contigs.fa.metabat-bins1500_cleaned
fi

# concoct purify
if [ ! -d "concoct_output_cleaned" ]; then
    cd concoct_output
    source ~/.bashrc
    conda activate python3
    echo "magpurify phylo-markers \$1 purify/\$1" > magpurify.sh
    echo "magpurify clade-markers \$1 purify/\$1" >> magpurify.sh
    echo "magpurify conspecific \$1 purify/\$1 /home/comp/zmzhang/Uhgg/Uhgg.msh" >> magpurify.sh
    echo "magpurify tetra-freq \$1 purify/\$1" >> magpurify.sh
    echo "magpurify gc-content \$1 purify/\$1" >> magpurify.sh
    echo "magpurify known-contam \$1 purify/\$1" >> magpurify.sh
    echo "magpurify coverage \$1 purify/\$1 ../align.megahit.bam" >> magpurify.sh
    echo "magpurify clean-bin \$1 purify/\$1 purify/\$1/\$1" >> magpurify.sh
    for file in *.fa; do
        nohup bash magpurify.sh $file > $file.log 2>&1 &
    done
    wait
    conda deactivate
    cd ../; mkdir concoct_output_cleaned
    cp concoct_output/purify/*/*.fa concoct_output_cleaned
fi

if [ ! -d "solidbin_cleaned" ]; then
    cd solidbin
    source ~/.bashrc
    conda activate python3
    echo "magpurify phylo-markers \$1 purify/\$1" > magpurify.sh
    echo "magpurify clade-markers \$1 purify/\$1" >> magpurify.sh
    echo "magpurify conspecific \$1 purify/\$1 /home/comp/zmzhang/Uhgg/Uhgg.msh" >> magpurify.sh
    echo "magpurify tetra-freq \$1 purify/\$1" >> magpurify.sh
    echo "magpurify gc-content \$1 purify/\$1" >> magpurify.sh
    echo "magpurify known-contam \$1 purify/\$1" >> magpurify.sh
    echo "magpurify coverage \$1 purify/\$1 ../align.megahit.bam" >> magpurify.sh
    echo "magpurify clean-bin \$1 purify/\$1 purify/\$1/\$1" >> magpurify.sh
    for file in *.fa; do
        nohup bash magpurify.sh $file > $file.log 2>&1 &
    done
    wait
    conda deactivate
    cd ../; mkdir solidbin_cleaned
    cp solidbin/purify/*/*.fa solidbin_cleaned
fi

if [ ! -d "mycc_cleaned" ]; then
    cd mycc
    source ~/.bashrc
    conda activate python3
    echo "magpurify phylo-markers \$1 purify/\$1" > magpurify.sh
    echo "magpurify clade-markers \$1 purify/\$1" >> magpurify.sh
    echo "magpurify conspecific \$1 purify/\$1 /home/comp/zmzhang/Uhgg/Uhgg.msh" >> magpurify.sh
    echo "magpurify tetra-freq \$1 purify/\$1" >> magpurify.sh
    echo "magpurify gc-content \$1 purify/\$1" >> magpurify.sh
    echo "magpurify known-contam \$1 purify/\$1" >> magpurify.sh
    echo "magpurify coverage \$1 purify/\$1 ../align.megahit.bam" >> magpurify.sh
    echo "magpurify clean-bin \$1 purify/\$1 purify/\$1/\$1" >> magpurify.sh
    for file in *.fa; do
        nohup bash magpurify.sh $file > $file.log 2>&1 &
    done
    wait
    conda deactivate
    cd ../; mkdir mycc_cleaned
    cp mycc/purify/*/*.fa mycc_cleaned
fi

python ~/software/GraphBin_New/support/prepResult.py --binned maxbin --assembler megahit --output maxbin
python ~/software/GraphBin_New/support/prepResult.py --binned maxbin_cleaned --assembler megahit --output maxbin_cleaned
python ~/software/GraphBin_New/support/prepResult.py --binned final.contigs.fa.metabat-bins1500 --assembler megahit --output final.contigs.fa.metabat-bins1500
python ~/software/GraphBin_New/support/prepResult.py --binned final.contigs.fa.metabat-bins1500_cleaned --assembler megahit --output final.contigs.fa.metabat-bins1500_cleaned
python ~/software/GraphBin_New/support/prepResult.py --binned concoct_output --assembler megahit --output concoct_output
python ~/software/GraphBin_New/support/prepResult.py --binned concoct_output_cleaned --assembler megahit --output concoct_output_cleaned
python ~/software/GraphBin_New/support/prepResult.py --binned solidbin --assembler megahit --output solidbin
python ~/software/GraphBin_New/support/prepResult.py --binned solidbin_cleaned --assembler megahit --output solidbin_cleaned
python ~/software/GraphBin_New/support/prepResult.py --binned mycc --assembler megahit --output mycc
python ~/software/GraphBin_New/support/prepResult.py --binned mycc_cleaned --assembler megahit --output mycc_cleaned

# checkm
# if [ ! -f "maxbin.txt" ]; then
#     checkm lineage_wf -f maxbin.txt -t 100 -x fasta maxbin checkm_maxbin > checkm_maxbin.log 2>&1 &
# fi
# if [ ! -f "maxbin_cleaned.txt" ]; then
#     checkm lineage_wf -f maxbin_cleaned.txt -t 100 -x fasta maxbin_cleaned checkm_maxbin_cleaned > checkm_maxbin_cleaned.log 2>&1 &
# fi
# wait

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
if [ ! -d "assembler" ]; then
    mkdir assembler
fi
cd assembler
start=`date +%s`
if [ ! -f "assembler.log" ]; then
    ~/code/assembler/assembler ngs_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv ../align.megahit.bam ../maxbin_cleaned/initial_contig_bins.csv $4 $k_mer $6 $7 $8> assembler.log 2>&1
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

~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv ../maxbin/initial_contig_bins.csv $k_mer > evaluate_megahit.maxbin
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv edge_graph_graphbin/graphbin_output.csv $k_mer > evaluate_megahit.edge_graph
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv final/final.cluster $k_mer > evaluate_megahit.amgl

python ~/code/assembler/parse_result.py --assembler megahit --cluster edge_graph_graphbin/graphbin_output.csv --output edge_graph_graphbin --contigs ../final.contigs.fa
python ~/code/assembler/parse_result.py --assembler megahit --cluster final/final.cluster --output final --contigs ../final.contigs.fa

# if [ ! -f "edge_graph_graphbin.txt" ]; then
#     checkm lineage_wf -f edge_graph_graphbin.txt -t 100 -x fasta edge_graph_graphbin checkm_edge_graph_graphbin > checkm_edge_graph_graphbin.log 2>&1 &
# fi
# if [ ! -f "final.txt" ]; then
#     checkm lineage_wf -f final.txt -t 100 -x fasta final checkm_final > checkm_final.log 2>&1 &
# fi
# wait

# maxbin with purify
if [ ! -d "edge_graph_graphbin_magpurify" ]; then
    ~/software/GraphBin_New/graphbin --binned ../maxbin_cleaned/initial_contig_bins.csv --graph ../final.contigs.gfa --output edge_graph_graphbin_magpurify --assembler megahit --contigs ../final.contigs.fa > edge_graph_graphbin_magpurify.log 2>&1
fi

if [ ! -d "final_magpurify" ]; then
    python ~/code/assembler/AMGL-SEMI/AMGL-SEMI.py --remove True --assembly_graph edge_graph.txt --PE_graph PE_graph.txt --binned ../maxbin_cleaned/initial_contig_bins.csv --max_iter 100 --thresh 0.00000001 > amgl_magpurify.log 2>&1
    mkdir final_magpurify; mv final.cluster final_magpurify
fi

~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv ../maxbin_cleaned/initial_contig_bins.csv $k_mer > evaluate_megahit_magpurify.maxbin
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv edge_graph_graphbin_magpurify/graphbin_output.csv $k_mer > evaluate_megahit_magpurify.edge_graph
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv final_magpurify/final.cluster $k_mer > evaluate_megahit_magpurify.amgl

python ~/code/assembler/parse_result.py --assembler megahit --cluster edge_graph_graphbin_magpurify/graphbin_output.csv --output edge_graph_graphbin_magpurify --contigs ../final.contigs.fa
python ~/code/assembler/parse_result.py --assembler megahit --cluster final_magpurify/final.cluster --output final_magpurify --contigs ../final.contigs.fa

# if [ ! -f "edge_graph_graphbin_magpurify.txt" ]; then
#     checkm lineage_wf -f edge_graph_graphbin_magpurify.txt -t 100 -x fasta edge_graph_graphbin_magpurify checkm_edge_graph_graphbin_magpurify > checkm_edge_graph_graphbin_magpurify.log 2>&1 &
# fi
# if [ ! -f "final_magpurify.txt" ]; then
#     checkm lineage_wf -f final_magpurify.txt -t 100 -x fasta final_magpurify checkm_final_magpurify > checkm_final_magpurify.log 2>&1 &
# fi
# wait

# amgl without reomve
if [ ! -d "final_noremove" ]; then
    python ~/code/assembler/AMGL-SEMI/AMGL-SEMI.py --remove False --assembly_graph edge_graph.txt --PE_graph PE_graph.txt --binned ../maxbin/initial_contig_bins.csv --max_iter 100 --thresh 0.00000001 > amgl_noremove.log 2>&1
    mkdir final_noremove; mv final.cluster final_noremove
fi
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv final_noremove/final.cluster $k_mer > evaluate_megahit_noremove.amgl
if [ ! -d "final_noremove_magpurify" ]; then
    python ~/code/assembler/AMGL-SEMI/AMGL-SEMI.py --remove False --assembly_graph edge_graph.txt --PE_graph PE_graph.txt --binned ../maxbin_cleaned/initial_contig_bins.csv --max_iter 100 --thresh 0.00000001 > amgl_noremove_magpurify.log 2>&1
    mkdir final_noremove_magpurify; mv final.cluster final_noremove_magpurify
fi
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv final_noremove_magpurify/final.cluster $k_mer > evaluate_megahit_magpurify_noremove.amgl

python ~/code/assembler/parse_result.py --assembler megahit --cluster final_noremove/final.cluster --output final_noremove --contigs ../final.contigs.fa
python ~/code/assembler/parse_result.py --assembler megahit --cluster final_noremove_magpurify/final.cluster --output final_noremove_magpurify --contigs ../final.contigs.fa

# if [ ! -f "final_noremove.txt" ]; then
#     checkm lineage_wf -f final_noremove.txt -t 100 -x fasta final_noremove checkm_final_noremove > checkm_final_noremove.log 2>&1 &
# fi
# if [ ! -f "final_noremove_magpurify.txt" ]; then
#     checkm lineage_wf -f final_noremove_magpurify.txt -t 100 -x fasta final_noremove_magpurify checkm_final_noremove_magpurify > checkm_final_noremove_magpurify 2>&1 &
# fi

if [ ! -d "final_noremove_magpurify_cleaned" ]; then
    cd final_noremove_magpurify
    source ~/.bashrc
    conda activate python3
    echo "magpurify phylo-markers \$1 purify/\$1" > magpurify.sh
    echo "magpurify clade-markers \$1 purify/\$1" >> magpurify.sh
    echo "magpurify conspecific \$1 purify/\$1 /home/comp/zmzhang/Uhgg/Uhgg.msh" >> magpurify.sh
    echo "magpurify tetra-freq \$1 purify/\$1" >> magpurify.sh
    echo "magpurify gc-content \$1 purify/\$1" >> magpurify.sh
    echo "magpurify known-contam \$1 purify/\$1" >> magpurify.sh
    echo "magpurify coverage \$1 purify/\$1 ../../align.megahit.bam" >> magpurify.sh
    echo "magpurify clean-bin \$1 purify/\$1 purify/\$1/\$1" >> magpurify.sh
    for file in *.fasta; do
        nohup bash magpurify.sh $file > $file.log 2>&1 &
    done
    wait
    cd ../; mkdir final_noremove_magpurify_cleaned
    conda deactivate
    cp final_noremove_magpurify/purify/*/*.fasta final_noremove_magpurify_cleaned
fi

python ~/software/GraphBin_New/support/prepResult.py --binned final_noremove_magpurify_cleaned --assembler megahit --output final_noremove_magpurify_cleaned
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv final_noremove_magpurify_cleaned/initial_contig_bins.csv $k_mer > evaluate_megahit_magpurify_noremove_cleaned.amgl
# if [ ! -f "final_noremove_magpurify_cleaned.txt" ]; then
#     checkm lineage_wf -f final_noremove_magpurify_cleaned.txt -t 100 -x fasta final_noremove_magpurify_cleaned checkm_final_noremove_magpurify_cleaned > checkm_final_noremove_magpurify_cleaned.log 2>&1
# fi

# for metabat
if [ ! -d "../assembler_metabat" ]; then
    mkdir ../assembler_metabat
fi
cd ../assembler_metabat
cp ../assembler/PE_graph.txt ../assembler/edge_graph.txt ./

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

~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv ../final.contigs.fa.metabat-bins1500/initial_contig_bins.csv $k_mer > evaluate_megahit.maxbin
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv edge_graph_graphbin/graphbin_output.csv $k_mer > evaluate_megahit.edge_graph
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv final/final.cluster $k_mer > evaluate_megahit.amgl

python ~/code/assembler/parse_result.py --assembler megahit --cluster edge_graph_graphbin/graphbin_output.csv --output edge_graph_graphbin --contigs ../final.contigs.fa
python ~/code/assembler/parse_result.py --assembler megahit --cluster final/final.cluster --output final --contigs ../final.contigs.fa

# if [ ! -f "edge_graph_graphbin.txt" ]; then
#     checkm lineage_wf -f edge_graph_graphbin.txt -t 100 -x fasta edge_graph_graphbin checkm_edge_graph_graphbin > checkm_edge_graph_graphbin.log 2>&1 &
# fi
# if [ ! -f "final.txt" ]; then
#     checkm lineage_wf -f final.txt -t 100 -x fasta final checkm_final > checkm_final.log 2>&1 &
# fi
# wait

# maxbin with purify
if [ ! -d "edge_graph_graphbin_magpurify" ]; then
    ~/software/GraphBin_New/graphbin --binned ../final.contigs.fa.metabat-bins1500_cleaned/initial_contig_bins.csv --graph ../final.contigs.gfa --output edge_graph_graphbin_magpurify --assembler megahit --contigs ../final.contigs.fa > edge_graph_graphbin_magpurify.log 2>&1
fi

if [ ! -d "final_magpurify" ]; then
    python ~/code/assembler/AMGL-SEMI/AMGL-SEMI.py --remove True --assembly_graph edge_graph.txt --PE_graph PE_graph.txt --binned ../final.contigs.fa.metabat-bins1500_cleaned/initial_contig_bins.csv --max_iter 100 --thresh 0.00000001 > amgl_magpurify.log 2>&1
    mkdir final_magpurify; mv final.cluster final_magpurify
fi

~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv ../final.contigs.fa.metabat-bins1500_cleaned/initial_contig_bins.csv $k_mer > evaluate_megahit_magpurify.maxbin
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv edge_graph_graphbin_magpurify/graphbin_output.csv $k_mer > evaluate_megahit_magpurify.edge_graph
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv final_magpurify/final.cluster $k_mer > evaluate_megahit_magpurify.amgl

python ~/code/assembler/parse_result.py --assembler megahit --cluster edge_graph_graphbin_magpurify/graphbin_output.csv --output edge_graph_graphbin_magpurify --contigs ../final.contigs.fa
python ~/code/assembler/parse_result.py --assembler megahit --cluster final_magpurify/final.cluster --output final_magpurify --contigs ../final.contigs.fa

# if [ ! -f "edge_graph_graphbin_magpurify.txt" ]; then
#     checkm lineage_wf -f edge_graph_graphbin_magpurify.txt -t 100 -x fasta edge_graph_graphbin_magpurify checkm_edge_graph_graphbin_magpurify > checkm_edge_graph_graphbin_magpurify.log 2>&1 &
# fi
# if [ ! -f "final_magpurify.txt" ]; then
#     checkm lineage_wf -f final_magpurify.txt -t 100 -x fasta final_magpurify checkm_final_magpurify > checkm_final_magpurify.log 2>&1 &
# fi
# wait

# amgl without reomve
if [ ! -d "final_noremove" ]; then
    python ~/code/assembler/AMGL-SEMI/AMGL-SEMI.py --remove False --assembly_graph edge_graph.txt --PE_graph PE_graph.txt --binned ../final.contigs.fa.metabat-bins1500/initial_contig_bins.csv --max_iter 100 --thresh 0.00000001 > amgl_noremove.log 2>&1
    mkdir final_noremove; mv final.cluster final_noremove
fi
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv final_noremove/final.cluster $k_mer > evaluate_megahit_noremove.amgl
if [ ! -d "final_noremove_magpurify" ]; then
    python ~/code/assembler/AMGL-SEMI/AMGL-SEMI.py --remove False --assembly_graph edge_graph.txt --PE_graph PE_graph.txt --binned ../final.contigs.fa.metabat-bins1500_cleaned/initial_contig_bins.csv --max_iter 100 --thresh 0.00000001 > amgl_noremove_magpurify.log 2>&1
    mkdir final_noremove_magpurify; mv final.cluster final_noremove_magpurify
fi
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv final_noremove_magpurify/final.cluster $k_mer > evaluate_megahit_magpurify_noremove.amgl

python ~/code/assembler/parse_result.py --assembler megahit --cluster final_noremove/final.cluster --output final_noremove --contigs ../final.contigs.fa
python ~/code/assembler/parse_result.py --assembler megahit --cluster final_noremove_magpurify/final.cluster --output final_noremove_magpurify --contigs ../final.contigs.fa

# if [ ! -f "final_noremove.txt" ]; then
#     checkm lineage_wf -f final_noremove.txt -t 100 -x fasta final_noremove checkm_final_noremove > checkm_final_noremove.log 2>&1 &
# fi
# if [ ! -f "final_noremove_magpurify.txt" ]; then
#     checkm lineage_wf -f final_noremove_magpurify.txt -t 100 -x fasta final_noremove_magpurify checkm_final_noremove_magpurify > checkm_final_noremove_magpurify 2>&1 &
# fi

if [ ! -d "final_noremove_magpurify_cleaned" ]; then
    cd final_noremove_magpurify
    source ~/.bashrc
    conda activate python3
    echo "magpurify phylo-markers \$1 purify/\$1" > magpurify.sh
    echo "magpurify clade-markers \$1 purify/\$1" >> magpurify.sh
    echo "magpurify conspecific \$1 purify/\$1 /home/comp/zmzhang/Uhgg/Uhgg.msh" >> magpurify.sh
    echo "magpurify tetra-freq \$1 purify/\$1" >> magpurify.sh
    echo "magpurify gc-content \$1 purify/\$1" >> magpurify.sh
    echo "magpurify known-contam \$1 purify/\$1" >> magpurify.sh
    echo "magpurify coverage \$1 purify/\$1 ../../align.megahit.bam" >> magpurify.sh
    echo "magpurify clean-bin \$1 purify/\$1 purify/\$1/\$1" >> magpurify.sh
    for file in *.fasta; do
        nohup bash magpurify.sh $file > $file.log 2>&1 &
    done
    wait
    cd ../; mkdir final_noremove_magpurify_cleaned
    conda deactivate
    cp final_noremove_magpurify/purify/*/*.fasta final_noremove_magpurify_cleaned
fi

python ~/software/GraphBin_New/support/prepResult.py --binned final_noremove_magpurify_cleaned --assembler megahit --output final_noremove_magpurify_cleaned
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv final_noremove_magpurify_cleaned/initial_contig_bins.csv $k_mer > evaluate_megahit_magpurify_noremove_cleaned.amgl
# if [ ! -f "final_noremove_magpurify_cleaned.txt" ]; then
#     checkm lineage_wf -f final_noremove_magpurify_cleaned.txt -t 100 -x fasta final_noremove_magpurify_cleaned checkm_final_noremove_magpurify_cleaned > checkm_final_noremove_magpurify_cleaned.log 2>&1
# fi

# for concoct
if [ ! -d "../assembler_concoct" ]; then
    mkdir ../assembler_concoct
fi
cd ../assembler_concoct
cp ../assembler/PE_graph.txt ../assembler/edge_graph.txt ./

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

~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv ../concoct_output/initial_contig_bins.csv $k_mer > evaluate_megahit.maxbin
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv edge_graph_graphbin/graphbin_output.csv $k_mer > evaluate_megahit.edge_graph
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv final/final.cluster $k_mer > evaluate_megahit.amgl

python ~/code/assembler/parse_result.py --assembler megahit --cluster edge_graph_graphbin/graphbin_output.csv --output edge_graph_graphbin --contigs ../final.contigs.fa
python ~/code/assembler/parse_result.py --assembler megahit --cluster final/final.cluster --output final --contigs ../final.contigs.fa

# if [ ! -f "edge_graph_graphbin.txt" ]; then
#     checkm lineage_wf -f edge_graph_graphbin.txt -t 100 -x fasta edge_graph_graphbin checkm_edge_graph_graphbin > checkm_edge_graph_graphbin.log 2>&1 &
# fi
# if [ ! -f "final.txt" ]; then
#     checkm lineage_wf -f final.txt -t 100 -x fasta final checkm_final > checkm_final.log 2>&1 &
# fi
# wait

# maxbin with purify
if [ ! -d "edge_graph_graphbin_magpurify" ]; then
    ~/software/GraphBin_New/graphbin --binned ../concoct_output_cleaned/initial_contig_bins.csv --graph ../final.contigs.gfa --output edge_graph_graphbin_magpurify --assembler megahit --contigs ../final.contigs.fa > edge_graph_graphbin_magpurify.log 2>&1
fi

if [ ! -d "final_magpurify" ]; then
    python ~/code/assembler/AMGL-SEMI/AMGL-SEMI.py --remove True --assembly_graph edge_graph.txt --PE_graph PE_graph.txt --binned ../concoct_output_cleaned/initial_contig_bins.csv --max_iter 100 --thresh 0.00000001 > amgl_magpurify.log 2>&1
    mkdir final_magpurify; mv final.cluster final_magpurify
fi

~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv ../concoct_output_cleaned/initial_contig_bins.csv $k_mer > evaluate_megahit_magpurify.maxbin
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv edge_graph_graphbin_magpurify/graphbin_output.csv $k_mer > evaluate_megahit_magpurify.edge_graph
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv final_magpurify/final.cluster $k_mer > evaluate_megahit_magpurify.amgl

python ~/code/assembler/parse_result.py --assembler megahit --cluster edge_graph_graphbin_magpurify/graphbin_output.csv --output edge_graph_graphbin_magpurify --contigs ../final.contigs.fa
python ~/code/assembler/parse_result.py --assembler megahit --cluster final_magpurify/final.cluster --output final_magpurify --contigs ../final.contigs.fa

# if [ ! -f "edge_graph_graphbin_magpurify.txt" ]; then
#     checkm lineage_wf -f edge_graph_graphbin_magpurify.txt -t 100 -x fasta edge_graph_graphbin_magpurify checkm_edge_graph_graphbin_magpurify > checkm_edge_graph_graphbin_magpurify.log 2>&1 &
# fi
# if [ ! -f "final_magpurify.txt" ]; then
#     checkm lineage_wf -f final_magpurify.txt -t 100 -x fasta final_magpurify checkm_final_magpurify > checkm_final_magpurify.log 2>&1 &
# fi
# wait

# amgl without reomve
if [ ! -d "final_noremove" ]; then
    python ~/code/assembler/AMGL-SEMI/AMGL-SEMI.py --remove False --assembly_graph edge_graph.txt --PE_graph PE_graph.txt --binned ../concoct_output/initial_contig_bins.csv --max_iter 100 --thresh 0.00000001 > amgl_noremove.log 2>&1
    mkdir final_noremove; mv final.cluster final_noremove
fi
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv final_noremove/final.cluster $k_mer > evaluate_megahit_noremove.amgl
if [ ! -d "final_noremove_magpurify" ]; then
    python ~/code/assembler/AMGL-SEMI/AMGL-SEMI.py --remove False --assembly_graph edge_graph.txt --PE_graph PE_graph.txt --binned ../concoct_output_cleaned/initial_contig_bins.csv --max_iter 100 --thresh 0.00000001 > amgl_noremove_magpurify.log 2>&1
    mkdir final_noremove_magpurify; mv final.cluster final_noremove_magpurify
fi
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv final_noremove_magpurify/final.cluster $k_mer > evaluate_megahit_magpurify_noremove.amgl

python ~/code/assembler/parse_result.py --assembler megahit --cluster final_noremove/final.cluster --output final_noremove --contigs ../final.contigs.fa
python ~/code/assembler/parse_result.py --assembler megahit --cluster final_noremove_magpurify/final.cluster --output final_noremove_magpurify --contigs ../final.contigs.fa

# if [ ! -f "final_noremove.txt" ]; then
#     checkm lineage_wf -f final_noremove.txt -t 100 -x fasta final_noremove checkm_final_noremove > checkm_final_noremove.log 2>&1 &
# fi
# if [ ! -f "final_noremove_magpurify.txt" ]; then
#     checkm lineage_wf -f final_noremove_magpurify.txt -t 100 -x fasta final_noremove_magpurify checkm_final_noremove_magpurify > checkm_final_noremove_magpurify 2>&1 &
# fi

if [ ! -d "final_noremove_magpurify_cleaned" ]; then
    cd final_noremove_magpurify
    source ~/.bashrc
    conda activate python3
    echo "magpurify phylo-markers \$1 purify/\$1" > magpurify.sh
    echo "magpurify clade-markers \$1 purify/\$1" >> magpurify.sh
    echo "magpurify conspecific \$1 purify/\$1 /home/comp/zmzhang/Uhgg/Uhgg.msh" >> magpurify.sh
    echo "magpurify tetra-freq \$1 purify/\$1" >> magpurify.sh
    echo "magpurify gc-content \$1 purify/\$1" >> magpurify.sh
    echo "magpurify known-contam \$1 purify/\$1" >> magpurify.sh
    echo "magpurify coverage \$1 purify/\$1 ../../align.megahit.bam" >> magpurify.sh
    echo "magpurify clean-bin \$1 purify/\$1 purify/\$1/\$1" >> magpurify.sh
    for file in *.fasta; do
        nohup bash magpurify.sh $file > $file.log 2>&1 &
    done
    wait
    cd ../; mkdir final_noremove_magpurify_cleaned
    conda deactivate
    cp final_noremove_magpurify/purify/*/*.fasta final_noremove_magpurify_cleaned
fi

python ~/software/GraphBin_New/support/prepResult.py --binned final_noremove_magpurify_cleaned --assembler megahit --output final_noremove_magpurify_cleaned
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv final_noremove_magpurify_cleaned/initial_contig_bins.csv $k_mer > evaluate_megahit_magpurify_noremove_cleaned.amgl
# if [ ! -f "final_noremove_magpurify_cleaned.txt" ]; then
#     checkm lineage_wf -f final_noremove_magpurify_cleaned.txt -t 100 -x fasta final_noremove_magpurify_cleaned checkm_final_noremove_magpurify_cleaned > checkm_final_noremove_magpurify_cleaned.log 2>&1
# fi

# for solidbin
if [ ! -d "../assembler_solidbin" ]; then
    mkdir ../assembler_solidbin
fi
cd ../assembler_solidbin
cp ../assembler/PE_graph.txt ../assembler/edge_graph.txt ./

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

~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv ../solidbin/initial_contig_bins.csv $k_mer > evaluate_megahit.maxbin
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv edge_graph_graphbin/graphbin_output.csv $k_mer > evaluate_megahit.edge_graph
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv final/final.cluster $k_mer > evaluate_megahit.amgl

python ~/code/assembler/parse_result.py --assembler megahit --cluster edge_graph_graphbin/graphbin_output.csv --output edge_graph_graphbin --contigs ../final.contigs.fa
python ~/code/assembler/parse_result.py --assembler megahit --cluster final/final.cluster --output final --contigs ../final.contigs.fa

# if [ ! -f "edge_graph_graphbin.txt" ]; then
#     checkm lineage_wf -f edge_graph_graphbin.txt -t 100 -x fasta edge_graph_graphbin checkm_edge_graph_graphbin > checkm_edge_graph_graphbin.log 2>&1 &
# fi
# if [ ! -f "final.txt" ]; then
#     checkm lineage_wf -f final.txt -t 100 -x fasta final checkm_final > checkm_final.log 2>&1 &
# fi
# wait

# maxbin with purify
if [ ! -d "edge_graph_graphbin_magpurify" ]; then
    ~/software/GraphBin_New/graphbin --binned ../solidbin_cleaned/initial_contig_bins.csv --graph ../final.contigs.gfa --output edge_graph_graphbin_magpurify --assembler megahit --contigs ../final.contigs.fa > edge_graph_graphbin_magpurify.log 2>&1
fi

if [ ! -d "final_magpurify" ]; then
    python ~/code/assembler/AMGL-SEMI/AMGL-SEMI.py --remove True --assembly_graph edge_graph.txt --PE_graph PE_graph.txt --binned ../solidbin_cleaned/initial_contig_bins.csv --max_iter 100 --thresh 0.00000001 > amgl_magpurify.log 2>&1
    mkdir final_magpurify; mv final.cluster final_magpurify
fi

~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv ../solidbin_cleaned/initial_contig_bins.csv $k_mer > evaluate_megahit_magpurify.maxbin
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv edge_graph_graphbin_magpurify/graphbin_output.csv $k_mer > evaluate_megahit_magpurify.edge_graph
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv final_magpurify/final.cluster $k_mer > evaluate_megahit_magpurify.amgl

python ~/code/assembler/parse_result.py --assembler megahit --cluster edge_graph_graphbin_magpurify/graphbin_output.csv --output edge_graph_graphbin_magpurify --contigs ../final.contigs.fa
python ~/code/assembler/parse_result.py --assembler megahit --cluster final_magpurify/final.cluster --output final_magpurify --contigs ../final.contigs.fa

# if [ ! -f "edge_graph_graphbin_magpurify.txt" ]; then
#     checkm lineage_wf -f edge_graph_graphbin_magpurify.txt -t 100 -x fasta edge_graph_graphbin_magpurify checkm_edge_graph_graphbin_magpurify > checkm_edge_graph_graphbin_magpurify.log 2>&1 &
# fi
# if [ ! -f "final_magpurify.txt" ]; then
#     checkm lineage_wf -f final_magpurify.txt -t 100 -x fasta final_magpurify checkm_final_magpurify > checkm_final_magpurify.log 2>&1 &
# fi
# wait

# amgl without reomve
if [ ! -d "final_noremove" ]; then
    python ~/code/assembler/AMGL-SEMI/AMGL-SEMI.py --remove False --assembly_graph edge_graph.txt --PE_graph PE_graph.txt --binned ../solidbin/initial_contig_bins.csv --max_iter 100 --thresh 0.00000001 > amgl_noremove.log 2>&1
    mkdir final_noremove; mv final.cluster final_noremove
fi
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv final_noremove/final.cluster $k_mer > evaluate_megahit_noremove.amgl
if [ ! -d "final_noremove_magpurify" ]; then
    python ~/code/assembler/AMGL-SEMI/AMGL-SEMI.py --remove False --assembly_graph edge_graph.txt --PE_graph PE_graph.txt --binned ../solidbin_cleaned/initial_contig_bins.csv --max_iter 100 --thresh 0.00000001 > amgl_noremove_magpurify.log 2>&1
    mkdir final_noremove_magpurify; mv final.cluster final_noremove_magpurify
fi
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv final_noremove_magpurify/final.cluster $k_mer > evaluate_megahit_magpurify_noremove.amgl

python ~/code/assembler/parse_result.py --assembler megahit --cluster final_noremove/final.cluster --output final_noremove --contigs ../final.contigs.fa
python ~/code/assembler/parse_result.py --assembler megahit --cluster final_noremove_magpurify/final.cluster --output final_noremove_magpurify --contigs ../final.contigs.fa

# if [ ! -f "final_noremove.txt" ]; then
#     checkm lineage_wf -f final_noremove.txt -t 100 -x fasta final_noremove checkm_final_noremove > checkm_final_noremove.log 2>&1 &
# fi
# if [ ! -f "final_noremove_magpurify.txt" ]; then
#     checkm lineage_wf -f final_noremove_magpurify.txt -t 100 -x fasta final_noremove_magpurify checkm_final_noremove_magpurify > checkm_final_noremove_magpurify 2>&1 &
# fi

if [ ! -d "final_noremove_magpurify_cleaned" ]; then
    cd final_noremove_magpurify
    source ~/.bashrc
    conda activate python3
    echo "magpurify phylo-markers \$1 purify/\$1" > magpurify.sh
    echo "magpurify clade-markers \$1 purify/\$1" >> magpurify.sh
    echo "magpurify conspecific \$1 purify/\$1 /home/comp/zmzhang/Uhgg/Uhgg.msh" >> magpurify.sh
    echo "magpurify tetra-freq \$1 purify/\$1" >> magpurify.sh
    echo "magpurify gc-content \$1 purify/\$1" >> magpurify.sh
    echo "magpurify known-contam \$1 purify/\$1" >> magpurify.sh
    echo "magpurify coverage \$1 purify/\$1 ../../align.megahit.bam" >> magpurify.sh
    echo "magpurify clean-bin \$1 purify/\$1 purify/\$1/\$1" >> magpurify.sh
    for file in *.fasta; do
        nohup bash magpurify.sh $file > $file.log 2>&1 &
    done
    wait
    cd ../; mkdir final_noremove_magpurify_cleaned
    conda deactivate
    cp final_noremove_magpurify/purify/*/*.fasta final_noremove_magpurify_cleaned
fi

python ~/software/GraphBin_New/support/prepResult.py --binned final_noremove_magpurify_cleaned --assembler megahit --output final_noremove_magpurify_cleaned
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv final_noremove_magpurify_cleaned/initial_contig_bins.csv $k_mer > evaluate_megahit_magpurify_noremove_cleaned.amgl
# if [ ! -f "final_noremove_magpurify_cleaned.txt" ]; then
#     checkm lineage_wf -f final_noremove_magpurify_cleaned.txt -t 100 -x fasta final_noremove_magpurify_cleaned checkm_final_noremove_magpurify_cleaned > checkm_final_noremove_magpurify_cleaned.log 2>&1
# fi

# for mycc
if [ ! -d "../assembler_mycc" ]; then
    mkdir ../assembler_mycc
fi
cd ../assembler_mycc
cp ../assembler/PE_graph.txt ../assembler/edge_graph.txt ./

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

~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv ../mycc/initial_contig_bins.csv $k_mer > evaluate_megahit.maxbin
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv edge_graph_graphbin/graphbin_output.csv $k_mer > evaluate_megahit.edge_graph
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv final/final.cluster $k_mer > evaluate_megahit.amgl

python ~/code/assembler/parse_result.py --assembler megahit --cluster edge_graph_graphbin/graphbin_output.csv --output edge_graph_graphbin --contigs ../final.contigs.fa
python ~/code/assembler/parse_result.py --assembler megahit --cluster final/final.cluster --output final --contigs ../final.contigs.fa

# if [ ! -f "edge_graph_graphbin.txt" ]; then
#     checkm lineage_wf -f edge_graph_graphbin.txt -t 100 -x fasta edge_graph_graphbin checkm_edge_graph_graphbin > checkm_edge_graph_graphbin.log 2>&1 &
# fi
# if [ ! -f "final.txt" ]; then
#     checkm lineage_wf -f final.txt -t 100 -x fasta final checkm_final > checkm_final.log 2>&1 &
# fi
# wait

# maxbin with purify
if [ ! -d "edge_graph_graphbin_magpurify" ]; then
    ~/software/GraphBin_New/graphbin --binned ../mycc_cleaned/initial_contig_bins.csv --graph ../final.contigs.gfa --output edge_graph_graphbin_magpurify --assembler megahit --contigs ../final.contigs.fa > edge_graph_graphbin_magpurify.log 2>&1
fi

if [ ! -d "final_magpurify" ]; then
    python ~/code/assembler/AMGL-SEMI/AMGL-SEMI.py --remove True --assembly_graph edge_graph.txt --PE_graph PE_graph.txt --binned ../mycc_cleaned/initial_contig_bins.csv --max_iter 100 --thresh 0.00000001 > amgl_magpurify.log 2>&1
    mkdir final_magpurify; mv final.cluster final_magpurify
fi

~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv ../mycc_cleaned/initial_contig_bins.csv $k_mer > evaluate_megahit_magpurify.maxbin
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv edge_graph_graphbin_magpurify/graphbin_output.csv $k_mer > evaluate_megahit_magpurify.edge_graph
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv final_magpurify/final.cluster $k_mer > evaluate_megahit_magpurify.amgl

python ~/code/assembler/parse_result.py --assembler megahit --cluster edge_graph_graphbin_magpurify/graphbin_output.csv --output edge_graph_graphbin_magpurify --contigs ../final.contigs.fa
python ~/code/assembler/parse_result.py --assembler megahit --cluster final_magpurify/final.cluster --output final_magpurify --contigs ../final.contigs.fa

# if [ ! -f "edge_graph_graphbin_magpurify.txt" ]; then
#     checkm lineage_wf -f edge_graph_graphbin_magpurify.txt -t 100 -x fasta edge_graph_graphbin_magpurify checkm_edge_graph_graphbin_magpurify > checkm_edge_graph_graphbin_magpurify.log 2>&1 &
# fi
# if [ ! -f "final_magpurify.txt" ]; then
#     checkm lineage_wf -f final_magpurify.txt -t 100 -x fasta final_magpurify checkm_final_magpurify > checkm_final_magpurify.log 2>&1 &
# fi
# wait

# amgl without reomve
if [ ! -d "final_noremove" ]; then
    python ~/code/assembler/AMGL-SEMI/AMGL-SEMI.py --remove False --assembly_graph edge_graph.txt --PE_graph PE_graph.txt --binned ../mycc/initial_contig_bins.csv --max_iter 100 --thresh 0.00000001 > amgl_noremove.log 2>&1
    mkdir final_noremove; mv final.cluster final_noremove
fi
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv final_noremove/final.cluster $k_mer > evaluate_megahit_noremove.amgl
if [ ! -d "final_noremove_magpurify" ]; then
    python ~/code/assembler/AMGL-SEMI/AMGL-SEMI.py --remove False --assembly_graph edge_graph.txt --PE_graph PE_graph.txt --binned ../mycc_cleaned/initial_contig_bins.csv --max_iter 100 --thresh 0.00000001 > amgl_noremove_magpurify.log 2>&1
    mkdir final_noremove_magpurify; mv final.cluster final_noremove_magpurify
fi
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv final_noremove_magpurify/final.cluster $k_mer > evaluate_megahit_magpurify_noremove.amgl

python ~/code/assembler/parse_result.py --assembler megahit --cluster final_noremove/final.cluster --output final_noremove --contigs ../final.contigs.fa
python ~/code/assembler/parse_result.py --assembler megahit --cluster final_noremove_magpurify/final.cluster --output final_noremove_magpurify --contigs ../final.contigs.fa

# if [ ! -f "final_noremove.txt" ]; then
#     checkm lineage_wf -f final_noremove.txt -t 100 -x fasta final_noremove checkm_final_noremove > checkm_final_noremove.log 2>&1 &
# fi
# if [ ! -f "final_noremove_magpurify.txt" ]; then
#     checkm lineage_wf -f final_noremove_magpurify.txt -t 100 -x fasta final_noremove_magpurify checkm_final_noremove_magpurify > checkm_final_noremove_magpurify 2>&1 &
# fi

if [ ! -d "final_noremove_magpurify_cleaned" ]; then
    cd final_noremove_magpurify
    source ~/.bashrc
    conda activate python3
    echo "magpurify phylo-markers \$1 purify/\$1" > magpurify.sh
    echo "magpurify clade-markers \$1 purify/\$1" >> magpurify.sh
    echo "magpurify conspecific \$1 purify/\$1 /home/comp/zmzhang/Uhgg/Uhgg.msh" >> magpurify.sh
    echo "magpurify tetra-freq \$1 purify/\$1" >> magpurify.sh
    echo "magpurify gc-content \$1 purify/\$1" >> magpurify.sh
    echo "magpurify known-contam \$1 purify/\$1" >> magpurify.sh
    echo "magpurify coverage \$1 purify/\$1 ../../align.megahit.bam" >> magpurify.sh
    echo "magpurify clean-bin \$1 purify/\$1 purify/\$1/\$1" >> magpurify.sh
    for file in *.fasta; do
        nohup bash magpurify.sh $file > $file.log 2>&1 &
    done
    wait
    cd ../; mkdir final_noremove_magpurify_cleaned
    conda deactivate
    cp final_noremove_magpurify/purify/*/*.fasta final_noremove_magpurify_cleaned
fi

python ~/software/GraphBin_New/support/prepResult.py --binned final_noremove_magpurify_cleaned --assembler megahit --output final_noremove_magpurify_cleaned
~/code/assembler/assembler evaluate_megahit ../final.contigs.fa ../final.contigs.fastg ../kraken2_0.9.tsv final_noremove_magpurify_cleaned/initial_contig_bins.csv $k_mer > evaluate_megahit_magpurify_noremove_cleaned.amgl
# if [ ! -f "final_noremove_magpurify_cleaned.txt" ]; then
#     checkm lineage_wf -f final_noremove_magpurify_cleaned.txt -t 100 -x fasta final_noremove_magpurify_cleaned checkm_final_noremove_magpurify_cleaned > checkm_final_noremove_magpurify_cleaned.log 2>&1
# fi