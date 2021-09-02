# METAMVGL

## Dependencies
- GCC with C++11, [HTSlib](https://github.com/samtools/htslib).
- Python3 with [numpy](https://numpy.org/install), [scipy](https://www.scipy.org/install.html), [networkx](http://networkx.github.io) and [Bio](https://biopython.org/wiki/Getting_Started).
- Initial binning tools [MaxBin2](https://sourceforge.net/projects/maxbin2) and [MetaBat2](https://bitbucket.org/berkeleylab/metabat/src/master).
- Alignment tools [bwa](https://github.com/lh3/bwa) and [samtools](http://www.htslib.org).

## Environment variables
- HTSlib: After compiling, there would be ```include``` and ```lib``` under your specified folder. You need append ```/path/to/include``` to ```$CPLUS_INCLUDE_PATH```, and ```/path/to/lib``` to both ```$LIBRARY_PATH``` and ```$LD_LIBRARY_PATH```.
- MaxBin2 : ```run_MaxBin.pl``` should be avalible.
- MetaBat2: ```jgi_summarize_bam_contig_depths``` and ```metabat2``` should be avalible.

## Install METAMVGL
Download and compile:
```
git clone https://github.com/ZhangZhenmiao/METAMVGL.git
cd METAMVGL && make && chmod +x *.py
```
Add components of METAMVGL to $PATH:
```
export PATH=/path/to/METAMVGL:$PATH
```

## Assembly Graph
For metaSPAdes, the assembly graph (assembly_graph.fastg) is already in the output folder.

For MEGAHIT, the assembly graph is derived from final.contigs.fa:
```
megahit_toolkit contig2fastg k_mer final.contigs.fa > final.contigs.fastg
```

## Run METAMVGL

### The Wrapper

The simplest way to run METAMVGL is the wrapper, that runs initial binning (MetaBat2 or MaxBin2), paired-end graph extraction, and METAMVGL in one command:
```
usage: run_METAMVGL.py [-h] -a ASSEMBLER [-b INITIAL_BINNING] -r1 READS1 -r2
                       READS2 -c CONTIGS [-p PATHS] -g ASSEMBLY_GRAPH
                       [-m MAPPING_QUALITY] [-s IDENTITY] [-i INSERT_SIZE]
                       [-n PE] [-u MAX_ITER] [-d THRESHOLD] [-t THREADS] -o
                       OUTPUT

optional arguments:
  -h, --help            show this help message and exit
  -a ASSEMBLER, --assembler ASSEMBLER
                        the assembler that was used to generate contigs
                        (choose from metaSPAdes and MEGAHIT)
  -b INITIAL_BINNING, --initial_binning INITIAL_BINNING
                        the initial binning tool to create initial binning
                        results (choose from MetaBat2 and MaxBin2, default
                        MetaBat2)
  -r1 READS1, --reads1 READS1
                        the path to reads1
  -r2 READS2, --reads2 READS2
                        the path to reads2
  -c CONTIGS, --contigs CONTIGS
                        the path to the contigs
  -p PATHS, --paths PATHS
                        the path to the contigs.paths file under metaSPAdes
                        output folder, needed for metaSPAdes
  -g ASSEMBLY_GRAPH, --assembly_graph ASSEMBLY_GRAPH
                        the path to the assembly_graph.fastg file under
                        metaSPAdes output folder (metaSPAdes), or fastg format
                        file ceretad by megahit_toolkit (MEGAHIT)
  -m MAPPING_QUALITY, --mapping_quality MAPPING_QUALITY
                        the threshold of mapping quality for reads alignment
                        (default 10)
  -s IDENTITY, --identity IDENTITY
                        the threshold of alignment identity for reads
                        alignment (default 0.95)
  -i INSERT_SIZE, --insert_size INSERT_SIZE
                        the insert size of paired-end reads (default 270)
  -n PE, --pe PE        the minimum number of paired-end reads to support a PE
                        link (default 3)
  -u MAX_ITER, --max_iter MAX_ITER
                        the maximun number of iteration of label propagation
                        (default 100)
  -d THRESHOLD, --threshold THRESHOLD
                        the threshold to stop iteration of label propagation
                        (default 0.00000001)
  -t THREADS, --threads THREADS
                        the number of threads for initial binning tools
                        (default 16)
  -o OUTPUT, --output OUTPUT
                        output dir
```
Example 1:
To run METAMVGL for metaSPAdes assembly, and use MetaBat2 to initial bin:
```
run_METAMVGL.py -a metaspades -r1 /path/to/reads1.fq.gz -r2 /path/to/reads2.fq.gz -c metaspades/contigs.fasta -p metaspades/contigs.paths -g metaspades/assembly_graph.fastg -o METAMVGL_output -t 100
```
Example2:
To run METAMVGL for MEGAHIT assembly, and use MaxBin2 to initial bin:
```
run_METAMVGL.py -a megahit -b maxbin2 -r1 /path/to/reads1.fq.gz -r2 /path/to/reads2.fq.gz -c megahit/final.contigs.fa -g megahit/final.contigs.fastg -o METAMVGL_output -t 100
```
The results would be in ```METAMVGL_output/METAMVGL```.

### Run METAMVGL Step by Step

#### Assembly

We support [metaSPAdes](https://github.com/ablab/spades) and [MEGAHIT](https://github.com/voutcn/megahit).

#### Initial binning

Initial binning can use any binning tools. To convert the initial binning results to the input format of METAMVGL, we modified [prepResult.py](https://github.com/Vini2/GraphBin/tree/master/support), and the usage is:
```
prepResult.py --binned /path/to/initial_binning_result --assembler assembler_type_(SPAdes/MEGAIHT) --output /path/to/output_folder
```
It will create ```initial_contig_bins.csv``` in ```/path/to/output_folder```. It is the input to METAMVGL.py (--binned).

#### Prepare graphs

We generate the assembly graph (.ag) and PE graph (.pe) by ```prep_graph```:
```
usage: prep_graph --assembler=string --assembly-graph=string --bam=string --output=string [options] ...
options:
  -a, --assembler          the assembler used to produce contigs, currently support metaSPAdes and MEGAHIT (string)
  -c, --contigs            the path to the contigs, only needed for MEGAHIT (string [=final.contigs.fa])
  -p, --paths              the path to the .paths file, only needed for metaSPAdes (string [=contigs.paths])
  -g, --assembly-graph     the path to the assembly graph in fastg (string)
  -b, --bam                the path to the alignment bam file (string)
  -m, --mapping-quality    the threshold of mapping quality (double [=10])
  -i, --identity           the threshold of identity (double [=0.95])
  -s, --insert-size        the insert size of paired-end reads (int [=270])
  -n, --pe                 the minimum number of paired-end reads to support a link (int [=3])
  -o, --output             the prefix to output (string)
  -?, --help               print this message
```

#### Multi-view graph-based binning
We create the binning results by ```METAMVGL.py```:
```
usage: METAMVGL.py [-h] --contigs CONTIGS --assembler ASSEMBLER
                   --assembly_graph ASSEMBLY_GRAPH --PE_graph PE_GRAPH
                   --binned BINNED [--max_iter MAX_ITER] [--thresh THRESH]
                   --output OUTPUT

optional arguments:
  -h, --help            show this help message and exit
  --contigs CONTIGS     path to contigs file
  --assembler ASSEMBLER
                        assembler used (metaSPAdes or MEGAHIT)
  --assembly_graph ASSEMBLY_GRAPH
                        path to the .ag file
  --PE_graph PE_GRAPH   path to the .pe file
  --binned BINNED       path to the .csv file as initial binning
  --max_iter MAX_ITER   max iteration (default 100)
  --thresh THRESH       stop threshold (default 0.00000001)
  --output OUTPUT       output folder
```

## Output
In the output folder, we provide two types of binning results:
- ```binning_result.csv```, each line is contig_name, cluster_id
- ```cluster.*.fasta```, the contigs in fasta format of each cluster

## Time and Memory comparing with GraphBin (v1.3)
- The comparison results can be accessed [here](https://drive.google.com/drive/folders/11U4YwiLLrcTCwpWy7Vax9n5Pk99E_8WL?usp=sharing).
- The machine used for comparison is CentOS 8.2 (64-bit), with Dual 26-core Intel Xeon Gold 6230R 2.10GHz CPU and 768GB RAM.
- The measured time and memory include GraphBin/METAMVGL binning on MaxBin2/MetaBAT2 initial binning results from metaSPAdes/MEGAHIT assembly on BMock12, SYNTH64 and Sharon datasets.
- The `time_memory/README.md` has the commands for binning, the evaluation results are in `time_memory/*/*/*.time` and generated by `time_memory/run_compare.sh`.
