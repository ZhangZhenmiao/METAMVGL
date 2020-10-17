# METAMVGL

## Dependencies
- gcc with C++11 support, [HTSlib](https://github.com/samtools/htslib) required
- python3 with [numpy](https://numpy.org/install/), [scipy](https://www.scipy.org/install.html) and [networkx](http://networkx.github.io/)

## Install
Download METAMVGL, ```cd METAMVGL``` and:
```
make
```
The commnad will generate ```prep_graph```. We need add ```prep_graph``` to environmental variables:
```
export PATH=`pwd`:$PATH
```

## Usage

### Assembly

Currently we support [metaSPAdes](https://github.com/ablab/spades) and [MEGAHIT](https://github.com/voutcn/megahit).

### Initial binning

We accept any initial binning tools. To convert the binning result to the input format of METAMVGL, we suggest to use [prepResult.py](https://github.com/Vini2/GraphBin/tree/master/support):
```
python prepResult.py --binned /path/to/folder_with_binning_result --assembler assembler_type_(SPAdes/MEGAIHT) --output /path/to/output_folder
```
This command will create a file named ```initial_contig_bins.csv``` in ```/path/to/output_folder```.

### Prepare graphs
We generate the assembly graph (.ag) and PE graph (.pe) by ```prep_graph```:
```
usage: ./prep_graph --assembler=string --assembly-graph=string --bam=string --output=string [options] ...
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

### Multi-view graph-based binning
We create the binning result of METAMVGL by ```METAMVGL.py```:
```
usage: METAMVGL.py [-h] --contigs CONTIGS --assembler ASSEMBLER
                   --assembly_graph ASSEMBLY_GRAPH --PE_graph PE_GRAPH
                   --binned BINNED --max_iter MAX_ITER --thresh THRESH
                   --output OUTPUT

optional arguments:
  -h, --help            show this help message and exit
  --contigs CONTIGS     path to contigs file
  --assembler ASSEMBLER
                        assembler used
  --assembly_graph ASSEMBLY_GRAPH
                        path to the .ag file
  --PE_graph PE_GRAPH   path to the .pe file
  --binned BINNED       path to the .csv file as initial binning
  --max_iter MAX_ITER   max iteration
  --thresh THRESH       stop threshold
  --output OUTPUT       output folder
```
In the OUTPUT folder, we provide two types of binning result:
- ```binning_result.csv```, each line is contig_name, cluster_id
- ```cluster.*.fasta```, the contigs in fasta format of each cluster
