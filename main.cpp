#include <string>
#include "graph.h"
#include "cmdline.h"
int main(int argc, char* argv[]) {
    cmdline::parser argParser;
    argParser.add<std::string>("assembler", 'a',
                               "the assembler used to produce contigs, currently support metaSPAdes and MEGAHIT",  true, "");
    argParser.add<std::string>("contigs", 'c', "the path to the contigs, only needed for MEGAHIT", false, "final.contigs.fa");
    argParser.add<std::string>("paths", 'p', "the path to the .paths file, only needed for metaSPAdes", false, "contigs.paths");
    argParser.add<std::string>("assembly-graph", 'g', "the path to the assembly graph in fastg", true, "");
    argParser.add<std::string>("bam", 'b', "the path to the alignment bam file", true, "");
    argParser.add<double>("mapping-quality", 'm', "the threshold of mapping quality", false, 10);
    argParser.add<double>("identity", 'i', "the threshold of identity", false, 0.95);
    argParser.add<int>("insert-size", 's', "the insert size of paired-end reads", false, 270);
    argParser.add<int>("pe", 'n', "the minimum number of paired-end reads to support a link", false, 3);
    argParser.add<std::string>("output", 'o', "the prefix to output", true, "graph_output");

    argParser.parse_check(argc, argv);
    std::string assembler = argParser.get<std::string>("assembler");
    std::string contigs = argParser.get<std::string>("contigs");
    std::string paths = argParser.get<std::string>("paths");
    std::string assembly = argParser.get<std::string>("assembly-graph");
    std::string bam = argParser.get<std::string>("bam");
    double mq = argParser.get<double>("mapping-quality");
    double idt = argParser.get<double>("identity");
    int insert = argParser.get<int>("insert-size");
    int pe = argParser.get<int>("pe");
    std::string output = argParser.get<std::string>("output");

    Graph graph;
    if (!assembler.compare("MEGAHIT")) {
        graph.load_from_megahit_fastg(contigs, assembly);
        graph.load_ngs_alignment_megahit(bam, mq, idt);
        graph.construct_PE_overlap_athena(insert, pe);
        graph.breadth_first_search_ngs(pe, output);
    }
    else if (!assembler.compare("metaSPAdes")) {
        graph.load_from_spades_fastg(assembly);
        graph.load_spades_paths(paths);
        graph.load_ngs_alignment_spades(bam, mq, idt);
        graph.construct_spades_links_by_edge();
        graph.construct_PE_overlap_athena(insert, pe);
        graph.breadth_first_search_ngs(pe, output);
    }
    else std::cout << "unsupported assembler " + assembler << std::endl;
    return 0;
}
