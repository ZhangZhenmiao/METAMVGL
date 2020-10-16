#include "graph.h"
int main(int argc, char* argv[]) {
    Graph graph;
    std::string graphType = argv[1];
    if (!graphType.compare("cloud_spades")) {
        graph.load_from_spades_fastg("/mnt/osf1/user/zhanglu/zhanglu/hybrid/hlj/hlj_stlfr/hlj_stlfr_assembly/assembly_graph.fastg");
        graph.load_spades_paths("/mnt/osf1/user/zhanglu/zhanglu/hybrid/hlj/hlj_stlfr/hlj_stlfr_assembly/contigs.paths", argv[2]);
        graph.load_spades_alignment_cloud("/mnt/osf2/user/zhangzhm/stlfr/Reads/stlfr_hlj/contigs.#_#_#.bam");
        graph.load_graphbin_result("/mnt/osf2/user/zhangzhm/stlfr/Reads/stlfr_hlj/contigs.graphbin/graphbin_output.csv");
        graph.construct_spades_links_by_edge();
        graph.construct_PE_overlap();
        graph.construct_barcode_overlap();
        graph.breadth_first_search_cloud();
    }
    else if (!graphType.compare("cloud_megahit")) {
        graph.load_from_megahit_fastg(
            "/mnt/osf2/user/zhangzhm/stlfr/megahit_10x_atcc/megahit/k141.contigs.fa",
            "/mnt/osf2/user/zhangzhm/stlfr/megahit_10x_atcc/megahit/k141.fastg",
            "/mnt/osf2/USERS/zhangzhm/stlfr/assembler/assembler_10x_mock/align_ref/answer.txt",
            "141"
        );
        graph.load_megahit_alignment_cloud("/mnt/osf2/USERS/zhangzhm/stlfr/assembler/assembler_10x_mock/10x_mock_pe/align.megahit.bam");
        graph.load_graphbin_result("/mnt/osf2/USERS/zhangzhm/stlfr/maxbin_10x_mock/initial_contig_bins.csv");
        graph.construct_PE_overlap();
        graph.construct_barcode_overlap();
        graph.breadth_first_search_cloud();
    }
    else if (!graphType.compare("ngs_megahit")) {
        // .fa, .fastg, answer, k_mer
        graph.load_from_megahit_fastg(argv[2], argv[3], argv[4], argv[8]);
        // .bam, mq, idt
        graph.load_ngs_alignment_megahit(argv[5], argv[9], argv[10]);
        graph.load_graphbin_result(argv[6]);
        // insert size, num of pes
        graph.construct_PE_overlap_athena(argv[7], argv[11]);
        graph.breadth_first_search_ngs(argv[11]);
    }
    else if (!graphType.compare("ngs_spades")) {
        graph.load_from_spades_fastg(argv[2]);
        graph.load_spades_paths(argv[3], argv[4]);
        // .bam, mq, idt
        graph.load_ngs_alignment_spades(argv[5], argv[8], argv[9]);
        graph.load_graphbin_result(argv[6]);
        graph.construct_spades_links_by_edge();
        // insert size, num of pes
        graph.construct_PE_overlap_athena(argv[7], argv[10]);
        graph.breadth_first_search_ngs(argv[10]);
    }
    else if (!graphType.compare("metacarvel")) {
        graph.load_from_megahit_fastg(
            "/mnt/osf2/user/zhangzhm/stlfr/megahit_10x_atcc/megahit/k141.contigs.fa",
            "/mnt/osf2/user/zhangzhm/stlfr/megahit_10x_atcc/megahit/k141.fastg",
            "/home/comp/zmzhang/12Mock/megahit/align_ref/answer.txt",
            "141"
        );
        graph.load_from_megahit_metacarvel("/mnt/osf2/USERS/zhangzhm/stlfr/metacarvel_10x_mock/contig_links");
        graph.load_graphbin_result(argv[2]);
        graph.breadth_first_search_metacarvel();
    }
    else if (!graphType.compare("evaluate_megahit")) {
        // .fa, .fastg, answer, k_mer
        graph.load_from_megahit_fastg(argv[2], argv[3], argv[4], argv[6]);
        graph.load_graphbin_result(argv[5]);
        graph.bin_evaluate();
    }
    else if (!graphType.compare("evaluate_spades")) {
        graph.load_from_spades_fastg(argv[2]);
        graph.load_spades_paths(argv[3], argv[4]);
        graph.load_graphbin_result(argv[5]);
        graph.construct_spades_links_by_edge();
        graph.bin_evaluate();
    }
    else if (!graphType.compare("test_atcc")) {
        // .fa, .fastg, answer, k_mer
        graph.load_from_megahit_fastg(argv[2], argv[3], argv[4], argv[8]);
        graph.breadth_first_search_ngs("0");
    }
    return 0;
}
