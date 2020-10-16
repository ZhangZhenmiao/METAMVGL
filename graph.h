#ifndef GUARD_graph_h
#define GUARD_graph_h
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>

class Edge {
public:
    int length;
    double depth;
    std::string forwardSequence, reverseSequence;
    Edge();
    Edge(std::string, bool, double, int);
    Edge(bool, std::string, double, int);
    std::string reverseComplementary(std::string);
    Edge& operator=(const Edge&);
};

class AlignedRead {
public:
    int startPos, endPos;
    int alignedBases, mappingQuality;
    double blastIdentity;
    std::string sequence, barcode;
    AlignedRead();
    AlignedRead& operator=(const AlignedRead&);
};

class Path {
public:
    int length, bin;
    double depth;
    std::vector<std::string> annotation;
    std::vector<std::vector<std::string>> pathEdge;
    std::unordered_map<std::string, std::vector<AlignedRead>> reads1vec;
    std::unordered_map<std::string, std::vector<AlignedRead>> reads2vec;
    std::unordered_map<std::string, std::unordered_set<std::string>> barcode2reads;
    std::unordered_set<std::string> linksByEdge;
    std::unordered_map<std::string, int> barcodeOverlap;
    std::unordered_map<std::string, int> PEOverlap;
    Path();
    Path(std::vector<std::vector<std::string>>, double, int);
    Path& operator=(const Path&);
};

class Graph {
public:
    std::unordered_map<std::string, Edge> edges;
    std::unordered_map<std::string, std::unordered_set<std::string>> forwardLinks, reverseLinks;
    std::unordered_map<std::string, Path> paths;
    std::unordered_map<std::string, std::unordered_set<std::string>> edge2paths;
    std::unordered_map<std::string, std::unordered_set<std::string>> barcode2paths;
    std::unordered_map<std::string, std::unordered_set<std::string>> reads12paths;
    std::unordered_map<std::string, std::unordered_set<std::string>> reads22paths;

    void load_from_spades_fastg(std::string);
    // .fa, .fastg, answer, k_mer
    void load_from_megahit_fastg(std::string, std::string, std::string, std::string);
    void load_from_megahit_metacarvel(std::string);

    void load_spades_paths(std::string, std::string);
    void construct_spades_links_by_edge();

    void load_spades_alignment_cloud(std::string);
    void load_megahit_alignment_cloud(std::string);
    void load_ngs_alignment_megahit(std::string, std::string, std::string);
    void load_ngs_alignment_spades(std::string, std::string, std::string);

    void load_graphbin_result(std::string);

    void construct_barcode_overlap();
    void construct_PE_overlap();
    // insert size
    void construct_PE_overlap_athena(std::string, std::string);

    void breadth_first_search_cloud();
    void breadth_first_search_ngs(std::string);
    void breadth_first_search_metacarvel();

    void bin_evaluate();
};

#endif
