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
    Edge();
    Edge(int, double);
    Edge& operator=(const Edge&);
};

class AlignedRead {
public:
    int startPos, endPos;
    int alignedBases, mappingQuality;
    double blastIdentity;
    AlignedRead();
    AlignedRead& operator=(const AlignedRead&);
};

class Path {
public:
    int length;
    double depth;
    std::vector<std::vector<std::string>> pathEdge;
    std::unordered_map<std::string, std::vector<AlignedRead>> reads1vec;
    std::unordered_map<std::string, std::vector<AlignedRead>> reads2vec;
    std::unordered_set<std::string> linksByEdge;
    std::unordered_map<std::string, int> PEOverlap;
    Path();
    Path& operator=(const Path&);
};

class Graph {
public:
    std::unordered_map<std::string, Edge> edges;
    std::unordered_map<std::string, std::unordered_set<std::string>> forwardLinks, reverseLinks;
    std::unordered_map<std::string, Path> paths;
    std::unordered_map<std::string, std::unordered_set<std::string>> edge2paths;
    std::unordered_map<std::string, std::unordered_set<std::string>> reads12paths;
    std::unordered_map<std::string, std::unordered_set<std::string>> reads22paths;

    void load_from_spades_fastg(std::string);
    void load_spades_paths(std::string);
    void construct_spades_links_by_edge();
    void load_ngs_alignment_spades(std::string, double, double);

    void load_from_megahit_fastg(std::string, std::string);
    void load_ngs_alignment_megahit(std::string, double, double);

    void construct_PE_overlap_athena(int, int);
    void breadth_first_search_ngs(int, std::string);
};

#endif
