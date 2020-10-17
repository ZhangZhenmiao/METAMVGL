#include "graph.h"
#include <map>
#include <set>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <htslib/sam.h>
#include <unordered_map>
#include <unordered_set>

Edge::Edge() {
    this->length = 0;
    this->depth = 0;
}
Edge::Edge(int len, double cov) {
    this->length = len;
    this->depth = cov;
}
Edge& Edge::operator=(const Edge& edge) {
    if (this != &edge) {
        this->length = edge.length;
        this->depth = edge.depth;
    }
    return *this;
}
AlignedRead::AlignedRead() {
    this->startPos = 0;
    this->endPos = 0;
    this->alignedBases = 0;
    this->mappingQuality = 0;
    this->blastIdentity = 0;
}
AlignedRead& AlignedRead::operator=(const AlignedRead& read) {
    if (this != &read) {
        this->startPos = read.startPos;
        this->endPos = read.endPos;
        this->alignedBases = read.alignedBases;
        this->mappingQuality = read.mappingQuality;
        this->blastIdentity = read.blastIdentity;
    }
    return *this;
}
Path::Path() {
    this->depth = 0;
    this->length = 0;
}
Path& Path::operator=(const Path& path) {
    if (this != &path) {
        this->length = path.length;
        this->depth = path.depth;
        this->pathEdge = path.pathEdge;
        this->reads1vec = path.reads1vec;
        this->reads2vec = path.reads2vec;
        this->linksByEdge = path.linksByEdge;
        this->PEOverlap = path.PEOverlap;
    }
    return *this;
}

// load assembly graph from spades fastg
// construct edges (unsingned), forwardLinks (signed), reverseLinks (signed)
void Graph::load_from_spades_fastg(std::string filename) {
    std::cout << "Start loading fastg ..." << std::endl;
    bool isPostive = false;
    int length = 0;
    double depth = 0;
    std::string line, header;
    std::fstream fastgFile(filename.c_str(), std::fstream::in);
    while (getline(fastgFile, line)) {
        // line is empty
        if (line.empty()) continue;
        // line is header
        if (line.at(0) == '>') {
            int flag = 0;
            std::size_t headerLeft;
            std::vector<std::string> link;
            std::vector<bool> linkPostive;
            // update edges
            if (!header.empty()) {
                Edge tmpEdge(length, depth);
                this->edges[header] = tmpEdge;
            }
            // get header and links
            for (std::size_t i = 0; i < line.size(); ++i) {
                // header start
                if (line.at(i) == '_' && flag == 0) {
                    headerLeft = i + 1;
                    flag = 1;
                }
                // header end
                else if (line.at(i) == '_' && flag == 1) {
                    header = line.substr(headerLeft, i - headerLeft);
                    flag = 2;
                }
                // link start after header
                else if (line.at(i) == ':' && flag == 2) {
                    if (line.at(i - 1) == '\'') isPostive = false;
                    else isPostive = true;
                    std::size_t depthIndex = line.find_last_of('_', i - 1);
                    depth = atof(line.substr(depthIndex + 1, i - depthIndex - 1).c_str());
                    std::size_t lengthIndex = line.find_last_of('_', depthIndex - 5);
                    length = atoi(line.substr(lengthIndex + 1, depthIndex - lengthIndex - 5).c_str());
                    flag = 3;
                }
                // link start after link
                else if (line.at(i) == ',' && flag == 5) {
                    if (line.at(i - 1) == '\'') linkPostive.push_back(false);
                    else linkPostive.push_back(true);
                    flag = 3;
                }
                // link header start
                else if (line.at(i) == '_' && flag == 3) {
                    headerLeft = i + 1;
                    flag = 4;
                }
                // link header end
                else if (line.at(i) == '_' && flag == 4) {
                    link.push_back(line.substr(headerLeft, i - headerLeft));
                    flag = 5;
                }
                // end line after header
                else if (line.at(i) == ';' && flag == 2) {
                    if (line.at(i - 1) == '\'') isPostive = false;
                    else isPostive = true;
                    std::size_t depthIndex = line.find_last_of('_', i - 1);
                    depth = atof(line.substr(depthIndex + 1, i - depthIndex - 1).c_str());
                    std::size_t lengthIndex = line.find_last_of('_', depthIndex - 5);
                    length = atoi(line.substr(lengthIndex + 1, depthIndex - lengthIndex - 5).c_str());
                }
                // end line after link
                else if (line.at(i) == ';' && flag == 5) {
                    if (line.at(i - 1) == '\'') linkPostive.push_back(false);
                    else linkPostive.push_back(true);
                }
            }
            // update links
            std::string signedHeader = header;
            std::string signedHeaderReverse = header;
            if (isPostive) {
                signedHeader += '+';
                signedHeaderReverse += '-';
            }
            else {
                signedHeader += '-';
                signedHeaderReverse += '+';
            }
            for (std::size_t i = 0; i < link.size(); ++i) {
                if (linkPostive.at(i)) {
                    this->forwardLinks[signedHeader].insert(link.at(i) + '+');
                    this->reverseLinks[link.at(i) + '+'].insert(signedHeader);
                    this->forwardLinks[link.at(i) + '-'].insert(signedHeaderReverse);
                    this->reverseLinks[signedHeaderReverse].insert(link.at(i) + '-');
                }
                else {
                    this->forwardLinks[signedHeader].insert(link.at(i) + '-');
                    this->reverseLinks[link.at(i) + '-'].insert(signedHeader);
                    this->forwardLinks[link.at(i) + '+'].insert(signedHeaderReverse);
                    this->reverseLinks[signedHeaderReverse].insert(link.at(i) + '+');
                }
            }
        }
    }
    // update remaining edges
    if (!header.empty()) {
        Edge tmpEdge(length, depth);
        this->edges[header] = tmpEdge;
    }
    fastgFile.close();
    std::cout << "Loaded " << this->edges.size() << " edges, ";
    int cnt_f = 0;
    for (auto i : this->forwardLinks) cnt_f += i.second.size();
    std::cout << cnt_f << " links.\nDone." << std::endl;
}

// load spades contigs.paths
// construct length, depth and pathEdge (signed) of paths, edge2paths (signed)
void Graph::load_spades_paths(std::string filename) {
    std::cout << "Start loading spades paths ..." << std::endl;
    std::string line, node;
    std::fstream pathFile(filename.c_str(), std::fstream::in);
    while (getline(pathFile, line)) {
        // line is empty
        if (line.empty()) continue;
        // line is NODE
        if (line.at(0) == 'N') {
            if (line.at(line.size() - 1) == '\'')
                node.clear();
            else {
                node = "NODE_";
                std::size_t i = 5;
                for (; line.at(i) != '_'; ++i)
                    node += line.at(i);
                std::string len;
                for (i += 8; line.at(i) != '_'; ++i)
                    len += line.at(i);
                this->paths[node].length = atoi(len.c_str());
                std::size_t depthIndex = line.find_last_of('_');
                this->paths[node].depth = atof(line.substr(depthIndex + 1, line.size() - depthIndex - 1).c_str());
            }
        }
        // line is path
        else if (!node.empty()) {
            std::string edge;
            std::vector<std::string> tmpPath;
            if (line.at(line.size() - 1) == ';')
                line.erase(line.end() - 1);
            for (std::size_t i = 0; i < line.size(); ++i) {
                if (line.at(i) == ',') {
                    this->edge2paths[edge].insert(node);
                    tmpPath.push_back(edge);
                    edge.clear();
                }
                else {
                    edge += line.at(i);
                }
            }
            this->edge2paths[edge].insert(node);
            tmpPath.push_back(edge);
            this->paths[node].pathEdge.push_back(tmpPath);
        }
    }
    pathFile.close();
    std::cout << "Loaded " << this->paths.size() << " paths.\nDone." << std::endl;
}

// load assembly graph from megahit fastg
// construct length, depth, linksByEdge
void Graph::load_from_megahit_fastg(std::string contigs, std::string filename) {
    // map fastg name to contig name
    std::fstream contigFile(contigs.c_str(), std::fstream::in);
    std::unordered_map<std::string, std::string> nameMap;
    int cnt = 0;
    std::string line;
    while (getline(contigFile, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            std::string node = "NODE_";
            nameMap[node + std::to_string(++cnt)] = line.substr(1, line.find(' ') - 1);
        }
    }
    contigFile.close();

    std::cout << "Start loading fastg ..." << std::endl;
    bool isPostive = false;
    int length = 0;
    double depth = 0;
    std::string header;
    std::fstream fastgFile(filename.c_str(), std::fstream::in);
    while (getline(fastgFile, line)) {
        // line is empty
        if (line.empty()) continue;
        // line is header
        if (line.at(0) == '>') {
            int flag = 0;
            std::size_t headerLeft;
            std::vector<std::string> link;
            std::vector<bool> linkPostive;
            // update nodes
            if (!header.empty()) {
                this->paths[nameMap["NODE_" + header]].length = length;
                this->paths[nameMap["NODE_" + header]].depth = depth;
            }
            // get header and links
            for (std::size_t i = 0; i < line.size(); ++i) {
                // header start
                if (line.at(i) == '_' && flag == 0) {
                    headerLeft = i + 1;
                    flag = 1;
                }
                // header end
                else if (line.at(i) == '_' && flag == 1) {
                    header = line.substr(headerLeft, i - headerLeft);
                    flag = 2;
                }
                // link start after header
                else if (line.at(i) == ':' && flag == 2) {
                    if (line.at(i - 1) == '\'') isPostive = false;
                    else isPostive = true;
                    std::size_t idIndex = line.find_last_of('_', i - 1);
                    std::size_t depthIndex = line.find_last_of('_', idIndex - 4);
                    depth = atof(line.substr(depthIndex + 1, idIndex - depthIndex - 4).c_str());
                    std::size_t lengthIndex = line.find_last_of('_', depthIndex - 5);
                    length = atoi(line.substr(lengthIndex + 1, depthIndex - lengthIndex - 5).c_str());
                    flag = 3;
                }
                // link start after link
                else if (line.at(i) == ',' && flag == 5) {
                    if (line.at(i - 1) == '\'') linkPostive.push_back(false);
                    else linkPostive.push_back(true);
                    flag = 3;
                }
                // link header start
                else if (line.at(i) == '_' && flag == 3) {
                    headerLeft = i + 1;
                    flag = 4;
                }
                // link header end
                else if (line.at(i) == '_' && flag == 4) {
                    link.push_back(line.substr(headerLeft, i - headerLeft));
                    flag = 5;
                }
                // end line after header
                else if (line.at(i) == ';' && flag == 2) {
                    if (line.at(i - 1) == '\'') isPostive = false;
                    else isPostive = true;
                    std::size_t idIndex = line.find_last_of('_', i - 1);
                    std::size_t depthIndex = line.find_last_of('_', idIndex - 4);
                    depth = atof(line.substr(depthIndex + 1, idIndex - depthIndex - 4).c_str());
                    std::size_t lengthIndex = line.find_last_of('_', depthIndex - 5);
                    length = atoi(line.substr(lengthIndex + 1, depthIndex - lengthIndex - 5).c_str());
                }
                // end line after link
                else if (line.at(i) == ';' && flag == 5) {
                    if (line.at(i - 1) == '\'') linkPostive.push_back(false);
                    else linkPostive.push_back(true);
                }
            }
            // update links
            for (std::size_t i = 0; i < link.size(); ++i) {
                if (header.compare(link.at(i))) {
                    this->paths[nameMap["NODE_" + header]].linksByEdge.insert(nameMap["NODE_" + link.at(i)]);
                    this->paths[nameMap["NODE_" + link.at(i)]].linksByEdge.insert(nameMap["NODE_" + header]);
                }
            }
        }
    }
    // update remaining nodes
    if (!header.empty()) {
        this->paths[nameMap["NODE_" + header]].length = length;
        this->paths[nameMap["NODE_" + header]].depth = depth;
    }
    fastgFile.close();
    std::cout << "Loaded " << this->paths.size() << " nodes, ";
    int cnt_f = 0;
    for (auto i : this->paths)
        cnt_f += i.second.linksByEdge.size();
    std::cout << cnt_f << " links.\nDone." << std::endl;
}

// load bam from spades alignment without barcode
// construct reads1vec and reads2vec of paths, reads12paths, reads22paths
void Graph::load_ngs_alignment_spades(std::string filename, double mq, double idt) {
    const int supplementaryAlignment = 0x00000800, secondaryAlignment = 0x00000100;
    const int Read1 = 0x00000040, Read2 = 0x00000080;
    samFile *bamFile = hts_open(filename.c_str(), "r");
    bam_hdr_t *bamHeader = sam_hdr_read(bamFile);
    bam1_t *bamAln = bam_init1();
    int cntLine = 0;
    std::cout << "Start parsing alignments ..." << std::endl;
    while (sam_read1(bamFile, bamHeader, bamAln) >= 0) {
        int flag = bamAln->core.flag;
        bool isSupplementary = flag & supplementaryAlignment, isSecondary = flag & secondaryAlignment;
        bool isRead1 = flag & Read1, isRead2 = flag & Read2;
        if (isSupplementary || isSecondary) continue;

        int pos = bamAln->core.pos + 1, endPos = bam_endpos(bamAln) + 1;
        if (endPos == pos + 1) continue;

        std::string readName, contigName;
        if (bam_get_qname(bamAln) != NULL) readName = bam_get_qname(bamAln);
        if (readName.size() == 0) continue;
        if (sam_hdr_tid2name(bamHeader, bamAln->core.tid) != NULL) contigName = bamHeader->target_name[bamAln->core.tid];
        if (contigName.size() == 0) continue;

        std::string node = "NODE_";
        for (std::size_t i = 5; contigName.at(i) != '_'; ++i) node += contigName.at(i);
        
        int alignmentColumns = 0, NM = 0;
        uint32_t *cigarPointer = bam_get_cigar(bamAln);
        for (int i = 0; i < bamAln->core.n_cigar; ++i) {
            char cigarOperator = bam_cigar_opchr(cigarPointer[i]);
            if (cigarOperator == 'M' || cigarOperator == 'I' || cigarOperator == 'D')
                alignmentColumns += bam_cigar_oplen(cigarPointer[i]);
        }
        if (alignmentColumns == 0) continue;
        uint8_t *tmpNM = bam_aux_get(bamAln, "NM");
        if (tmpNM != NULL) NM = bam_aux2i(tmpNM);
        double blastIdentity = 1.0 * (alignmentColumns - NM) / alignmentColumns;
        if (blastIdentity < idt) continue;

        int mappingQuality = bamAln->core.qual;
        if (mappingQuality < mq) continue;

        // update paths
        AlignedRead tmpRead;
        tmpRead.startPos = pos;
        tmpRead.endPos = endPos;
        tmpRead.alignedBases = alignmentColumns;
        tmpRead.mappingQuality = mappingQuality;
        tmpRead.blastIdentity = blastIdentity;
        if (isRead1) {
            this->paths[node].reads1vec[readName].push_back(tmpRead);
            this->reads12paths[readName].insert(node);
        }
        if (isRead2) {
            this->paths[node].reads2vec[readName].push_back(tmpRead);
            this->reads22paths[readName].insert(node);
        }
        // logging
        if (++cntLine % 10000000 == 0)
            std::cout << "Parsed " << cntLine << " lines ..." << std::endl;
    }
    if (cntLine % 10000000)
        std::cout << "Parsed " << cntLine << " lines ...\nDone." << std::endl;
    else
        std::cout << "Done." << std::endl;
    bam_hdr_destroy(bamHeader);
    bam_destroy1(bamAln);
    sam_close(bamFile);
}

// load bam from megahit alignment without barcode
// construct reads1vec and reads2vec of paths, reads12paths, reads22paths
void Graph::load_ngs_alignment_megahit(std::string filename, double mq, double idt) {
    const int supplementaryAlignment = 0x00000800, secondaryAlignment = 0x00000100;
    const int Read1 = 0x00000040, Read2 = 0x00000080;
    samFile *bamFile = hts_open(filename.c_str(), "r");
    bam_hdr_t *bamHeader = sam_hdr_read(bamFile);
    bam1_t *bamAln = bam_init1();
    int cntLine = 0;
    std::cout << "Start parsing alignments ..." << std::endl;
    while (sam_read1(bamFile, bamHeader, bamAln) >= 0) {
        int flag = bamAln->core.flag;
        bool isSupplementary = flag & supplementaryAlignment, isSecondary = flag & secondaryAlignment;
        bool isRead1 = flag & Read1, isRead2 = flag & Read2;
        if (isSupplementary || isSecondary) continue;

        int pos = bamAln->core.pos + 1, endPos = bam_endpos(bamAln) + 1;
        if (endPos == pos + 1) continue;

        std::string readName, contigName;
        if (bam_get_qname(bamAln) != NULL) readName = bam_get_qname(bamAln);
        if (readName.size() == 0) continue;
        if (sam_hdr_tid2name(bamHeader, bamAln->core.tid) != NULL) contigName = bamHeader->target_name[bamAln->core.tid];
        if (contigName.size() == 0) continue;
        std::string node = contigName;

        int alignmentColumns = 0, NM = 0;
        uint32_t *cigarPointer = bam_get_cigar(bamAln);
        for (int i = 0; i < bamAln->core.n_cigar; ++i) {
            char cigarOperator = bam_cigar_opchr(cigarPointer[i]);
            if (cigarOperator == 'M' || cigarOperator == 'I' || cigarOperator == 'D')
                alignmentColumns += bam_cigar_oplen(cigarPointer[i]);
        }
        if (alignmentColumns == 0) continue;
        uint8_t *tmpNM = bam_aux_get(bamAln, "NM");
        if (tmpNM != NULL) NM = bam_aux2i(tmpNM);
        double blastIdentity = 1.0 * (alignmentColumns - NM) / alignmentColumns;
        if (blastIdentity < idt) continue;

        int mappingQuality = bamAln->core.qual;
        if (mappingQuality < mq) continue;

        // update paths
        AlignedRead tmpRead;
        tmpRead.startPos = pos;
        tmpRead.endPos = endPos;
        tmpRead.alignedBases = alignmentColumns;
        tmpRead.mappingQuality = mappingQuality;
        tmpRead.blastIdentity = blastIdentity;
        if (isRead1) {
            this->paths[node].reads1vec[readName].push_back(tmpRead);
            this->reads12paths[readName].insert(node);
        }
        if (isRead2) {
            this->paths[node].reads2vec[readName].push_back(tmpRead);
            this->reads22paths[readName].insert(node);
        }
        // logging
        if (++cntLine % 10000000 == 0)
            std::cout << "Parsed " << cntLine << " lines ..." << std::endl;
    }
    if (cntLine % 10000000)
        std::cout << "Parsed " << cntLine << " lines ...\nDone." << std::endl;
    else
        std::cout << "Done." << std::endl;
    bam_hdr_destroy(bamHeader);
    bam_destroy1(bamAln);
    sam_close(bamFile);
}

// construct linksByEdge of paths from spades edges
void Graph::construct_spades_links_by_edge() {
    std::cout << "Start constructing links between paths ..." << std::endl;
    int cnt = 0, sum = 0;
    for (auto &i : this->paths) {
        for (auto j : i.second.pathEdge) {
            for (auto edge : j) {
                std::string unsignedEdge = edge.substr(0, edge.size() - 1);
                std::unordered_set<std::string> tmpEdge1 = this->forwardLinks[unsignedEdge + '+'];
                std::unordered_set<std::string> tmpEdge2 = this->reverseLinks[unsignedEdge + '+'];
                std::unordered_set<std::string> tmpEdge3 = this->forwardLinks[unsignedEdge + '-'];
                std::unordered_set<std::string> tmpEdge4 = this->reverseLinks[unsignedEdge + '-'];
                // tmpEdge1.insert(unsignedEdge + '+');
                // tmpEdge1.insert(unsignedEdge + '-');
                tmpEdge1.insert(tmpEdge2.begin(), tmpEdge2.end());
                tmpEdge1.insert(tmpEdge3.begin(), tmpEdge3.end());
                tmpEdge1.insert(tmpEdge4.begin(), tmpEdge4.end());
                for (auto neighbour : tmpEdge1) {
                    for (auto path : this->edge2paths[neighbour])
                        if (i.first.compare(path)) i.second.linksByEdge.insert(path);
                }
            }
        }
        sum += i.second.linksByEdge.size();
        if (++cnt % 10000 == 0) {
            std::cout << "Constructed links of " << cnt << " paths, new links " << sum << " ..." << std::endl;
            sum = 0;
        }
    }
    if (cnt % 10000)
        std::cout << "Constructed links of " << cnt << " paths, new links " << sum << " ...\nDone." << std::endl;
    else
        std::cout << "Done." << std::endl;
}

// construct PEOverlap of paths using athena method
void Graph::construct_PE_overlap_athena(int insert_size, int pe_num) {
    std::cout << "Start constructing PE overlap ..." << std::endl;
    int cnt_l = 0, sum = 0;
    for (auto& i : this->paths) {
        if (!i.second.reads1vec.size() && !i.second.reads2vec.size()) continue;
        std::unordered_map<std::string, std::vector<AlignedRead>> reads_i, reads_j;
        for (auto ar : i.second.reads1vec) {
            for (auto j : this->reads22paths[ar.first]) {
                if (i.first.compare(j)) {
                    reads_i[j].insert(reads_i[j].end(), ar.second.begin(), ar.second.end());
                    reads_j[j].insert(reads_j[j].end(), this->paths[j].reads2vec[ar.first].begin(), this->paths[j].reads2vec[ar.first].end());
                }
            }
        }
        for (auto ar : i.second.reads2vec) {
            for (auto j : this->reads12paths[ar.first]) {
                if (i.first.compare(j)) {
                    reads_i[j].insert(reads_i[j].end(), ar.second.begin(), ar.second.end());
                    reads_j[j].insert(reads_j[j].end(), this->paths[j].reads1vec[ar.first].begin(), this->paths[j].reads1vec[ar.first].end());
                }
            }
        }
        for (auto ar : reads_i) {
            std::map<int, int> count_i, count_j;
            for (std::size_t r1 = 0; r1 < ar.second.size(); ++r1) {
                for (std::size_t r2 = 0; r2 < ar.second.size(); ++r2) {
                    if (ar.second[r2].startPos - insert_size < ar.second[r1].endPos && ar.second[r1].startPos < ar.second[r2].endPos + insert_size)
                        count_i[r2] += 1;
                }
            }
            for (std::size_t r1 = 0; r1 < reads_j[ar.first].size(); ++r1) {
                for (std::size_t r2 = 0; r2 < reads_j[ar.first].size(); ++r2) {
                    if (reads_j[ar.first][r2].startPos - insert_size < reads_j[ar.first][r1].endPos && reads_j[ar.first][r1].startPos < reads_j[ar.first][r2].endPos + insert_size)
                        count_j[r2] += 1;
                }
            }
            int overl1 = 0, overl2 = 0;
            for (auto cnt : count_i) {
                if (cnt.second > overl1) overl1 = cnt.second;
            }
            for (auto cnt : count_j) {
                if (cnt.second > overl2) overl2 = cnt.second;
            }
            if (overl1 >= pe_num && overl2 >= pe_num && overl1 * 2 >= count_i.size() && overl2 * 2 >= count_j.size()) {
                i.second.PEOverlap[ar.first] = overl2;
                this->paths[ar.first].PEOverlap[i.first] = overl1;
                sum += overl2;
            }
            cnt_l += 1;
            if (cnt_l % 1000 == 0) {
                std::cout << "Constructed PE overlap of " << cnt_l << " paths, new overlaps " << sum << " ..." << std::endl;
                sum = 0;
            }
        }
    }
    if (cnt_l % 1000)
        std::cout << "Constructed PE overlap of " << cnt_l << " paths, new overlaps " << sum << " ...\nDone." << std::endl;
    else
        std::cout << "Done." << std::endl;
}

// output the two graphs
void Graph::breadth_first_search_ngs(int pe_num, std::string output) {
    std::cout << "Start traversing assembly graph ..." << std::endl;
    int cnt = 0;
    std::fstream edgeGraph(output + ".ag", std::fstream::out);
    for (auto i : this->paths) {
        edgeGraph << i.first << " depth" << i.second.depth << " length" << i.second.length << ":\n";
        for (auto j : i.second.linksByEdge) {
            edgeGraph << '\t' << j << " depth" << this->paths[j].depth << " length" << this->paths[j].length << " overlapPE" << i.second.PEOverlap[j] << ";\n";
        }
        if (++cnt % 10000 == 0)
            std::cout << "Traversed edge neighbours of " << cnt << " paths ..." << std::endl;
    }
    if (cnt % 10000)
        std::cout << "Traversed edge neighbours of " << cnt << " paths ...\nDone." << std::endl;
    edgeGraph.close();

    std::cout << "Start traversing PE graph ..." << std::endl;
    cnt = 0;
    std::fstream PEGraph(output + ".pe", std::fstream::out);
    for (auto i : this->paths) {
        PEGraph << i.first << " depth" << i.second.depth << " length" << i.second.length  << ":\n";
        for (auto j : i.second.PEOverlap) {
            if (j.second < pe_num) continue;
            PEGraph << '\t' << j.first << " depth" << this->paths[j.first].depth << " length" << this->paths[j.first].length  << " overlapPE" << j.second << ";\n";
        }
        if (++cnt % 10000 == 0)
            std::cout << "Traversed PE neighbours of " << cnt << " paths ..." << std::endl;
    }
    if (cnt % 10000)
        std::cout << "Traversed PE neighbours of " << cnt << " paths ..." << std::endl;
    std::cout << "Done." << std::endl;
    PEGraph.close();
}