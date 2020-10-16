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

// Edge constructer
Edge::Edge() {
    this->length = 0;
    this->depth = 0;
}
Edge::Edge(std::string forward, bool reverse, double cover, int len) {
    this->forwardSequence = forward;
    this->reverseSequence = this->reverseComplementary(forward);
    this->depth = cover;
    this->length = len;
}
Edge::Edge(bool forward, std::string reverse, double cover, int len) {
    this->reverseSequence = reverse;
    this->forwardSequence = this->reverseComplementary(reverse);
    this->depth = cover;
    this->length = len;
}
std::string Edge::reverseComplementary(std::string s) {
    std::string rc;
    for (int i = s.size() - 1; i >= 0; --i) {
        if (s.at(i) == 'A') rc += 'T';
        else if (s.at(i) == 'T') rc += 'A';
        else if (s.at(i) == 'C') rc += 'G';
        else if (s.at(i) == 'G') rc += 'C';
        else if (s.at(i) == 'N') rc += 'N';
        else rc += 'N';
    }
    return rc;
}
Edge& Edge::operator=(const Edge& edge) {
    if (this != &edge) {
        this->forwardSequence = edge.forwardSequence;
        this->reverseSequence = edge.reverseSequence;
        this->depth = edge.depth;
        this->length = edge.length;
    }
    return *this;
}

// Path constructer
Path::Path() {
    this->depth = 0;
    this->length = 0;
    this->bin = 0;
}
Path::Path(std::vector<std::vector<std::string>> edge, double cover, int len) {
    this->pathEdge = edge;
    this->depth = cover;
    this->length = len;
}
Path& Path::operator=(const Path& path) {
    if (this != &path) {
        this->pathEdge = path.pathEdge;
        this->depth = path.depth;
        this->length = path.length;
        this->bin = path.bin;
        this->reads1vec = path.reads1vec;
        this->reads2vec = path.reads2vec;
        this->barcode2reads = path.barcode2reads;
        this->linksByEdge = path.linksByEdge;
        this->barcodeOverlap = path.barcodeOverlap;
        this->PEOverlap = path.PEOverlap;
        this->annotation = path.annotation;
    }
    return *this;
}

// AlignedRead constructer
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
        this->sequence = read.sequence;
        this->barcode = read.barcode;
        this->alignedBases = read.alignedBases;
        this->mappingQuality = read.mappingQuality;
        this->blastIdentity = read.blastIdentity;
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
    std::string line, header, sequence;
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
                if (isPostive) {
                    Edge tmpEdge(sequence, false, depth, length);
                    this->edges[header] = tmpEdge;
                }
                else {
                    Edge tmpEdge(false, sequence, depth, length);
                    this->edges[header] = tmpEdge;
                }
                sequence.clear();
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
        // line is sequence
        else {
            sequence += line;
        }
    }
    // update remaining edges
    if (!header.empty()) {
        if (isPostive) {
            Edge tmpEdge(sequence, false, depth, length);
            this->edges[header] = tmpEdge;
        }
        else {
            Edge tmpEdge(false, sequence, depth, length);
            this->edges[header] = tmpEdge;
        }
        sequence.clear();
    }
    fastgFile.close();
    std::cout << "Loaded " << this->edges.size() << " edges, ";
    int cnt_f = 0;
    for (auto i : this->forwardLinks) cnt_f += i.second.size();
    std::cout << cnt_f << " links.\nDone." << std::endl;
}

// load spades contigs.paths
// construct length, depth and pathEdge (signed) of paths, edge2paths (signed)
void Graph::load_spades_paths(std::string filename, std::string answer) {
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

    std::cout << "Start parsing annotation ..." << std::endl;
    std::fstream answerFile(answer.c_str(), std::fstream::in);
    if (answer.find("kraken") == std::string::npos) {
        while (getline(answerFile, line)) {
            if (line.empty()) continue;
            if (line.at(0) != 'N' || line.at(1) != 'O' || line.at(2) != 'D' || line.at(3) != 'E') continue;
            node = "NODE_";
            std::size_t i = 5;
            for (; line.at(i) != '_'; ++i)
                node += line.at(i);
            while (line.at(i) != '\t') ++i;
            while (i < line.size() && line.at(i) == '\t') {
                std::string anno;
                for (++i; i < line.size() && line.at(i) != '\t'; ++i) anno += line.at(i);
                this->paths[node].annotation.push_back(anno);
            }
        }
    }
    else {
        std::map<std::string, std::string> taxid2species;
        std::fstream taxid("/home/comp/zmzhang/software/kraken2_install/parse_inspect.txt", std::fstream::in);
        while (getline(taxid, line)) {
            if (line.empty()) continue;
            std::size_t i = 0;
            std::string tax;
            while (line.at(i) != '\t') tax += line.at(i++);
            std::string species;
            for (++i; i < line.size(); ++i) species += line.at(i);
            taxid2species[tax] = species;
        }
        taxid.close();
        std::set<std::string> all_annos;
        while (getline(answerFile, line)) {
            if (line.empty()) continue;
            if (line.at(0) != 'C') continue;
            std::set<std::string> annotation;
            std::string node;
            std::string maintax;
            std::size_t i = 7;
            node = "NODE_";
            for (; line.at(i) != '_'; ++i) node += line.at(i);
            while (line.at(i) != '\t') i++;
            for (++i; line.at(i) != '\t'; ++i) maintax += line.at(i);
            for (++i; line.at(i) != '\t'; ++i);
            while (i < line.size() && (line.at(i) == '\t' || line.at(i) == ' ')) {
                std::string anno;
                for (++i; line.at(i) != ':'; ++i) anno += line.at(i);
                while (i < line.size() && line.at(i) != ' ') ++i;
                if (anno.compare("0")) {
                    anno = taxid2species[anno];
                    if (anno.compare("0") && anno.compare("9606")) annotation.insert(anno);
                }
            }
            std::cout << node;
            if (maintax.compare("0")) {
                maintax = taxid2species[maintax];
                if (maintax.compare("0") && maintax.compare("9606")) {
                    this->paths[node].annotation.push_back(maintax);
                    all_annos.insert(maintax);
                    std::cout << ' ' << maintax;
                }
            }
            std::cout << std::endl;
        }
        std::cout << "All annotations: " << all_annos.size() << std::endl;
    }
    answerFile.close();
    std::cout << "Done." << std::endl;
} // Graph load spades paths

// load assembly graph from megahit fastg
// construct length, depth, linksByEdge and annotation (kraken filename must contain kraken) of paths
void Graph::load_from_megahit_fastg(std::string contigs, std::string filename, std::string answer, std::string k_mer) {
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
    std::string header, sequence;
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
                sequence.clear();
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
        // line is sequence
        else {
            sequence += line;
        }
    }
    // update remaining nodes
    if (!header.empty()) {
        this->paths[nameMap["NODE_" + header]].length = length;
        this->paths[nameMap["NODE_" + header]].depth = depth;
        sequence.clear();
    }
    fastgFile.close();
    std::cout << "Loaded " << this->paths.size() << " nodes, ";
    int cnt_f = 0;
    for (auto i : this->paths)
        cnt_f += i.second.linksByEdge.size();
    std::cout << cnt_f << " links.\nDone." << std::endl;

    std::cout << "Start parsing annotation ..." << std::endl;
    std::fstream answerFile(answer.c_str(), std::fstream::in);
    if (answer.find("kraken") == std::string::npos) {
        while (getline(answerFile, line)) {
            if (line.empty()) continue;
            bool is_contig = true;
            if (line.at(0) != 'k') continue;
            for (std::size_t it_k = 0; it_k < k_mer.size(); ++it_k) {
                if (is_contig && line.at(it_k + 1) != k_mer.at(it_k)) is_contig = false;
            }
            if (!is_contig) continue;
            std::string node = "k" + k_mer + "_";
            std::size_t i = k_mer.size() + 2;
            for (; line.at(i) != '\t'; ++i) node += line.at(i);
            while (i < line.size() && line.at(i) == '\t') {
                std::string anno;
                for (++i; i < line.size() && line.at(i) != '\t'; ++i) anno += line.at(i);
                this->paths[node].annotation.push_back(anno);
            }
        }
    }
    else {
        std::map<std::string, std::string> taxid2species;
        std::fstream taxid("/home/comp/zmzhang/software/kraken2_install/parse_inspect.txt", std::fstream::in);
        while (getline(taxid, line)) {
            if (line.empty()) continue;
            std::size_t i = 0;
            std::string tax;
            while (line.at(i) != '\t') tax += line.at(i++);
            std::string species;
            for (++i; i < line.size(); ++i) species += line.at(i);
            taxid2species[tax] = species;
        }
        taxid.close();
        std::set<std::string> all_annos;
        while (getline(answerFile, line)) {
            if (line.empty()) continue;
            if (line.at(0) != 'C') continue;
            std::set<std::string> annotation;
            std::string node;
            std::string maintax;
            std::size_t i = 2;
            while (line.at(i) != '\t') node += line.at(i++);
            for (++i; line.at(i) != '\t'; ++i) maintax += line.at(i);
            for (++i; line.at(i) != '\t'; ++i);
            while (i < line.size() && (line.at(i) == '\t' || line.at(i) == ' ')) {
                std::string anno;
                for (++i; line.at(i) != ':'; ++i) anno += line.at(i);
                while (i < line.size() && line.at(i) != ' ') ++i;
                if (anno.compare("0")) {
                    anno = taxid2species[anno];
                    if (anno.compare("0") && anno.compare("9606")) annotation.insert(anno);
                }
            }
            std::cout << node;
            if (maintax.compare("0")) {
                maintax = taxid2species[maintax];
                if (maintax.compare("0") && maintax.compare("9606")) {
                    this->paths[node].annotation.push_back(maintax);
                    all_annos.insert(maintax);
                    std::cout << ' ' << maintax;
                }
            }
            std::cout << std::endl;
        }
        std::cout << "All annotations: " << all_annos.size() << std::endl;
    }
    answerFile.close();
    std::cout << "Done." << std::endl;
}

// load PE gfa from metacarvel contig_links (megahit assembly)
// construct PEOverlap of paths
void Graph::load_from_megahit_metacarvel(std::string filename) {
    std::string line;
    std::fstream metacarvelFile(filename.c_str(), std::fstream::in);
    while (getline(metacarvelFile, line)) {
        if (line.empty()) continue;
        std::size_t i = 0;
        std::string node1, node2;
        for (; line.at(i) != '\t'; ++i) node1 += line.at(i);
        for (i += 3; line.at(i) != '\t'; ++i) node2 += line.at(i);
        this->paths[node1].PEOverlap[node2] = 1;
        this->paths[node2].PEOverlap[node1] = 1;
    }
    metacarvelFile.close();
}

// load bam from spades alignment without barcode
// construct reads1vec and reads2vec of paths, reads12paths, reads22paths
void Graph::load_ngs_alignment_spades(std::string filename, std::string mq, std::string idt) {
    int mappingq = atoi(mq.c_str());
    double identity = atof(idt.c_str());
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
        int mappingQuality = bamAln->core.qual;
        std::string readName = bam_get_qname(bamAln);
        if (readName.size() == 0 || readName.compare("*") == 0) continue;
        std::string contigName = bamHeader->target_name[bamAln->core.tid];
        if (contigName.size() == 0 || contigName.compare("*") == 0) continue;
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
        if (blastIdentity < identity) continue;
        if (mappingQuality < mappingq) continue;
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

// load bam from spades alignment with barcode
// construct reads1vec, reads2vec and barcode2reads of paths, reads12paths, reads22paths, barcode2paths
void Graph::load_spades_alignment_cloud(std::string filename) {
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
        int mappingQuality = bamAln->core.qual;
        std::string readName = bam_get_qname(bamAln);
        if (readName.size() == 0 || readName.compare("*") == 0) continue;
        std::string contigName = bamHeader->target_name[bamAln->core.tid];
        if (contigName.size() == 0 || contigName.compare("*") == 0) continue;
        std::string node = "NODE_";
        for (std::size_t i = 5; contigName.at(i) != '_'; ++i) node += contigName.at(i);
        // int contigLength = bamHeader->target_len[bamAln->core.tid];
        /*std::string readSeq;
        uint8_t *readPointer = bam_get_seq(bamAln);
        for (int i = 0; i < bamAln->core.l_qseq; ++i) {
            int base = bam_seqi(readPointer, i);
            if (base == 1) readSeq += 'A';
            else if (base == 2) readSeq += 'C';
            else if (base ==  4) readSeq += 'G';
            else if (base == 8) readSeq += 'T';
            else if (base == 15) readSeq += 'N';
            else std::cout<<"Read contains non-ATCGN characters."<<std::endl;
        }*/
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
        if (blastIdentity < 0.98) continue;
        if (mappingQuality < 20) continue;
        std::string barcode;
        uint8_t *tmpBC = bam_aux_get(bamAln, "BX");
        if (tmpBC != NULL) barcode = bam_aux2Z(tmpBC);
        // update paths
        AlignedRead tmpRead;
        tmpRead.startPos = pos;
        tmpRead.endPos = endPos;
        // tmpRead.sequence = readSeq;
        tmpRead.barcode = barcode;
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
        this->paths[node].barcode2reads[barcode].insert(readName);
        this->barcode2paths[barcode].insert(node);
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
void Graph::load_ngs_alignment_megahit(std::string filename, std::string mq, std::string idt) {
    int mappingq = atoi(mq.c_str());
    double identity = atof(idt.c_str());
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
        int mappingQuality = bamAln->core.qual;
        std::string readName = bam_get_qname(bamAln);
        if (readName.size() == 0 || readName.compare("*") == 0) continue;
        std::string contigName = bamHeader->target_name[bamAln->core.tid];
        if (contigName.size() == 0 || contigName.compare("*") == 0) continue;
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
        if (blastIdentity < identity) continue;
        if (mappingQuality < mappingq) continue;
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

// load bam from megahit alignment with barcode
// construct reads1vec, reads2vec and barcode2reads of paths, reads12paths, reads22paths, barcode2paths
void Graph::load_megahit_alignment_cloud(std::string filename) {
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
        int mappingQuality = bamAln->core.qual;
        std::string readName = bam_get_qname(bamAln);
        if (readName.size() == 0 || readName.compare("*") == 0) continue;
        std::string contigName = bamHeader->target_name[bamAln->core.tid];
        if (contigName.size() == 0 || contigName.compare("*") == 0) continue;
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
        if (blastIdentity < 0.98) continue;
        if (mappingQuality < 20) continue;
        std::string barcode;
        uint8_t *tmpBC = bam_aux_get(bamAln, "BX");
        if (tmpBC != NULL) barcode = bam_aux2Z(tmpBC);
        // update paths
        AlignedRead tmpRead;
        tmpRead.startPos = pos;
        tmpRead.endPos = endPos;
        tmpRead.barcode = barcode;
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
        this->paths[node].barcode2reads[barcode].insert(readName);
        this->barcode2paths[barcode].insert(node);
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

// construct barcodeOverlap of paths, requiring supporting reads >= 5
void Graph::construct_barcode_overlap() {
    std::cout << "Start constructing barcode overlap ..." << std::endl;
    int cnt = 0, sum = 0;
    for (auto& i : this->paths) {
        for (auto j : i.second.barcode2reads) {
            if (!j.first.empty() && j.first.compare("0_0_0") && j.second.size() >= 5) {
                for (auto path : this->barcode2paths[j.first])
                    if (i.first.compare(path) && this->paths[path].barcode2reads[j.first].size() >= 5) i.second.barcodeOverlap[path] += 1;
            }
        }
        sum += i.second.barcodeOverlap.size();
        if (++cnt % 10000 == 0) {
            std::cout << "Constructed barcode overlap of " << cnt << " paths, new overlaps " << sum << " ..." << std::endl;
            sum = 0;
        }
    }
    if (cnt % 10000)
        std::cout << "Constructed barcode overlap of " << cnt << " paths, new overlaps " << sum << " ...\nDone." << std::endl;
    else
        std::cout << "Done." << std::endl;
}

// construct PEOverlap of paths
void Graph::construct_PE_overlap() {
    std::cout << "Start constructing PE overlap ..." << std::endl;
    int cnt = 0, sum = 0;
    for (auto& i : this->paths) {
        for (auto j : i.second.reads1vec) {
            for (auto path : this->reads22paths[j.first])
                if (i.first.compare(path)) i.second.PEOverlap[path] += 1;
        }
        for (auto j : i.second.reads2vec) {
            for (auto path : this->reads12paths[j.first])
                if (i.first.compare(path)) i.second.PEOverlap[path] += 1;
        }
        sum += i.second.PEOverlap.size();
        if (++cnt % 10000 == 0) {
            std::cout << "Constructed PE overlap of " << cnt << " paths, new overlaps " << sum << " ..." << std::endl;
            sum = 0;
        }
    }
    if (cnt % 10000)
        std::cout << "Constructed PE overlap of " << cnt << " paths, new overlaps " << sum << " ...\nDone." << std::endl;
    else
        std::cout << "Done." << std::endl;
}

// construct PEOverlap of paths using athena method
void Graph::construct_PE_overlap_athena(std::string insert, std::string num) {
    int insert_size = atoi(insert.c_str());
    int pe_num = atoi(num.c_str());
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
            // std::cout << i.first << "\t" << ar.first << "\t" << overl1 << "/" << count_i.size() << "\t" << overl2 << "/" << count_j.size() << std::endl;
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

// load binning results
void Graph::load_graphbin_result(std::string filename) {
    std::cout << "Start loading graphbin results ..." << std::endl;
    std::string line;
    std::fstream graphbin(filename.c_str(), std::fstream::in);
    while (getline(graphbin, line)) {
        std::string node;
        std::size_t i = 0;
        for (; line.at(i) != ','; ++i)
            node += line.at(i);
        std::string tmpBin;
        for (++i; i < line.size(); ++i)
            tmpBin += line.at(i);
        this->paths[node].bin = atoi(tmpBin.c_str());
    }
    graphbin.close();
    std::cout << "Done." << std::endl;
}

// edge_graph.txt (with Size bin0 -> # !bin0)
// barcode_graph.txt (with Deadend)
// PE_graph.txt (with Deadend)
// without annotation
void Graph::breadth_first_search_cloud() {
    std::cout << "Start traversing edge graph ..." << std::endl;
    int cnt = 0;
    std::fstream edgeGraph("edge_graph.txt", std::fstream::out);
    for (auto i : this->paths) {
        edgeGraph << i.first << " bin" << i.second.bin << " depth" << i.second.depth << " length" << i.second.length << ":\n";
        std::unordered_set<int> tmpBin;
        for (auto j : i.second.linksByEdge) {
            int binNum = this->paths[j].bin;
            if (i.second.bin == 0 && binNum != 0) tmpBin.insert(binNum);
            edgeGraph << '\t' << j << " bin" << binNum << " depth" << this->paths[j].depth << " length" << this->paths[j].length << " overlapPE" << i.second.PEOverlap[j] << " overlapBC" << i.second.barcodeOverlap[j] << ";\n";
        }
        if (i.second.bin == 0 && tmpBin.size() >= 2) edgeGraph << "\tSize " << tmpBin.size() << ";\n";
        if (++cnt % 10000 == 0)
            std::cout << "Searched edge neighbours of " << cnt << " paths ..." << std::endl;
    }
    if (cnt % 10000)
        std::cout << "Searched edge neighbours of " << cnt << " paths ...\nDone." << std::endl;
    else
        std::cout << "Done." << std::endl;
    edgeGraph.close();
    std::cout << "Start traversing barcode graph ..." << std::endl;
    cnt = 0;
    std::fstream barcodeGraph("barcode_graph.txt", std::fstream::out);
    for (auto i : this->paths) {
        barcodeGraph << i.first << " bin" << i.second.bin << " depth" << i.second.depth << " length" << i.second.length  << ":\n";
        int neighbourNum = i.second.linksByEdge.size();
        for (auto j : i.second.barcodeOverlap) {
            barcodeGraph << '\t' << j.first << " bin" << this->paths[j.first].bin << " depth" << this->paths[j.first].depth << " length" << this->paths[j.first].length  << " overlapPE" << i.second.PEOverlap[j.first] << " overlapBC" << j.second;
            int flag1 = 0, flag2 = 0;
            if (neighbourNum == 0) ++flag1;
            if (this->paths[j.first].linksByEdge.size() == 0) ++flag2;
            if (flag1 && !flag2)
                barcodeGraph << " Deadend_left;\n";
            else if (!flag1 && flag2)
                barcodeGraph << " Deadend_right;\n";
            else if (flag1 && flag2)
                barcodeGraph << " Deadend_both;\n";
            else
                barcodeGraph << ";\n";
        }
        if (++cnt % 10000 == 0)
            std::cout << "Searched barcode neighbours of " << cnt << " paths ..." << std::endl;
    }
    if (cnt % 10000)
        std::cout << "Searched barcode neighbours of " << cnt << " paths ...\nDone." << std::endl;
    else
        std::cout << "Done." << std::endl;
    barcodeGraph.close();
    std::cout << "Start traversing PE graph ..." << std::endl;
    cnt = 0;
    std::fstream PEGraph("PE_graph.txt", std::fstream::out);
    for (auto i : this->paths) {
        PEGraph << i.first << " bin" << i.second.bin << " depth" << i.second.depth << " length" << i.second.length  << ":\n";
        int neighbourNum = i.second.linksByEdge.size();
        for (auto j : i.second.PEOverlap) {
            PEGraph << '\t' << j.first << " bin" << this->paths[j.first].bin << " depth" << this->paths[j.first].depth << " length" << this->paths[j.first].length  << " overlapPE" << j.second << " overlapBC" << i.second.barcodeOverlap[j.first];
            int flag1 = 0, flag2 = 0;
            if (neighbourNum == 0) ++flag1;
            if (this->paths[j.first].linksByEdge.size() == 0) ++flag2;
            if (flag1 && !flag2)
                PEGraph << " Deadend_left;\n";
            else if (!flag1 && flag2)
                PEGraph << " Deadend_right;\n";
            else if (flag1 && flag2)
                PEGraph << " Deadend_both;\n";
            else
                PEGraph << ";\n";
        }
        if (++cnt % 10000 == 0)
            std::cout << "Searched PE neighbours of " << cnt << " paths ..." << std::endl;
    }
    if (cnt % 10000)
        std::cout << "Searched PE neighbours of " << cnt << " paths ...\nDone." << std::endl;
    else
        std::cout << "Done." << std::endl;
    PEGraph.close();
}

// dominant annotation for each bin, dominant bin for each annotation
// edge_graph.txt (with Size bin0 -> # !bin0, annotation)
// PE_graph.txt (PE >=5, with Size !bin0 -> #!bin0 & !equal, Deadend, annotation)
// Dominant line and #, not-Dominant line and #
// If bins are true, # all PEs, # correct PEs and # wrong PEs
// According to annotation, wrong lines of edge links, # wrong (with line)/right/ambigous PE links, and # in future might be corrected deadends and non-deadends (with line)
void Graph::breadth_first_search_ngs(std::string num) {
    int pe_num = atoi(num.c_str());
    // dominant annotation for each bin, dominant bin for each annotation
    std::map<int, std::unordered_map<std::string, int>> bin2anno;
    std::map<std::string, std::unordered_map<int, int>> anno2bin;
    std::map<int, int> bin2max;
    std::map<std::string, int> anno2max;
    std::map<int, std::string> bin2dominant;
    std::map<std::string, int> anno2dominant;
    for (auto i : this->paths) {
        for (auto a : i.second.annotation) {
            bin2anno[i.second.bin][a] += 1;
            if (bin2anno[i.second.bin][a] > bin2max[i.second.bin]) {
                bin2dominant[i.second.bin] = a;
                bin2max[i.second.bin] = bin2anno[i.second.bin][a];
            }
            if (i.second.bin != 0) {
                anno2bin[a][i.second.bin] += 1;
                if (anno2bin[a][i.second.bin] > anno2max[a]) {
                    anno2dominant[a] = i.second.bin;
                    anno2max[a] = anno2bin[a][i.second.bin];
                }
            }
        }
    }
    std::cout << "Dominant annotation for each bin:" << std::endl;
    for (auto i : bin2dominant) if (i.first) std::cout << i.first << "\t" << i.second << std::endl;
    std::cout << "Dominant bin for each annotation:" << std::endl;
    for (auto i : anno2dominant) std::cout << i.first << "\t" << i.second << std::endl;

    int cnt_line = 0, not_domi = 0, domi = 0;
    std::cout << "Start traversing edge graph ..." << std::endl;
    int cnt = 0;
    std::fstream edgeGraph("edge_graph.txt", std::fstream::out);
    for (auto i : this->paths) {
        cnt_line += 1;
        edgeGraph << i.first << " bin" << i.second.bin << " annotation[ ";
        bool flag_domi = false;
        for (auto a : i.second.annotation) {
            edgeGraph << a << ' ';
            if (i.second.bin && !a.compare(bin2dominant[i.second.bin])) flag_domi = true;
        }
        if (!flag_domi && i.second.annotation.size() && i.second.bin) {
            not_domi += 1;
            std::cout << "Not Dominant Annotation: line " << cnt_line << std::endl;
        }
        else if (flag_domi) {
            domi += 1;
        }
        edgeGraph << "] depth" << i.second.depth << " length" << i.second.length << ":\n";
        std::unordered_set<int> tmpBin;
        for (auto j : i.second.linksByEdge) {
            int binNum = this->paths[j].bin;
            if (i.second.bin == 0 && binNum != 0) tmpBin.insert(binNum);
            edgeGraph << '\t' << j << " bin" << binNum << " annotation[ ";
            for (auto a : paths[j].annotation) edgeGraph << a << ' ';
            edgeGraph << "] depth" << this->paths[j].depth << " length" << this->paths[j].length << " overlapPE" << i.second.PEOverlap[j] << ";\n";
            cnt_line += 1;
            bool flag = false;
            for (auto a1 : i.second.annotation) {
                for (auto a2 : this->paths[j].annotation)
                    if (!a1.compare(a2)) flag = true;
            }
            if (i.second.annotation.size() && paths[j].annotation.size() && !flag) {
                std::cout << "Edge is wrong: line " << cnt_line << std::endl;
            }
        }
        if (i.second.bin == 0 && tmpBin.size() >= 2) {
            cnt_line += 1;
            edgeGraph << "\tSize " << tmpBin.size() << ";\n";
        }
        if (++cnt % 10000 == 0)
            std::cout << "Searched edge neighbours of " << cnt << " paths ..." << std::endl;
    }
    if (cnt % 10000)
        std::cout << "Searched edge neighbours of " << cnt << " paths ..." << std::endl;
    std::cout << "Dominant annotation: " << domi << std::endl;
    std::cout << "Not dominant annotation: " << not_domi << std::endl;
    std::cout << "Done." << std::endl;
    edgeGraph.close();
    std::cout << "Start traversing PE graph ..." << std::endl;
    cnt = 0;
    cnt_line = 0;
    not_domi = 0;
    domi = 0;
    int cnt_pe_ok = 0, cnt_pe_wa = 0, cnt_pe_all = 0;
    int cnt_wrongPE = 0, cnt_rightPE = 0, cnt_Unknown = 0, cnt_deadendPE = 0, cnt_deadendPE_co = 0, cnt_nodeadPE = 0, cnt_nodeadPE_co = 0;
    std::fstream PEGraph("PE_graph.txt", std::fstream::out);
    std::vector<int> correct;
    std::vector<int> wrong;
    for (auto i : this->paths) {
        cnt_line += 1;
        std::unordered_set<std::string> tmpBin;
        PEGraph << i.first << " bin" << i.second.bin << " annotation[ ";
        bool flag_domi = false;
        for (auto a : i.second.annotation) {
            PEGraph << a << ' ';
            if (i.second.bin && !a.compare(bin2dominant[i.second.bin])) flag_domi = true;
        }
        if (!flag_domi && i.second.annotation.size() && i.second.bin) {
            not_domi += 1;
            std::cout << "Not Dominant Annotation: line " << cnt_line << std::endl;
        }
        else if (flag_domi) {
            domi += 1;
        }
        PEGraph << "] depth" << i.second.depth << " length" << i.second.length  << ":\n";
        for (auto j : i.second.PEOverlap) {
            if (j.second < pe_num) continue;
            cnt_line += 1;
            int neighbourNum = i.second.linksByEdge.size();
            int binNum = this->paths[j.first].bin;
            if (i.second.bin != 0 && i.second.bin != binNum && binNum != 0 && j.second != 0) {
                tmpBin.insert(j.first);
                cnt_pe_wa += j.second;
                wrong.push_back(j.second);
            }
            if (i.second.bin != 0 && i.second.bin == binNum && j.second != 0) {
                cnt_pe_ok += j.second;
                correct.push_back(j.second);
            }
            cnt_pe_all += j.second;
            bool flag = false;
            for (auto a1 : i.second.annotation) {
                for (auto a2 : this->paths[j.first].annotation)
                    if (!a1.compare(a2)) flag = true;
            }
            if (!i.second.annotation.size() || !paths[j.first].annotation.size()) cnt_Unknown += 1;
            else if (!flag) {
                cnt_wrongPE += 1;
                std::cout << "PE is wrong: line " << cnt_line << std::endl;
            }
            else {
                cnt_rightPE += 1;
                if (!neighbourNum || !this->paths[j.first].linksByEdge.size()) {
                    cnt_deadendPE += 1;
                    if (i.second.bin != binNum) {
                        cnt_deadendPE_co += 1;
                        std::cout << "Error correction of deadend: line " << cnt_line << std::endl;
                    }
                }
                else {
                    cnt_nodeadPE += 1;
                    if (i.second.bin != binNum) {
                        cnt_nodeadPE_co += 1;
                        std::cout << "Error correction of non-deadend: line " << cnt_line << std::endl;
                    }
                }
            }
            PEGraph << '\t' << j.first << " bin" << binNum << " annotation[ ";
            bool isconnectdomi = false;
            for (auto a : paths[j.first].annotation) {
                PEGraph << a << ' ';
                if (binNum == i.second.bin && !a.compare(bin2dominant[i.second.bin])) isconnectdomi = true;
            }
            if (!flag_domi && i.second.annotation.size() && i.second.bin && isconnectdomi) PEGraph << "ToDomin ";
            PEGraph << "] depth" << this->paths[j.first].depth << " length" << this->paths[j.first].length  << " overlapPE" << j.second;
            int flag1 = 0, flag2 = 0;
            if (neighbourNum == 0) ++flag1;
            if (this->paths[j.first].linksByEdge.size() == 0) ++flag2;
            if (flag1 && !flag2)
                PEGraph << " Deadend_left;\n";
            else if (!flag1 && flag2)
                PEGraph << " Deadend_right;\n";
            else if (flag1 && flag2)
                PEGraph << " Deadend_both;\n";
            else
                PEGraph << ";\n";
        }
        if (i.second.bin != 0 && tmpBin.size() >= 1) {
            PEGraph << "\tSize " << tmpBin.size() << ";\n";
            cnt_line += 1;
        }
        if (++cnt % 10000 == 0)
            std::cout << "Searched PE neighbours of " << cnt << " paths ..." << std::endl;
    }
    if (cnt % 10000)
        std::cout << "Searched PE neighbours of " << cnt << " paths ..." << std::endl;
    std::cout << "Dominant annotation: " << domi << std::endl;
    std::cout << "Not dominant annotation: " << not_domi << std::endl;
    std::cout << "All PE: " << cnt_pe_all << ", correct PE: " << cnt_pe_ok << ", wrong PE: " << cnt_pe_wa << std::endl;
    std::cout << "Wrong PE edges: " << cnt_wrongPE << ", Right PE edges: " << cnt_rightPE << ", Ambigous PE edges: " << cnt_Unknown << ", Error correction (deadend): " << cnt_deadendPE_co << '/' << cnt_deadendPE << ", Error correction (not-deadend): " << cnt_nodeadPE_co << '/' << cnt_nodeadPE << std::endl;
    std::cout << "Done." << std::endl;
    PEGraph.close();
}

// PE_graph.txt (with Size !bin0 -> #!bin0 & !equal, Deadend, annotation)
// If bins are true, # all PEs, # correct PEs and # wrong PEs
// According to annotation, # wrong / right / ambigous PE links, and # in future might be corrected deadendsand non - deadends(with line)
void Graph::breadth_first_search_metacarvel() {
    std::cout << "Start traversing PE graph ..." << std::endl;
    int cnt = 0, cnt_pe_ok = 0, cnt_pe_wa = 0, cnt_pe_all = 0, cnt_line = 0;
    int cnt_wrongPE = 0, cnt_rightPE = 0, cnt_Unknown = 0, cnt_deadendPE = 0, cnt_deadendPE_co = 0, cnt_nodeadPE = 0, cnt_nodeadPE_co = 0;
    std::fstream PEGraph("PE_graph.txt", std::fstream::out);
    std::vector<int> correct;
    std::vector<int> wrong;
    for (auto i : this->paths) {
        std::unordered_set<std::string> tmpBin;
        PEGraph << i.first << " bin" << i.second.bin << " annotation[ ";
        for (auto a : i.second.annotation) PEGraph << a << ' ';
        PEGraph << "]:\n";
        cnt_line += 1;
        for (auto j : i.second.PEOverlap) {
            cnt_line += 1;
            int neighbourNum = i.second.linksByEdge.size();
            int binNum = this->paths[j.first].bin;
            if (i.second.bin != 0 && i.second.bin != binNum && binNum != 0 && j.second != 0) {
                tmpBin.insert(j.first);
                cnt_pe_wa += j.second;
                wrong.push_back(j.second);
            }
            if (i.second.bin != 0 && i.second.bin == binNum && j.second != 0) {
                cnt_pe_ok += j.second;
                correct.push_back(j.second);
            }
            cnt_pe_all += j.second;
            bool flag = false;
            for (auto a1 : i.second.annotation) {
                for (auto a2 : this->paths[j.first].annotation)
                    if (!a1.compare(a2)) flag = true;
            }
            if (!i.second.annotation.size() || !paths[j.first].annotation.size()) cnt_Unknown += 1;
            else if (!flag) cnt_wrongPE += 1;
            else {
                cnt_rightPE += 1;
                if (!neighbourNum || !this->paths[j.first].linksByEdge.size()) {
                    cnt_deadendPE += 1;
                    if (i.second.bin != binNum) {
                        cnt_deadendPE_co += 1;
                        std::cout << "Error correction of deadend: line " << cnt_line << std::endl;
                    }
                }
                else {
                    cnt_nodeadPE += 1;
                    if (i.second.bin != binNum) {
                        cnt_nodeadPE_co += 1;
                        std::cout << "Error correction of non-deadend: line " << cnt_line << std::endl;
                    }
                }
            }
            PEGraph << '\t' << j.first << " bin" << binNum << " annotation[ ";
            for (auto a : paths[j.first].annotation) PEGraph << a << ' ';
            PEGraph << "]";
            int flag1 = 0, flag2 = 0;
            if (neighbourNum == 0) ++flag1;
            if (this->paths[j.first].linksByEdge.size() == 0) ++flag2;
            if (flag1 && !flag2)
                PEGraph << " Deadend_left;\n";
            else if (!flag1 && flag2)
                PEGraph << " Deadend_right;\n";
            else if (flag1 && flag2)
                PEGraph << " Deadend_both;\n";
            else
                PEGraph << ";\n";
        }
        if (i.second.bin != 0 && tmpBin.size() >= 1) {
            PEGraph << "\tSize " << tmpBin.size() << ";\n";
            cnt_line += 1;
        }
        if (++cnt % 10000 == 0)
            std::cout << "Searched PE neighbours of " << cnt << " paths ..." << std::endl;
    }
    if (cnt % 10000)
        std::cout << "Searched PE neighbours of " << cnt << " paths ..." << std::endl;
    std::cout << "All PE: " << cnt_pe_all << ", correct PE: " << cnt_pe_ok << ", wrong PE: " << cnt_pe_wa << std::endl;
    std::cout << "Wrong PE edges: " << cnt_wrongPE << ", Right PE edges: " << cnt_rightPE << ", Ambigous PE edges: " << cnt_Unknown << ", Error correction (deadend): " << cnt_deadendPE_co << '/' << cnt_deadendPE << ", Error correction (not-deadend): " << cnt_nodeadPE_co << '/' << cnt_nodeadPE << std::endl;
    std::cout << "Done." << std::endl;
    PEGraph.close();
} // Graph breadth first search

// output preceision, recall, f-1, ARI
void::Graph::bin_evaluate() {
    int all_binned = 0, all_annoed = 0, all_binned_annoed = 0;
    std::map<int, std::unordered_map<std::string, int>> bin2anno;
    std::map<int, int> bin2max;
    std::map<int, std::string> bin2dominant;
    std::map<std::string, std::unordered_map<int, int>> anno2bin;
    std::map<std::string, int> anno2max;
    std::map<std::string, int> anno2dominant;

    for (auto &i : this->paths) {
        // if (i.second.annotation.size() > 1) i.second.annotation.clear();
        if (i.second.bin) all_binned += 1;
        if (i.second.annotation.size()) all_annoed += 1;
        if (i.second.bin && i.second.annotation.size()) {
            all_binned_annoed += 1;
            if (i.second.annotation.size() == 1)
                std::cout << i.second.bin << '\t' << i.second.annotation[0] << std::endl;
            for (auto a : i.second.annotation) {
                bin2anno[i.second.bin][a] += 1;
                if (bin2anno[i.second.bin][a] > bin2max[i.second.bin]) {
                    bin2dominant[i.second.bin] = a;
                    bin2max[i.second.bin] = bin2anno[i.second.bin][a];
                }
                anno2bin[a][i.second.bin] += 1;
                if (anno2bin[a][i.second.bin] > anno2max[a]) {
                    anno2dominant[a] = i.second.bin;
                    anno2max[a] = anno2bin[a][i.second.bin];
                }
            }
        }
    }
    // std::cout << "Dominant annotation for each bin:" << std::endl;
    // for (auto i : bin2dominant) std::cout << i.first << "\t" << i.second << std::endl;
    // std::cout << "Dominant bin for each annotation:" << std::endl;
    // for (auto i : anno2dominant) std::cout << i.first << "\t" << i.second << std::endl;

    int prec_up = 0, rec_up = 0;
    for (auto i : this->paths) {
        if (i.second.bin == 0) continue;
        bool flag_prec = false, flag_rec = false;
        for (auto a : i.second.annotation) {
            if (!a.compare(bin2dominant[i.second.bin])) flag_prec = true;
            if (i.second.bin == anno2dominant[a]) flag_rec = true;
        }
        if (flag_prec) prec_up += 1;
        if (flag_rec) rec_up += 1;
    }

    double precision = 1.0 * prec_up / all_binned_annoed;
    double recall = 1.0 * rec_up / all_annoed;
    double f1 = 2.0 * precision * recall / (precision + recall);

    int t0 = 0, t1 = 0, t2 = 0, N = 0;
    for (auto i : bin2anno) {
        int sum_t1 = 0;
        for (auto j : i.second) {
            t0 = t0 + j.second * (j.second - 1) / 2;
            sum_t1 += j.second;
        }
        t1 = t1 + sum_t1 * (sum_t1 - 1) / 2;
        N += sum_t1;
    }
    for (auto i : anno2bin) {
        int sum_t2 = 0;
        for (auto j : i.second) sum_t2 += j.second;
        t2 = t2 + sum_t2 * (sum_t2 - 1) / 2;
    }
    double t3 = 2.0 * t1 * t2 / (N * (N - 1));
    double ARI = 0;
    if (t0 - t3 == 0 && 0.5 * (t1 + t2) - t3 == 0) ARI = 1;
    else ARI = 1.0 * (t0 - t3) / (0.5 * (t1 + t2) - t3);
    std::cout << "t1: " << t1 << " t2: " << t2 << " t3: " << t3 << std::endl;
    std::cout << "Precision: " << precision << '\t' << prec_up << '\t' << all_binned_annoed << std::endl;
    std::cout << "Recall: " << recall << '\t' << rec_up << '\t' << all_annoed << std::endl;
    std::cout << "F1: " << f1 << std::endl;
    std::cout << "ARI: " << ARI << std::endl;
}