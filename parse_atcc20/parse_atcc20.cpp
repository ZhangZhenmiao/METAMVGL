#include <unordered_map>
#include <unordered_set>
#include <string>
#include <htslib/sam.h>
#include <iostream>
void load_alignment(std::string filename, std::string mq) {
    int mappingq = atoi(mq.c_str());
    const int Read1 = 0x00000040, Read2 = 0x00000080;
    const int supplementaryAlignment = 0x00000800, secondaryAlignment = 0x00000100;
    samFile *bamFile = hts_open(filename.c_str(), "r");
    bam_hdr_t *bamHeader = sam_hdr_read(bamFile);
    bam1_t *bamAln = bam_init1();
    int cntLine = 0;
    int reads1Len = 127, reads2Len = 150;
    std::unordered_map<std::string, std::string> seq2ref({
        {"CP000521.1", "Acinetobacter_baumannii"}, {"CP000522.1", "Acinetobacter_baumannii"}, {"CP000523.1", "Acinetobacter_baumannii"}, {"DS264586.1", "Actinomyces_odontolyticus"}, {"DS264585.1", "Actinomyces_odontolyticus"}, {"AE017194.1", "Bacillus_cereus"}, {"AE017195.1", "Bacillus_cereus"}, {"CP000139.1", "Bacteroides_vulgatus"}, {"AP009256.1", "Bifidobacterium_adolescentis"}, {"CP006777.1", "Clostridium_beijerinckii"}, {"AE000513.1", "Deinococcus_radiodurans"}, {"AE001825.1", "Deinococcus_radiodurans"}, {"AE001827.1", "Deinococcus_radiodurans"}, {"AE001826.1", "Deinococcus_radiodurans"}, {"CP002621.1", "Enterococcus_faecalis"}, {"U00096.3", "Escherichia_coli"}, {"AE000511.1", "Helicobacter_pylori"}, {"CP000413.1", "Lactobacillus_gasseri"}, {"AE002098.2", "Neisseria_meningitidis"}, {"AP009380.1", "Porphyromonas_gingivalis"}, {"CP003084.1", "Propionibacterium_acnes"}, {"PDLX01000010.1", "Pseudomonas_aeruginosa"}, {"PDLX01000011.1", "Pseudomonas_aeruginosa"}, {"PDLX01000012.1", "Pseudomonas_aeruginosa"}, {"PDLX01000013.1", "Pseudomonas_aeruginosa"}, {"PDLX01000014.1", "Pseudomonas_aeruginosa"}, {"PDLX01000015.1", "Pseudomonas_aeruginosa"}, {"PDLX01000016.1", "Pseudomonas_aeruginosa"}, {"PDLX01000017.1", "Pseudomonas_aeruginosa"}, {"PDLX01000018.1", "Pseudomonas_aeruginosa"}, {"PDLX01000019.1", "Pseudomonas_aeruginosa"}, {"PDLX01000001.1", "Pseudomonas_aeruginosa"}, {"PDLX01000020.1", "Pseudomonas_aeruginosa"}, {"PDLX01000021.1", "Pseudomonas_aeruginosa"}, {"PDLX01000022.1", "Pseudomonas_aeruginosa"}, {"PDLX01000023.1", "Pseudomonas_aeruginosa"}, {"PDLX01000024.1", "Pseudomonas_aeruginosa"}, {"PDLX01000025.1", "Pseudomonas_aeruginosa"}, {"PDLX01000026.1", "Pseudomonas_aeruginosa"}, {"PDLX01000027.1", "Pseudomonas_aeruginosa"}, {"PDLX01000028.1", "Pseudomonas_aeruginosa"}, {"PDLX01000029.1", "Pseudomonas_aeruginosa"}, {"PDLX01000002.1", "Pseudomonas_aeruginosa"}, {"PDLX01000030.1", "Pseudomonas_aeruginosa"}, {"PDLX01000031.1", "Pseudomonas_aeruginosa"}, {"PDLX01000032.1", "Pseudomonas_aeruginosa"}, {"PDLX01000033.1", "Pseudomonas_aeruginosa"}, {"PDLX01000034.1", "Pseudomonas_aeruginosa"}, {"PDLX01000035.1", "Pseudomonas_aeruginosa"}, {"PDLX01000036.1", "Pseudomonas_aeruginosa"}, {"PDLX01000037.1", "Pseudomonas_aeruginosa"}, {"PDLX01000038.1", "Pseudomonas_aeruginosa"}, {"PDLX01000039.1", "Pseudomonas_aeruginosa"}, {"PDLX01000003.1", "Pseudomonas_aeruginosa"}, {"PDLX01000040.1", "Pseudomonas_aeruginosa"}, {"PDLX01000041.1", "Pseudomonas_aeruginosa"}, {"PDLX01000042.1", "Pseudomonas_aeruginosa"}, {"PDLX01000043.1", "Pseudomonas_aeruginosa"}, {"PDLX01000044.1", "Pseudomonas_aeruginosa"}, {"PDLX01000045.1", "Pseudomonas_aeruginosa"}, {"PDLX01000046.1", "Pseudomonas_aeruginosa"}, {"PDLX01000047.1", "Pseudomonas_aeruginosa"}, {"PDLX01000048.1", "Pseudomonas_aeruginosa"}, {"PDLX01000049.1", "Pseudomonas_aeruginosa"}, {"PDLX01000004.1", "Pseudomonas_aeruginosa"}, {"PDLX01000050.1", "Pseudomonas_aeruginosa"}, {"PDLX01000051.1", "Pseudomonas_aeruginosa"}, {"PDLX01000052.1", "Pseudomonas_aeruginosa"}, {"PDLX01000053.1", "Pseudomonas_aeruginosa"}, {"PDLX01000054.1", "Pseudomonas_aeruginosa"}, {"PDLX01000055.1", "Pseudomonas_aeruginosa"}, {"PDLX01000056.1", "Pseudomonas_aeruginosa"}, {"PDLX01000057.1", "Pseudomonas_aeruginosa"}, {"PDLX01000058.1", "Pseudomonas_aeruginosa"}, {"PDLX01000059.1", "Pseudomonas_aeruginosa"}, {"PDLX01000005.1", "Pseudomonas_aeruginosa"}, {"PDLX01000006.1", "Pseudomonas_aeruginosa"}, {"PDLX01000007.1", "Pseudomonas_aeruginosa"}, {"PDLX01000008.1", "Pseudomonas_aeruginosa"}, {"PDLX01000009.1", "Pseudomonas_aeruginosa"}, {"CP000577.1", "Rhodobacter_sphaeroides"}, {"CP000578.1", "Rhodobacter_sphaeroides"}, {"CP000579.1", "Rhodobacter_sphaeroides"}, {"CP000255.1", "Staphylococcus_aureus"}, {"CP000256.1", "Staphylococcus_aureus"}, {"CP000257.1", "Staphylococcus_aureus"}, {"CP000258.1", "Staphylococcus_aureus"}, {"AE015929.1", "Staphylococcus_epidermidis"}, {"AE015930.1", "Staphylococcus_epidermidis"}, {"AE015931.1", "Staphylococcus_epidermidis"}, {"AE015932.1", "Staphylococcus_epidermidis"}, {"AE015933.1", "Staphylococcus_epidermidis"}, {"AE015934.1", "Staphylococcus_epidermidis"}, {"AE015935.1", "Staphylococcus_epidermidis"}, {"NC_004116.1", "Streptococcus_agalactiae"}, {"AE014133.2", "Streptococcus_mutans"}
    });
    std::unordered_set<std::string> seqs;
    std::unordered_map<std::string, long long> ref2len;
    std::unordered_map<std::string, long long> ref2reads1;
    std::unordered_map<std::string, long long> ref2reads2;
    std::cout << "Start parsing alignments ..." << std::endl;
    while (sam_read1(bamFile, bamHeader, bamAln) >= 0) {
        try {
            int flag = bamAln->core.flag;
            bool isRead1 = flag & Read1, isRead2 = flag & Read2;
            bool isSupplementary = flag & supplementaryAlignment, isSecondary = flag & secondaryAlignment;
            std::string contigName;
            if (sam_hdr_tid2name(bamHeader, bamAln->core.tid) != NULL)
                contigName = sam_hdr_tid2name(bamHeader, bamAln->core.tid);
            if (contigName.size() == 0 || contigName.compare("*") == 0) continue;
            int contigLen = sam_hdr_tid2len(bamHeader, bamAln->core.tid);
            if (contigLen == 0) continue;
            if (seqs.find(contigName) == seqs.end()) {
                seqs.insert(contigName);
                ref2len[seq2ref[contigName]] += contigLen;
                std::cout << "Reference length " << seq2ref[contigName] << " updated to " << ref2len[seq2ref[contigName]] << " ..." << std::endl;
            }
            int mappingQuality = bamAln->core.qual;
            if (mappingQuality < mappingq) continue;
            if (isRead1) ++ref2reads1[seq2ref[contigName]];
            else if (isRead2) ++ref2reads2[seq2ref[contigName]];
            // logging
            if (++cntLine % 10000000 == 0) {
                std::cout << "Parsed " << cntLine << " lines ..." << std::endl;
                for (auto ref : ref2len) {
                    std::cout << ref.first << ": length " << ref.second << ", reads1 " << ref2reads1[ref.first] << ", reads2 " << ref2reads2[ref.first] << ", coverage " << (1.0 * reads1Len * ref2reads1[ref.first] + 1.0 * reads2Len * ref2reads2[ref.first]) / ref.second << std::endl;
                }
            }
        }
        catch (...) {
            std::cout << "An error occured: line " << cntLine << " ..." << std::endl;
        }
    }
    if (cntLine % 10000000)
        std::cout << "Parsed " << cntLine << " lines ...\nDone." << std::endl;
    else
        std::cout << "Done." << std::endl;
    for (auto ref : ref2len) {
        std::cout << ref.first << ": length " << ref.second << ", reads1 " << ref2reads1[ref.first] << ", reads2 " << ref2reads2[ref.first] << ", coverage " << (1.0 * reads1Len * ref2reads1[ref.first] + 1.0 * reads2Len * ref2reads2[ref.first]) / ref.second << std::endl;
    }
    bam_hdr_destroy(bamHeader);
    bam_destroy1(bamAln);
    sam_close(bamFile);
}
int main(int argc, char* argv[]) {
    load_alignment(argv[1], argv[2]);
    return 0;
}