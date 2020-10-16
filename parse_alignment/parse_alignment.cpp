#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <htslib/sam.h>

std::unordered_map<std::string, std::unordered_set<std::string>> bc2reads_1, bc2reads_2;
std::unordered_set<std::string> reads1_1, reads2_1, reads1_2, reads2_2;

void load_alignment_1(std::string filename) {
	const int supplementaryAlignment = 0x00000800, secondaryAlignment = 0x00000100;
	const int Read1 = 0x00000040, Read2 = 0x00000080;
	samFile* bamFile = hts_open(filename.c_str(), "r");
	bam_hdr_t* bamHeader = sam_hdr_read(bamFile);
	bam1_t* bamAln = bam_init1();
	int cntLine = 0;
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
		int contigLength = bamHeader->target_len[bamAln->core.tid];
		int readLength = bamAln->core.l_qseq;
		int alignmentColumns = 0, NM = 0;
		uint32_t* cigarPointer = bam_get_cigar(bamAln);
		for (int i = 0; i < bamAln->core.n_cigar; ++i) {
			char cigarOperator = bam_cigar_opchr(cigarPointer[i]);
			if (cigarOperator == 'M' || cigarOperator == 'I' || cigarOperator == 'D')
				alignmentColumns += bam_cigar_oplen(cigarPointer[i]);
		}
		if (alignmentColumns == 0) continue;
		uint8_t* tmpNM = bam_aux_get(bamAln, "NM");
		if (tmpNM != NULL) NM = bam_aux2i(tmpNM);
		double blastIdentity = 1.0 * (alignmentColumns - NM) / alignmentColumns;
		std::string barcode;
		uint8_t *tmpBC = bam_aux_get(bamAln, "BX");
		if (tmpBC != NULL) barcode = bam_aux2Z(tmpBC);
		if (blastIdentity < 0.98) continue;
		if (mappingQuality < 20) continue;
		if(isRead1) reads1_1.insert(readName);
		else if(isRead2) reads2_1.insert(readName);
		else std::cout << "abnormal reads" << std::endl;
		if(!barcode.empty()) bc2reads_1[barcode].insert(readName);
		std::cout << filename << '\t' << readName << '\t' << barcode << '\t' << isRead1 << '\t' << isRead2 << '\t' << mappingQuality << '\t' << blastIdentity << '\t' << alignmentColumns << '\t' << pos << '\t' << endPos << '\t' << contigLength << std::endl;
	}
	bam_hdr_destroy(bamHeader);
	bam_destroy1(bamAln);
	sam_close(bamFile);
}

void load_alignment_2(std::string filename) {
        const int supplementaryAlignment = 0x00000800, secondaryAlignment = 0x00000100;
        const int Read1 = 0x00000040, Read2 = 0x00000080;
        samFile* bamFile = hts_open(filename.c_str(), "r");
        bam_hdr_t* bamHeader = sam_hdr_read(bamFile);
        bam1_t* bamAln = bam_init1();
        int cntLine = 0;
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
                int contigLength = bamHeader->target_len[bamAln->core.tid];
                int readLength = bamAln->core.l_qseq;
                int alignmentColumns = 0, NM = 0;
                uint32_t* cigarPointer = bam_get_cigar(bamAln);
                for (int i = 0; i < bamAln->core.n_cigar; ++i) {
                        char cigarOperator = bam_cigar_opchr(cigarPointer[i]);
                        if (cigarOperator == 'M' || cigarOperator == 'I' || cigarOperator == 'D')
                                alignmentColumns += bam_cigar_oplen(cigarPointer[i]);
                }
                if (alignmentColumns == 0) continue;
                uint8_t* tmpNM = bam_aux_get(bamAln, "NM");
                if (tmpNM != NULL) NM = bam_aux2i(tmpNM);
                double blastIdentity = 1.0 * (alignmentColumns - NM) / alignmentColumns;
                std::string barcode;
                uint8_t *tmpBC = bam_aux_get(bamAln, "BX");
                if (tmpBC != NULL) barcode = bam_aux2Z(tmpBC);
                if (blastIdentity < 0.98) continue;
                if (mappingQuality < 20) continue;
                if(isRead1) reads1_2.insert(readName);
                else if(isRead2) reads2_2.insert(readName);
                else std::cout << "abnormal reads" << std::endl;
                if(!barcode.empty()) bc2reads_2[barcode].insert(readName);
                std::cout << filename << '\t' << readName << '\t' << barcode << '\t' << isRead1 << '\t' << isRead2 << '\t' << mappingQuality << '\t' << blastIdentity << '\t' << alignmentColumns << '\t' << pos << '\t' << endPos << '\t' << contigLength << std::endl;
        }
        bam_hdr_destroy(bamHeader);
        bam_destroy1(bamAln);
        sam_close(bamFile);
}

int main(int args, char** argv) {
	load_alignment_1(argv[1]);
	load_alignment_2(argv[2]);
	std::cout << "Results:" << std::endl;
	std::cout << "Reads1 of " << argv[1] << ": " << reads1_1.size() << std::endl;
	for (auto i : reads1_1) std::cout << i << std::endl;
	std::cout << "Reads2 of " << argv[1] << ": " << reads2_1.size() << std::endl;
	for (auto i : reads2_1) std::cout << i << std::endl;
	std::cout << "Reads1 of " << argv[2] << ": " << reads1_2.size() << std::endl;
	for (auto i : reads1_2) std::cout << i << std::endl;
	std::cout << "Reads2 of " << argv[2] << ": " << reads2_2.size() << std::endl;
	for (auto i : reads2_2) std::cout << i << std::endl;
	std::cout << "Barcodes of " << argv[1] << ": " << bc2reads_1.size() << std::endl;
	for (auto i : bc2reads_1){
		if (i.second.size() < 5) continue;
		std::cout << i.first << ' ' << i.second.size();
		for (auto j : i.second) std::cout << '\t' << j;
		std::cout << std::endl;
	}
	std::cout << "Barcodes of " << argv[2] << ": " << bc2reads_2.size() << std::endl;
	for (auto i : bc2reads_2){
                if (i.second.size() < 5) continue;
                std::cout << i.first << ' ' << i.second.size();
                for (auto j : i.second) std::cout << '\t' << j;
                std::cout << std::endl;
        }
	std::cout << "reads1 overlap reads2:" << std::endl;
	int cnt = 0;
	for (auto i : reads1_1) { if (reads2_2.find(i) != reads2_2.end()) {std::cout << i << std::endl; cnt += 1;} }
	std::cout << cnt << std::endl;
	std::cout << "reads2 overlap reads1:" << std::endl;
	cnt = 0;
	for (auto i : reads2_1) { if (reads1_2.find(i) != reads1_2.end()) {std::cout << i << std::endl; cnt += 1;} }
	std::cout << cnt << std::endl;
	std::cout << "barcode overlap:" << std::endl;
	cnt = 0;
	for (auto i : bc2reads_1){
		if (i.second.size() >= 5 && bc2reads_2.find(i.first) != bc2reads_2.end() && bc2reads_2[i.first].size() >= 5){
			std::cout << i.first << std::endl;
			cnt += 1;
		}
	}
	std::cout << cnt << std::endl;
	return 0;
}
