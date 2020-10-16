#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <htslib/sam.h>

std::unordered_map<std::string, std::unordered_set<std::string>> map;

void load_alignment(std::string filename) {
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
		// if (isSupplementary || isSecondary) continue;
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
		if (blastIdentity < 0.95) continue;
		// if (mappingQuality < 20) continue;
		map[readName].insert(filename);
		std::cout << filename << '\t' << readName << '\t' << pos << '\t' << endPos << '\t' << mappingQuality << '\t' << blastIdentity << '\t' << alignmentColumns << '\t' << readLength << '\t' << node << '\t' << contigLength << std::endl;
	}
	bam_hdr_destroy(bamHeader);
	bam_destroy1(bamAln);
	sam_close(bamFile);
}

int main() {
	load_alignment("1.bam");
	load_alignment("2.bam");
	load_alignment("3.bam");
	load_alignment("4.bam");
	load_alignment("5.bam");
	load_alignment("6.bam");
	load_alignment("7.bam");
	load_alignment("8.bam");
	load_alignment("9.bam");
	load_alignment("10.bam");
	load_alignment("11.bam");
	load_alignment("12.bam");
	load_alignment("13.bam");
	load_alignment("14.bam");
	load_alignment("15.bam");
	load_alignment("16.bam");
	load_alignment("17.bam");
	load_alignment("18.bam");
	load_alignment("19.bam");
	load_alignment("20.bam");
	load_alignment("21.bam");
	load_alignment("22.bam");
	load_alignment("23.bam");
	load_alignment("24.bam");
	load_alignment("25.bam");
	load_alignment("26.bam");
	load_alignment("27.bam");
	load_alignment("28.bam");
	load_alignment("29.bam");
	load_alignment("30.bam");
	load_alignment("31.bam");
	load_alignment("32.bam");
	load_alignment("33.bam");
	load_alignment("34.bam");
	load_alignment("35.bam");
	load_alignment("36.bam");
	load_alignment("37.bam");
	load_alignment("38.bam");
	load_alignment("39.bam");
	load_alignment("40.bam");
	load_alignment("41.bam");
	load_alignment("42.bam");
	load_alignment("43.bam");
	load_alignment("44.bam");
	load_alignment("45.bam");
	load_alignment("46.bam");
	load_alignment("47.bam");
	load_alignment("48.bam");
	load_alignment("49.bam");
	load_alignment("50.bam");
	load_alignment("51.bam");
	load_alignment("52.bam");
	load_alignment("53.bam");
	load_alignment("54.bam");
	load_alignment("55.bam");
	load_alignment("56.bam");
	load_alignment("57.bam");
	load_alignment("58.bam");
	load_alignment("59.bam");
	load_alignment("60.bam");
	load_alignment("61.bam");
	load_alignment("62.bam");
	load_alignment("63.bam");
	load_alignment("64.bam");
	std::cout << "Start final output ..." << std::endl;
	for (auto i : map) {
		std::cout << i.first;
		for (auto j : i.second) {
			std::cout << '\t' << j;
		}
		std::cout << std::endl;
	}
	return 0;
}
