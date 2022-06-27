#include <iostream>
#include <array>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <string>

#include <boost/serialization/unordered_map.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#include <signal.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

const int kNumOfTexts = 2000;
const char* fname = "serialization.txt";
// const char* fname = "/mnt/cryptoLab/serialization.txt";

std::string timeNow()
{
    const boost::posix_time::ptime now = 
        boost::posix_time::microsec_clock::local_time();

    const boost::posix_time::time_duration td = now.time_of_day();

    const long hours        = td.hours();
    const long minutes      = td.minutes();
    const long seconds      = td.seconds();
    const long milliseconds = td.total_milliseconds() -
                              ((hours * 3600 + minutes * 60 + seconds) * 1000);
    char buf[40];
    sprintf(buf, "%02ld:%02ld:%02ld.%03ld", 
        hours, minutes, seconds, milliseconds);

    return buf;
}

using DpList = std::unordered_map<uint16_t, int>;
using DpTable = std::unordered_map<uint16_t, DpList>;
DpTable precalc;

void findStats() {
	std::cout << __FUNCTION__ << ": " << timeNow() << '\n';
	int min = 0xffff;
	int max = 0;
	int sum = 0;
	int count = 0;
	double avg = 0;
	for (const auto& [a, aList]:precalc) {
		for (const auto& [b, freq]: aList) {
			if (freq > max) max = freq;
			if (freq < min) min = freq;
			sum += freq;
		}
		count += aList.size();
	}
	avg = (double)sum/count;
	std::cout << "min " << min << ", max " << max
		<< ", avg " << avg << '\n';
	std::cout << __FUNCTION__ << ": " << timeNow() << '\n';
}

void readFromFile() {
	std::cout << __FUNCTION__ << ": " << timeNow() << '\n';
	std::cout << "Reading...\n";
	std::ifstream ifs(fname);
	boost::archive::text_iarchive ia(ifs);
	ia >> precalc;
	std::cout << "Reading complete\n";
	// findStats();
	std::cout << __FUNCTION__ << ": " << timeNow() << '\n';
}

void saveToFile() {
	std::cout << __FUNCTION__ << ": " << timeNow() << '\n';
	std::cout << "Saving...\n";
	std::ofstream ofs{fname};
	boost::archive::text_oarchive oa(ofs);
	oa << precalc;
	std::cout << "Saving complete\n";
	findStats();
	std::cout << __FUNCTION__ << ": " << timeNow() << '\n';
}

void addSignalHandler() {
	struct sigaction sigIntHandler;

	sigIntHandler.sa_handler = [](int s) {
		saveToFile();
		exit(1);
	};
	sigemptyset(&sigIntHandler.sa_mask);
	sigIntHandler.sa_flags = 0;

	sigaction(SIGINT, &sigIntHandler, NULL);
}

void outputBytes(const auto& vec) {
	auto toAscii = [](uint8_t c) -> char {
		return c < 10 ? c + 48 : c + 55;
	};
	for (auto el : vec) {
		std::cout << toAscii(el>>4) << toAscii(el&0xf) << ' ';
	}
	std::cout << '\n';
}

uint16_t heys_round(const uint16_t x, const uint16_t k) {
	static uint8_t s[16] =  {4, 0xB, 1, 0xF, 9, 2, 0xE, 0xC, 6, 0xA, 8, 7, 3, 5, 0, 0xD};

	uint16_t y = x^k;

	uint8_t z[4] = {
		s[y      &0xf],
		s[(y>>4) &0xf],
		s[(y>>8) &0xf],
		s[(y>>12)],
	};

	uint16_t res = 0;
	for(int j = 0; j < 4; ++j) {
		for (int i = 0; i < 4; ++i) {
			res |= ((z[j] >> i) & 0x1) << (j + i*4);
		}
	}

	return res;
}

std::vector<uint8_t> heys(const std::vector<uint8_t>& x, const std::array<uint8_t, 14>& k) {
	if (x.size() < 2) {
		std::cout << "Size of input should be at least 2";
		exit(1);
	}

	int size_oddified = x.size() - x.size()%2;
	std::vector<uint8_t> res(size_oddified, 0);

	for (int i = 0; i < size_oddified; i += 2) {
		uint16_t tmp = (x[i]&0xff) | (x[i+1]<<8);
		for (int j = 0; j < 6; ++j) {
			// outputBytes(heys_round(tmp, {k[j*2], k[j*2+1]}));
			uint16_t roundKey = (k[j*2]&0xff) | (k[j*2+1]<<8);
			tmp = heys_round(tmp, roundKey);
		}
		res[i] = tmp ^ k[12];
		res[i+1] = (tmp>>8) ^ k[13];
	}

	return res;
}

uint16_t invL(uint16_t val) {

}

uint16_t invS(uint16_t val) {

}

// todo: precalc: only non-trivial differentials! is only (0,0) trivial?
void performPrecalculation(const uint16_t k) {
	std::cout << __FUNCTION__ << ": " << timeNow() << '\n';
	addSignalHandler();
	// 65534
	const double incr = 1/kNumOfTexts;
	for (uint16_t x1{0}; x1<=kNumOfTexts; ++x1) {
		for (uint16_t a{0}; a<=65534; ++a) {
			uint16_t x2 = x1 + a;
			uint16_t tmp1 = heys_round(x1, k);
			uint16_t tmp2 = heys_round(x2, k);
			uint16_t b = tmp2 - tmp1;

			if (a==0 && b == 0) continue;

			if(auto it = precalc.find(a); it !=precalc.end()) {
				auto& listPerA = it->second;
				if (auto it2 = listPerA.find(b); it2 != listPerA.end()) {
					auto& freq = it2->second;
					freq += 1;
					// freq += incr;
				} else {
					listPerA.emplace(b, 1);
				}
			} else {
				precalc.emplace(a, DpList{{b, 1}});
			}
		}
		
		if (x1%1000 == 0) std::cout <<  '\n' << x1 << '\n';
	}

	// won't reach here if too many iterations
	saveToFile();
}


DpList branchAndBound(uint16_t a) {
	DpList gPrev{{a, 1.0}};
	for (int t=1; t<=5; ++t) {
		DpList gCurr;
		for (const auto& prevState: gPrev) {
			for (const auto& nextState: precalc[prevState.first]) {
				// auto& gama = precalc.first;

				if (auto it = gCurr.find(nextState.first); it!=gCurr.end()) {
					it->second += 
				}
			}
		}
	}
}




int main(int argc, char* argv[]) {
	

	// std::ifstream filePt{argv[1], std::ios::binary | std::ios::ate};
	// int filePtSize = filePt.tellg();
	// std::vector<uint8_t> pt(filePtSize, 0);
	// filePt.seekg(0, filePt.beg);
	// filePt.read((char*)pt.data(), filePtSize);

	// std::ifstream fileKey{argv[3], std::ios::binary | std::ios::ate};
	// int fileKeySize = fileKey.tellg();
	// if (fileKeySize > 14) {
	// 	std::cout << "Key should be 14 byte long(if shorter, then padded with 0s)";
	// 	exit(1);
	// }
	// std::array<uint8_t, 14> k{0};
	// fileKey.seekg(0, fileKey.beg);
	// fileKey.read((char*)k.data(), fileKeySize);

	// auto res = heys(pt, k);
	// // outputBytes(pt);
	// outputBytes(res);

	// std::ofstream fileCt{argv[2], std::ios::binary};
	// fileCt.write((char*)res.data(), res.size());
	
	// performPrecalculation(1234U);

	readFromFile();
	// findStats();
}