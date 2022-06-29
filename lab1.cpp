#include <iostream>
#include <array>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <string>
#include <functional>
#include <set>

#include <boost/date_time/posix_time/posix_time.hpp>

#include <signal.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

const int kNumOfMaxs = 10;

const int kNumOfTexts = 10000;
const char* fname = "serialization.txt";
// don't need it
const double kQStar = 0.0001;
// works fast anyway
const std::array<double, 5> pStars = {0.0001, 0.0001, 0.00001, 0.00001, 0.00001};
uint8_t S[16] =  {4, 0xB, 1, 0xF, 9, 2, 0xE, 0xC, 6, 0xA, 8, 7, 3, 5, 0, 0xD};

// const int kNumOfTextsForBrute = 10000;
const int kNumOfTextsForBrute = 10000;
// the more, the better
const int kMaxKeyCandidatesPerDiff = 30;

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

using DpList = std::unordered_map<uint16_t, double>;
using DpTable = std::unordered_map<uint16_t, DpList>;
DpTable precalc;

auto findStats() {
	// std::cout << __FUNCTION__ << ": " << timeNow() << '\n';
	double min = 1.0;
	std::multimap<double, std::pair<uint16_t, uint16_t>> max{};
	double sum = 0;
	int count = 0;
	for (const auto& [a, aList]:precalc) {
		for (const auto [b, prob]: aList) {
			// if (prob > max) { max = prob; maxA=a; maxB=b; }
			if(max.size()>0 && prob > max.begin()->first) {
				if (max.size() >= kNumOfMaxs) {
					// *max.begin() = prob;
					max.erase(max.begin());
				} 
				max.emplace(prob, std::make_pair(a, b));
			} else if (max.size() == 0) {
				max.emplace(prob, std::make_pair(a, b));
			}
			if (prob < min) min = prob;
			sum += prob;
		}
		count += aList.size();
	}
	double avg = (double)sum/count;
	std::string maxs = "[";
	for (const auto [p, dif]:max) {maxs += "(" + std::to_string(dif.first) + "," + std::to_string(dif.second) + ")=" + std::to_string(p) + ", ";}
	maxs += "]";
	std::cout << "min " << min << ", max " << maxs
		<< ", avg " << avg << '\n';
	// std::cout << __FUNCTION__ << ": " << timeNow() << '\n';
	return max;
}

void findStatsList(uint16_t a, const DpList& l, int round) {
	// std::cout << __FUNCTION__ << ": " << timeNow() << '\n';
	double min = 1.0;
	std::multimap<double, uint16_t> max{};
	double sum = 0;
	for (const auto [b, prob]: l) {
		if(max.size()>0 && prob > max.begin()->first) {
			if (max.size() >= kNumOfMaxs) {
				// *max.begin() = prob;
				max.erase(max.begin());
			} 
			max.emplace(prob, b);
		} else if (max.size() == 0) {
			max.emplace(prob, b);
		}
		if (prob < min) min = prob;
		sum += prob;
	}
	double avg = (double)sum/l.size();
	std::string maxs = "[";
	for (const auto [p, b]:max) {maxs += "(" + std::to_string(a) + "," + std::to_string(b) + ")=" + std::to_string(p) + ", ";}
	maxs += "]";
	std::cout << "min " << min << ", max " << round << "-round diff:" << maxs
		<< ", avg " << avg << '\n';
	if (round == 5) {
		maxs = "";
		for (auto it=max.rbegin();it!=max.rend();++it) {maxs += "{" + std::to_string(a) + "," + std::to_string(it->second) + "}, ";}
		std::cout << "max = " << maxs << '\n';
	}
	// std::cout << __FUNCTION__ << ": " << timeNow() << '\n';
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
	uint8_t S[16] =  {4, 0xB, 1, 0xF, 9, 2, 0xE, 0xC, 6, 0xA, 8, 7, 3, 5, 0, 0xD};

	uint16_t y = x^k;

	uint8_t z[4] = {
		S[y      &0xf],
		S[(y>>4) &0xf],
		S[(y>>8) &0xf],
		S[(y>>12)],
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
	uint16_t res = 0;
	for(int j = 0; j < 4; ++j) {
		for (int i = 0; i < 4; ++i) {
			res |= ((val >> (i+j*4)) & 0x1) << (j + i*4);
		}
	}
	return res;
}

uint16_t invS(uint16_t val) {
	static const std::array<uint8_t,16> inverseS = [](){
		std::array<uint8_t,16> res{0};
		for (int i=0; i<16; ++i) {
			res[S[i]] = i;
		}
		return res;
	}();
	uint16_t res = inverseS[(val)&0xf] |
		(inverseS[(val>> 4)&0xf] <<  4) |
		(inverseS[(val>> 8)&0xf] <<  8) |
		(inverseS[(val>>12)&0xf] << 12);
	return res;
}

auto performPrecalculation(const uint16_t k) {
	std::cout << "Precalculation using key " << std::hex << k << std::dec<< '\n';
	const double incr = 1.0/kNumOfTexts;
	for (uint16_t x1{0}; x1<kNumOfTexts; ++x1) {
		for (int aExt{0}; aExt<=65535; ++aExt) {
			uint16_t a = (uint16_t)aExt;
			uint16_t x2 = x1 ^ a;
			uint16_t tmp1 = heys_round(x1, k);
			uint16_t tmp2 = heys_round(x2, k);
			uint16_t b = tmp2 ^ tmp1;

			if (a==0 && b == 0) continue;

			if(auto it = precalc.find(a); it !=precalc.end()) {
				auto& listPerA = it->second;
				listPerA[b] += incr;
			} else {
				precalc.emplace(a, DpList{{b, incr}});
			}
		}
		
		if (x1%1000 == 999) std::cout <<  '\n' << x1 << '\n';
	}

	findStats();

	int countRemoved = 0;
	for (auto aListIt = precalc.begin();aListIt!=precalc.end();++aListIt) {
		auto& listPerA = aListIt->second;
		for (auto bListIt = listPerA.begin();bListIt!=listPerA.end();) {
			if (bListIt->second <= kQStar) {
				++countRemoved;
				bListIt = listPerA.erase(bListIt);
			} else {
				++bListIt;
			}
		}
	}

	std::cout << "Removed (prob less then q*) diffs " << countRemoved << '\n';

	auto max = findStats();

	std::cout << __FUNCTION__ << ": " << timeNow() << "\n\n";
	return max;
}


DpList branchAndBound(uint16_t a) {
	DpList gCurr{{a, 1.0}};
	DpList gNext;
	for (int t=1; t<=5; ++t) {
		gNext.clear();
		// branch
		for (const auto& currState: gCurr) {
			for (const auto& nextState: precalc[currState.first]) {
				gNext[nextState.first] += currState.second * nextState.second;
			}
		}

		findStatsList(a, gNext, t);
		// if (t == 5) findStatsList(a, gNext, t);

		//bound
		for (auto it=gNext.begin();it!=gNext.end();) {
			if (it->second <= pStars[t-1]) {
				it = gNext.erase(it);
			} else {
				++it;
			}
		}
		gCurr = std::move(gNext);
	}
	std::cout << __FUNCTION__ << ": " << timeNow() << "\n\n";
	return gCurr;
}

// 1 text = 2 bytes
auto ctGenMyHeys(int numOfTexts, uint16_t a, const std::array<uint8_t, 14>& k) {
    std::vector<uint8_t> pt1(2 * numOfTexts, 0);
    std::vector<uint8_t> pt2(2 * numOfTexts, 0);
    for (int i{0}; i < numOfTexts; ++i) {
        pt1[2*i] = i&0xff;
        pt1[2*i+1] = (i>>8)&0xff;
        uint16_t xXORa = (uint16_t)i^a;
        pt2[2*i] = (xXORa)&0xff;
        pt2[2*i+1] = ((xXORa)>>8)&0xff;
    }

    std::vector<uint8_t> ct1 = heys(pt1, k);
    std::vector<uint8_t> ct2 = heys(pt2, k);

    std::cout << __FUNCTION__ << ": " << timeNow() << '\n';
    return std::make_tuple(std::move(ct1), std::move(ct2));
}

void ptWriteNotMyHeys(int numOfTexts, uint16_t a) {
    std::vector<uint8_t> pt1(2 * numOfTexts, 0);
    std::vector<uint8_t> pt2(2 * numOfTexts, 0);
    for (int i{0}; i < numOfTexts; ++i) {
        pt1[2*i] = i&0xff;
        pt1[2*i+1] = (i>>8)&0xff;
        uint16_t xXORa = (uint16_t)i^a;
        pt2[2*i] = (xXORa)&0xff;
        pt2[2*i+1] = ((xXORa)>>8)&0xff;
    }

    // they're gonna be heys.bin's input PT
    std::ofstream fileInputPt1{"largePt1.txt", std::ios::binary};
	fileInputPt1.write((char*)pt1.data(), pt1.size());

	std::ofstream fileInputPt2{"largePt2.txt", std::ios::binary};
	fileInputPt2.write((char*)pt2.data(), pt2.size());

    std::cout << __FUNCTION__ << ": " << timeNow() << '\n';
}

auto ctReadNotMyHeys(int numOfTexts) {
    std::vector<uint8_t> ct1(2 * numOfTexts, 0);
    std::vector<uint8_t> ct2(2 * numOfTexts, 0);
    
    // they're gonna be heys.bin's output CT
    std::ifstream fileOutputCt1{"largeCt1.txt", std::ios::binary};
	fileOutputCt1.read((char*)ct1.data(), ct1.size());

	std::ifstream fileOutputCt2{"largeCt2.txt", std::ios::binary};
	fileOutputCt2.read((char*)ct2.data(), ct2.size());

    std::cout << __FUNCTION__ << ": " << timeNow() << '\n';
    return std::make_tuple(std::move(ct1), std::move(ct2));
}

// 1 text = 2 bytes
void bruteK7(int numOfTexts, uint16_t a, uint16_t b, const std::vector<uint8_t>& ct1, const std::vector<uint8_t>& ct2) {
    std::unordered_map<uint16_t, int> perKeySuccessCount;
    for (uint16_t i{0}; i < numOfTexts; ++i) {
        for (int k7Ext{0}; k7Ext < 65536; ++k7Ext) {
            uint16_t k7 = (uint16_t)k7Ext;
            uint16_t x71 = ct1[2 * i] | ((uint16_t)ct1[2 * i + 1] << 8);
            uint16_t x72 = ct2[2 * i] | ((uint16_t)ct2[2 * i + 1] << 8);
            uint16_t x61 = x71 ^ k7;
            uint16_t x62 = x72 ^ k7;
            uint16_t x51XORk6 = invS(invL(x61));
            uint16_t x52XORk6 = invS(invL(x62));
            uint16_t potentialB = x51XORk6 ^ x52XORk6;
            if (potentialB == b) {
                perKeySuccessCount[k7] += 1;
            }
        }
    }

    std::multimap<int, uint16_t> max{};
    for (const auto [k, freq] : perKeySuccessCount) {
        if (max.size() > 0 && freq > max.begin()->first) {
            if (max.size() >= kMaxKeyCandidatesPerDiff) {
                max.erase(max.begin());
            }
            max.emplace(freq, k);
        } else if (max.size() == 0) {
            max.emplace(freq, k);
        }
    }

    std::cout << "Most probable key for (" << a << ',' << b << ") in (key, frequency) format = ";
    for (const auto [freq, k] : max) {
        std::cout << "(" << std::hex << k << std::dec<< "," << freq << "), ";
    }
    std::cout << '\n';
    std::cout << __FUNCTION__ << ": " << timeNow() << '\n';
}

int main(int argc, char* argv[]) {
	std::cout << __FUNCTION__ << ": " << timeNow() << '\n';

	// --------------------------------------- SEARCHING FOR MOST PROBABLE DIFFS --------------------------------------

	// auto maxAB = performPrecalculation(0x0000);
	// // // taking As with most probable difs in precalc
	// std::set<uint16_t> aUsed;
	// for (auto it=maxAB.rbegin();it!=maxAB.rend();++it) {
	// 	if (aUsed.find(it->second.first) != aUsed.end()) {continue;}
	// 	std::cout << "Performing branch and bound on a=" << it->second.first << '\n';
	// 	branchAndBound(it->second.first);
	// 	aUsed.insert(it->second.first);
	// }

	// --------------------------------------- TRYING TO BRUTEFORCE KEY SUPPLIED BY ME --------------------------------------

    // auto mostProbable5RoundDiffs = std::vector<std::pair<uint16_t, uint16_t>>{
    // 	// {40960,1092}, {40960,52428}, {40960,4}, {40960,8736}, {40960,34952}, {40960,17472},
    // 	{53248,34952}, {53248,17472}, {53248,4}, {53248,2056}, {53248,1092}, {53248,64}, {53248,8736}, {53248,16384}, {53248,17476}, {53248,52428},
    // 	// {53248,8736},
    // };
    // std::vector<uint8_t> ct1, ct2;
    // uint16_t currA = 0;
    // // std::array<uint8_t, 14> testKey{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0x11, 0x11};
    // // std::array<uint8_t, 14> testKey{0x15, 0x57, 0x57, 0x67, 0x57, 0x11, 0x57, 0x05, 0x54, 0x1, 0x2, 0x3, 0x81, 0x54};
    // std::array<uint8_t, 14> testKey{0x15, 1, 5, 3, 4, 100, 67, 6, 3, 7, 0x2, 0x3, 0x0f, 0x0f};
    // for (const auto [a, b] : mostProbable5RoundDiffs) {
    // 	if (currA != a) {
    //         std::tie(ct1, ct2) = ctGenMyHeys(kNumOfTextsForBrute, a, testKey);
    //         currA = a;
    // 	}
    //     bruteK7(kNumOfTextsForBrute, a, b, ct1, ct2);
    // }

 //    // --------------------------------------- TRYING TO BRUTEFORCE KEY EMBEDDED IN HEYS.BIN -> FAILS UNFORTENETLY -----------------------------
    
    ptWriteNotMyHeys(kNumOfTextsForBrute, 53248);
    // // ptWriteNotMyHeys(kNumOfTextsForBrute, 20480);

	// can use only the same a for all diff when bruteforcing embedded key
    auto mostProbable5RoundDiffs = std::vector<std::pair<uint16_t, uint16_t>>{
    	{53248,8736}, {53248,17476}, {53248,17472}, {53248,34952}, {53248,4}, {53248,2056}, {53248,1092}, {53248,64}, {53248,16384}, {53248,52428},
    	// {20480,34824}, {20480,32776}, {20480,4369}, {20480,8738}, {20480,4097}, {20480,4353}, {20480,257}, {20480,17}, {20480,514}, {20480,8706}
    };
    std::vector<uint8_t> ct1, ct2;
    uint16_t currA = 0;
    for (const auto [a, b] : mostProbable5RoundDiffs) {
    	if (currA != a) {
            std::tie(ct1, ct2) = ctReadNotMyHeys(kNumOfTextsForBrute);
            currA = a;
    	}
        bruteK7(kNumOfTextsForBrute, a, b, ct1, ct2);
    }
}