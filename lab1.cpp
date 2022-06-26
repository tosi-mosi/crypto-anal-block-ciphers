#include <iostream>
#include <array>
#include <vector>
#include <fstream>

void outputBytes(const auto& vec) {
	auto toAscii = [](uint8_t c) -> char {
		return c < 10 ? c + 48 : c + 55;
	};
	for (auto el : vec) {
		std::cout << toAscii(el>>4) << toAscii(el&0xf) << ' ';
	}
	std::cout << '\n';
}

// arr has most significan bytes on lowest index
uint16_t arrToUint16 (const auto& vec) {
	return vec[0]<<8 | vec[1]&&0xff;
}

std::array<uint8_t, 2> uint16ToArr (uint16_t val) {
	return {val>>8, val&&0xff};
}

void precalculation() {
	std::array<uint8_t, 2> k{0x10, 0x11};
	for (uint16_t a=0; a<=65535; ++a) {
		for (uint16_t b=0; b<=65535; ++b) {
			for (uint16_t x=0; x<=65535; ++x){
				auto tmp1 = heys_round(uint16ToArr(x), k);
				uint8_t tmp1u = (tmp1[0] << 8) + tmp1[1] & 0xff;
				if( 
					== heys_round());
			}
		}
	}
}

std::array<uint8_t, 2> heys_round(const std::array<uint8_t, 2>& x, const std::array<uint8_t, 2>& k) {
	static uint8_t s[16] =  {4, 0xB, 1, 0xF, 9, 2, 0xE, 0xC, 6, 0xA, 8, 7, 3, 5, 0, 0xD};

	uint8_t y[2] = {
		x[0] ^ k[0],
		x[1] ^ k[1],
	};

	uint8_t z[4] = {
		s[y[0]&0xf],
		s[y[0]>>4],
		s[y[1]&0xf],
		s[y[1]>>4],
	};

	std::array<uint8_t, 2> res{0};
	for(int j = 0; j < 4; ++j) {
		for (int i = 0; i < 4; ++i) {
			res[i/2] |= ((z[j] >> i) & 0x1) << (i%2==0 ? j : j+4);
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
		std::array<uint8_t, 2> tmp{x[i], x[i+1]};
		for (int j = 0; j < 6; ++j) {
			// std::cout << "r" << j+1 << ": ";
			// outputBytes(heys_round(tmp, {k[j*2], k[j*2+1]}));
			tmp = heys_round(tmp, {k[j*2], k[j*2+1]});
		}
		res[i] = tmp[0] ^ k[12];
		res[i+1] = tmp[1] ^ k[13];
	}

	return res;
}

// a is fixed
using DpList = std::vector<std::pair<uint16_t, double>>;
// a is not fixed
// no need to use std::pair, since vector's index will be a, cuz we iterate
//  through all possible a
using DpTable = std::vector<DpList>;

std::vector<std::pair<uint16_t, double>> branchAndBound(
	uint16_t a, roundDpTable) {

}

// precalc: only non-trivial differentials!

int main() {
	std::ifstream filePt{"pt.txt", std::ios::binary | std::ios::ate};
	int filePtSize = filePt.tellg();
	std::vector<uint8_t> pt(filePtSize, 0);
	filePt.seekg(0, filePt.beg);
	filePt.read((char*)pt.data(), filePtSize);

	std::ifstream fileKey{"k.txt", std::ios::binary | std::ios::ate};
	int fileKeySize = fileKey.tellg();
	if (fileKeySize > 14) {
		std::cout << "Key should be 14 byte long(if shorter, then padded with 0s)";
		exit(1);
	}
	std::array<uint8_t, 14> k{0};
	fileKey.seekg(0, fileKey.beg);
	fileKey.read((char*)k.data(), fileKeySize);

	auto res = heys(pt, k);
	// outputBytes(pt);
	outputBytes(res);

}