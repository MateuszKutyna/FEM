#include "pch.h"

void GlobalData::readFromFile(std::string file_name) {
	std::fstream file;
	file.open(file_name, std::ios::in);
	if (!file.is_open()) {
		std::cerr << "There was a problem opening the input file!\n";
		exit(1);
	}
	std::string line;
	std::array<double, 13> tab;
	for (size_t i = 0; i <13; i++) {
		getline(file, line);
		tab[i] = stod(line);
	}
	H  = tab[0];
	W  = tab[1];
	nH = tab[2];
	nW = tab[3];
	nN = nH * nW;
	nE = (nH - 1)*(nW - 1);
	integralPoints = tab[4];
	k = tab[5];
	alfa=tab[6];
	cw = tab[7];
	ro = tab[8];
	t8 = tab[9];
	dt = tab[10];
	sim_time = tab[11];
	for (int i = 0; i < nN; i++) {
		t0.emplace_back(tab[12]);
	}

}