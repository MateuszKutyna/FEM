#include"pch.h"

void FEM::fillNodesPosition() {
	double dx, dy;
	dx = data.W / (data.nW - 1);
	dy = data.H / (data.nH - 1);
	int j = 0, k = 0;
	for (size_t i = 0; i < data.nN; i++) {
		if (i%data.nH == 0 && i != 0) {
			j++;
			k = 0;
		}

		nodes.emplace_back(dx * j, dy * k);
		k++;
	}
}
void FEM::fillElementsWithNodes() {
	std::array<int, 4> tab;
	for (size_t i = 0, j = 1; i < data.nE; i++, j++) {
		
		if (i%(data.nH-1) == 0 && i != 0) {
			j++;
		}
		tab[0] = j;
		tab[1] = tab[0] + data.nH;
		tab[2] = tab[1] + 1;
		tab[3] = tab[0] + 1;
		elements[i].assignNodes(nodes, tab);
	}
}


void FEM::printElements() {
	fillNodesPosition();
	fillElementsWithNodes();
	for (size_t i = 0; i < data.nE; i++) {
		std::cout << elements[i].ID << ":";
		for (size_t j = 0; j < 4; j++) {
			std::cout <<" "<< elements[i].element[j].t <<"("<<elements[i].element[j].x<<"," << elements[i].element[j].y <<"):";
		}
		std::cout << std::endl;
	}
}