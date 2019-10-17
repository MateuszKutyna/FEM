#include"pch.h"

int Element::licznik = 0;

void Element::assignNodes(const std::vector<Node>& ele, const std::array<int, 4>& tab_id) {
	for (size_t i = 0; i < 4; i++) {
		element[i] = ele[tab_id[i] - 1];
	}
}