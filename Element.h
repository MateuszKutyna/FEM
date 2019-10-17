#pragma once
#include"Node.h"

struct Element
{
	std::array<Node, 4> element;
	int ID;
	static int licznik;
	Element() {
		licznik++;
		ID = licznik;
	}
	void assignNodes(const std::vector<Node>&, const std::array<int, 4>&);
};