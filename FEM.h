#pragma once

struct FEM
{
	std::vector<Node> nodes;
	std::vector<Element> elements;
	GlobalData data;
	FEM(const GlobalData &data1) :data(data1) {
		nodes.reserve(data.nN);
		for (std::size_t i = 0; i < data.nN; ++i)
			elements.emplace_back();
	}

	void fillNodesPosition();
	void fillElementsWithNodes();
	void printElements();
};