#pragma once

struct FEM
{
	std::vector<Node> nodes;
	std::vector<Element> elements;
	GlobalData data;
	FEM() = default;
	FEM(const GlobalData &data1) :data(data1) {
		nodes.reserve(data.nN);
		elements.resize(data.nE);
		//for (std::size_t i = 0; i < data.nN; ++i)
		//	elements.emplace_back();
		fillNodesPosition();
		fillElementsWithNodes();
	}

	void fillNodesPosition();
	void fillElementsWithNodes();
	void printElements();
};