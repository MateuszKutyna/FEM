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
		fillNodesPosition();
		fillElementsWithNodes();
	}

	void fillNodesPosition();
	void fillElementsWithNodes();
	void printElements();
};