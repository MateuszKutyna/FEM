#pragma once
#include"pch.h"
struct Grid
{
	std::vector<Node> nodes;
	std::vector<Element> elements;
	GlobalData data;
	Grid() = default;
	Grid(const GlobalData &data1) :data(data1) {
		nodes.reserve(data.nN);
		elements.resize(data.nE);
		fillNodesPosition();
		fillElementsWithNodes();
	}

	void fillNodesPosition();
	void fillElementsWithNodes();
	void printElements();
};