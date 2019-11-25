#pragma once
#include"pch.h"
struct GlobalData{
	int nH, nW, integralPoints;
	double H, W, nN, nE, k, alfa, ro, cw;
	void readFromFile(std::string);
};