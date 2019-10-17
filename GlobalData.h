#pragma once
#include"pch.h"
struct GlobalData{
	int nH, nW;
	double H, W, nN, nE;
	void readFromFile(std::string);
};