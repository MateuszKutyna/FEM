#pragma once
#include"pch.h"
struct GlobalData{
	int nH, nW, integralPoints;
	double H,
		W,
		nN,
		nE,
		k,//conductivity
		alfa,
		ro,//density
		cw,//specific heat
		t8,//ambient temperature
		dt,//temperature step
		sim_time;//simulation time
	std::vector<double> t0;//initial temperature
	void readFromFile(std::string);
};