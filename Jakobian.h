#pragma once
#include "pch.h"
class Jakobian
{
public:
	std::vector<double> dx_dE;
	std::vector<double> dx_dN;
	std::vector<double> dy_dE;
	std::vector<double> dy_dN;
	std::vector<double> det_Jakobian;
	std::vector<std::vector<double>> C;

	std::vector<std::vector<double>> J1;

	std::vector<std::vector<std::vector<double>>> dNdx_4_4;
	std::vector<std::vector<std::vector<double>>> dNdy_4_4;

	std::vector<std::vector<std::vector<double>>> Hpc;
	std::vector<std::vector<double>> H;
	std::vector<std::vector<double>> Hbc;

	std::vector<std::vector<double>> dNi_dx;
	std::vector<std::vector<double>> dNi_dy;

	GlobalData data;
	universalElement *uni_ele;
	FEM *siatka;
	Jakobian() = default;
	Jakobian(const GlobalData& _data) :data(_data) {
		dx_dE.resize(data.integralPoints);
		dx_dN.resize(data.integralPoints);
		dy_dE.resize(data.integralPoints);
		dy_dN.resize(data.integralPoints);
		C.resize(4, std::vector<double>(4)); //2D
		det_Jakobian.resize(data.integralPoints);
		J1.resize(data.integralPoints, std::vector<double>(4)); //2D
		dNi_dx.resize(data.integralPoints, std::vector<double>(4)); //2D
		dNi_dy.resize(data.integralPoints, std::vector<double>(4)); //2D
		Hpc.resize(data.integralPoints, std::vector<std::vector<double> >(4, std::vector<double>(4))); //3D
		dNdx_4_4.resize(data.integralPoints, std::vector<std::vector<double> >(4, std::vector<double>(4))); //3D
		dNdy_4_4.resize(data.integralPoints, std::vector<std::vector<double> >(4, std::vector<double>(4))); //3D
		H.resize(4, std::vector<double>(4));
		Hbc.resize(4, std::vector<double>(4));
		uni_ele = new universalElement(data);
		siatka = new FEM(data);
		
	}

	void calculate_H(int elementId);
	void calculate_Hbc(int elementId);
	void calculate_C(int elementId);
	void print_2D(std::vector<std::vector<double>> tab);
};
