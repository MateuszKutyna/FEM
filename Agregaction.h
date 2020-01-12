#pragma once
class Agregaction {
public:
	std::vector<std::vector<double>> Hg;
	std::vector<std::vector<double>> Cg;
	std::vector<double> Pg;
	std::vector<double> t1;
	
	std::vector<Jakobian> elements;
	GlobalData data;
	FEM siatka;
	
	Agregaction(const GlobalData& _data,const FEM& _siatka) :data(_data),siatka(_siatka) {
		Hg.resize(data.nN, std::vector<double>(data.nN+1));
		Cg.resize(data.nN, std::vector<double>(data.nN));
		Pg.resize(data.nN);
		t1.resize(data.nN);
		elements.resize(data.nE);
	}
	void print2D(const std::vector<std::vector<double>>&);
	void agregate();
	void calculateTemperature();
	bool gauss(const int&, std::vector<std::vector<double>>,  std::vector<double>&);
	void gaussElimination(std::vector<double>&,int time_step);
	double searchMin(std::vector<double>&);
	double searchMax(std::vector<double>&);
};