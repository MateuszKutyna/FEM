#include "pch.h"
void Agregaction::agregate() {
	for (int i = 0; i < data.nE; i++) {
		Jakobian jak(data,siatka);
		jak.calculate_H(i);
		jak.calculate_Hbc(i);
		jak.calculate_C(i);
		jak.calculate_P(i);
		elements[i] = jak;
	}
	
	std::array<Node,4> ID;
	for (int i = 0; i < data.nE;i++) {
		ID = elements[i].siatka.elements[i].element;
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				Hg[ID[j].t-1][ID[k].t-1] += elements[i].Hbc[k][j];	
				Cg[ID[j].t-1][ID[k].t-1] += elements[i].C[k][j];	
			}
		}
	}
	for (int i = 0; i < data.nE; i++) {
		ID = elements[i].siatka.elements[i].element;
		for (int j = 0; j < 4; j++) {
				Pg[ID[j].t - 1] += elements[i].P[j];
		}
	}
	
	

}

void Agregaction::print2D(const std::vector<std::vector<double>>& tab) {
	for (int i = 0; i < tab.size(); i++) {
		for (int j = 0; j < tab.size(); j++) {
			std::cout << tab[i][j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

void Agregaction::calculateTemperature() {
	agregate();
	//[C]/dt
	for (int i = 0; i < data.nN; i++) {
		for (int j = 0; j < data.nN; j++) {
			Cg[i][j] /= data.dt;
		}
	}
	//[H]+[C]/dT
	for (int i = 0; i < data.nN; i++) {
		for (int j = 0; j < data.nN; j++) {
			Hg[i][j] += Cg[i][j];
		}
	}
	std::vector<double> temp;
	temp.resize(data.nN);
	std::cout << "Time[s]\t" << "MinTemp[C]\t" << "MaxTemp[C]" << std::endl;
	for (int i = 0; i < data.sim_time; i += data.dt){	
		//[C]/dT*{t0}+{P}
		for (int i = 0; i < temp.size(); i++)
			temp[i] = 0;

		for (int i = 0; i < data.nN; i++) {
			for (int j = 0; j < data.nN; j++) {
				temp[i] += (Cg[i][j]) * data.t0[j];
			}
		}
		for (int i = 0; i < data.nN; i++) {
			temp[i] += Pg[i];
		}
		//[Hg]{t1]={Pg}
		gaussElimination(temp,i);
	}
	
}
bool Agregaction::gauss(const int& n, std::vector<std::vector<double>> AB, std::vector<double>& X) {
	int i, j, k;
	double m, s;
	const double eps = 1e-12; 
	for (i = 0; i < n - 1; i++)
	{
		for (j = i + 1; j < n; j++)
		{
			if (fabs(AB[i][i]) < eps) return false;
			m = -AB[j][i] / AB[i][i];
			for (k = i + 1; k <= n; k++)
				AB[j][k] += m * AB[i][k];
		}
	}
	for (i = n - 1; i >= 0; i--)
	{
		s = AB[i][n];
		for (j = n - 1; j >= i + 1; j--)
			s -= AB[i][j] * X[j];
		if (fabs(AB[i][i]) < eps) return false;
		X[i] = s / AB[i][i];
	}
	return true;
}
void Agregaction::gaussElimination(std::vector<double>& temp,int time_step) {
	for (int i = 0; i < data.nN; i++) {
		Hg[i][data.nN] = temp[i];
	}
	bool ifZero = false;
	ifZero=gauss(t1.size(), Hg, t1);
	if (!ifZero) exit(1);
	
	double max = searchMax(t1);
	double min = searchMin(t1);
	std::cout << time_step+1 << "\t" << min << "\t\t" << max<<std::endl;
	
	data.t0 = t1;
}
double Agregaction::searchMax(std::vector<double>& tab) {
	double max = -100000000000000.0;
	for (int i = 0; i < tab.size(); i++) {
		if (max < tab[i])max = tab[i];
	}
	return max;
}
double Agregaction::searchMin(std::vector<double>& tab) {
	double min = 100000000000000.0;
	for (int i = 0; i < tab.size(); i++) {
		if (min > tab[i])min = tab[i];
	}
	return min;
}