#pragma once
#include"pch.h"
class universalElement {
public:
	GlobalData data;
	std::vector<double> p;
	std::vector<double> w;
	std::vector<std::vector<double>> NiS;
	std::vector<std::vector<double>> NiV;
	std::vector<std::vector<double>> dNi_dE;
	std::vector<std::vector<double>> dNi_dN;
		 
	


	universalElement() = default;
	universalElement(const GlobalData& _data) :data(_data){
		//Shape function and their derivatives
		NiS.resize(8, std::vector<double>(4));
		NiV.resize(8, std::vector<double>(4));
		dNi_dE.resize(data.integralPoints, std::vector<double>(4));
		dNi_dN.resize(data.integralPoints, std::vector<double>(4));

		//Integral pontis and their wages
		p.resize(data.integralPoints/sqrt(data.integralPoints));
		w.resize(data.integralPoints/sqrt(data.integralPoints));
		//Integral points cases 
		if (data.integralPoints == 4) {
			p[0] = -1.0 / sqrt(3);
			p[1] = 1.0 / sqrt(3);
			w[0] = 1;//Sum = 2
			w[1] = 1;//
		}
		else if(data.integralPoints==9){
			p[0] = -0.77;
			p[1] = 0.0;
			p[2] = 0.77;
			w[0] = 5.0 / 9.0;//
			w[1] = 8.0 / 9.0;//Sum = 2
			w[2] = 5.0 / 9.0;//
		}
		Complete();
	}
	//To use class metods as a argument of function we need template
	template<typename Fun>
	void Fill(const Fun&, std::vector<std::vector<double>>&,int x );
	template<typename FunN>
	void FillN(const FunN&, std::vector<std::vector<double>>&,int x );

	void Complete();
	//Shape functions
	double N1(double E, double N);
	double N2(double E, double N);
	double N3(double E, double N);
	double N4(double E, double N);
	//Derivative of a shape function (Xi)
	double dN1_dE(double E, double N);
	double dN2_dE(double E, double N);
	double dN3_dE(double E, double N);
	double dN4_dE(double E, double N);
	//Derivative of a shape function (Eta)
	double dN1_dN(double E, double N);
	double dN2_dN(double E, double N);
	double dN3_dN(double E, double N);
	double dN4_dN(double E, double N);
	
};
//Calculating for given points
template<typename Fun>

void universalElement::Fill(const Fun& shapeFunction, std::vector<std::vector<double>> &tab2D,int iterationNumber) {
	if (data.integralPoints == 4) {
			tab2D[0][iterationNumber] = shapeFunction(p[0], p[0]);
			tab2D[1][iterationNumber] = shapeFunction(p[0], p[1]);
			tab2D[2][iterationNumber] = shapeFunction(p[1], p[1]);
			tab2D[3][iterationNumber] = shapeFunction(p[1], p[0]);

	}
	else if (data.integralPoints == 9) {
		
			tab2D[0][iterationNumber] = shapeFunction(p[2], p[2]);
			tab2D[1][iterationNumber] = shapeFunction(p[2], p[1]);
			tab2D[2][iterationNumber] = shapeFunction(p[2], p[0]);
			tab2D[3][iterationNumber] = shapeFunction(p[1], p[0]);
			tab2D[4][iterationNumber] = shapeFunction(p[0], p[0]);
			tab2D[5][iterationNumber] = shapeFunction(p[0], p[1]);
			tab2D[6][iterationNumber] = shapeFunction(p[0], p[2]);
			tab2D[7][iterationNumber] = shapeFunction(p[1], p[2]);
			tab2D[8][iterationNumber] = shapeFunction(p[1], p[1]);
	}

}
template<typename FunN>
void universalElement::FillN(const FunN& shapeFunction, std::vector<std::vector<double>> &tab2D, int iterationNumber) {
	if (data.integralPoints == 4) {
		tab2D[0][iterationNumber] = shapeFunction(-1 / sqrt(3), -1);
		tab2D[1][iterationNumber] = shapeFunction(1 / sqrt(3), -1);
		tab2D[2][iterationNumber] = shapeFunction(-1 / sqrt(3), 1);
		tab2D[3][iterationNumber] = shapeFunction(1 / sqrt(3), 1);
		tab2D[4][iterationNumber] = shapeFunction(-1, -1 / sqrt(3));
		tab2D[5][iterationNumber] = shapeFunction(-1, 1 / sqrt(3));
		tab2D[6][iterationNumber] = shapeFunction(1, -1 / sqrt(3));
		tab2D[7][iterationNumber] = shapeFunction(1, 1 / sqrt(3));
	}
	else if (data.integralPoints == 9) {}
}