#include"pch.h"

void universalElement::Complete() {
	
	auto f1 = [=](double E, double N) { return this->N1(E, N); };
	auto f2 = [=](double E, double N) { return this->N2(E, N); };
	auto f3 = [=](double E, double N) { return this->N3(E, N); };
	auto f4 = [=](double E, double N) { return this->N4(E, N); };

	auto f1_dE = [=](double E, double N) { return this->dN1_dE(E, N); };
	auto f2_dE = [=](double E, double N) { return this->dN2_dE(E, N); };
	auto f3_dE = [=](double E, double N) { return this->dN3_dE(E, N); };
	auto f4_dE = [=](double E, double N) { return this->dN4_dE(E, N); };

	auto f1_dN = [=](double E, double N) { return this->dN1_dN(E, N); };
	auto f2_dN = [=](double E, double N) { return this->dN2_dN(E, N); };
	auto f3_dN = [=](double E, double N) { return this->dN3_dN(E, N); };
	auto f4_dN = [=](double E, double N) { return this->dN4_dN(E, N); };

	Fill(f1, Ni, 0);
	Fill(f2, Ni, 1);
	Fill(f3, Ni, 2);
	Fill(f4, Ni, 3);

	Fill(f1_dE, dNi_dE, 0);
	Fill(f2_dE, dNi_dE, 1);
	Fill(f3_dE, dNi_dE, 2);
	Fill(f4_dE, dNi_dE, 3);

	Fill(f1_dN, dNi_dN, 0);
	Fill(f2_dN, dNi_dN, 1);
	Fill(f3_dN, dNi_dN, 2);
	Fill(f4_dN, dNi_dN, 3);
}
// Shape functions
double universalElement::N1(double E, double N){
	return 1.0 / 4.0 * ((1.0 - E)*(1.0 - N));
}
double universalElement::N2(double E, double N){
	return 1.0 / 4.0 * ((1.0 + E)*(1.0 - N));
}
double universalElement::N3(double E, double N){
	return 1.0 / 4.0 * ((1.0 + E)*(1.0 + N));
}
double universalElement::N4(double E, double N){
	return 1.0 / 4.0 * ((1.0 - E)*(1.0 + N));
}
//Derivative of a shape function (Xi)
double universalElement::dN1_dE(double E,double N) {
	return -1.0 / 4.0 * (1.0 - N);
}
double universalElement::dN2_dE(double E, double N) {
	return 1.0 / 4.0 * (1.0 - N);
}
double universalElement::dN3_dE(double E, double N) {
	return 1.0 / 4.0 * (1.0 + N);
}
double universalElement::dN4_dE(double E, double N) {
	return -1.0 / 4.0 * (1.0 + N);
}
//Derivative of a shape function (Eta)
double universalElement::dN1_dN(double E, double N) {
	return -1.0 / 4.0 * (1.0 - E);
}
double universalElement::dN2_dN(double E, double N) {
	return -1.0 / 4.0 * (1.0 + E);
}
double universalElement::dN3_dN(double E, double N) {
	return 1.0 / 4.0 * (1.0 + E);
}
double universalElement::dN4_dN(double E, double N) {
	return 1.0 / 4.0 * (1.0 - E);
}