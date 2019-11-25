#include "pch.h"

int main()
{
	GlobalData data;
	data.readFromFile("data.txt");
	// fem(data);
	//universalElement test(data);
	//test.Complete();
	Jakobian J(data);
	J.siatka->printElements();
	std::cout << std::endl;
	J.calculate_H(0);
	std::cout << std::endl;
	J.calculate_Hbc(0);
	std::cout << std::endl;
	J.calculate_C(0);
	

	//fem.printElements();

}
