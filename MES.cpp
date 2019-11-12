#include "pch.h"

int main()
{
	GlobalData data;
	data.readFromFile("data.txt");
	FEM fem(data);
	universalElement test(data);
	test.Complete();
	for (int i = 0; i < data.integralPoints; i++) {
		for (int j = 0; j < 4; j++) {
			std::cout << test.Ni[i][j]<<" " ;
		}
		std::cout << std::endl;
	}
	//fem.printElements();

}
