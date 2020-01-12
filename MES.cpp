#include "pch.h"

int main()
{
	GlobalData data;
	data.readFromFile("data1.txt");
	FEM fem(data);

	Agregaction agre(data,fem);
	agre.calculateTemperature();
	
}
