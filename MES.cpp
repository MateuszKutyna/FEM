#include "pch.h"

int main()
{
	GlobalData data;
	data.readFromFile("data1.txt");
	Grid fem(data);

	Aggregation agre(data,fem);
	agre.calculateTemperature();
	
}
