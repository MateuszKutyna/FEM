#include "pch.h"

int main()
{
	GlobalData data;
	data.readFromFile("data.txt");
	FEM fem(data);
	fem.printElements();
}
