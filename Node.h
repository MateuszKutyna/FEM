#pragma once
#include"GlobalData.h"
struct Node
{
	Node() = default;

	Node(double _x, double _y,bool _warunek) : x(_x), y(_y) {

		licznik++;
		t = licznik;
		warunek_brzegowy = _warunek;
	}
	int t;
	double x = 0.0, y = 0.0;
	static int licznik;
	bool warunek_brzegowy;
};