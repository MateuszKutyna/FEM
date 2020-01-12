#include"pch.h"
void Jakobian::calculate_H(int elementId) {
		
		for (int i = 0; i < data.integralPoints; i++) {
			for (int j = 0; j < 4; j++) {
				dx_dE[i] += uni_ele->dNi_dE[i][j] * siatka.elements[elementId].element[j].x;
				dx_dN[i] += uni_ele->dNi_dN[i][j] * siatka.elements[elementId].element[j].x;
				dy_dE[i] += uni_ele->dNi_dE[i][j] * siatka.elements[elementId].element[j].y;
				dy_dN[i] += uni_ele->dNi_dN[i][j] * siatka.elements[elementId].element[j].y;
			}
		}
		for (int i = 0; i < data.integralPoints; i++) {
				det_Jakobian[i] = dx_dE[i] * dy_dN[i] - dy_dE[i] * dx_dN[i];
		}
		for (int i = 0; i < data.integralPoints; i++) {
			J1[i][0] = (1 / det_Jakobian[i]) * dy_dN[i];
			J1[i][1] = (1 / det_Jakobian[i]) * (-dy_dE[i]);
			J1[i][2] = (1 / det_Jakobian[i]) * (-dx_dN[i]);
			J1[i][3] = (1 / det_Jakobian[i]) * dx_dE[i];

		}
		for(int i=0;i<data.integralPoints;i++){
			for (int j = 0; j < 4; j++) {
				dNi_dx[i][j] += J1[i][0] * uni_ele->dNi_dE[i][j] + J1[i][1] * uni_ele->dNi_dN[i][j];
				dNi_dy[i][j] += J1[i][2] * uni_ele->dNi_dE[i][j] + J1[1][3] * uni_ele->dNi_dN[i][j];
			}
		}
	
		for (int pc = 0; pc < data.integralPoints; pc++){
			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					dNdx_4_4[pc][i][j] = dNi_dx[pc][i] * dNi_dx[pc][j];
					dNdy_4_4[pc][i][j] = dNi_dy[pc][i] * dNi_dy[pc][j];
				}
			}
		}
		for (int pc = 0; pc < data.integralPoints; pc++) {
			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					Hpc[pc][i][j] = data.k*(dNdx_4_4[pc][i][j] + dNdy_4_4[pc][i][j]);
				}
			}
		}
		if(data.integralPoints==4){
			for (int pc = 0; pc < data.integralPoints; pc++) {
				for (int i = 0; i < 4; i++) {
					for (int j = 0; j < 4; j++) {
						H[i][j] += Hpc[pc][i][j] * uni_ele->w[0] *uni_ele->w[1]*det_Jakobian[pc];
					}
				}
			}
		}
		else if(data.integralPoints==9){
			for (int pc = 0; pc < data.integralPoints; pc++) {
				for (int i = 0; i < 4; i++) {
					for (int j = 0; j < 4; j++) {
						H[i][j] += Hpc[pc][i][j] * uni_ele->w[0] * uni_ele->w[1] * uni_ele->w[2] * det_Jakobian[pc];
					}
				}
			}
		}
		
}
void Jakobian::calculate_Hbc(int elementId) {
	bool indexy[4];
	for (int i = 0; i < 4; i++)
		indexy[i]=0;
	for (int i = 0; i < 4; i++)
		if (siatka.elements[elementId].element[i].warunek_brzegowy == true)
			indexy[i] = true;
	std::vector<std::vector<double>> N;
	std::vector<std::vector<double>> N1;
	N.resize(4, std::vector<double>(4)); //2D
	N1.resize(4, std::vector<double>(4)); //2D
	double detJx;
	double detJy;
	detJx = (siatka.elements[elementId].element[1].x - siatka.elements[elementId].element[0].x) / 2;
	detJy = (siatka.elements[elementId].element[3].y - siatka.elements[elementId].element[0].y) / 2;
	if (indexy[0] && indexy[1] && !indexy[2] && indexy[3]) {
		
		for (int i = 0; i < 4; i++) {

			for (int j = 0; j < 4; j++) {
				N[i][j] += uni_ele->NiS[0][i] * uni_ele->NiS[0][j] * uni_ele->w[0];
				N[i][j] += uni_ele->NiS[1][i] * uni_ele->NiS[1][j] * uni_ele->w[1];
				N1[i][j] += uni_ele->NiS[6][i] * uni_ele->NiS[6][j] * uni_ele->w[0];
				N1[i][j] += uni_ele->NiS[7][i] * uni_ele->NiS[7][j] * uni_ele->w[1];
			}
		}
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				N[i][j] *= detJx;
				N1[i][j] *= detJy;
			}
		}
	}
	else if (indexy[0] &&!indexy[1]&& !indexy[2] && indexy[3]) {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				N[i][j] += uni_ele->NiS[6][i] * uni_ele->NiS[6][j] * uni_ele->w[0];
				N[i][j] += uni_ele->NiS[7][i] * uni_ele->NiS[7][j] * uni_ele->w[1];
			}
		}
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				N[i][j] *= detJx;
			}
		}
	}
	else if (indexy[0]&& !indexy[1] && indexy[2] && indexy[3]) {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				N[i][j] += uni_ele->NiS[4][i] * uni_ele->NiS[4][j] * uni_ele->w[0];
				N[i][j] += uni_ele->NiS[5][i] * uni_ele->NiS[5][j] * uni_ele->w[1];
				N1[i][j] += uni_ele->NiS[6][i] * uni_ele->NiS[6][j] * uni_ele->w[0];
				N1[i][j] += uni_ele->NiS[7][i] * uni_ele->NiS[7][j] * uni_ele->w[1];
			}
		}
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				N[i][j] *= detJx;
				N1[i][j] *= detJy;
			}
		}
	}
	else if (indexy[0] && indexy[1] && !indexy[2] && !indexy[3]) {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				N[i][j] += uni_ele->NiS[0][i] * uni_ele->NiS[0][j] * uni_ele->w[0];
				N[i][j] += uni_ele->NiS[1][i] * uni_ele->NiS[1][j] * uni_ele->w[1];
			}
		}
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				N[i][j] *= detJx;
			}
		}
	}
	else if (!indexy[0] && !indexy[1]&&indexy[2] && indexy[3]) {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				N[i][j] += uni_ele->NiS[4][i] * uni_ele->NiS[4][j] * uni_ele->w[0];
				N[i][j] += uni_ele->NiS[5][i] * uni_ele->NiS[5][j] * uni_ele->w[1];
			}
		}
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				N[i][j] *= detJy;
			}
		}
	}
	else if (indexy[0] && indexy[1] && indexy[2] && !indexy[3]) {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				N[i][j] += uni_ele->NiS[0][i] * uni_ele->NiS[0][j] * uni_ele->w[0];
				N[i][j] += uni_ele->NiS[1][i] * uni_ele->NiS[1][j] * uni_ele->w[1];
				N1[i][j] += uni_ele->NiS[2][i] * uni_ele->NiS[2][j] * uni_ele->w[0];
				N1[i][j] += uni_ele->NiS[3][i] * uni_ele->NiS[3][j] * uni_ele->w[1];
			}
		}
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				N[i][j] *= detJx;
				N1[i][j] *= detJy;
			}
		}
	}
	else if (!indexy[0]&&indexy[1] && indexy[2]&&!indexy[3]) {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				N[i][j] += uni_ele->NiS[2][i] * uni_ele->NiS[2][j] * uni_ele->w[0];
				N[i][j] += uni_ele->NiS[3][i] * uni_ele->NiS[3][j] * uni_ele->w[1];
			}
		}
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				N[i][j] *= detJy;
			}
		}
	}
	else if (!indexy[0] && indexy[1] && indexy[2] &&indexy[3]) {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				N[i][j] += uni_ele->NiS[2][i] * uni_ele->NiS[2][j] * uni_ele->w[0];
				N[i][j] += uni_ele->NiS[3][i] * uni_ele->NiS[3][j] * uni_ele->w[1];
				N1[i][j] += uni_ele->NiS[4][i] * uni_ele->NiS[4][j] * uni_ele->w[0];
				N1[i][j] += uni_ele->NiS[5][i] * uni_ele->NiS[5][j] * uni_ele->w[1];
			}
		}
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				N[i][j] *= detJx;
				N1[i][j] *= detJy;
			}
		}
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			N[i][j] *= data.alfa;
			N1[i][j] *= data.alfa;
		}

	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			N[i][j] += N1[i][j];
		}
	}
	
	
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			Hbc[i][j] = N[i][j]+H[i][j];
		}
	}
	
	
}
void Jakobian::calculate_C(int elementId) {
	
	for (int pc = 0; pc < data.integralPoints; pc++) {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				C[i][j] += uni_ele->NiV[pc][i] * uni_ele->NiV[pc][j] * det_Jakobian[pc] * uni_ele->w[0] * uni_ele->w[1]*data.ro*data.cw;
			}
		}
	}
	
}
void Jakobian::calculate_P(int elementId) {
	bool indexy[4];
	for (int i = 0; i < 4; i++)
		indexy[i] = 0;
	for (int i = 0; i < 4; i++)
		if (siatka.elements[elementId].element[i].warunek_brzegowy == true)
			indexy[i] = true;
	double detJx;
	double detJy;
	detJx = (siatka.elements[elementId].element[1].x - siatka.elements[elementId].element[0].x) / 2;
	detJy = (siatka.elements[elementId].element[3].y - siatka.elements[elementId].element[0].y) / 2;
	std::vector<double> P1;
	std::vector<double>P2;
	P1.resize(4); //2D
	P2.resize(4); //2D
	if (indexy[0] && indexy[1] && !indexy[2] && indexy[3]) {
		for (int i = 0; i < 4; i++) {
			P1[i] += uni_ele->NiS[0][i] * data.t8*data.alfa*uni_ele->w[0];
			P1[i] += uni_ele->NiS[1][i] * data.t8*data.alfa*uni_ele->w[1];
			P2[i] += uni_ele->NiS[6][i] * data.t8*data.alfa*uni_ele->w[0];
			P2[i] += uni_ele->NiS[7][i] * data.t8*data.alfa*uni_ele->w[1];
		}
		for (int i = 0; i < 4; i++) {
			P1[i] *= detJx;
			P2[i] *= detJy;
		}
	}
	else if (indexy[0] && !indexy[1] && !indexy[2] && indexy[3]) {
		for (int i = 0; i < 4; i++) {
			P2[i] += uni_ele->NiS[6][i] * data.t8*data.alfa*uni_ele->w[0];
			P2[i] += uni_ele->NiS[7][i] * data.t8*data.alfa*uni_ele->w[1];
		}
		for (int i = 0; i < 4; i++) {
			P2[i] *= detJy;
		}

	}
	else if (indexy[0] && !indexy[1] && indexy[2] && indexy[3]) {
		for (int i = 0; i < 4; i++) {
			P1[i] += uni_ele->NiS[6][i] * data.t8*data.alfa*uni_ele->w[0];
			P1[i] += uni_ele->NiS[7][i] * data.t8*data.alfa*uni_ele->w[1];
			P2[i] += uni_ele->NiS[4][i] * data.t8*data.alfa*uni_ele->w[0];
			P2[i] += uni_ele->NiS[5][i] * data.t8*data.alfa*uni_ele->w[1];
		}
		for (int i = 0; i < 4; i++) {
			P1[i] *= detJx;
			P2[i] *= detJy;
		}
	}
	else if (indexy[0] && indexy[1] && !indexy[2] && !indexy[3]) {
		for (int i = 0; i < 4; i++) {
			P1[i] += uni_ele->NiS[0][i] * data.t8*data.alfa*uni_ele->w[0];
			P1[i] += uni_ele->NiS[1][i] * data.t8*data.alfa*uni_ele->w[1];
		}
		for (int i = 0; i < 4; i++) {
			P1[i] *= detJx;
		}
	}
	else if (!indexy[0] && !indexy[1] && indexy[2] && indexy[3]) {
		for (int i = 0; i < 4; i++) {
			P1[i] += uni_ele->NiS[4][i] * data.t8*data.alfa*uni_ele->w[0];
			P1[i] += uni_ele->NiS[5][i] * data.t8*data.alfa*uni_ele->w[1];
		}
		for (int i = 0; i < 4; i++) {
			P1[i] *= detJx;
		}
	}
	else if (indexy[0] && indexy[1] && indexy[2] && !indexy[3]) {
		for (int i = 0; i < 4; i++) {
			P1[i] += uni_ele->NiS[0][i] * data.t8*data.alfa*uni_ele->w[0];
			P1[i] += uni_ele->NiS[1][i] * data.t8*data.alfa*uni_ele->w[1];
			P2[i] += uni_ele->NiS[2][i] * data.t8*data.alfa*uni_ele->w[0];
			P2[i] += uni_ele->NiS[3][i] * data.t8*data.alfa*uni_ele->w[1];
		}
		for (int i = 0; i < 4; i++) {
			P1[i] *= detJx;
			P2[i] *= detJy;
		}
	}
	else if (!indexy[0] && indexy[1] && indexy[2] && !indexy[3]) {
		for (int i = 0; i < 4; i++) {
			P2[i] += uni_ele->NiS[2][i] * data.t8*data.alfa*uni_ele->w[0];
			P2[i] += uni_ele->NiS[3][i] * data.t8*data.alfa*uni_ele->w[1];
		}
		for (int i = 0; i < 4; i++) {
			P2[i] *= detJy;
		}
	}
	else if (!indexy[0] && indexy[1] && indexy[2] && indexy[3]) {
		for (int i = 0; i < 4; i++) {
			P1[i] += uni_ele->NiS[4][i] * data.t8*data.alfa*uni_ele->w[0];
			P1[i] += uni_ele->NiS[5][i] * data.t8*data.alfa*uni_ele->w[1];
			P2[i] += uni_ele->NiS[2][i] * data.t8*data.alfa*uni_ele->w[0];
			P2[i] += uni_ele->NiS[3][i] * data.t8*data.alfa*uni_ele->w[1];
		}
		for (int i = 0; i < 4; i++) {
			P1[i] *= detJx;
			P2[i] *= detJy;
		}
	}
	
	for (int i = 0; i < 4; i++) {
		P[i]=P1[i] + P2[i];
	}
	
}
void Jakobian::print_2D(std::vector<std::vector<double>> tab) {
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			std::cout << tab[i][j] << " ";
		}
		std::cout << std::endl;
	}
}
