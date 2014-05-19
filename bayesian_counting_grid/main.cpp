#include "GeneralHeader.h"

int main()
{
	const int NO_GAMMAS = 50000;
	map<float, float> gammaLookUp;

	std::cout << "Ciao!";


	for (int g = 0; g < NO_GAMMAS; g++)
	{
		float newGamma = lgammaf(g + BASE_PRIOR);
		gammaLookUp.insert(std::pair<float, float>(float(g) + BASE_PRIOR, newGamma));
	}
	std::cout << "Gamma lookup initialized";
	system("Pause");

	return 0;

}