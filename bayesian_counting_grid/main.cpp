#include "GeneralHeader.h"
#include "CountingGrid.h"
#include "Datapoint.h"
#include "DataReader.h"

int main()
{
	const int NO_GAMMAS = 50000;
	map<float, float> gammaLookUp;

	std::cout << "Ciao!" << endl;


	for (int g = 0; g < NO_GAMMAS; g++)
	{
		float newGamma = lgammaf(g + BASE_PRIOR);
		gammaLookUp.insert(std::pair<float, float>(float(g) + BASE_PRIOR, newGamma));
	}
	std::cout << "Gamma lookup initialized" << endl;
	float gm = 10 + BASE_PRIOR;
	cout << "Gamma of " << gm << " equal to " << gammaLookUp[gm] << endl;
	system("Pause");

	//CountingGrid cgProva = CountingGrid( &gammaLookUp );
	cout << "Minchia se funziona..." << endl;
        
        DataReader* dr = new DataReader("/home/mzanotto/projects/code/Gibbs_CG/science.txt");
        dr->loadData();

	return 0;

}