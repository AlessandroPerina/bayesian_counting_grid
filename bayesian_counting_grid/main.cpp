#include "GeneralHeader.h"
#include "CountingGrid.h"
#include "Datapoint.h"
#include "DataReader.h"

int main()
{
	const int NO_GAMMAS = 5000;
	const float base_prior = BASE_PRIOR*WD_ROWS*WD_COLS;
	map<int, float> gammaLookUp;
	

	std::cout << "Ciao!" << endl;
        std::map<int, Datapoint*>* localData;
        Datapoint* d;

//	for (int g = 0; g < NO_GAMMAS; g++)
//	{
//		int key = g;
//		float newGamma = lgammaf(g + base_prior);
//		gammaLookUp.insert(std::pair<int, float>(key, newGamma));
//	}
//	std::cout << "Gamma lookup initialized" << endl;
//	float gm = 10 + base_prior;
//	cout << "Gamma of " << gm << " equal to " << gammaLookUp[gm] << endl;
//	system("Pause");
//
//	CountingGrid cgProva = CountingGrid( &gammaLookUp );
//	cgProva.get_a().print("a: ");
//	cgProva.get_Aw().print("Aw: ");
//	cgProva.get_logG().print("logG: ");

	cout << "Minchia se funziona..." << endl;
        
        DataReader* dr = new DataReader("/home/mzanotto/projects/code/Gibbs_CG/science_sub.txt");
        dr->loadData();
        
        localData = dr->getData();
        d = localData->at(3);
        cout<<d->getWords()[0]<<" "<<d->getWords()[1]<<" "<<d->getWords()[2]<<endl;
        cout<<d->getCountsArray()[d->getWords()[0]]<<" "<<d->getCountsArray()[d->getWords()[1]]<<" "<<d->getCountsArray()[d->getWords()[2]]<<endl;

	return 0;

}