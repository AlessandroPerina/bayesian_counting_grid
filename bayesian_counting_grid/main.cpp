#include "GeneralHeader.h"
#include "CountingGrid.h"
#include "Datapoint.h"
#include "DataReader.h"

int main()
{
	const int NO_GAMMAS = 5000;
	const float base_prior = BASE_PRIOR*WD_ROWS*WD_COLS;
	map<int, float> gammaLookUp;
	
        boost::mt19937 rng;
        rng.seed(123);

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

        //Gibbs Sampling
        
        int gibbsIter;
        std::vector<int> keysVec;
        int currKey;
        int asgnIdx;
        
        int rowAsgn, colAsgn;
        
        CountingGrid cg = CountingGrid( &gammaLookUp );
        fcolvec locPost;
        
        //Load all map keys in a vector to simplify shuffling
        for(std::map<int,Datapoint*>::iterator it = localData->begin(); it != localData->end(); ++it)
            keysVec.push_back(it->first);
        
        gibbsIter = 10;
        
        for (int iterId = 0; iterId < gibbsIter; iterId++)
        {
            //Shuffle keys vector to iterate over points in random order
            std::random_shuffle(keysVec.begin(), keysVec.end());
            
            for (std::vector<int>::iterator it = keysVec.begin(); it != keysVec.end(); it++)
            {
                currKey = *it;
                if (localData->at(currKey)->checkAsgn() == 1)
                {
                    //Remove Datapoint form Counting Grid
                    cg.removeDatapoint(localData->at(currKey));
                    //Reinitialise Datapoint location assignment
                    //Reinitialise Datapoint word assignment to empty
                }
                
                //Compute Datapoint location posterior 
                locPost = cg.locationPosterior(localData->at(currKey));
                
                //Sample location
                localData->at(currKey)->sampleLocation(locPost, &rng);

                //Map Datapoint to grid (distribute token)
                
                //Add Datapoint to Counting Grid and update counts
                cg.addDatapoint(localData->at(currKey));
            }
            
        }
        
        
        
	return 0;

}