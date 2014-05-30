#include "GeneralHeader.h"
#include "CountingGrid.h"
#include "Datapoint.h"
#include "DataReader.h"

int main()
{
	//const double NO_GAMMAS = MAX_GAMMA;
	const float base_prior = BASE_PRIOR*WD_ROWS*WD_COLS;
	map<unsigned long, float> gammaLookUp;
	
        boost::mt19937 rng;
        rng.seed(123);

        std::map<unsigned long, Datapoint*>* localData;
        Datapoint* d;

	for (unsigned long g = 0; g < MAX_GAMMA; g++)
	{
		//int key = g;
		//float newGamma = lgammaf(g + base_prior);
		//gammaLookUp.insert(std::pair<int, float>(key, newGamma));
                float newGamma = lgammaf((float)g);
		//gammaLookUp.insert(std::pair<int, float>(g, newGamma));
                gammaLookUp.insert(std::pair<unsigned long, float>(g, newGamma));
                //if (g<100)
                //    cout<<g<<endl;
	}
	std::cout << "Gamma lookup initialized [" << gammaLookUp.size()<<"]"<<endl;
	//float gm = 10 + base_prior;

//	system("Pause");
//
//	CountingGrid cgProva = CountingGrid( &gammaLookUp );
//	cgProva.get_a().print("a: ");
//	cgProva.get_Aw().print("Aw: ");
//	cgProva.get_logG().print("logG: ");

 //	cout << "Minchia se funziona..." << endl;
        
        //DataReader* dr = new DataReader("C:\\Users\\APerina\\Documents\\DataCG\\science_nips_reduced.txt");
        //DataReader* dr = new DataReader("/home/mzanotto/projects/code/Gibbs_CG/science_nips.txt");
        //DataReader* dr = new DataReader("/home/mzanotto/projects/code/Gibbs_CG/trees.txt");
        //DataReader* dr = new DataReader("/home/mzanotto/projects/code/Gibbs_CG/nips12.txt");
        DataReader* dr = new DataReader(FILENAME);
        dr->loadData();
        cout << "Data Loaded" << endl;
        localData = dr->getData();
        //d = localData->at(3);

        //Gibbs Sampling
        
        int gibbsIter;
        std::vector<unsigned long> keysVec;
        int currKey;
        int asgnIdx;
        
        clock_t begin, end;
        double time_spent;
        
        int rowAsgn, colAsgn;
        
        CountingGrid cg = CountingGrid( &gammaLookUp );
        
        cout<<"Counting Grid initialised"<<endl;
        fcolvec locPost;
        
        
        //Load all map keys in a vector to simplify shuffling
        for(std::map<unsigned long ,Datapoint*>::iterator it = localData->begin(); it != localData->end(); ++it)
            keysVec.push_back(it->first);
        
        gibbsIter = 5000;
        
        for (int iterId = 0; iterId < gibbsIter; iterId++)
        {

            //Shuffle keys vector to iterate over points in random order
            std::random_shuffle(keysVec.begin(), keysVec.end());
            //int t = 0;
            begin = clock();
            for (std::vector<unsigned long>::iterator it = keysVec.begin(); it != keysVec.end(); it++)
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
                //begin = clock();
                //locPost = cg.locationPosterior(localData->at(currKey));
                //locPost = cg.locationPosteriorLoop(localData->at(currKey));
                //t++;
                //cout<<t<<endl;
                locPost = cg.locationPosteriorLoopFast(localData->at(currKey));
                //locPost = cg.locationPosteriorLoopPar(localData->at(currKey));
                //locPost.t().print("locPost");
                //cout<<endl<<endl<<sum(locPost)<<endl;
                
                //cout<<localData->at(currKey)->getCountsArray()<<endl;
                //if (t == 2)
                    //break;
                //end = clock();
                //time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
                //cout<<time_spent<<endl;
                
                //Sample location
                //begin = clock();
                localData->at(currKey)->sampleLocation(locPost, &rng);
                //end = clock();
                //time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
                //cout<<time_spent<<endl;

                //Map Datapoint to grid (distribute token)
                //begin = clock();
		localData->at(currKey)->sampleTokenLocation(&cg, &rng);
                //end = clock();
                //time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
                //cout<<time_spent<<endl;

                //Add Datapoint to Counting Grid and update counts
                //begin = clock();
                cg.addDatapoint(localData->at(currKey));
                //end = clock();
                //time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
                //cout<<time_spent<<endl;
                //cout<<"Datapoint "<<currKey<<" [it: "<<iterId <<"]"<<endl;
            }
        end = clock();
        time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
        cout << "Iteration " << iterId << " completed in "<< time_spent<<" seconds"<<endl;
        if (iterId%500 == 0)
        {
            string filename;
            filename = (string)"cg" + to_string(iterId) + (string)".cg";
            cg.get_a().save(filename, raw_ascii);
        }
        //break;
        }
        
        
	cg.get_a().save("CountingGrid.cg", raw_ascii);
	return 0;

}