#include "CountingGrid.h"
// #include "Datapoint.h"
// #include "GeneralHeader.h"

CountingGrid::CountingGrid(map<int, float>* gl)
{
	this->gammaLookUp = gl;

	this->a = fcube(CG_ROWS, CG_COLS, Z);
	//this->a = floor( 5 * arma::randu<fcube>(CG_ROWS, CG_COLS, Z)) + BASE_PRIOR;
	this->a.fill(BASE_PRIOR);

	this->Aw = fcube(CG_ROWS, CG_COLS, Z);
	this->Aw.fill(0);

	this->logG = fcube(CG_ROWS, CG_COLS, Z);
	this->logG.fill(0);

	this->logGsum = fmat(CG_ROWS, CG_COLS);
	this->logGsum.fill(0);

	/*
	fcube tmp = fcube(CG_ROWS+WD_ROWS-1, CG_COLS+WD_COLS-1, Z);
	//tmp.fill(1);
        tmp.fill(BASE_PRIOR);
	tmp(span(0, CG_ROWS - 1), span(0, CG_COLS - 1), span::all ) = this->a(span(0, CG_ROWS - 1), span(0, CG_COLS - 1), span::all);
	tmp(span(0, CG_ROWS - 1), span(CG_COLS, CG_COLS + WD_COLS - 2), span::all) = this->a(span(0, CG_ROWS - 1), span(0, WD_COLS - 2), span::all);
	tmp(span(CG_ROWS, CG_ROWS + WD_ROWS - 2), span(0, CG_COLS - 1), span::all) = this->a(span(0, WD_ROWS - 2), span(0, CG_COLS - 1), span::all);
	tmp( span(CG_ROWS, CG_ROWS + WD_ROWS - 2), span(CG_COLS, CG_COLS + WD_COLS - 2), span::all) = this->a(span(0, WD_ROWS - 2), span(0,WD_COLS - 2), span::all);

	for (int r = 0; r < CG_ROWS; r++)
	{
		for (int c = 0; c < CG_COLS; c++)
		{
			fcube tmp2 = reshape(tmp( span(r, r + WD_ROWS - 1), span(c, c + WD_COLS - 1),span::all), WD_ROWS*WD_COLS, Z, 1);
			fvec tmp3 = sum(tmp2.slice(0),0);
			this->Aw(span(r, r), span(c, c), span::all) = tmp3;
		}
	}
	*/

	for (int z = 0; z < Z; z++)
	{
		for (int c = 0; c < CG_COLS; c++)
		{
			for (int r = 0; r < CG_ROWS; r++)
			{
				for (int wRow = 0; wRow < WD_ROWS; wRow++)
				{
					for (int wCol = 0; wCol < WD_COLS; wCol++)
					{
						//int rIdx = (int)(r + wRow) % WD_ROWS;
						//int cIdx = (int)(c + wCol) % WD_COLS;
                                                int rIdx = (int)(r + wRow) % CG_ROWS;
						int cIdx = (int)(c + wCol) % CG_COLS;
						this->Aw(r, c, z) += this->a(rIdx, cIdx, z);
					}
				}
			}
		}
	}

	// Inefficiente ma va fatto una volta sola
	for (int r = 0; r < CG_ROWS; r++)
	{
		for (int c = 0; c < CG_COLS; c++)
		{
			// Calcolo logGamma
			for (int z = 0; z < Z; z++)
			{

				map<int, float>::iterator it;
				float tmpA = this->Aw(r, c, z);
				float costa = WD_ROWS*WD_COLS*BASE_PRIOR;
				int key = (int)round(this->Aw(r, c, z) - WD_ROWS*WD_COLS*BASE_PRIOR);
				it = this->gammaLookUp->find( key );
				if (it == this->gammaLookUp->end()){
					// Calcola la nuova gamma
					float newGamma = lgammaf(this->Aw(r, c, z));
					// Aggiungila alla lookup table
					this->gammaLookUp->insert(std::pair<int, float>(key, newGamma));
					// Aggiorna la Counting Grid
					this->logG(r, c, z) = newGamma;
                                        cout<<newGamma<<" inserted"<<endl;
				}
				else{
					// Aggiorna la Counting Grid
					this->logG(r, c, z) = it->second;
				}
			}

			// Calcolo logGammaSum
			map<int, float>::iterator it;
			fvec tmp1 = this->Aw(span(r, r), span(c, c), span::all);
			tmp1 = this->Aw(span(r, r), span(c, c), span::all) - WD_ROWS*WD_COLS*BASE_PRIOR;
			float sumTmp1 = sum(fvec(this->Aw(span(r, r), span(c, c), span::all)));

			int keySum = (int)round(sum(fvec(this->Aw(span(r, r), span(c, c), span::all)) - WD_ROWS*WD_COLS*BASE_PRIOR));
			float trueValue = sum(fvec(this->Aw(span(r, r), span(c, c), span::all)));

			it = this->gammaLookUp->find(keySum);
			if (it == this->gammaLookUp->end()){
				// Calcola la nuova gamma
				float newGamma = lgammaf(trueValue);
				// Aggiungila alla lookup table
				this->gammaLookUp->insert(std::pair<int, float>(keySum, newGamma));
				// Aggiorna la Counting Grid
				this->logGsum(r, c ) = newGamma;
			}
			else{
				// Aggiorna la Counting Grid
				this->logGsum(r, c) = it->second;
			}

		}
	}
}

CountingGrid::~CountingGrid()
{
}

int CountingGrid::sumAllWindows()
{

	this->Aw = fcube(CG_ROWS, CG_COLS, Z);
	this->Aw.fill(0);

	fcube padded = fcube(CG_ROWS + WD_ROWS - 1, CG_COLS + WD_COLS - 1, Z);
	padded(span(0, CG_ROWS - 1), span(0, CG_COLS - 1), span::all) = this->a(span(0, CG_ROWS - 1), span(0, CG_COLS - 1), span::all);
	padded(span(0, CG_ROWS - 1), span(CG_COLS, CG_COLS + WD_COLS - 2), span::all) = this->a(span(0, CG_ROWS - 1), span(0, WD_COLS - 2), span::all);
	padded(span(CG_ROWS, CG_ROWS + WD_ROWS - 2), span(0, CG_COLS - 1), span::all) = this->a(span(0, WD_ROWS - 2), span(0, CG_COLS - 1), span::all);
	padded(span(CG_ROWS, CG_ROWS + WD_ROWS - 2), span(CG_COLS, CG_COLS + WD_COLS - 2), span::all) = this->a(span(0, WD_ROWS - 2), span(0, WD_COLS - 2), span::all);

	for (int r = 0; r < CG_ROWS; r++)
	{
		for (int c = 0; c < CG_COLS; c++)
		{
			fcube tmp2 = reshape(padded(span(r, r + WD_ROWS - 1), span(c, c + WD_COLS - 1), span::all), WD_ROWS*WD_COLS, Z, 1);
			this->Aw(span(r, r), span(c, c), span::all).fill((float)arma::accu(tmp2.slice(0)));
		}
	}

	return 0;
}

int CountingGrid::sumAllWindowsLoop()
{
    int rIdx, cIdx;
	//this->Aw.fill(0);
        //loop over nonzero 
        for (int z = 0;  z < Z; z++)
        {
            for (int r = 0; r < CG_ROWS; r++)
            {
                for (int c = 0; c < CG_COLS; c++)
                {
                    for (int wRow = 0; wRow < WD_ROWS; wRow++)
                    {

                        for (int wCol = 0; wCol < WD_COLS; wCol++)
                        {
                            //rIdx = (int)(r + wRow)%WD_ROWS;
                            //cIdx = (int)(c + wCol)%WD_COLS;
                            rIdx = (int)(r + wRow)%CG_ROWS;
                            cIdx = (int)(c + wCol)%CG_COLS;
                            this->Aw(r, c, z) += this->a(rIdx,cIdx,z);
                        }
                    }
                }
            }
        }
	return 0;
}

int CountingGrid::updateAw(Datapoint* dp)
{
	// changes only slices where the Datapoint has data
    urowvec localWords;
    urowvec::iterator wordsIt;
    int rIdx, cIdx;
    int r, c, wRow, wCol;

    localWords = dp->getWords();
    //localWords.print("Words:");
        for( wordsIt = localWords.begin(); wordsIt != localWords.end(); wordsIt++)
        {
            //cout<<"Word: "<<*wordsIt<<endl;
            this->Aw.slice(*wordsIt).fill(0);
            
            for (c = 0; c < CG_COLS; c++)
            {
                for (r = 0; r < CG_ROWS; r++)
                {
                    for (wRow = 0; wRow < WD_ROWS; wRow++)
                    {

                        for (wCol = 0; wCol < WD_COLS; wCol++)
                        {
                            //rIdx = (int)(r + wRow)%WD_ROWS;
                            //cIdx = (int)(c + wCol)%WD_COLS;
                            rIdx = (int)(r + wRow)%CG_ROWS;
                            cIdx = (int)(c + wCol)%CG_COLS;
                            this->Aw(r, c, *wordsIt) += this->a(rIdx,cIdx, *wordsIt);
                        }
                    }
                }
            }
        }
	return 0;
}

int CountingGrid::addDatapoint(Datapoint* dp )
{
	map<int, arma::sp_fmat> tmpMap = dp->getTokenLoc();

        //cout<<"update loop"<<endl;
	for (map<int, arma::sp_fmat>::iterator it = tmpMap.begin(); it != tmpMap.end(); it++)
	{
		this->a.slice(it->first) += it->second;
	}
	// Aggiorno la variabile somma
        //cout<<"updating sum"<<endl;
	//this->sumAllWindowsLoop();
        this->updateAw(dp);
        //cout<<"computing logGamma"<<endl;
	this->computeLogGammaCG(dp->getTokenLoc());
	//cout<<"done"<<endl;
	return 0;

}

int CountingGrid::removeDatapoint( Datapoint* dp ) 
{
	map<int, arma::sp_fmat> tmpMap = dp->getTokenLoc();

	for (map<int, arma::sp_fmat>::iterator it = tmpMap.begin(); it != tmpMap.end(); it++)
	{
		this->a.slice(it->first) -= it->second;
	}
	// Aggiorno la variabile somma
	//this->sumAllWindowsLoop();
        this->updateAw(dp);
	this->computeLogGammaCG(dp->getTokenLoc());

	return 0;
}


fcolvec CountingGrid::locationPosterior(Datapoint* dp)
{

	fcube tmp = reshape(this->logG, CG_ROWS*CG_COLS, Z, 1);
	fmat T3 = reshape(sum(tmp.slice(0), 1), CG_ROWS, CG_COLS);
	T3 = arma::accu(T3) - T3;

	//fmat T1 = fmat(CG_ROWS, CG_COLS, fill::zeros);
	fmat T1 = reshape(sum(tmp.slice(0), 1), CG_ROWS, CG_COLS);

	fmat T2 = fmat(CG_ROWS, CG_COLS, fill::zeros);

	fmat T4 = fmat(CG_ROWS, CG_COLS, fill::zeros);
	T4 = arma::accu(this->logGsum) - this->logGsum;


        mapdata tmpMap = dp->getCountsDict();

	
	//	fmat(CG_COLS, CG_ROWS, fill::zeros);
	int accumT2Key = 0;
	float accumT2 = 0;
	// itera solo sulle features che ha il campione corrente
        clock_t begin, end;
        double time_spent;
        begin = clock();
        
	for (int c = 0; c < CG_COLS; c++)
	{
		for (int r = 0; r < CG_ROWS; r++)
		{

			
			accumT2Key = (int) round(sum(fvec(Aw(span(r, r), span(c, c), span::all) - WD_ROWS*WD_COLS*BASE_PRIOR)));
			accumT2 = sum(fvec(Aw(span(r, r), span(c, c), span::all) ));
			for (mapdata::iterator featureIt = tmpMap.begin(); featureIt != tmpMap.end(); featureIt++)
			{

				accumT2Key += featureIt->second;
				accumT2 += float(featureIt->second);
				// sn qui
				map<int, float>::iterator itLogGamma;

				int keyTmpCount = (int)round(Aw(r, c, featureIt->first) - WD_ROWS*WD_COLS*BASE_PRIOR + featureIt->second );
				float trueCount = Aw(r, c, featureIt->first) + float(featureIt->second);
				itLogGamma = this->gammaLookUp->find(keyTmpCount);
				if (itLogGamma == this->gammaLookUp->end()){
					float newGamma = lgammaf(trueCount);
					this->gammaLookUp->insert(std::pair<int, float>(keyTmpCount, newGamma));
					T1(r, c) = T1(r, c) - this->logG(r, c, featureIt->first) + newGamma;
                                        cout<<newGamma<<" not found"<<endl;
				}
				else{
					T1(r, c) = T1(r, c) - this->logG(r, c, featureIt->first) + itLogGamma->second;
				}

			}

			map<int, float>::iterator itLogGamma;
			itLogGamma = this->gammaLookUp->find(accumT2Key);
			if (itLogGamma == this->gammaLookUp->end()){
				float newGamma = lgammaf(accumT2);
				this->gammaLookUp->insert(std::pair<int, float>(accumT2Key, newGamma));
				T2(r, c) = newGamma;
                                cout<<newGamma<<" not found"<<endl;
			}
			else{
				T2(r, c) = itLogGamma->second;
			}

		}
	}
	
	fmat likelihoodLocation = (T3 - T4 + T1 - T2);
	likelihoodLocation.reshape(CG_ROWS*CG_COLS, 1);
        
	fcolvec posteriorTmp = arma::conv_to<fcolvec>::from(likelihoodLocation);
        end = clock();
        time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
        cout<<"reshaping takes: "<<time_spent<<endl;
	//fmat posterior = reshape(exp(posteriorTmp - posteriorTmp.max() - log(sum(exp(posteriorTmp - posteriorTmp.max())))), CG_ROWS, CG_COLS);
        //return posterior;
        return exp(posteriorTmp - posteriorTmp.max() - log(sum(exp(posteriorTmp - posteriorTmp.max()))));
}

fcolvec CountingGrid::locationPosteriorLoop2(Datapoint* dp)
{
    fmat lg_A = fmat(CG_ROWS, CG_COLS, fill::zeros);
    fmat lg_AN = fmat(CG_ROWS, CG_COLS, fill::zeros);
    fmat slg_a = fmat(CG_ROWS, CG_COLS, fill::zeros);
    fmat slg_an = fmat(CG_ROWS, CG_COLS, fill::zeros);
    
    int r, c, z;
    float Aw_supp;
    float accA, accAN;
    
    frowvec n = dp->getCountsArray();
    
    for (r = 0; r < CG_ROWS; r++)
        for (c = 0; c < CG_COLS; c++)
        {   
            accA = 0;
            accAN = 0;
            
            for (z = 0; z < Z; z++)
            {
                Aw_supp = this->Aw(r,c,z);
                accA += Aw_supp;
                accAN += Aw_supp + n(z);
                
            }
            
            
        }
}

fcolvec CountingGrid::locationPosteriorLoop(Datapoint* dp)
{
        //clock_t begin, end;
        //double time_spent;
        
        //begin = clock();
        
	//fcube tmp = reshape(this->logG, CG_ROWS*CG_COLS, Z, 1);
	//fmat T3 = reshape(sum(tmp.slice(0), 1), CG_ROWS, CG_COLS);
	//T3 = arma::accu(T3) - T3;
        fmat T3 = fmat(CG_ROWS, CG_COLS, fill::zeros);

	fmat T1 = fmat(CG_ROWS, CG_COLS, fill::zeros);
	//fmat T1 = reshape(sum(tmp.slice(0), 1), CG_ROWS, CG_COLS);

	fmat T2 = fmat(CG_ROWS, CG_COLS, fill::zeros);

	fmat T4 = fmat(CG_ROWS, CG_COLS, fill::zeros);
	//T4 = arma::accu(this->logGsum) - this->logGsum;

	int accumT2Key;
        int accumT1Key;
	float accumT2;
        float accumT1;
        float accumT4;
	
        float T3sum = 0;
        float logG_temp = 0;
        float Aw_temp;
        float T1arg;
        frowvec countsArray = dp->getCountsArray();
        
        map<int, float>::iterator itLogGamma;
        
	for (int c = 0; c < CG_COLS; c++)
	{
            for (int r = 0; r < CG_ROWS; r++)
            {
                //accumT2Key = (int) round(sum(fvec(Aw(span(r, r), span(c, c), span::all) - WD_ROWS*WD_COLS*BASE_PRIOR)));
                accumT2 = 0;
                accumT1 = 0;
                
                for (int z = 0; z < Z; z++)
                {
                    logG_temp = this->logG(r,c,z);
                    Aw_temp = this->Aw(r,c,z);
                    
                    T1arg = (Aw_temp + countsArray(z));
                    itLogGamma = this->gammaLookUp->find(T1arg - WD_ROWS*WD_COLS*BASE_PRIOR);
                    
                    if (itLogGamma == this->gammaLookUp->end()){
                            float newGamma = lgammaf(T1arg);
                            this->gammaLookUp->insert(std::pair<int, float>((int) round(accumT1 - WD_ROWS*WD_COLS*BASE_PRIOR), newGamma));
                            accumT1 += newGamma;
                            cout<<T1arg<<" not found"<<endl;
                    }
                    else{
                            accumT1 += itLogGamma->second;
                    }
                    //accumT1 += (Aw_temp + countsArray(z));
                    
                    
                    T3(r,c) += logG_temp;
                    //T3(r,c) = T1(r,c);
                    T3sum += logG_temp; //logGsum and T3sum look like they are the same thing
                    accumT2 += T1arg;
                }
                
                accumT4 += this->logGsum(r,c);
                
                //cout<<"done accumulating"<<endl;
                accumT2Key = (int) round(accumT2 - WD_ROWS*WD_COLS*BASE_PRIOR);
                accumT1Key = (int) round(accumT1 - WD_ROWS*WD_COLS*BASE_PRIOR);
                
                /*
                try{
                    T1(r, c) = this->gammaLookUp->at(accumT1Key);
                }
                catch(const std::out_of_range& oor){
                    float newGamma = lgammaf(accumT1);
                    this->gammaLookUp->insert(std::pair<int, float>(accumT1, newGamma));
                    T1(r, c) = newGamma;
                    cout<<accumT1<<" not found"<<endl;
                }
                //cout<<"T1 assigned"<<endl;
                try{
                    T2(r, c) = this->gammaLookUp->at(accumT2Key);
                }
                catch(const std::out_of_range& oor){
                    float newGamma = lgammaf(accumT2);
                    this->gammaLookUp->insert(std::pair<int, float>(accumT2, newGamma));
                    T2(r, c) = newGamma;
                    cout<<accumT2<<" not found"<<endl;
                }
                //cout<<"T2 assigned"<<endl;
                */
                /*
                for (mapdata::iterator featureIt = tmpMap.begin(); featureIt != tmpMap.end(); featureIt++)
                {

                        accumT2Key += featureIt->first;
                        accumT2 += float(featureIt->first);
                        // sn qui
                        map<int, float>::iterator itLogGamma;

                        int keyTmpCount = (int)round(Aw(r, c, featureIt->first) - WD_ROWS*WD_COLS*BASE_PRIOR + featureIt->second );
                        float trueCount = Aw(r, c, featureIt->first) + float(featureIt->second);
                        itLogGamma = this->gammaLookUp->find(keyTmpCount);
                        if (itLogGamma == this->gammaLookUp->end()){
                                float newGamma = lgammaf(trueCount);
                                this->gammaLookUp->insert(std::pair<int, float>(keyTmpCount, newGamma));
                                T1(r, c) = T1(r, c) - this->logG(r, c, featureIt->first) + newGamma;
                                cout<<keyTmpCount<<" not found"<<endl;
                        }
                        else{
                                T1(r, c) = T1(r, c) - this->logG(r, c, featureIt->first) + itLogGamma->second;
                        }

                }
                */

                /*
                map<int, float>::iterator itLogGamma;
                itLogGamma = this->gammaLookUp->find(accumT1Key);
                if (itLogGamma == this->gammaLookUp->end()){
                        float newGamma = lgammaf(accumT1);
                        this->gammaLookUp->insert(std::pair<int, float>(accumT1Key, newGamma));
                        T1(r, c) = newGamma;
                        cout<<accumT1Key<<" not found"<<endl;
                }
                else{
                        T1(r, c) = itLogGamma->second;
                }
                */
                T1(r,c) = accumT1;
                
                //map<int, float>::iterator itLogGamma;
                itLogGamma = this->gammaLookUp->find(accumT2Key);
                if (itLogGamma == this->gammaLookUp->end()){
                        float newGamma = lgammaf(accumT2);
                        this->gammaLookUp->insert(std::pair<int, float>(accumT2Key, newGamma));
                        T2(r, c) = newGamma;
                        cout<<accumT2Key<<" not found"<<endl;
                }
                else{
                        T2(r, c) = itLogGamma->second;
                }
                
            }
	}
        
        //end = clock();
        //time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
        //cout<<"reshaping and summing takes: "<<time_spent<<endl;
	
        T3 = T3sum - T3;
        T4 = accumT4 - this->logGsum;
        
	fmat likelihoodLocation = (T3 - T4 + T1 - T2);
	likelihoodLocation.reshape(CG_ROWS*CG_COLS, 1);
        
	fcolvec posteriorTmp = arma::conv_to<fcolvec>::from(likelihoodLocation);
        
	//fmat posterior = reshape(exp(posteriorTmp - posteriorTmp.max() - log(sum(exp(posteriorTmp - posteriorTmp.max())))), CG_ROWS, CG_COLS);
        //return posterior;
        return exp(posteriorTmp - posteriorTmp.max() - log(sum(exp(posteriorTmp - posteriorTmp.max()))));
}



double CountingGrid::computeEnergy(Datapoint* dp)
{
	fmat posterior = locationPosterior(dp);
	return posterior(dp->getRow(), dp->getCol());
}


int CountingGrid::computeLogGammaCG(map<int, arma::sp_fmat> tokenLoc)
{
	// togli e aggiorna log G
	for (map<int, arma::sp_fmat>::iterator featureIt = tokenLoc.begin(); featureIt != tokenLoc.end(); featureIt++)
	{
		// Trick per convertire da matrice sparsa. Faccio questo perch� la conversione non funziona.
		fmat tmp = fmat(CG_ROWS, CG_COLS);
		tmp.fill(0);
		tmp += featureIt->second; // tmp � la matrice della slice z
                //tmp.print("tmp");

		uvec nonZeroIds = arma::find(tmp);
		for (uvec::iterator itv = nonZeroIds.begin(); itv != nonZeroIds.end(); itv++)
		{
			// Iteratore su una slice, devo tenere conto che sto analizzando la z-sima feature.
			/* linear index -> Non so se funziona
			int idInCubeCoordinates = *itv + (Z - 1)*featureIt->first; // BUG_AlERT
			*/

			int tmpRow = *itv % CG_ROWS;
			int tmpCol = *itv / CG_ROWS;

			map<int, float>::iterator itLogGamma;

			int key = (int) round( this->Aw(tmpRow, tmpCol, featureIt->first) - BASE_PRIOR*WD_ROWS*WD_COLS);
			itLogGamma = this->gammaLookUp->find(key);
			if (itLogGamma == this->gammaLookUp->end()){
                            cout<<"cannot find "<<key<<endl;
				float newGamma = lgammaf(this->Aw(tmpRow, tmpCol, featureIt->first));
				this->gammaLookUp->insert(std::pair<int, float>(key, newGamma));
				this->logG(tmpRow,tmpCol,featureIt->first) = newGamma;
			}
			else{
				this->logG(tmpRow, tmpCol, featureIt->first) = itLogGamma->second;
			}

			// Vado a modificare la variabile logGsum.
			itLogGamma = this->gammaLookUp->begin(); // Risposto all'inizio l'iteratore
			int keySumA = (int)round(sum(fvec(this->Aw(span(tmpRow, tmpRow), span(tmpCol, tmpCol), span::all)) - BASE_PRIOR*CG_COLS*CG_ROWS));
			int sumA = sum(fvec(this->Aw(span(tmpRow, tmpRow), span(tmpCol, tmpCol), span::all)));

			itLogGamma = this->gammaLookUp->find(keySumA);
			if (itLogGamma == this->gammaLookUp->end()){
				// Calcola la nuova gamma
				float newGamma = lgammaf(sumA);
				// Aggiungila alla lookup table
				this->gammaLookUp->insert(std::pair<int, float>(keySumA, newGamma));
				// Aggiorna la Counting Grid
				this->logGsum(tmpRow, tmpCol) = newGamma;
			}
			else{
				// Aggiorna la Counting Grid
				this->logGsum(tmpRow, tmpCol) = itLogGamma->second;
			}

		}
	}
	return 0;
}


int CountingGrid::printCg(int sl)
{
	(this->a.slice(sl)).print("--> pi: ");
	(this->Aw.slice(sl)).print("--> h: ");
	(this->logG.slice(sl)).print("--> logGamma( h ): ");
	return 0;
}

int CountingGrid::saveCg(string path)
{
	this->a.save(path + "\\pi.mat", arma_ascii);
	return 0;
}


fcube CountingGrid::get_a()
{
	return this->a;
}

fcube CountingGrid::get_Aw()
{
	return this->Aw;
}

fcube CountingGrid::get_logG()
{
	return this->logG;
}

