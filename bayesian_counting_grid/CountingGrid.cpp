#include "CountingGrid.h"
// #include "Datapoint.h"
// #include "GeneralHeader.h"

CountingGrid::CountingGrid(map<int, float>* gl)
{
	this->gammaLookUp = gl;

	this->a = fcube(CG_ROWS, CG_COLS, Z);
	this->a = floor( 5 * arma::randu<fcube>(CG_ROWS, CG_COLS, Z)) + BASE_PRIOR;
	// this->a.fill(BASE_PRIOR);

	this->Aw = fcube(CG_ROWS, CG_COLS, Z);
	this->Aw.fill(0);

	this->logG = fcube(CG_ROWS, CG_COLS, Z);
	this->logG.fill(0);

	this->logGsum = fmat(CG_ROWS, CG_COLS);
	this->logGsum.fill(0);

	// pad a
	fcube tmp = fcube(CG_ROWS+WD_ROWS-1, CG_COLS+WD_COLS-1, Z);
	tmp.fill(1);
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


int CountingGrid::addDatapoint(Datapoint* dp )
{
	for (map<int, arma::sp_fmat>::iterator it = dp->getTokenLoc().begin(); it != dp->getTokenLoc().end(); it++)
	{
		this->a.slice(it->first) += it->second;
	}
	// Aggiorno la variabile somma
	this->sumAllWindows();
	this->computeLogGammaCG(dp->getTokenLoc());
	
	return 0;

}

int CountingGrid::removeDatapoint( Datapoint* dp ) 
{
	for (map<int, arma::sp_fmat>::iterator it = dp->getTokenLoc().begin(); it != dp->getTokenLoc().end(); it++)
	{
		this->a.slice(it->first) -= it->second;
	}
	// Aggiorno la variabile somma
	this->sumAllWindows();
	this->computeLogGammaCG(dp->getTokenLoc());

	return 0;
}


fcolvec CountingGrid::locationPosterior(Datapoint* dp)
{
	fmat locationLikelihood = fmat(CG_ROWS, CG_COLS, fill::zeros);

	fmat T1 = fmat(CG_ROWS, CG_COLS, fill::zeros);
	T1 = this->logGsum;
	fmat T2 = fmat(CG_ROWS, CG_COLS, fill::zeros);

	fmat T4 = arma::accu(this->logGsum) - this->logGsum;

	fcube tmp = reshape(this->logG, CG_ROWS*CG_COLS, Z, 1);
	fmat T3 = reshape(sum( tmp.slice(0),1), CG_ROWS, CG_COLS);
	T3 = arma::accu(T3) - T3;

	//	fmat(CG_COLS, CG_ROWS, fill::zeros);
	int accumT2Key = 0;
	float accumT2 = 0;
	// itera solo sulle features che ha il campione corrente
	for (int r = 0; r < CG_ROWS; r++)
	{
		for (int c = 0; c < CG_COLS; c++)
		{

			accumT2Key = (int)round(sum(fvec(Aw(span(r, r), span(c, c), span::all) - WD_ROWS*WD_COLS*BASE_PRIOR)));
			accumT2 = sum(fvec(Aw(span(r, r), span(c, c), span::all) ));
			for (urowvec::iterator featureIt = dp->getWords().begin(); featureIt != dp->getWords().end(); featureIt++)
			{		
				accumT2Key += dp->getSingleCountsDict(*featureIt);
				accumT2 += float(dp->getSingleCountsDict(*featureIt));
				// sn qui
				map<int, float>::iterator itLogGamma;
				int keyTmpCount = (int)round(Aw(r, c, *featureIt) - WD_ROWS*WD_COLS*BASE_PRIOR + dp->getSingleCountsDict(*featureIt));
				float trueCount = Aw(r, c, *featureIt) + float(dp->getSingleCountsDict(*featureIt));
				itLogGamma = this->gammaLookUp->find(keyTmpCount);
				if (itLogGamma == this->gammaLookUp->end()){
					float newGamma = lgammaf(trueCount);
					this->gammaLookUp->insert(std::pair<float, float>(keyTmpCount, newGamma));
					T1(r, c) = T1(r, c) - this->logG(r, c, *featureIt) + newGamma;
				}
				else{
					T1(r, c) = T1(r, c) - this->logG(r, c, *featureIt) + itLogGamma->second;
				}

			}

			map<int, float>::iterator itLogGamma;
			itLogGamma = this->gammaLookUp->find(accumT2Key);
			if (itLogGamma == this->gammaLookUp->end()){
				float newGamma = lgammaf(accumT2);
				this->gammaLookUp->insert(std::pair<int, float>(accumT2Key, newGamma));
				T2(r, c) = newGamma;
			}
			else{
				T2(r, c) = itLogGamma->second;
			}

		}
	}
	
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
		fmat tmp;
		tmp.fill(0);
		tmp += featureIt->second; // tmp � la matrice della slice z

		uvec nonZeroIds = arma::find(tmp, 0);
		for (uvec::iterator itv = nonZeroIds.begin(); itv != nonZeroIds.end(); itv++)
		{
			// Iteratore su una slice, devo tenere conto che sto analizzando la z-sima feature.
			/* linear index -> Non so se funziona
			int idInCubeCoordinates = *itv + (Z - 1)*featureIt->first; // BUG_AlERT
			*/

			int tmpRow = *itv % CG_ROWS;
			int tmpCol = (int)floor(*itv / CG_ROWS) % CG_COLS;

			map<int, float>::iterator itLogGamma;
			int key = (int) this->Aw(tmpRow, tmpCol, featureIt->first) - BASE_PRIOR;
			itLogGamma = this->gammaLookUp->find(key);
			if (itLogGamma == this->gammaLookUp->end()){
				float newGamma = lgammaf(this->Aw(tmpRow, tmpCol, featureIt->first));
				this->gammaLookUp->insert(std::pair<int, float>(key, newGamma));
				this->logG(tmpRow,tmpCol,featureIt->first) = newGamma;
			}
			else{
				this->logG(tmpRow, tmpCol, featureIt->first) = itLogGamma->second;
			}

			// Vado a modificare la variabile logGsum.
			itLogGamma = this->gammaLookUp->begin(); // Risposto all'inizio l'iteratore
			int keySumA = (int)sum(fvec(this->Aw(span(tmpRow, tmpRow), span(tmpCol, tmpCol), span::all)) - BASE_PRIOR);
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

int CountingGrid::computeLogGammaCG()
{
	for (int r = 0; r < CG_ROWS; r++)
	{
		map<int, float>::iterator itLogGamma;
		for (int c = 0; c < CG_COLS; c++)
		{
			for (int z = 0; z < Z; z++)
			{
				itLogGamma = this->gammaLookUp->begin();
				int key = (int) this->Aw(r, c, z) - BASE_PRIOR;
				itLogGamma = this->gammaLookUp->find(key);
				if (itLogGamma == this->gammaLookUp->end()){

					// Calcola la nuova gamma
					float newGamma = lgammaf(this->Aw(r, c, z));
					// Aggiungila alla lookup table
					this->gammaLookUp->insert(std::pair<float, float>(key, newGamma));
					// Aggiorna la Counting Grid
					this->logG(r, c, z) = newGamma;
				}
				else{
					// Aggiorna la Counting Grid
					this->logG(r, c, z) = itLogGamma->second;
				}
			}

			// Calcolo logGammaSum
			itLogGamma = this->gammaLookUp->begin();
			int keySumA = (int) sum(fvec(this->Aw(span(r, r), span(c, c), span::all)) - BASE_PRIOR );
			float sumA = sum(fvec(this->Aw(span(r, r), span(c, c), span::all)));

			itLogGamma = this->gammaLookUp->find(keySumA);
			if (itLogGamma == this->gammaLookUp->end()){
				// Calcola la nuova gamma
				float newGamma = lgammaf(sumA);
				// Aggiungila alla lookup table
				this->gammaLookUp->insert(std::pair<int, float>(keySumA, newGamma));
				// Aggiorna la Counting Grid
				this->logGsum(r, c) = newGamma;
			}
			else{
				// Aggiorna la Counting Grid
				this->logGsum(r, c) = itLogGamma->second;
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
	this->a.save(path + "\pi.mat", arma_ascii);
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

