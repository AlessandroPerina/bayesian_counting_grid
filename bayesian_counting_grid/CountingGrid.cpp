#include "CountingGrid.h"
// #include "Datapoint.h"
// #include "GeneralHeader.h"

CountingGrid::CountingGrid(map<float, float>* gl)
{
	this->gammaLookUp = gl;
	this->a = fcube(CG_ROWS, CG_COLS, Z);
	this->a.fill(BASE_PRIOR);

	this->Aw = fcube(CG_ROWS, CG_COLS, Z);
	this->Aw.fill(0);

	// pad a
	fcube tmp = fcube(CG_ROWS+WD_ROWS-1, CG_COLS+WD_COLS-1, Z);
	tmp.fill(0);
	tmp.tube(0, CG_ROWS - 1, 0, CG_COLS - 1) = a;
	tmp.tube(0, CG_ROWS - 1, CG_COLS, CG_COLS + WD_COLS - 2) = a.tube(0, CG_ROWS-1, 0, WD_COLS - 2);
	tmp.tube(CG_ROWS, CG_ROWS + WD_ROWS - 2, 0, CG_COLS - 1) = a.tube(0, WD_ROWS - 2, 0, CG_COLS - 1);
	tmp.tube(CG_ROWS, CG_ROWS + WD_ROWS - 2, CG_COLS, CG_COLS + WD_COLS - 2) = a.tube(0, WD_ROWS - 2, WD_COLS - 2);

	for (int r = 0; r < CG_ROWS; r++)
	{
		for (int c = 0; c < CG_COLS; c++)
		{
			fcube tmp2 = reshape(tmp.tube(r, r + WD_ROWS - 1, c, c + WD_ROWS - 1), WD_ROWS*WD_COLS, Z, 1);
			this->Aw.tube(r, r, c, c) = arma::accu(tmp2.slice(0));
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
				map<float, float>::iterator it;
				it = this->gammaLookUp->find(this->Aw(r, c, z));
				if (it == this->gammaLookUp->end()){
					// Calcola la nuova gamma
					float newGamma = lgammaf(this->Aw(r, c, z));
					// Aggiungila alla lookup table
					this->gammaLookUp->insert(std::pair<float, float>(this->Aw(r, c, z), newGamma));
					// Aggiorna la Counting Grid
					this->logG.at(r, c, z) = newGamma;
				}
				else{
					// Aggiorna la Counting Grid
					this->logG.at(r, c, z) = it->second;
				}
			}

			// Calcolo logGammaSum
			map<float, float>::iterator it;
			float sumA = sum(fvec(this->Aw.tube(r, r, c, c)));
			it = this->gammaLookUp->find(sumA);
			if (it == this->gammaLookUp->end()){
				// Calcola la nuova gamma
				float newGamma = lgammaf(sumA);
				// Aggiungila alla lookup table
				this->gammaLookUp->insert(std::pair<float, float>(sumA, newGamma));
				// Aggiorna la Counting Grid
				this->logGsum.at(r, c ) = newGamma;
			}
			else{
				// Aggiorna la Counting Grid
				this->logGsum.at(r, c) = it->second;
			}

		}
	}


}

CountingGrid::CountingGrid( map<float,float>* gl, ucube prior)
{
	this->gammaLookUp = gl;
	this->a = arma::conv_to<fcube>::from(prior);

	this->Aw = fcube(CG_ROWS, CG_COLS, Z);
	this->Aw.fill(0);

	// pad a
	fcube tmp = fcube(CG_ROWS + WD_ROWS - 1, CG_COLS + WD_COLS - 1, Z);
	tmp.fill(0);
	tmp.tube(0, CG_ROWS - 1, 0, CG_COLS - 1) = a;
	tmp.tube(0, CG_ROWS - 1, CG_COLS, CG_COLS + WD_COLS - 2) = a.tube(0, CG_ROWS - 1, 0, WD_COLS - 2);
	tmp.tube(CG_ROWS, CG_ROWS + WD_ROWS - 2, 0, CG_COLS - 1) = a.tube(0, WD_ROWS - 2, 0, CG_COLS - 1);
	tmp.tube(CG_ROWS, CG_ROWS + WD_ROWS - 2, CG_COLS, CG_COLS + WD_COLS - 2) = a.tube(0, WD_ROWS - 2, WD_COLS - 2);

	for (int r = 0; r < CG_ROWS; r++)
	{
		for (int c = 0; c < CG_COLS; c++)
		{
			fcube tmp2 = reshape(tmp.tube(r, r + WD_ROWS - 1, c, c + WD_ROWS - 1), WD_ROWS*WD_COLS, Z, 1);
			this->Aw.tube(r, r, c, c) = arma::accu(tmp2.slice(0));
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
				map<float, float>::iterator it;
				it = this->gammaLookUp->find(this->Aw(r, c, z));
				if (it == this->gammaLookUp->end()){
					// Calcola la nuova gamma
					float newGamma = lgammaf(this->Aw(r, c, z));
					// Aggiungila alla lookup table
					this->gammaLookUp->insert(std::pair<float, float>(this->Aw(r, c, z), newGamma));
					// Aggiorna la Counting Grid
					this->logG.at(r, c, z) = newGamma;
				}
				else{
					// Aggiorna la Counting Grid
					this->logG.at(r, c, z) = it->second;
				}
			}

			// Calcolo logGammaSum
			map<float, float>::iterator it;
			float sumA = sum(fvec(this->Aw.tube(r, r, c, c)));
			it = this->gammaLookUp->find(sumA);
			if (it == this->gammaLookUp->end()){
				// Calcola la nuova gamma
				float newGamma = lgammaf(sumA);
				// Aggiungila alla lookup table
				this->gammaLookUp->insert(std::pair<float, float>(sumA, newGamma));
				// Aggiorna la Counting Grid
				this->logGsum.at(r, c) = newGamma;
			}
			else{
				// Aggiorna la Counting Grid
				this->logGsum.at(r, c) = it->second;
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
	padded.tube(0, CG_ROWS - 1, 0, CG_COLS - 1) = a;
	padded.tube(0, CG_ROWS - 1, CG_COLS, CG_COLS + WD_COLS - 1) = a.tube(0, CG_ROWS - 1, 0, WD_COLS - 2);
	padded.tube(CG_ROWS, CG_ROWS + WD_ROWS - 1, 0, CG_COLS - 1) = a.tube(0, WD_ROWS - 2, 0, CG_COLS - 1);
	padded.tube(CG_ROWS, CG_ROWS + WD_ROWS - 1, CG_COLS, CG_COLS + WD_COLS - 1) = a.tube(0, WD_ROWS - 2, WD_COLS - 2);

	for (int r = 0; r < CG_ROWS; r++)
	{
		for (int c = 0; c < CG_COLS; c++)
		{
			fcube tmp2 = reshape(padded.tube(r, r + WD_ROWS - 1, c, c + WD_ROWS - 1), WD_ROWS*WD_COLS, Z, 1);
			this->Aw.tube(r, r, c, c) = sum(tmp2.slice(0), 0);
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


fmat CountingGrid::locationPosterior(Datapoint* dp)
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
	float accumT2 = 0;
	// itera solo sulle features che ha il campione corrente
	for (int r = 0; r < CG_ROWS; r++)
	{
		for (int c = 0; c < CG_COLS; c++)
		{

			accumT2 = sum(fvec(Aw.tube(r, r, c, c)));
			for (urowvec::iterator featureIt = dp->getWords().begin(); featureIt != dp->getWords().end(); featureIt++)
			{		
				accumT2 += float(dp->getSingleCountsDict(*featureIt));

				map<float, float>::iterator itLogGamma;
				float tmpCount = Aw[r, c, *featureIt] + float(dp->getSingleCountsDict(*featureIt));
				itLogGamma = this->gammaLookUp->find(tmpCount);
				if (itLogGamma == this->gammaLookUp->end()){
					float newGamma = lgammaf(tmpCount);
					this->gammaLookUp->insert(std::pair<float, float>(tmpCount, newGamma));
					T1(r, c) = T1(r, c) - this->logG(r, c, *featureIt) + newGamma;
				}
				else{
					T1(r, c) = T1(r, c) - this->logG(r, c, *featureIt) + itLogGamma->second;
				}

			}

			map<float, float>::iterator itLogGamma;
			itLogGamma = this->gammaLookUp->find(accumT2);
			if (itLogGamma == this->gammaLookUp->end()){
				float newGamma = lgammaf(accumT2);
				this->gammaLookUp->insert(std::pair<float, float>(accumT2, newGamma));
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
	fmat posterior = reshape(exp(posteriorTmp - posteriorTmp.max() - log(sum(exp(posteriorTmp - posteriorTmp.max())))), CG_ROWS, CG_COLS);
	return posterior;
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
		// Trick per convertire da matrice sparsa. Faccio questo perchè la conversione non funziona.
		fmat tmp;
		tmp.fill(0);
		tmp += featureIt->second; // tmp è la matrice della slice z

		uvec nonZeroIds = arma::find(tmp, 0);
		for (uvec::iterator itv = nonZeroIds.begin(); itv != nonZeroIds.end(); itv++)
		{
			// Iteratore su una slice, devo tenere conto che sto analizzando la z-sima feature.
			/* linear index -> Non so se funziona
			int idInCubeCoordinates = *itv + (Z - 1)*featureIt->first; // BUG_AlERT
			*/

			int tmpRow = *itv % CG_ROWS;
			int tmpCol = (int)floor(*itv / CG_ROWS) % CG_COLS;


			map<float, float>::iterator itLogGamma;
			itLogGamma = this->gammaLookUp->find(Aw(tmpRow, tmpCol, featureIt->first));
			if (itLogGamma == this->gammaLookUp->end()){
				float newGamma = lgammaf(Aw(tmpRow, tmpCol, featureIt->first));
				this->gammaLookUp->insert(std::pair<float, float>(Aw(tmpRow, tmpCol, featureIt->first), newGamma));
				this->logG.at(tmpRow,tmpCol,featureIt->first) = newGamma;
			}
			else{
				this->logG.at(tmpRow, tmpCol, featureIt->first) = itLogGamma->second;
			}

			// Vado a modificare la variabile logGsum.
			itLogGamma = this->gammaLookUp->begin(); // Risposto all'inizio l'iteratore
			float sumA = sum(fvec(Aw.tube(tmpRow, tmpRow, tmpCol, tmpCol)));
			itLogGamma = this->gammaLookUp->find(sumA);
			if (itLogGamma == this->gammaLookUp->end()){
				// Calcola la nuova gamma
				float newGamma = lgammaf(sumA);
				// Aggiungila alla lookup table
				this->gammaLookUp->insert(std::pair<float, float>(sumA, newGamma));
				// Aggiorna la Counting Grid
				this->logGsum.at(tmpRow, tmpCol) = newGamma;
			}
			else{
				// Aggiorna la Counting Grid
				this->logGsum.at(tmpRow, tmpCol) = itLogGamma->second;
			}

		}
	}
	return 0;
}

int CountingGrid::computeLogGammaCG()
{
	for (int r = 0; r < CG_ROWS; r++)
	{
		for (int c = 0; c < CG_COLS; c++)
		{
			for (int z = 0; z < Z; z++)
			{
				map<float, float>::iterator it;
				it = this->gammaLookUp->find(this->Aw(r, c, z));
				if (it == this->gammaLookUp->end()){

					// Calcola la nuova gamma
					float newGamma = lgammaf(this->Aw(r, c, z));
					// Aggiungila alla lookup table
					this->gammaLookUp->insert(std::pair<float, float>(this->Aw(r, c, z), newGamma));
					// Aggiorna la Counting Grid
					this->logG.at(r, c, z) = newGamma;
				}
				else{
					// Aggiorna la Counting Grid
					this->logG.at(r, c, z) = it->second;
				}
			}

			// Calcolo logGammaSum
			map<float, float>::iterator it;
			float sumA = sum(fvec(this->Aw.tube(r, r, c, c)));
			it = this->gammaLookUp->find(sumA);
			if (it == this->gammaLookUp->end()){
				// Calcola la nuova gamma
				float newGamma = lgammaf(sumA);
				// Aggiungila alla lookup table
				this->gammaLookUp->insert(std::pair<float, float>(sumA, newGamma));
				// Aggiorna la Counting Grid
				this->logGsum.at(r, c) = newGamma;
			}
			else{
				// Aggiorna la Counting Grid
				this->logGsum.at(r, c) = it->second;
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


