#include "CountingGrid.h"
#include "Datapoint.h"

CountingGrid::CountingGrid()
{
	this->a = fcube(CG_ROWS, CG_COLS, Z);
	this->a.fill(BASE_PRIOR);

	this->Aw = fcube(CG_ROWS, CG_COLS, Z);
	this->Aw.fill(0);

	// pad a
	fcube tmp = fcube(CG_ROWS+WD_ROWS-1, CG_COLS+WD_COLS-1, Z);
	tmp.fill(0);
	tmp.tube(0, CG_ROWS - 1, 0, CG_COLS - 1) = a;
	tmp.tube(0, CG_ROWS - 1, CG_COLS, CG_COLS + WD_COLS - 1) = a.tube(0, CG_ROWS-1, 0, WD_COLS - 2);
	tmp.tube(CG_ROWS, CG_ROWS + WD_ROWS - 1, 0, CG_COLS - 1) = a.tube(0, WD_ROWS - 2, 0, CG_COLS - 1);
	tmp.tube(CG_ROWS, CG_ROWS + WD_ROWS - 1, CG_COLS, CG_COLS + WD_COLS - 1) = a.tube(0, WD_ROWS - 2, WD_COLS - 2);

	for (int r = 0; r < CG_ROWS; r++)
	{
		for (int c = 0; c < CG_COLS; c++)
		{		
			this->Aw[r, c] = arma::accu( tmp.tube(r, r + WD_ROWS - 1, c, c + WD_ROWS - 1));
		}
	}

	// Inefficiente ma va fatto una volta sola
	for (int r = 0; r < CG_ROWS; r++)
	{
		for (int c = 0; c < CG_COLS; c++)
		{
			for (int z = 0; z < Z; z++)
			{
				map<int, float>::iterator it;
				int deprioredValue = int(Aw[r, c, z] - BASE_PRIOR);
				it = gammaLookUp.find(deprioredValue);
				if (it == gammaLookUp.end()){
					gammaLookUp.insert(std::pair<int, float>(deprioredValue, lgammaf(deprioredValue)));
				}
				else{
					this->logG.at(r, c, z) = it->second;
				}
			}

		}
	}


}

CountingGrid::CountingGrid( ucube prior)
{
	this->a = arma::conv_to<fcube>::from(prior);

	this->Aw = fcube(CG_ROWS, CG_COLS, Z);
	this->Aw.fill(0);

	// pad a
	fcube tmp = fcube(CG_ROWS + WD_ROWS - 1, CG_COLS + WD_COLS - 1, Z);
	tmp.tube(0, CG_ROWS - 1, 0, CG_COLS - 1) = a;
	tmp.tube(0, CG_ROWS - 1, CG_COLS, CG_COLS + WD_COLS - 1) = a.tube(0, CG_ROWS - 1, 0, WD_COLS - 2);
	tmp.tube(CG_ROWS, CG_ROWS + WD_ROWS - 1, 0, CG_COLS - 1) = a.tube(0, WD_ROWS - 2, 0, CG_COLS - 1);
	tmp.tube(CG_ROWS, CG_ROWS + WD_ROWS - 1, CG_COLS, CG_COLS + WD_COLS - 1) = a.tube(0, WD_ROWS - 2, WD_COLS - 2);

	for (int r = 0; r < CG_ROWS; r++)
	{
		for (int c = 0; c < CG_COLS; c++)
		{
			this->Aw[r, c] = arma::accu(tmp.tube(r, r + WD_ROWS - 1, c, c + WD_ROWS - 1));
		}
	}

	// Inefficiente ma va fatto una volta sola
	for (int r = 0; r < CG_ROWS; r++)
	{
		for (int c = 0; c < CG_COLS; c++)
		{
			for (int z = 0; z < Z; z++)
			{
				map<int, float>::iterator it;
				int deprioredValue = int(Aw[r, c, z] - BASE_PRIOR);
				it = gammaLookUp.find( deprioredValue );
				if (it == gammaLookUp.end()){
					gammaLookUp.insert(std::pair<int, float>(deprioredValue, lgammaf( deprioredValue ) ));
				}
				else{
					this->logG.at(r, c, z) = it->second;
				}
			}

		}
	}

}

CountingGrid::~CountingGrid()
{
}

fcube CountingGrid::sumAllWindows()
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
			this->Aw[r, c] = arma::accu(padded.tube(r, r + WD_ROWS - 1, c, c + WD_ROWS - 1));
		}
	}
	return Aw;
}


void CountingGrid::addDatapoint(Datapoint* dp )
{
	for (map<int, arma::sp_fmat>::iterator it = dp->getTokenLoc().begin(); it != dp->getTokenLoc().end(); it++)
	{
		this->a.slice(it->first) += it->second;
	}
	// Aggiorno la variabile somma
	this->sumAllWindows();

}

void CountingGrid::removeDatapoint( Datapoint* dp ) 
{
	for (map<int, arma::sp_fmat>::iterator it = dp->getTokenLoc().begin(); it != dp->getTokenLoc().end(); it++)
	{
		this->a.slice(it->first) -= it->second;
	}
	// Aggiorno la variabile somma
	this->sumAllWindows();

	// togli e aggiorna log G
	for (map<int, arma::sp_fmat>::iterator it = dp->getTokenLoc().begin(); it != dp->getTokenLoc().end(); it++)
	{	
		uvec nonZero = find( it->second, 0 );
		for ( uvec::iterator itv = nonZero.begin(); itv != itv.end(); it++ )
		{
			log[2]
		}
	}


}


double CountingGrid::locationPosterior(Datapoint* dp)
{
	for (int r = 0; r < CG_ROWS; r++)
	{
		for (int c = 0; c < length; i++)
		{
			// aggiungi cz
			for (map<int, arma::sp_fmat>::iterator it = dp->getTokenLoc().begin(); it != dp->getTokenLoc().end(); it++)
			{
				T1 = reshape(cube, n_rows, n_cols, n_slices)
			}
		}

	}

}

double CountingGrid::computeEnergy(Datapoint* dp)
{

}