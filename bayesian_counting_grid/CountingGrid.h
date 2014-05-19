#ifndef COUNTINGGRID_H
#define COUNTINGGRID_H

#include "GeneralHeader.h"
#include "Datapoint.h"

class CountingGrid
{
public:
	// Costruttori
	CountingGrid(map<float, float>* ); // Usa il prior uguale per tutte le locazioni definito in generalHeader
	CountingGrid(map<float, float>*, ucube); // Prende in ingresso un prior predefinito nel main
	~CountingGrid();

	int addDatapoint( Datapoint* );
	int removeDatapoint( Datapoint* );
	int sumAllWindows();
	int computeLogGammaCG();
	int computeLogGammaCG( map<int, arma::sp_fmat> );

	fmat locationPosterior( Datapoint* );
	double computeEnergy( Datapoint* );

	int printCg(int);
	int saveCg(string);

private:
	fcube a;
	fcube Aw;
	fcube logG;
	fmat logGsum;
	map<float, float>* gammaLookUp;
};

#endif