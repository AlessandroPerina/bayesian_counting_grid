#ifndef COUNTINGGRID_H
#define COUNTINGGRID_H

#include "GeneralHeader.h" // Da fare

class CountingGrid
{
public:
	// Costruttori
	CountingGrid(); // Usa il prior uguale per tutte le locazioni definito in generalHeader
	CountingGrid( ucube ); // Prende in ingresso un prior predefinito nel main
	~CountingGrid();

	int addDatapoint( Datapoint* );
	int removeDatapoint( Datapoint* );
	int sumAllWindows();
	int computeLogGammaCG();
	int computeLogGammaCG( map<int, arma::sp_fmat> );

	fmat locationPosterior( Datapoint* );
	double computeEnergy( Datapoint* );

private:
	fcube a;
	fcube Aw;
	fcube logG;
	fmat logGsum;
};

#endif