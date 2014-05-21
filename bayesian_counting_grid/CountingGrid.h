#ifndef COUNTINGGRID_H
#define COUNTINGGRID_H

#include "GeneralHeader.h"
#include "Datapoint.h"

class CountingGrid
{
public:
	// Costruttori
	CountingGrid(map<int, float>* ); // Usa il prior uguale per tutte le locazioni definito in generalHeader
	// CountingGrid(map<int, float>*, ucube); // Prende in ingresso un prior predefinito nel main
	~CountingGrid();

	int addDatapoint( Datapoint* );
	int removeDatapoint( Datapoint* );
	int sumAllWindows();
	int computeLogGammaCG();
	int computeLogGammaCG( map<int, arma::sp_fmat> );

	//fmat locationPosterior( Datapoint* );
        fcolvec locationPosterior( Datapoint* );
	double computeEnergy( Datapoint* );

	int printCg(int);
	int saveCg(string);

	fcube get_a();
	fcube get_Aw();
	fcube get_logG();

private:
	fcube a;
	fcube Aw;
	fcube logG;
	fmat logGsum;
	map<int, float>* gammaLookUp;
};

#endif