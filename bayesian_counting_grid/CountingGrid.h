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
        int sumAllWindowsLoop();
        int updateAw(Datapoint*);
	int computeLogGammaCG();
	int computeLogGammaCG( map<int, arma::sp_fmat> );

	//fmat locationPosterior( Datapoint* );
        fcolvec locationPosterior( Datapoint* );
        fcolvec locationPosteriorLoop( Datapoint* );
        fcolvec locationPosteriorLoopFast( Datapoint* );
        fcolvec locationPosteriorLoopPar( Datapoint* );
        std::map<int, arma::sp_fmat> findUpdateLoc(map<int, arma::sp_fmat>);
	double computeEnergy( Datapoint* );

	int printCg(int);
	int saveCg(string);

	fcube get_a();
	fcube get_Aw();
	fcube get_logG();

private:
	fcube a;
	fcube Aw;
	fmat Aw_sum;
        fcube logG;
	fmat logGsum;
        fmat sum_logG;
	map<int, float>* gammaLookUp;
};

#endif