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

	void addDatapoint( Datapoint* );
	void removeDatapoint( Datapoint* );
	fcube sumAllWindows();

	double locationPosterior( Datapoint* );
	double computeEnergy( Datapoint* );

private:
	fcube a;
	fcube Aw;
	fcube logG;

};

#endif