#ifndef DATAPOINT_H
#define DATAPOINT_H

#include "GeneralHeader.h" // Da fare

class Datapoint {

public:
	// Costruttori
	Datapoint( mapdata );
	Datapoint( mapdata, string );
	virtual ~Datapoint();
	void setLocation( int, int );
	void setTokenLocations( map<int, sp_fmat > );
	string getName();
	int getnUniWords();
	int getRow();
	int getCol();
	urowvec getWords();
	mapdata getCountsDict();
	int getSingleCountsDict(int);
	frowvec getCountsArray();
	map<int, sp_fmat> getTokenLoc();
        int sampleLocation(fcolvec, boost::mt19937*);
        int checkAsgn();
        
private:
	int isAssigned;
	string pointName;
	int rowMap;
	int colMap;
	int nUniWords; // numero di parole presenti nel documento
	int nTokens;
	urowvec words;
	frowvec countsArray; 
	map<int, sp_fmat> tokenLoc;
	mapdata countsDict;

};
#endif