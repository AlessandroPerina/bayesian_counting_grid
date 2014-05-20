#include "Datapoint.h"


Datapoint::Datapoint( mapdata counts )
{   
	this->countsDict = counts;
	this->pointName = "no_name";
	this->isAssigned = -1;
	this->rowMap = -1;
	this->colMap = -1;
	this->nUniWords = this->countsDict.size();
        //cout<<counts.size()<<" - "<<this->nUniWords<<endl;
	// frowvec tmp = frowvec(this->nonZero); Controlla se funziona
	this->countsArray  = zeros<frowvec>(  Z );
	this->words = zeros<urowvec>( this->nUniWords );

	int k=0;
        cout<<"starting inner loop"<<endl; 
        
	for (mapdata::iterator it = this->countsDict.begin(); it != this->countsDict.end(); it++ )
	{
		this->countsArray[it->first] = float( it->second );
		this->words[k] = it->first;
		this->tokenLoc.insert(std::pair<int, arma::sp_fmat>(it->first, arma::sp_fmat(CG_ROWS, CG_COLS)));
		k++;
	}
        cout<<"now computing sum"<<endl;
	this->nTokens = int(arma::sum( this->countsArray ));
        cout<<"Datapoint is initialised"<<endl;
}


Datapoint::Datapoint(mapdata counts, string name)
{
	this->countsDict = counts;
	this->pointName = name;
	this->isAssigned = -1;
	this->rowMap = -1;
	this->colMap = -1;
	this->nUniWords = this->countsDict.size();

	// frowvec tmp = frowvec(this->nonZero); Controlla se funziona
	this->countsArray  = zeros<frowvec>(  Z );
	this->words = zeros<urowvec>( this->nUniWords );

	int k=0;
	for ( mapdata::iterator it = this->countsDict.begin(); it != this->countsDict.end(); it ++ )
	{
		this->countsArray[it->first] = float( it->second );
		this->words[k] = it->first;
		this->tokenLoc.insert(std::pair<int, arma::sp_fmat>(int(it->first), arma::sp_fmat(CG_ROWS, CG_COLS)));
		k++;
	}
	this->nTokens = int(arma::sum( this->countsArray ));
}

Datapoint::~Datapoint()
{}


void Datapoint::setLocation( int row, int col )
{
	this->rowMap = row;
	this->colMap = col;
}


void Datapoint::setTokenLocations( map<int, arma::sp_fmat> maptoken)
{
	for (map<int,arma::sp_fmat>::iterator it = maptoken.begin(); it != maptoken.end(); it++)
	{
		this->tokenLoc.insert(*it); // Mistero della fede... 
	}
}

string Datapoint::getName()
{
	return this->pointName;
}
 
int Datapoint::getnUniWords()
{
	return this->nUniWords;
}

int Datapoint::getRow()
{
	return this->rowMap;
}

int Datapoint::getCol()
{
	return this->colMap;
}

mapdata Datapoint::getCountsDict()
{
	return this->countsDict;
}

frowvec Datapoint::getCountsArray()
{
	return this->countsArray;
}

map<int, sp_fmat> Datapoint::getTokenLoc()
{
	return this->tokenLoc;
}

urowvec Datapoint::getWords()
{
	return this->words;
}

int Datapoint::getSingleCountsDict(int id)
{
	return this->countsDict[id];
}

