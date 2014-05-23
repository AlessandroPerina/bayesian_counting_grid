#include "Datapoint.h"
#include "CountingGrid.h"

Datapoint::Datapoint( mapdata counts )
{   
	this->countsDict = counts;
	this->pointName = "no_name";
	this->isAssigned = -1;
	this->rowMap = -1;
	this->colMap = -1;
	this->nUniWords = this->countsDict.size();
	// frowvec tmp = frowvec(this->nonZero); Controlla se funziona
	this->countsArray  = zeros<frowvec>(  Z );
	this->words = zeros<urowvec>( this->nUniWords );

	int k=0;
        
	for (mapdata::iterator it = this->countsDict.begin(); it != this->countsDict.end(); it++ )
	{
		this->countsArray[it->first] = float( it->second );
		this->words[k] = it->first;
		this->tokenLoc.insert(std::pair<int, arma::sp_fmat>(it->first, arma::sp_fmat(CG_ROWS, CG_COLS)));
		k++;
	}
        
	this->nTokens = int(arma::sum( this->countsArray ));
        
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
        this->isAssigned = 1;
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

int Datapoint::checkAsgn(){
    return this->isAssigned;
}

int Datapoint::sampleLocation(fcolvec p, boost::mt19937* rng){
    
    int asgnIdx;
    
    std::vector<double> z = conv_to< std::vector<double> >::from(p);
    boost::random::discrete_distribution<> dist(z.begin(),z.end());    
    asgnIdx = dist(*rng);
    // cout<<asgnIdx<<endl;
    //p.print("p");
    //cout<<asgnIdx<<endl;
    //cout<<(int) asgnIdx % CG_ROWS<<" -"<<(int) asgnIdx / CG_ROWS<<endl;
    //Convert linear index in valid location
    this->setLocation((int) asgnIdx % CG_ROWS, (int) asgnIdx / CG_ROWS);
    //cout<<asgnIdx<<" r: "<<(int) asgnIdx % CG_ROWS<<" c: "<<(int) asgnIdx / CG_ROWS<<endl;
    //this->colMap = (int)asgnIdx / CG_ROWS;
    //this->rowMap = (int)asgnIdx % CG_ROWS;
    
    return 0;
}

int Datapoint::sampleTokenLocation(CountingGrid* cg, boost::mt19937* rng)
{
        int rowTok, colTok, rowFinal, colFinal;
        int asgnIdx;
        fcolvec aTmp;
        sp_fmat tmp;
        
	if (ASSIGN_TOKEN)
	{
		fcube suba = cg->get_a().tube(span(this->rowMap, this->rowMap + WD_ROWS - 1), span(this->colMap, this->colMap + WD_COLS - 1));
		for (mapdata::iterator itFeature = this->countsDict.begin(); itFeature != this->countsDict.end(); itFeature++)
		{
			aTmp = reshape(suba.slice(itFeature->first), WD_ROWS*WD_COLS, 1);
			aTmp = aTmp / sum(aTmp);
                        //aTmp.t().print("aTmp: ");
                        //cout<<"Sum: "<<arma::sum(aTmp)<<endl;

			std::vector<double> z = conv_to< std::vector<double> >::from(aTmp);
			boost::random::discrete_distribution<> dist(z.begin(), z.end());
			this->tokenLoc[itFeature->first].zeros(CG_ROWS, CG_COLS);
			for (int i = 0; i < itFeature->second; i++)
			{
				asgnIdx = dist(*rng);
				//cout<<"asgnIdx: "<<asgnIdx<<endl;
                                rowTok = (int)asgnIdx % WD_ROWS;
				colTok = (int)asgnIdx / WD_ROWS;
                                
                                rowFinal = (int)(rowTok+this->rowMap) % CG_ROWS;
                                colFinal = (int)(colTok+this->colMap) % CG_COLS;
                                
                                //cout<<rowTok<<" "<<colTok<<endl;
                                //cout<<itFeature->first<<endl;
                                //this->tokenLoc[itFeature->first].print();
                                
                                //cout<<rowTok+this->rowMap<<" "<<colTok+this->colMap<<endl;
				this->tokenLoc[itFeature->first](rowFinal, colFinal) += 1;
			}
		}
	}
	else
	{
		fcube suba = cg->get_a().tube(span(this->rowMap, this->rowMap + WD_ROWS - 1), span(this->colMap, this->colMap + WD_COLS - 1));
		for (mapdata::iterator itFeature = this->countsDict.begin(); itFeature != this->countsDict.end(); itFeature++)
		{
			aTmp = reshape(suba.slice(itFeature->first), WD_ROWS*WD_COLS, 1);
			aTmp = aTmp / sum(aTmp);

			
			std::vector<double> z = conv_to< std::vector<double> >::from(aTmp);
			boost::random::discrete_distribution<> dist(z.begin(), z.end());
			this->tokenLoc[itFeature->first].clear();

			asgnIdx = dist(*rng);
			rowTok = (int)asgnIdx / WD_ROWS;
			colTok = (int)asgnIdx % WD_COLS;
                        
                        rowFinal = (int)(rowTok+this->rowMap) / CG_ROWS;
                        colFinal = (int)(colTok+this->colMap) % CG_ROWS;
                        
			//this->rowMap = (int)asgnIdx % CG_ROWS;
			this->tokenLoc[itFeature->first](rowFinal, colFinal) = itFeature->second;

		}
	}
	return 0;
}