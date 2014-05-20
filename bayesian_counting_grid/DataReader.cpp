#include "DataReader.h"
#include "Datapoint.h"

DataReader::DataReader(string data_filename) {
    
    this->data_filename = data_filename;
}

DataReader::~DataReader() {
}

int DataReader::loadData(){
    
    //FILE *fileptr;
    int length, count, word, n;
    int pointCounter;
    string dataLine;
    
    mapdata tmpMap;
    //string pointName;
    ostringstream convert;
    
    cout<<"reading data from"<<this->data_filename;
    
    ifstream dataFile ((char*)this->data_filename.c_str());
    //fileptr = fopen(data_filename, "r");
    
    pointCounter = 0;
    //while ((fscanf(fileptr, "%10d", &length) != EOF))
    while ( getline (dataFile, dataLine) )
    {
        for (n = 0; n < length; n++)
        {
            sscanf((char*)dataLine.c_str(), "%10d:%10d", &word, &count);
            
            tmpMap.insert(std::pair<int,int>(word, count));
        }
        
        convert << pointCounter;
        
        this->data.insert(std::pair<int, Datapoint*>(pointCounter, new Datapoint(tmpMap)));
        pointCounter++;
    }
    
    return 0;

}

std::map<int, Datapoint*>* DataReader::getData(){
    return &this->data;   
}
