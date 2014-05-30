#include "DataReader.h"
#include "Datapoint.h"

DataReader::DataReader(string data_filename) {
    
    this->data_filename = data_filename;
}

DataReader::~DataReader() {
}

int DataReader::loadData(){
    
    int length, count, word, n;
    int pointCounter;
    string dataLine;
    string s;
    
    mapdata tmpMap;
    
    ostringstream convert;
    
    ifstream dataFile ((char*)this->data_filename.c_str());
    
    pointCounter = 0;
    
    while ( getline (dataFile, dataLine) )
    {   
        //cout<<dataLine<<endl;
        
        tmpMap.clear();
    
        std::stringstream linestream(dataLine);
        while(getline(linestream, s, ' '))
        {
            //cout<<s<<endl;
            sscanf((char*)s.c_str(), "%10d:%10d", &word, &count);
            tmpMap.insert(std::pair<int,int>(word, count));
        }

        convert << pointCounter;
        
        //cout<<pointCounter<<endl;
        //cout<<"pre-point"<<endl;
        this->data.insert(std::pair<unsigned long, Datapoint*>(pointCounter, new Datapoint(tmpMap)));
        //cout<<"post-point"<<endl;
        pointCounter++;

    }
    
    return 0;

}

std::map<unsigned long, Datapoint*>* DataReader::getData(){
    return &this->data;   
}
