#include "DataReader.h"
#include "Datapoint.h"

DataReader::DataReader(string data_filename) {
    
    this->data_filename = data_filename;
}

DataReader::~DataReader() {
}

int DataReader::loadData(){
    
    cout<<"let's see"<<endl;
    int length, count, word, n;
    int pointCounter;
    string dataLine;
    string s;
    cout<<"standard types OK"<<endl;
    mapdata tmpMap;
    
    cout<<"map type OK"<<endl;
    ostringstream convert;
    cout<<"stream OK"<<endl;
    
    //cout<<"reading data from";//<<this->data_filename;
    
    ifstream dataFile ((char*)this->data_filename.c_str());
    cout<<"ifstream OK"<<endl;
    //fileptr = fopen(data_filename, "r");
    
    pointCounter = 0;
    //while ((fscanf(fileptr, "%10d", &length) != EOF))
    while ( getline (dataFile, dataLine) )
    {   
        tmpMap.clear();
        cout<<"map is clear"<<endl;
        std::stringstream linestream(dataLine);
        while(getline(linestream, s, ' '))
        {
            sscanf((char*)s.c_str(), "%10d:%10d", &word, &count);
            //cout<<s<<" :: "<<word<<" - "<<count<<endl;
            tmpMap.insert(std::pair<int,int>(word, count));
        }
//        for (n = 0; n < length; n++)
//        {
//            sscanf((char*)dataLine.c_str(), "%10d:%10d", &word, &count);
//            
//            tmpMap.insert(std::pair<int,int>(word, count));
//        }
        
        convert << pointCounter;
        
        cout<<"done with tmpMap of "<<pointCounter<<endl;
        
        this->data.insert(std::pair<int, Datapoint*>(pointCounter, new Datapoint(tmpMap)));
        pointCounter++;
        cout<<"point "<<pointCounter -1 <<" done"<<endl;
    }
    
    return 0;

}

std::map<int, Datapoint*>* DataReader::getData(){
    return &this->data;   
}
