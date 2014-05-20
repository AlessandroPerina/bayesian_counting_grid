#ifndef DATAREADER_H
#define	DATAREADER_H

#include "GeneralHeader.h"

class DataReader {
public:
    DataReader(string data_filename);
    virtual ~DataReader();
    int loadData();
    std::map<int, Datapoint*>* getData();

private:
    string data_filename;
    std::map<int, Datapoint*> data;
};

#endif	/* DATAREADER_H */
