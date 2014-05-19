ifndef DATAREADER_H
#define	DATAREADER_H

#include "GeneralHeader.h"

class DataReader {
public:
    DataReader(char* data_filename);
    int loadData();
    std::map<int, Datapoint*>* getData();

private:
    char* data_filename
    std::map<int, Datapoint*> data;
};

#endif	/* DATAPOINT_H */
