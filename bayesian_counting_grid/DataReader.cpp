#include "DataReader.h"

DataReader::DataReader(const char* data_filename) {
    
    this->data_filename = data_filename;
}

DataReader::~DataReader() {
}

int DataReader::loadData(){
    
    FILE *fileptr;
    int length, count, word, n, nd, nw, corpus_total = 0;
    
    printf("reading data from %s\n", data_filename);
    fileptr = fopen(data_filename, "r");
    
    while ((fscanf(fileptr, "%10d", &length) != EOF))
    {
        for (n = 0; n < length; n++)
        {
            fscanf(fileptr, "%10d:%10d", &word, &count);
        }
    }
    
    return 0;

}
