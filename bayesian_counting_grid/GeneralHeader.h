#ifndef GENERALHEADER_H
#define	GENERALHEADER_H

#define _USE_MATH_DEFINES

#define ARMA_NO_DEBUG

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <omp.h>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/gamma_distribution.hpp>

// #include <boost/python.hpp>
// #include <Python.h>
//#include <numpy/arrayobject.h>
//#include <numpy/ndarrayobject.h>

//#include <boost/math/special_functions/gamma.hpp>
//#include <boost/lexical_cast.hpp>
//#include <boost/random.hpp>
//#include <boost/math/special_functions/digamma.hpp>
//#include <omp.h>

#include <armadillo>
#endif	/* GENERALHEADER_H */
#include <assert.h>


//class Dataframe;
class Datapoint;
class DataReader;
class CountingGrid;


using namespace std;
using namespace arma;
typedef std::map<int,int> mapdata;

#define DIM 2
#define CG_ROWS 32
#define CG_COLS 32
#define WD_ROWS 5
#define WD_COLS 5
// #define Z 20
#define Z 24678
#define BASE_PRIOR 0.1
#define ASSIGN_TOKEN true


/*#ifndef DEFINE_CONST
map<float, float> gammaLookUp;
#endif
*/

// Define Counting Grid parameters
/*
const int DIM = 2;
const int CG_ROWS = 32;
const int CG_COLS = 32;
const int WD_ROWS = 5;
const int WD_COLS = 5;
const int Z = 2000;
const uword LOCS = CG_ROWS*CG_COLS*Z;
const double BASE_PRIOR = 0.1;
map<float, float> gammaLookUp; // Poi eliminalo
*/

// da aggiungere nel path in caso non salvi  C:\opencv\build\include;C:\Program Files\boost\boost_1_55_0\;C:\Python27\include;C:\Program Files (x86)\Inkscape\python\Lib\site-packages\numpy\core\include;C:\Users\APerina\Documents\Visual Studio 2010\Armadillo\include;
