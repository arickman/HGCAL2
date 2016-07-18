#ifndef parameterScan_hh
#define parameterScan_hh

#include <string>
#include <vector>
#include <map>

#include "TCanvas.h"
#include "TGraph.h"

class parameterScan {

public:

  parameterScan(){};
  ~parameterScan(){};

  unsigned nLayers;
  double maxThick;
  double minX0;
  double maxLambda;
  
  double stepSizeW;
  double stepSizePb;
  
  unsigned nSteps;
  
  double length;
  double X0tot;
  double L0tot;

  std::vector<unsigned> iS ;
  unsigned validModels;
  unsigned totModels;


  std::map<double,unsigned> modelMap;
  std::pair<std::map<double,unsigned>::iterator,bool> isInserted;

  unsigned nC;
  std::vector<TCanvas *> myc;
  std::vector<bool> first;
  std::vector<unsigned> counter;

  double x0w;
  double x0pb;
  double l0w;
  double l0pb;


  void process30layers();
  void process28layers();
  void process26layers();
  void process24layers();
  void process22layers();
  void process20layers();

  bool processModel(const unsigned nBlocks,
		    const std::vector<unsigned> & idx);

  void print();

};//class


#endif
