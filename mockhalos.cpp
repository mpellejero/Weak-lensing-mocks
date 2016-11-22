
/*
 * multiplane.cpp
 *
 *  Created on: Oct 25, 2016
 *      Author: bmetcalf & mpellejero
 *
 *      This program is for testing the adaptive griding on a random field by adapting to the
 *      high convergence regions.
 */

#include <slsimlib.h>
#include <standard.h>
#include <sstream>
#include <string.h>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <omp.h>
#include <thread>
#include <mutex>
//#include <geometry.h>
#include "elliptic.h"
#include "gridmap.h"
#include "lightcone_construction.h"


class TestClass{
  
  int i;
  
public:
  TestClass(int ii):i(ii){}
  double func(double x) const{ return x*i;}
  
};


using namespace std;

static std::mutex barrier;

int main(int arg,char **argv){

  /*{
    COSMOLOGY c(BigMultiDark);
    
    cout << c.rho_crit(0)*c.getOmega_matter() << " Msun/Mpc^3" <<endl;
    
    double sigma = 0.1,res=arcsecTOradians;
    Point_2d center;
    PixelMap map(center.x,1000,1000,res);
     Utilities::RandomNumbers_NR ran(1277);
    
    
    for(size_t i=0; i<map.size() ;++i){
      map[i] = sigma*ran.gauss();
    }
    
    std::vector<double> ps(10),ll(10);
    
    map.PowerSpectrum(ps,ll,1,false,true);
    for(int i=0 ; i < ps.size() ; ++i){
      cout << ll[i] << "  " << ps[i] << endl;
    }
    
    cout << 1.0/sqrt(map.size()) << endl;
    cout << res*res*sigma*sigma/pi/pi/4 << endl;
    //exit(1);
  }*/
  
  time_t t0,t1;
  time(&t0);
  long seed = -1827674;

  printf("initializing model\n");
  string paramfile;
  if(arg > 1) paramfile.assign(argv[1],strlen(argv[1]));
  else paramfile = "ParamFiles/paramfile_field";
  
  cout << "using parameter file: " << paramfile << endl;

  if(arg > 2) seed = stoi(argv[2]);

  cout << " seed = " << seed << endl;
  cout << "Create model" << endl;
  
  InputParams params(paramfile);

  // Make maps from Halo catalogs
  Lens lens(params,&seed,BigMultiDark,false);
  
  // source redhsifts
  //std::vector<PosType> zss = {2.297,2.119,1.955,1.802,1.66,1.527, 1.403, 1.287, 1.178, 1.075,0.9774, 0.8854, 0.7982, 0.7154,0.6365,0.5612, 0.4892,0.4201, 0.3538,0.2899, 0.2282,0.1686, 0.1108, 0.05465};

  std::vector<PosType> zss = {2.297,1.075,0.4892};

  PosType center[2] = {0,0};
  size_t NpixX = 512*2;
  
  //NpixX = 64;
  
  for(int i=0;i<zss.size();++i){
  //for(int i=0;i<1;++i){
        
    lens.ResetSourcePlane(zss[i],false);
    cout << "   making Grid for source plane " + std::to_string(i) << "...." << endl;
    GridMap grid(&lens,NpixX,center,3*1.3*degreesTOradians);
    
    cout << "   making fits images for source plane " + std::to_string(i)
    << " ...." << endl;
   
    std::string tag;
    params.get("outputfile",tag);
    tag = tag + std::to_string(zss[i]);
    tag = "Output_3deg_1Rvir_large/" + tag;
    
    //    grid.writePixelMapUniform(map,KAPPA);
    PixelMap map=grid.writePixelMapUniform(KAPPA);
    map.printFITS("!" + tag + ".kappa.fits");
    std::vector<PosType> pspectrum(40),multipole(40);
    map.PowerSpectrum(pspectrum,multipole,1,false,true);
    
    std::ofstream ps_file(tag + "PS" + ".csv");
    ps_file << "l,PS" << endl;
    
    for(int i=0;i<pspectrum.size();++i){
      ps_file << multipole[i] << "," << pspectrum[i] << endl;

      cout << multipole[i] << "   " << multipole[i]*multipole[i]*pspectrum[i] << endl;
    }
    ps_file.close();
    
    grid.writePixelMapUniform(map,GAMMA);
    map.printFITS("!" + tag + ".gamma.fits");
    grid.writePixelMapUniform(map,GAMMA1);
    map.printFITS("!" + tag + ".gamma1.fits");
    grid.writePixelMapUniform(map,GAMMA2);
    map.printFITS("!" + tag + ".gamma2.fits");
    grid.writePixelMapUniform(map,ALPHA);
    map.printFITS("!" + tag + ".alpha.fits");
    grid.writePixelMapUniform(map,ALPHA1);
    map.printFITS("!" + tag + ".alpha1.fits");
    grid.writePixelMapUniform(map,ALPHA2);
    map.printFITS("!" + tag + ".alpha2.fits");
  }
  
  cout << "   finished" << endl;
  time(&t1);
  cout << "   time: " << difftime(t1,t0) << endl;
  return 0;
}


