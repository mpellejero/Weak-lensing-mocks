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

using namespace std;

static std::mutex barrier;

int main(int arg,char **argv){
 
  /******************** test PixelMap::PowerSpectrum **************************
  {

    int n=0;
    std::vector<PosType> pspectrum(58),multipole(58);
    {
      PixelMap kappa_map("Vipers_w4_1/vipers_w4_0.kappa.fits");
      kappa_map.PowerSpectrum(pspectrum,multipole,true);
      ++n;
    }
    /*{
      PixelMap kappa_map("Vipers_w4_2/vipers_w4_0.kappa.fits");
      kappa_map.PowerSpectrum(pspectrum,multipole,false);
     ++n;
    }
    {
      PixelMap kappa_map("Vipers_w4_3/vipers_w4_0.kappa.fits");
      kappa_map.PowerSpectrum(pspectrum,multipole,false);
     ++n;
    }
    {
      PixelMap kappa_map("Vipers_w4_4/vipers_w4_0.kappa.fits");
      kappa_map.PowerSpectrum(pspectrum,multipole,false);
     ++n;
    }
    {
      PixelMap kappa_map("Vipers_w4_5/vipers_w4_0.kappa.fits");
      kappa_map.PowerSpectrum(pspectrum,multipole,false);
     ++n;
    }
    {
      PixelMap kappa_map("Vipers_w4_6/vipers_w4_0.kappa.fits");
      kappa_map.PowerSpectrum(pspectrum,multipole,false);
     ++n;
    }
    {
      PixelMap kappa_map("Vipers_w4_7/vipers_w4_0.kappa.fits");
      kappa_map.PowerSpectrum(pspectrum,multipole,false);
     ++n;
    }
    {
      PixelMap kappa_map("Vipers_w4_8/vipers_w4_0.kappa.fits");
      kappa_map.PowerSpectrum(pspectrum,multipole,false);
     ++n;
    }
    {
      PixelMap kappa_map("Vipers_w4_9/vipers_w4_0.kappa.fits");
      kappa_map.PowerSpectrum(pspectrum,multipole,false);
     ++n;
    }
    {
      PixelMap kappa_map("Vipers_w4_10/vipers_w4_0.kappa.fits");
      kappa_map.PowerSpectrum(pspectrum,multipole,false);
     ++n;
    }*/
    
/*    for(int i=0;i<pspectrum.size();++i){
      cout << multipole[i] << "   " << multipole[i]*multipole[i]*pspectrum[i]/n << endl;
    }

    exit(0);
    cout << endl << endl;

    PixelMap kappa_map2("/Users/bmetcalf/Desktop/kappa.fits");
    kappa_map2.PowerSpectrum(pspectrum,multipole);
    
    for(int i=0;i<pspectrum.size();++i){
      cout << multipole[i] << "   " << multipole[i]*multipole[i]*pspectrum[i] << endl;
    }
 
    
  exit(0);
  }
  //****************************************************************************/
  
  
  long seed = -1827674;

  printf("initializing model\n");
  string paramfile
  ;
  if(arg > 1) paramfile.assign(argv[1],strlen(argv[1]));
  else paramfile = "ParamFiles/paramfile_field";
  cout << "using parameter file: " << paramfile << endl;
  
  cout << "Create model" << endl;
  
  InputParams params(paramfile);

  /*
  
   // testing mass by integration of kappa and gamma for different ellipticities
  for(double f=1.0;f>0.2;f-=0.1){
    
    cout << "f = " << f << endl;
    LensHaloPowerLaw halo(1.0e13,0.7,0.5,1,1, f, 0, 0);
    //LensHaloNFW halo(1.0e13, 0.7, 0.5, 0.3, f, 0, 0);
    cout<< "   mass from 1D integral R=Rmax: " << halo.MassBy1DIntegation(halo.get_Rmax()*0.9999) << endl;
    cout<< "   mass from 1D integral R>Rmax: " << halo.MassBy1DIntegation(halo.get_Rmax()*2) << endl;
    cout<< "   mass from 2D integral R=Rmax: " << halo.MassBy2DIntegation(halo.get_Rmax()*0.9999) << endl;

  
  }

  {
   // test Elliptic class
    LensHaloPowerLaw halo(1.0e13,0.25,0.5 ,1 , 1.0, 1.0, 0, 0);
        
    PosType x[2] = {halo.get_Rmax()*1.2,0},alpha[2];
    KappaType kappa,gamma[2],phi;
    
    halo.force_halo(alpha,&kappa,gamma,&phi,x);
    cout << x[0]*pi*alpha[0] << "  " << alpha[1] << endl;
    
    double tmp = alpha[0];
    alpha[0] = alpha[1] = 0.0;
    Elliptic elliptic(&halo,1.0,0.0);
    elliptic.alpha(x,alpha);

    cout << x[0]*pi*alpha[0] << "  " << alpha[1] << "  " << alpha[0]/tmp << endl;
    exit(0);
  }
  
  //LensHaloPowerLaw halo(1.0e13,0.7,0.5,1,1.0, 1.0, 0, 0);
  //LensHaloNFW halo(1.0e13, 0.7, 0.5, 0.3, 1.0, 0, 0);
 
  Lens testlens(params,&seed);
  //testlens.insertMainHalo(&halo);
  //testlens.ResetSourcePlane(2.0,false);
  
  //cout << testlens.getZlens() << "  " << testlens.getSourceZ() << endl;
  LensHaloNFW *halo = testlens.getMainHalo<LensHaloNFW>(0);
  PosType center_tmp[2] = {0,0};
  PosType range=2.5*halo->get_Rmax()/testlens.getCosmo().angDist(testlens.getZlens());
  Grid testgrid(&testlens,512,center_tmp,range);
  
  testgrid.writeFits(center_tmp,512,range/512,KAPPA,"!testhalo");
  testgrid.writeFits(center_tmp,512,range/512,ALPHA,"!testhalo");
  testgrid.writeFits(center_tmp,512,range/512,ALPHA1,"!testhalo");
  testgrid.writeFits(center_tmp,512,range/512,ALPHA2,"!testhalo");
  testgrid.writeFits(center_tmp,512,range/512,GAMMA,"!testhalo");
  testgrid.writeFits(center_tmp,512,range/512,GAMMA1,"!testhalo");
  testgrid.writeFits(center_tmp,512,range/512,GAMMA2,"!testhalo");
  
  exit(0);
   /**/
  
  // Make maps from Halo catalogs
  Lens lens(params,&seed,WMAP5yr,true);
  
  std::vector<PosType> zss = {2.297,2.119,1.955,1.802,1.66,1.527, 1.403, 1.287, 1.178, 1.075,0.9774, 0.8854, 0.7982, 0.7154,0.6365,0.5612, 0.4892,0.4201, 0.3538,0.2899, 0.2282,0.1686, 0.1108, 0.05465};
  
//  PosType center[2] = {0,0},fieldrangeX = 8.7,fieldrangeY = 1.8;
   //PosType center[2] = {0,0},fieldrangeX = 1.6,fieldrangeY = 1.6;
  //size_t NpixX = 512;
  
  // w4
  PosType center[2] = {0,0},fieldrangeX = 5.5,fieldrangeY = 1.6;
  size_t NpixX = 3300,NpixY = 960;
  
  //size_t NpixX = 5220*1.6/5.5; //x1080
  //size_t NpixY = (size_t)(NpixX*fieldrangeY/fieldrangeX + 0.5);
  //NpixY=NpixX;
  
  PixelMap map(center,NpixX,NpixY,fieldrangeX*degreesTOradians/(NpixX-1));

    //for(int i=0;i<zss.size();++i){
        for(int i=0;i<1;++i){
        
    lens.ResetSourcePlane(zss[i],false);
    cout << "   making Grid for source plane " + std::to_string(i) << "...." << endl;
    GridMap grid(&lens,NpixX/10.0,center,fieldrangeX*degreesTOradians,fieldrangeY*degreesTOradians);
    
    /*PixelMap map(center,NpixX,NpixY,fieldrangeX*degreesTOradians/(NpixX-1));
    map.AddGrid(grid,KAPPA);
    map.printFITS("!test_vipers_w4_" + std::to_string(i));
    map.Clean();*/
    
    /*
    grid.writeFits(center,NpixX,NpixY,fieldrangeX*degreesTOradians/(NpixX-1), KAPPA,"!vipers_w4_"
                   + std::to_string(i));
    grid.writeFits(center,NpixX,NpixY,fieldrangeX*degreesTOradians/(NpixX-1), GAMMA1,"!vipers_w4_"
                   + std::to_string(i));
    grid.writeFits(center,NpixX,NpixY,fieldrangeX*degreesTOradians/(NpixX-1), GAMMA2,"!vipers_w4_"
                   + std::to_string(i));
    grid.writeFits(center,NpixX,NpixY,fieldrangeX*degreesTOradians/(NpixX-1), GAMMA,"!vipers_w4_"
                   + std::to_string(i));
    grid.writeFits(center,NpixX,NpixY,fieldrangeX*degreesTOradians/(NpixX-1), ALPHA,"!vipers_w4_"
                   + std::to_string(i));
    grid.writeFits(center,NpixX,NpixY,fieldrangeX*degreesTOradians/(NpixX-1), ALPHA1,"!vipers_w4_"
                   + std::to_string(i));
    grid.writeFits(center,NpixX,NpixY,fieldrangeX*degreesTOradians/(NpixX-1), ALPHA2,"!vipers_w4_"
                   + std::to_string(i));
    */
    
    cout << "   making fits images for source plane " + std::to_string(i) << "...." << endl;
    std::string tag;
    params.get("outputfile",tag);
    tag = "!" + tag + std::to_string(i);
            
    grid.writePixelMapUniform(map,KAPPA);
    map.printFITS(tag + ".kappa.fits");
    std::vector<PosType> pspectrum(58),multipole(58);
    map.PowerSpectrum(pspectrum,multipole);
          
    for(int i=0;i<pspectrum.size();++i){
        cout << multipole[i] << "   " << multipole[i]*multipole[i]*pspectrum[i] << endl;
    }
          
    grid.writePixelMapUniform(map,GAMMA);
    map.printFITS(tag + ".gamma.fits");
    grid.writePixelMapUniform(map,GAMMA1);
    map.printFITS(tag + ".gamma1.fits");
    grid.writePixelMapUniform(map,GAMMA2);
    map.printFITS(tag + ".gamma2.fits");
    grid.writePixelMapUniform(map,ALPHA);
    map.printFITS(tag + ".alpha.fits");
    grid.writePixelMapUniform(map,ALPHA1);
    map.printFITS(tag + ".alpha1.fits");
    grid.writePixelMapUniform(map,ALPHA2);
    map.printFITS(tag + ".alpha2.fits");
    
  }
  
  cout << "   finished" << endl;
  return 0;
}


