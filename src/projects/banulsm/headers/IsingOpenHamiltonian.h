#ifndef ISINGOPEN_H
#define ISINGOPEN_H

#include "Hamiltonian.h"
#include <vector>

class IsingOpenHamiltonian: 
public Hamiltonian{
  int L;
  double tau1,tau2,g,gamma1,gamma2;
  MPO hamil,Lindblad,Iden,L2,L5,L5_2,L6;;
  mwArray Z,Zz,z3,Z4,Z5,Z5_2,Z6; 
 
  //  std::vector<double> gis;

 public:
  IsingOpenHamiltonian(int L,double* tau1is,double* tau2is,double* gis,double* gamma1is,double* gamma2is);
   ~IsingOpenHamiltonian(){};
   const MPO& getHMPO() const {return hamil;}
   //MPO& getHMPO() const {return hamil;}
   MPO& getLMPO(){return Lindblad;}
   MPO& getIDMPO(){return Iden;}
   MPO& getL2MPO(){return L2;}
   MPO& getL5MPO(){return L5;}
   MPO& getL5_2MPO(){return L5_2;}
   MPO& getL6MPO(){return L6;}

  bool hasHMPO() const {return true;}

  int getBondDimension() const {return 4;}


 private:
  void initZ();
  void initHMPO();
  void initZz();
  void initLMPO();
  void initz3();
  void initIden();
  void initz4();
  void initL2MPO();
  void initz5();
  void initL5MPO();
  void initz5_2();
  void initL5_2MPO();
  void initz6();
  void initL6MPO();
  //void initHMPO(double* tis=0,double* Uis=0,double* muis=0,
  //	 double* Vis=0);

};


#endif // ISINGOPEN_H
