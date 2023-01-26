#include "FiniteVolume.h"
#include "DataFile.h"
#include "Mesh.h"
#include "Physics.h"

#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Sparse"

#include <iostream>
#include <cmath>
#include <algorithm>


//--------------------------------------------------------//
//---------------Classe mère flux numérique---------------//
//--------------------------------------------------------//
FiniteVolume::FiniteVolume()
{
}



FiniteVolume::FiniteVolume(DataFile* DF, Mesh* mesh, Physics* physics):
  _DF(DF), _mesh(mesh), _physics(physics), _fluxVector(_mesh->getNumberOfCells(), 2)
{
}



void FiniteVolume::Initialize(DataFile* DF, Mesh* mesh, Physics* physics)
{
  _DF = DF;
  _mesh = mesh;
  _physics = physics;
  _fluxVector.resize(_mesh->getNumberOfCells(), 2);
  std::cout << "Ini okay !" << std::endl;

}




// Minmod slope limiter
double FiniteVolume::minmod(double r) const
{
  if (r < 0)
    return 0.;
  else if (r > 1)
    return r;
  else
    return 1;
}


  Eigen::Matrix<double, Eigen::Dynamic, 2> FiniteVolume::buildGidemi2( double t, Eigen::Matrix<double, Eigen::Dynamic, 2> Sol)
  {
    _fluxVector.setZero();

    // Get mesh parameters
    int nCells(_mesh->getNumberOfCells());
    double dx(_mesh->getSpaceStep());

    // Get gravity
    double g(_DF->getGravityAcceleration());

    // Vectors to store the reconstruted values at the left and right of each interface
    Eigen::Matrix<double, Eigen::Dynamic, 2> SolD, SolG;
    SolD.resize(nCells + 1, 2);
    SolG.resize(nCells + 1, 2);

    Eigen::Matrix<double, Eigen::Dynamic, 2> slopes, limSlopes;
    slopes.resize(nCells + 1, 2);
    limSlopes.resize(nCells, 2);

      // Compute the slopes
      // Left boundary
      Eigen::Vector2d leftBoundarySol(_physics->leftBoundaryFunction(t + _DF->getTimeStep(), Sol));
      slopes(0,0) = (Sol(0,0) - leftBoundarySol(0)) / dx;
      slopes(0,1) = (Sol(0,1) - leftBoundarySol(1)) / dx;
      // Right boundary
      Eigen::Vector2d rightBoundarySol(_physics->rightBoundaryFunction(t+ _DF->getTimeStep(), Sol));
      slopes(nCells, 0) = (rightBoundarySol(0) - Sol(nCells - 1, 0)) / dx;
      slopes(nCells, 1) = (rightBoundarySol(1) - Sol(nCells - 1, 1)) / dx;
      // Interior edges
      for (int i(1) ; i < nCells ; ++i)
        {
          slopes.row(i) = (Sol.row(i) - Sol.row(i-1)) / dx;
        }

      // Limit the slopes
      for (int i(0) ; i < nCells - 1 ; ++i)
        {
          limSlopes(i,0) = minmod((slopes(i,0)/slopes(i+1,0)));
          limSlopes(i,1) = minmod((slopes(i,1)/slopes(i+1,1)));
        }

      // Reconstruct the values at each edge
      // Left boundary
      SolG.row(0) = leftBoundarySol;
      SolD.row(0) = Sol.row(0) - 0.5 * dx * limSlopes.row(0);
      // Right boundary
      SolG.row(nCells) = Sol.row(nCells - 1) + 0.5 * dx * limSlopes.row(nCells - 1);
      SolD.row(nCells) = rightBoundarySol;
      // Interior edges
      for (int i(1) ; i < nCells ; ++i)
        {
          SolG.row(i) = Sol.row(i-1) + 0.5 * dx * limSlopes.row(i-1);
          SolD.row(i) = Sol.row(i) - 0.5 * dx * limSlopes.row(i);
        }
        _fluxVector.row(0) += buildGidemi(SolG.row(0), SolD.row(0));
    // Interior fluxes contribution
    for (int i(1) ; i < nCells; ++i)
      {
        Eigen::Vector2d flux(buildGidemi(SolG.row(i), SolD.row(i)));
        _fluxVector.row(i-1) -= flux;
        _fluxVector.row(i) += flux;
      }
    // Right boundary contribution
    _fluxVector.row(nCells - 1) -= buildGidemi(SolG.row(nCells), SolD.row(nCells));
    return _fluxVector;
    std::cout << "Finite Volume methode !" << std::endl;

  }




  //-----------------------------------------------//
  //---------------Build G et Gidemi-------------------------//
  //-----------------------------------------------//

      //Construit Gidemi
      Eigen::Vector2d FiniteVolume::buildGidemi(const Eigen::Vector2d& SolD,const Eigen::Vector2d& SolG)
      {
        Eigen::Vector2d fluxidemi;
        double amin(0.);
        double amax(0.);
        _physics->computeWaveSpeed(SolG, SolD, &amin, &amax);
    //    double dt(_DF->getTimeStep());
    //    double dx(_mesh->getSpaceStep());
    //    Eigen::Vector2d PenteG();
    //.    Eigen::Vector2d PenteD;
    std::cout << amin << amax << '\n';
        if (amin>0 )
        {
          fluxidemi=(_physics->G(SolG));
        }
        else if (amin<0 && amax>0)
        {
          fluxidemi=(amax*_physics->G(SolG)-amin*_physics->G(SolD)-amin*amax*(SolD-SolG))/(amax-amin);
        }
        else if(amax<0)
        {
          fluxidemi=(_physics->G(SolD));
        }
        return fluxidemi;
        std::cout << "G !" << std::endl;

      }
