#include "Physics.h"
#include "DataFile.h"
#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Sparse"

#include <iostream>
#include <fstream>
#include <regex>
#include <cmath>
#include <algorithm>
#include <complex>



//------------------------------------------//
//---------------Constructors---------------//
//------------------------------------------//
Physics::Physics()
{
}



Physics::Physics(DataFile* DF, Mesh* mesh):
  _DF(DF), _mesh(mesh), _xmin(mesh->getxMin()), _xmax(mesh->getxMax()), _g(_DF->getGravityAcceleration()), _nCells(mesh->getNumberOfCells()), _i(0)
{
}



//--------------------------------------------//
//---------------Initialization---------------//
//--------------------------------------------//

void Physics::Initialize(DataFile* DF, Mesh* mesh)
{
  _DF = DF;
  _mesh = mesh;
  _xmin = mesh->getxMin();
  _xmax = mesh->getxMax();
  _g = DF->getGravityAcceleration();
  _i = 0;
  _nCells = mesh->getNumberOfCells();
  this->Initialize();
}



void Physics::Initialize()
{
  // Build
  buildTopography();
  buildInitialCondition();
}



//----------------------------------------------//
//---------------Build Topography---------------//
//----------------------------------------------//


void Physics::buildTopography()
{
  _topography.resize(_nCells);
  const Eigen::VectorXd& cellCenters(_mesh->getCellCenters());

  // Flat botttom
  if (_DF->getTopographyType() == "FlatBottom"||_DF->getTopographyType() == "RestingLake")
    {
      _topography.setZero();
    }
  // XingandShu test case topography
  else if (_DF->getTopographyType() == "XingandShu")
    {
      _topography.setZero();
      double xmin(_DF->getXmin()), xmax(_DF->getXmax()), L(xmax - xmin);
      double a(1.), h0(0.5);
      for (int i(0) ; i < _nCells ; ++i)
        {
          double x(cellCenters(i));
          _topography(i) = pow(sin(M_PI*x),2);
        }
    }
  // MovingBquad test case topography
  else if (_DF->getTopographyType() == "MovingBquad" || _DF->getTopographyType() == "MovingBPara")
    {
      _topography.setZero();
      double xmin(_DF->getXmin()), xmax(_DF->getXmax()), L(xmax - xmin);
      double a(1.), h0(0.5);
      for (int i(0) ; i < _nCells ; ++i)
        {
          double x(cellCenters(i));
          _topography(i) = 1-10*pow(x/3000,2);
        }
    }
  // MovingBquad test case topography
  else if (_DF->getTopographyType() == "SubcriticalFlow")
    {
      _topography.setZero();
      double xmin(_DF->getXmin()), xmax(_DF->getXmax()), L(xmax - xmin);
      double a(1.), h0(0.5);
      for (int i(0) ; i < _nCells ; ++i)
        {
          double x(cellCenters(i));
          _topography(i) = 0.8+0.25*exp(-135*pow((x-75)/150,2)/4);
        }
    }
  else
    {
      std::cout <<  "ERROR::TOPOGRAPHY : Case not implemented" << std::endl;
      std::cout <<  "====================================================================================================" << std::endl;
      exit(-1);
    }
}



//-----------------------------------------------------//
//---------------Build Initial Condition---------------//
//-----------------------------------------------------//
void Physics::buildInitialCondition()
{
  _Sol0.resize(_nCells, 2);
  const Eigen::VectorXd& cellCenters(_mesh->getCellCenters());
  if (_DF->getInitialCondition() == "UniformHeightAndDischarge")
    {
      double H0(_DF->getInitialHeight()), q0(_DF->getInitialDischarge());
      for (int i(0) ; i < _nCells ; ++i)
        {
          _Sol0(i,0) = std::max(H0 - _topography(i), 0.);
          _Sol0(i,1) = q0;
        }
    }

  else if (_DF->getInitialCondition() == "XingandShu")
    {
      _Sol0.resize(_nCells, 2);
      for (int i(0) ; i < _nCells ; ++i)
        {
          double x(cellCenters(i));
          _Sol0(i,0) = 5+exp(cos(2*M_PI*x));
          _Sol0(i,1) = sin(cos(2*M_PI*x));
        }
    }
  else if (_DF->getInitialCondition() == "MovingBquad"||_DF->getInitialCondition() == "MovingBquad")
    {
      _Sol0.col(1).setZero();
      _Sol0.col(0).fill(10.);
    }
  else if (_DF->getInitialCondition() == "SubcriticalFlow"||_DF->getInitialCondition() == "RestingLake")
    {
      _Sol0.col(1).setZero();
      _Sol0.col(0).setZero();
      }
  else
    {
      std::cout <<  "ERROR::INITIALCONDITION : Case not implemented" << std::endl;
      std::cout <<  "====================================================================================================" << std::endl;
      exit(-1);
    }
}





//--------------------------------------------------//
//---------------Build Exact Solution---------------//
//--------------------------------------------------//
void Physics::buildExactSolution(double t)
 {
   _exactSol.resize(_nCells, 2);
   const std::string& testCase(_DF->getTestCase());
   const Eigen::VectorXd& cellCenters(_mesh->getCellCenters());
   // Resting lake solutions
//   if (testCase == "RestingLake")
//     {
//       for (int i(0) ; i < _nCells ; ++i)
//         {
//           _exactSol(i,0) = std::max(_DF->getLeftBCImposedHeight() - _topography(i), 0.);
//           _exactSol(i,1) = 0.;
//         }
//     }
   if (testCase == "MovingBquad"||testCase == "MovingBPara")
     {
       double s(sqrt(8.*9.81*10./pow(3000.,2)-pow(0.001,2)));
       for (int i(0) ; i < _nCells ; ++i)
         {
           double x(cellCenters(i));
           _exactSol(i,0) = pow(3000*2/9.81,2)*exp(-0.001*t)*(-s*0.001*sin(2*s*t)+(pow(0.001,2)/4-pow(s,2))*cos(2*s*t))/(8*10)-9.81*exp(-0.001*t)-exp(-s*t/2)*(2*s*cos(s*t)+0.001*sin(s*t))*x/9.81;
           _exactSol(i,1) = 2*exp(-s*t/2)*sin(s*t)*_exactSol(i,0);
         }
     }
//   // Non hydrostatic stationnary solutions
//   else if (testCase == "SubcriticalFlow" || testCase == "TranscriticalFlowWithoutShock" || testCase == "TranscriticalFlowWithShock")
//     {
//       double epsilon(1.0 / _nCells);
//       double qIn(_DF->getLeftBCImposedDischarge());
//       double hOut(_DF->getRightBCImposedHeight());
//       double hMiddle(pow(pow(qIn,2)/_g, 1./3.)); // = critical height
//       double hMax(3.), zMax(0.2);
//       double a(0.), b(0.), c(0.), d(0.);
//       double p, q;
//       // Discharge must be constant in the whole domain
//       for (int i(0) ; i < _nCells ; ++i)
//         {
//           _exactSol(i,1) = qIn;
//         }
//       // Subcritical flow
//       if (testCase == "SubcriticalFlow")
//         {
//           for (int i(_nCells - 1) ; i >= 0  ; --i)
//             {
//               double z(_topography(i)), zEnd(_topography(_nCells - 1));
//               computeCoeffabcd(qIn, hOut, z, zEnd, &a, &b, &c, &d);
//               p = cardanP(a, b, c);
//               q = cardanQ(a, b, c, d);
//               double hnear;
//               if (i == _nCells - 1)
//                 hnear = hOut;
//               else
//                 hnear = _exactSol(i+1, 0);
//               _exactSol(i,0) = exactHeight(p, q, a, b, hnear, hMax);
//             }
//         }
//       else if (testCase == "TranscriticalFlowWithoutShock")
//         {
//           // Subcritical part (before the bump)
//           for (int i(2. * _nCells / 5. - 1) ; i >= 0 ; --i)
//             {
//               double z(_topography(i));
//               computeCoeffabcd(qIn, hMiddle, z, zMax, &a, &b, &c, &d);
//               p = cardanP(a, b, c);
//               q = cardanQ(a, b, c, d);
//               double hnear;
//               if (i == 2. * _nCells / 5. - 1)
//                 hnear = hMiddle;
//               else
//                 hnear = _exactSol(i+1, 0);
//               _exactSol(i,0) = exactHeight(p, q, a, b, hnear*(1+epsilon), hMax);
//             }
//           // Critical part (middle of the bump)
//           double z(_topography(2. * _nCells / 5.));
//           computeCoeffabcd(qIn, hMiddle, z, zMax, &a, &b, &c, &d);
//           p = cardanP(a, b, c); q = cardanQ(a, b, c, d);
//           _exactSol(2. * _nCells / 5., 0) = exactHeight(p, q, a, b, hMiddle, hMax);
//           // Supercritical part (after the bump)
//           for (int i(2. * _nCells / 5. + 1) ; i < _nCells ; ++i)
//             {
//               double z(_topography(i));
//               computeCoeffabcd(qIn, hMiddle, z, zMax, &a, &b, &c, &d);
//               p = cardanP(a, b, c); q = cardanQ(a, b, c, d);
//               _exactSol(i,0) = exactHeight(p, q, a, b, _exactSol(i-1, 0)*(1-epsilon), hMax);
//             }
//         }
//       else if (testCase == "TranscriticalFlowWithShock")
//         {
//           // Search for the limit between sub-super-sub
//           double test(100.);
//           double hplus(0.), hminus(0.);
//           int abslim(2 * _nCells / 5);
//           double epsi(10.0/_nCells);
//           double zEnd(_topography(_nCells - 1));
//           while(test > epsi && abslim < _nCells)
//             {
//               computeCoeffabcd(qIn, hOut, _topography(abslim), zEnd, &a, &b, &c, &d);
//               p = cardanP(a, b, c);
//               q = cardanQ(a, b, c, d);
//               hplus = exactHeight(p, q, a, b, hOut, hMax);
//               computeCoeffabcd(qIn, hMiddle, _topography(abslim), zMax, &a, &b, &c, &d);
//               p = cardanP(a, b, c);
//               q = cardanQ(a, b, c, d);
//               hminus = exactHeight(p, q, a, b, hMiddle, hMax);
//               test = RHJump(hplus, hminus, qIn);
//               ++abslim;
//             }
//           _exactSol(abslim, 0) = hminus;

//           for (int i(abslim - 1) ; i >= 0 ; --i)
//             {
//               computeCoeffabcd(qIn, hMiddle, _topography(i), zMax, &a, &b, &c, &d);
//               p = cardanP(a, b, c);
//               q = cardanQ(a, b, c, d);
//               _exactSol(i,0) = exactHeight(p, q, a, b, _exactSol(i+1,0) + epsilon, hMax);
//             }
//           for (int i(_nCells - 1) ; i > abslim ; --i)
//             {
//               computeCoeffabcd(qIn, hOut, _topography(i), zEnd, &a, &b, &c, &d);
//               p = cardanP(a, b, c);
//               q = cardanQ(a, b, c, d);
//               double hnear;
//               if (i == _nCells - 1)
//                 hnear = hOut;
//               else
//                 hnear = _exactSol(i+1,0);
//               _exactSol(i,0) = exactHeight(p, q, a, b, hnear, hMax);
//             }
//         }
//     }
 }


// Save the exact solution in a file
void Physics::saveExactSolution(std::string& fileName) const
{
#if VERBOSITY>0
  std::cout << "Saving exact solution" << std::endl;
#endif
  std::ofstream outputFile(fileName, std::ios::out);
  const Eigen::VectorXd& cellCenters(_mesh->getCellCenters());
  outputFile << "# x  H=h+z   h       u       q       Fr=|u|/sqrt(gh)" << std::endl;
  for (int i(0) ; i < _exactSol.rows() ; ++i)
    {
      outputFile << cellCenters(i) << " " <<
        _exactSol(i,0) + _topography(i) << " " <<
        _exactSol(i,0) << " " <<
        _exactSol(i,1)/_exactSol(i,0) << " " <<
        _exactSol(i,1) << " " <<
        abs(_exactSol(i,1)/_exactSol(i,0))/sqrt(_g * _exactSol(i,0)) << std::endl;
    }
}


//-------------------------------------------//
//---------------Physical flux---------------//
//-------------------------------------------//
Eigen::Vector2d Physics::G(const Eigen::Vector2d& Sol) const
{
  Eigen::Vector2d flux;
  double h(Sol(0)), qx(Sol(1));
  flux(0) = qx;
  flux(1) = qx*qx/h + 0.5*_g*h*h;
  return flux;
}


//----------------------------------------//
//---------------Wave speed---------------//
//----------------------------------------//
void Physics::computeWaveSpeed(const Eigen::Vector2d& SolG, const Eigen::Vector2d& SolD, double* amin, double* amax) const
{
  double hG(SolG(0)), hD(SolD(0));
  double uG(SolG(1)/hG), uD(SolD(1)/hD);
  if (hG < 1e-6)
    uG = 0.;
  if (hD < 1e-6)
    uD = 0.;
  *amin = std::min(uG - sqrt(_g * hG), uD - sqrt(_g * hD));
  *amax = std::max(uG + sqrt(_g * hG), uD + sqrt(_g * hD));
}

//------------------------------------------------------//
//---------------Left Boundary Conditions---------------//
//------------------------------------------------------//
Eigen::Vector2d Physics::leftBoundaryFunction(double t, const Eigen::Matrix<double, Eigen::Dynamic, 2>& Sol)
{
  Eigen::Vector2d SolG(0.,0.);

  // Calcul du nombre de Froude au bord
  double h(Sol(0,0)), q(Sol(0,1));
  double Fr(abs(q)/(h * sqrt(_g * h)));


  // Choix entre les differentes CL
  if (_DF->getLeftBC() == "Neumann")
    {
      SolG(0) = Sol(0,0);
      SolG(1) = Sol(0,1);
    }
  else if (_DF->getLeftBC() == "Wall")
    {
      SolG(0) = Sol(0,0);
      SolG(1) = 0.;
    }
  else if (_DF->getLeftBC() == "ImposedConstantDischarge")
    {
      // Entrée/sortie fluviale
      if (Fr < 1)
        {
          SolG(0) = Sol(0,0);
          SolG(1) = _DF->getLeftBCImposedDischarge();
        }
      // Sortie torrentielle (sortie libre, on n'impose rien)
      else if (Fr > 1 && q < 0)
        {
          SolG(0) = Sol(0,0);
          SolG(1) = Sol(0,1);
        }
      // Entrée torrentielle (on impose une hauteur et un debit)
      else if (Fr > 1 && q > 0)
        {
          SolG(0) = _DF->getLeftBCImposedHeight();
          SolG(1) = _DF->getLeftBCImposedDischarge();
        }
    }
  else if (_DF->getLeftBC() == "PeriodicWaves" || _DF->getLeftBC() == "DataFile" || _DF->getLeftBC() == "ImposedConstantHeight")
    {
      // Recupere la solution dans les mailles de centre x1 et x2 ainsi que dx et dt
      double h1(Sol(0,0)), h2(Sol(1,0));
      double u1(Sol(0,1)/h1), u2(Sol(1,1)/h2);
      double dx(_DF->getDx()), dt(_DF->getTimeStep());
      double x1(_DF->getXmin() + 0.5*dx);
      double a(pow(1 + dt/dx * (u2 - u1), 2));
      double b(2*dt*(u1 - x1/dx * (u2 - u1)) * (1 + dt/dx * (u2 - u1)) - dt*dt*_g*(h2 - h1)/dx);
      double c(pow(dt*u1 - dt/dx * x1 * (u2 - u1), 2) - dt*dt * _g * (h1 - x1/dx * (h2 - h1)));
      // std::cout << a << " " << b << " " << c << std::endl;
      double xe(FindRacine(a, b, c));
      // std::cout << xe << std::endl;
      double uXe(u1 + (xe - x1)*(u2 - u1)/dx);
      double hXe(h1 + (xe - x1)*(h2 - h1)/dx);
      if (_DF->getLeftBC() =="MovingBquad")
        {
          double s(sqrt(8.*9.81*10./pow(3000.,2)-pow(0.001,2)));
          SolG(0) = 3000*3000*exp(-0.001*t)*(-s*0.001*sin(2*s*t)+(0.001*0.001/4-s*s)*cos(2*s*t))/(2*9.81*9.81*10)-exp(-0.001*t)/9.81;
          SolG(1) = SolG(0)*2*exp(-0.0005*t)*sin(s*t);
        }
    }
  return SolG;
}



//-------------------------------------------------------//
//---------------Right Boundary Conditions---------------//
//-------------------------------------------------------//
Eigen::Vector2d Physics::rightBoundaryFunction(double t, const Eigen::Matrix<double, Eigen::Dynamic, 2>& Sol)
{
  Eigen::Vector2d SolD(0.,0.);

  // Calcul du nombre de Froude au bord
  double h(Sol(_nCells - 1,0)), q(Sol(_nCells - 1,1));
  double Fr(abs(q)/(h * sqrt(_g * h)));

  // Choix entre les differentes CL
  if (_DF->getRightBC() == "Neumann")
    {
      SolD(0) = Sol(_nCells - 1,0);
      SolD(1) = Sol(_nCells - 1,1);
    }
  else if (_DF->getRightBC() == "Wall")
    {
      SolD(0) = Sol(_nCells - 1,0);
      SolD(1) = 0.;
    }
  else if (_DF->getRightBC() == "ImposedConstantDischarge")
    {
      // Entrée/sortie fluviale
      if (Fr < 1)
        {
          SolD(0) = Sol(_nCells - 1,0);
          SolD(1) = _DF->getRightBCImposedDischarge();
        }
      // Sortie torrentielle (sortie libre, on n'impose rien)
      else if (Fr > 1 && q > 0)
        {
          SolD(0) = Sol(_nCells - 1,0);
          SolD(1) = Sol(_nCells - 1,1);
        }
      // Entrée torrentielle (on impose une hauteur et un debit)
      else if (Fr > 1 && q < 0)
        {
          SolD(0) = _DF->getRightBCImposedHeight() - _topography(_nCells - 1);
          SolD(1) = _DF->getRightBCImposedDischarge();
        }
    }
  else if (_DF->getRightBC() == "PeriodicWaves" || _DF->getRightBC() == "DataFile" || _DF->getRightBC() == "ImposedConstantHeight")
    {
      // Recupere la solution dans la maille de bord
      double h1(Sol(_nCells - 1,0)), u1(Sol(_nCells - 1,1)/h1);
      if (_DF->getRightBC() == "ImposedConstantHeight")
        {
          // Entrée/sortie fluviale
          if (Fr < 1)
            {
              SolD(0) = _DF->getRightBCImposedHeight() - _topography(_nCells - 1);
              SolD(1) = SolD(0) * (u1 + 2. * sqrt(_g * h1) - 2. * sqrt(_g * SolD(0)));
            }
          // Sortie torrentielle (sortie libre, on n'impose rien)
          else if (Fr > 1 && q > 0)
            {
              SolD(0) = Sol(_nCells - 1,0);
              SolD(1) = Sol(_nCells - 1,1);
            }
          // Entrée torrentielle (on impose une hauteur et un debit)
          else if (Fr > 1 && q < 0)
            {
              SolD(0) = _DF->getRightBCImposedHeight() - _topography(_nCells - 1);
              SolD(1) = _DF->getRightBCImposedDischarge();
            }
        }
      else if (_DF->getRightBC() == "PeriodicWaves")
        {
          SolD(0) = 3. + 0.1*sin(5 * M_PI * t);
          SolD(1) = SolD(0) * (u1 + 2. * sqrt(_g * h1) - 2. * sqrt(_g * SolD(0)));
        }
    }
  return SolD;
}



// Other
double Physics::FindRacine(double a, double b, double c)
{
  double delta;
  delta = b*b - 4*a*c;
  if (delta < 0)
    {
      std::cout << "Pas de racine reelle" << std::endl;
      exit(1);
    }
  else if (delta == 0)
    {
      double r(-b/(2*a));
      return r;
    }
  else
    {
      double r1, r2;
      r1 = (-b - sqrt(delta))/(2*a);
      r2 = (-b + sqrt(delta))/(2*a);
      return r2;
    }
}
