#include "TimeScheme.h"
#include "DataFile.h"
#include "Mesh.h"
#include "Physics.h"
#include "FiniteVolume.h"

#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Sparse"

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>



//----------------------------------------------------------//
//------------------Time Scheme base class------------------//
//----------------------------------------------------------//
TimeScheme::TimeScheme()
{
}



TimeScheme::TimeScheme(DataFile* DF, Mesh* mesh, Physics* physics, FiniteVolume* finVol):
  _DF(DF), _mesh(mesh), _physics(physics), _finVol(finVol), _Sol(_physics->getInitialCondition()), _timeStep(DF->getTimeStep()), _initialTime(DF->getInitialTime()), _finalTime(DF->getFinalTime()), _currentTime(_initialTime)
{
}



void TimeScheme::Initialize(DataFile* DF, Mesh* mesh, Physics* physics, FiniteVolume* finVol)
{
  _DF = DF;
  _mesh = mesh;
  _physics = physics;
  _finVol = finVol;
  _Sol = _physics->getInitialCondition();
  _timeStep = DF->getTimeStep();
  _initialTime = DF->getInitialTime();
  _finalTime = DF->getFinalTime();
  _currentTime = _initialTime;
}





void TimeScheme::saveCurrentSolution(std::string& fileName) const
{
  std::ofstream outputFile(fileName, std::ios::out);
  const Eigen::VectorXd& cellCenters(_mesh->getCellCenters());
  double g(_DF->getGravityAcceleration());
  // Gnuplot comments for the user
  outputFile << "# x  H=h+z   h       u       q       Fr=|u|/sqrt(gh)" << std::endl;
  for (int i(0) ; i < _Sol.rows() ; ++i)
    {
      outputFile << cellCenters(i) << " " <<
        _Sol(i,0) + _physics->getTopography()(i) << " " <<
        _Sol(i,0) << " " <<
        _Sol(i,1)/_Sol(i,0) << " " <<
        _Sol(i,1) << " " <<
        abs(_Sol(i,1)/_Sol(i,0))/sqrt(g * _Sol(i,0)) << std::endl;
    }
}




void TimeScheme::Solve()
{
  // Variables pratiques
  int n(0);
  double dx(_mesh->getSpaceStep());
  std::string resultsDir(_DF->getResultsDirectory());
  std::string fluxName(_finVol->getFluxName());

  // Sauvegarde la condition initiale
  std::string fileName(resultsDir + "/Solution_" + fluxName + "_" + std::to_string(n) + ".txt");
  saveCurrentSolution(fileName);

  // Sauvegarde la topographie
  std::string topoFileName(resultsDir + "/topography.txt");
  std::ofstream topoFile(topoFileName, std::ios::out);
  for (int i(0) ; i < _Sol.rows() ; ++i)
    {
      topoFile << _mesh->getCellCenters()(i) << " " << _physics->getTopography()(i) << std::endl;
    }

  // Boucle en temps
  while (_currentTime < _finalTime)
    {
      std::cout << (_currentTime < _finalTime) << "          " <<_currentTime << _finalTime<< '\n';
      oneStep(_currentTime);
      _timeStep=dx/(2*_amaxtot);
      n=n+1;
      _currentTime += _timeStep;
      // Save Solution at time t
      if (!_DF->isSaveFinalTimeOnly() &&  n % _DF->getSaveFrequency() == 0)
        {
          std::string fileName(resultsDir + "/Solution_" + fluxName + "_" + std::to_string(n/_DF->getSaveFrequency()) + ".txt");
          saveCurrentSolution(fileName);
        }
    }
  // End of time loop
  if (_DF->isSaveFinalTimeOnly())
    {
      std::string fileName(resultsDir + "/Solution_" + fluxName + "_" + std::to_string(n/_DF->getSaveFrequency()) + ".txt");
      saveCurrentSolution(fileName);
    }
  if (_DF->isTestCase())
    {
      _physics->buildExactSolution(_currentTime);
      std::string fileName(resultsDir + "/Solution_exacte.txt");
      _physics->saveExactSolution(fileName);
      Eigen::Vector2d L2error(computeL2Error());
      std::cout << "Error h  L2 = " << L2error(0) << " and error q L2 = " << L2error(1) << " for dx = " << _DF->getDx() << std::endl;
      Eigen::Vector2d L1error(computeL1Error());
      std::cout << "Error h  L1 = " << L1error(0) << " and error q L1 = " << L1error(1) << " for dx = " << _DF->getDx() << std::endl;
    }
}


Eigen::Vector2d TimeScheme::computeL2Error() const
{
  Eigen::Vector2d error(0., 0.);
  const Eigen::Matrix<double, Eigen::Dynamic, 2>& exactSol(_physics->getExactSolution());
  error(0) = (_Sol.col(0) - exactSol.col(0)).norm();
  error(1) = (_Sol.col(1) - exactSol.col(1)).norm();
  error *= _DF->getDx();
  return error;
}


Eigen::Vector2d TimeScheme::computeL1Error() const
{
  Eigen::Vector2d error(0., 0.);
  const Eigen::Matrix<double, Eigen::Dynamic, 2>& exactSol(_physics->getExactSolution());
  for (int i(0) ; i < _Sol.rows() ; ++i)
    {
      error(0) += abs(_Sol(i,0) - exactSol(i,0));
      error(1) += abs(_Sol(i,1) - exactSol(i,1));
    }
  error *= _DF->getDx();
  return error;
}




//-------------------------------------------------//
//------------------SV sans friction---------------//
//-------------------------------------------------//
SV::SV():
  TimeScheme()
{
}



SV::SV(DataFile* DF, Mesh* mesh, Physics* physics, FiniteVolume* finVol):
  TimeScheme(DF, mesh, physics, finVol)
{
}



void SV::Initialize(DataFile* DF, Mesh* mesh, Physics* physics, FiniteVolume* finVol)
{
  _DF = DF;
  _mesh = mesh;
  _physics = physics;
  _finVol = finVol;
  _Sol.resize(mesh->getNumberOfCells(), 2);
  _timeStep = DF->getTimeStep();
  _initialTime = DF->getInitialTime();
  _finalTime = DF->getFinalTime();
  _currentTime = _initialTime;
  _g = DF->getGravityAcceleration();
}



void SV::oneStep(double t)
{
  // Récupération des trucs importants
  double dt(_timeStep);
  double dx(_mesh->getSpaceStep());
  double aminG(0.);
  double amaxG(0.);
  double aminD(0.);
  double amaxD(0.);
  double amaxtot(0.);
  int i(0.);

  Eigen::Vector2d GidemiD;
  Eigen::Vector2d GidemiG;
  Eigen::Vector2d SolG;
  Eigen::Vector2d SolD;
  Eigen::Vector2d SolC;
  Eigen::Vector2d SolBC;

  for (int i(1) ; i < _Sol.rows()-2 ; ++i)
    {
      SolC(0)=_Sol(i,0);
      SolG(0)=_Sol(i-1,0);
      SolD(0)=_Sol(i+1,0);
      SolC(1)=_Sol(i,1);
      SolG(1)=_Sol(i-1,1);
      SolD(1)=_Sol(i+1,1);
      std::cout << _currentTime << '\n';
      GidemiG=_finVol->buildGidemi(SolC,SolG) ;
      GidemiD=_finVol->buildGidemi(SolD,SolC);
      _Sol(i,0) = SolC(0)-dt*(GidemiD(0)-GidemiG(0))/dx;
      _Sol(i,1) = SolC(1)-dt*(GidemiD(1)-GidemiG(1))/dx;
    }
}


//-------------------------------------------------//
//------------------SV avec friction---------------//
//-------------------------------------------------//
SV_f::SV_f():
  TimeScheme()
{
}



SV_f::SV_f(DataFile* DF, Mesh* mesh, Physics* physics, FiniteVolume* finVol):
  TimeScheme(DF, mesh, physics, finVol)
{
}



void SV_f::Initialize(DataFile* DF, Mesh* mesh, Physics* physics, FiniteVolume* finVol)
{
  _DF = DF;
  _mesh = mesh;
  _physics = physics;
  _finVol = finVol;
  _Sol.resize(mesh->getNumberOfCells(), 2);
  _timeStep = DF->getTimeStep();
  _initialTime = DF->getInitialTime();
  _finalTime = DF->getFinalTime();
  _currentTime = _initialTime;
  _eta=DF->getEta();
  _kappa=DF->getKappa();
  _g = DF->getGravityAcceleration();
}



void SV_f::oneStep(double t)
{
  // Récupération des trucs importants
  double dt(_timeStep);
  double dx(_mesh->getSpaceStep());
  double alphaG;
  double alphaD;
  double aminG(0.);
  double amaxG(0.);
  double aminD(0.);
  double amaxD(0.);
  double amaxtot(0.);
  int i(0.);

  Eigen::Vector2d GidemiD;
  Eigen::Vector2d GidemiG;
  Eigen::Vector2d SolG;
  Eigen::Vector2d SolD;
  Eigen::Vector2d SolC;
  Eigen::Vector2d SolBC;

  for (int i(1) ; i < _Sol.rows()-2 ; ++i)
    {
      SolC(0)=_Sol(i,0);
      SolG(0)=_Sol(i-1,0);
      SolD(0)=_Sol(i+1,0);
      SolC(1)=_Sol(i,1);
      SolG(1)=_Sol(i-1,1);
      SolD(1)=_Sol(i+1,1);
      std::cout << _currentTime << '\n';
      GidemiG=_finVol->buildGidemi(SolC,SolG) ;
      alphaG=pow(SolG(0)+SolC(0),_eta)*(amaxG-aminG)/(pow(SolG(0)+SolC(0),_eta)*(amaxG-aminG)+_kappa*dx);
      GidemiD=_finVol->buildGidemi(SolD,SolC);
      alphaD=pow(SolC(0)+SolD(0),_eta)*(amaxD-aminD)/(pow(SolC(0)+SolD(0),_eta)*(amaxD-aminD)+_kappa*dx);

      _Sol(i,0) = SolC(0)-dt*(GidemiD(0)-GidemiG(0))/dx;
      _Sol(i,1) = SolC(1)-dt*(alphaD*GidemiD(1)-alphaG*GidemiG(1))/dx-dt*((1-alphaG)*amaxG/dx+(1-alphaD)*aminG/dx)*abs(SolC(1))*SolC(1)+dt*(alphaG-alphaD)*(pow(SolC(1),2)/SolC(0)+_g*pow(SolC(0),2)/2)/dx;
    }
}
