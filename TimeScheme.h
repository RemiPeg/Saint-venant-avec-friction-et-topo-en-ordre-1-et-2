#ifndef TIME_SCHEME_H
#define TIME_SCHEME_H

#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Sparse"
#include "DataFile.h"
#include "Mesh.h"
#include "Physics.h"
#include "FiniteVolume.h"

#include <vector>



class TimeScheme
{
protected:
  // Pointeur vers les trucs importants
  DataFile* _DF;
  Mesh* _mesh;
  Physics* _physics;
  FiniteVolume* _finVol;

  // Vecteur solution
  Eigen::Matrix<double, Eigen::Dynamic, 2> _Sol;

  // Paramètres de temps
  double _timeStep;
  double _initialTime;
  double _finalTime;
  double _currentTime;

  double _g;
  double _eta;
  double _kappa;
  double _amaxtot;

public:
  // Constructeurs
  TimeScheme();
  TimeScheme(DataFile* DF, Mesh* mesh, Physics* physics, FiniteVolume* finVol);

  // Initialiseur
  void Initialize(DataFile* DF, Mesh* mesh, Physics* physics, FiniteVolume* finVol);
  // Destructeur
  virtual ~TimeScheme() = default;

  // Getters
  const Eigen::Matrix<double, Eigen::Dynamic, 2>& getSolution() const {return _Sol;};
  double getTimeStep() const {return _timeStep;};
  double getInitialTime() const {return _initialTime;};
  double getFinalTime() const {return _finalTime;};
  double getCurrentTime() const {return _currentTime;};
  double getamaxtot() const {return _amaxtot;};


  // Solve and save solution
  virtual void oneStep(double t) = 0;
  void saveCurrentSolution(std::string& fileName) const;
  void Solve();

  // Error
  Eigen::Vector2d computeL2Error() const;
  Eigen::Vector2d computeL1Error() const;
};




class SV: public TimeScheme
{
public:
  // Constructeurs
  SV();
  SV(DataFile* DF, Mesh* mesh, Physics* physics, FiniteVolume* finVol);

  // Initialiseur
  void Initialize(DataFile* DF, Mesh* mesh, Physics* physics, FiniteVolume* finVol);

  // One time step
  void oneStep(double t);
};



class SV_f: public TimeScheme
{
public:
  // Constructeurs
  SV_f();
  SV_f(DataFile* DF, Mesh* mesh, Physics* physics, FiniteVolume* finVol);

  // Initialiseur
  void Initialize(DataFile* DF, Mesh* mesh, Physics* physics, FiniteVolume* finVol);

  // One time step
  void oneStep(double t);
};



#endif // TIME_SCHEME_H
