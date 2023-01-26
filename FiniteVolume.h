#ifndef FINITE_VOLUME_H
#define FINITE_VOLUME_H

#include "DataFile.h"
#include "Mesh.h"
#include "Physics.h"

#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Sparse"



class FiniteVolume
{
protected:
  // Pointeurs vers les trucs importants
  DataFile* _DF;
  Mesh* _mesh;
  Physics* _physics;

  // Nom du flux numérique
  std::string _fluxName;

  // Vecteur des flux
  Eigen::Matrix<double, Eigen::Dynamic, 2> _fluxVector;

public:
  // Constructeurs
  FiniteVolume();
  FiniteVolume(DataFile* DF, Mesh* mesh, Physics* physics);

  // Destructeur
  virtual ~FiniteVolume() = default;

  // Initialisation
  void Initialize(DataFile* DF, Mesh* mesh, Physics* physics);

  // Getters
  const std::string& getFluxName() const {return _fluxName;};
  const Eigen::Matrix<double, Eigen::Dynamic, 2>& getFluxVector() const {return _fluxVector;};

  // Build the flux vector

    // Build flux vector
    Eigen::Vector2d buildGidemi(const Eigen::Vector2d& SolD,const Eigen::Vector2d& SolG);


  Eigen::Matrix<double, Eigen::Dynamic, 2> buildGidemi2(double t, Eigen::Matrix<double, Eigen::Dynamic, 2> Sol);


protected:
  // Minmod slope limiter for the 2nd order MUSCL schemes
  double minmod(double r) const;
};


#endif //FINITE_VOLUME_H
