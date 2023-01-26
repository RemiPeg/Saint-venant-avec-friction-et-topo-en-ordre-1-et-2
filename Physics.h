#ifndef PHYSICS_H
#define PHYSICS_H

#include "DataFile.h"
#include "Mesh.h"
#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Sparse"



class Physics
{
private:
  // Pointeur vers le fichier de paramètres pour récupérer les
  // conditions initiales, aux limites et le fichier de topographie (terme source).
  DataFile* _DF;
  Mesh* _mesh;

  // Variables pratiques
  double _xmin, _xmax;
  double _g;
  int _nCells;

  int _i;

  // Condition initiale
  Eigen::Matrix<double, Eigen::Dynamic, 2> _Sol0;

  // Topographie pour le terme de topographie.
  Eigen::VectorXd _topography;


  // Exact solution
  Eigen::Matrix<double, Eigen::Dynamic, 2> _exactSol;




public:
  // Constructeur
  Physics();
  Physics(DataFile* DF, Mesh* mesh);

  // Initialisation
  void Initialize();
  void Initialize(DataFile* DF, Mesh* mesh);

  // Getters
  const Eigen::Matrix<double, Eigen::Dynamic, 2>& getInitialCondition() const {return _Sol0;};
  const Eigen::VectorXd& getTopography() const {return _topography;};
  const Eigen::Matrix<double, Eigen::Dynamic, 2>& getExactSolution() const {return _exactSol;};

  // Construit/Sauvegarde la solution exacte
  void buildExactSolution(double t);
  void saveExactSolution(std::string& fileName) const;



      void computeWaveSpeed(const Eigen::Vector2d& SolG, const Eigen::Vector2d& SolD, double* amin, double* amax) const;
      Eigen::Vector2d G(const Eigen::Vector2d& Sol) const;

  // Conditions aux limites
  Eigen::Vector2d leftBoundaryFunction(double t, const Eigen::Matrix<double, Eigen::Dynamic, 2>& Sol);
  Eigen::Vector2d rightBoundaryFunction(double t, const Eigen::Matrix<double, Eigen::Dynamic, 2>& Sol);

protected:
  void buildTopography();
  void buildInitialCondition();
  // Boundary conditions

  // Resolution equation second ordre
  double FindRacine(double a, double b, double c);
  // Exact solution

};

#endif // PHYSICS_H
