#include "DataFile.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <regex>

DataFile::DataFile()
{
}

DataFile::DataFile(const std::string& fileName):
  _fileName(fileName), _initialCondition("none")
{
}

void DataFile::Initialize(const std::string& fileName)
{
  _fileName = fileName;
  _initialCondition = "none";
}

std::string DataFile::cleanLine(std::string &line)
{
  std::string res = line;
  // Remove everything after a possible #
  res = regex_replace(res, std::regex("#.*$"), std::string(""));
  // Replace tabulation(s) by space(s)
  res = regex_replace(res, std::regex("\t"), std::string(" "), std::regex_constants::match_any);
  // Replace multiple spaces by 1 space
  res = regex_replace(res, std::regex("\\s+"), std::string(" "), std::regex_constants::match_any);
  // Remove any leading spaces
  res = regex_replace(res, std::regex("^ *"), std::string(""));

  return res;
}

// Lit le fichier de paramètres ligne par ligne et affecte la valeur
// adéquate à chaque paramètre. Des vérifications sont faites et la valeur
// de certain paramètres peut changer selon celles d'autres paramètres. Dans ce
// cas là, un warning est affiché sur le terminal.
void DataFile::readDataFile()
{
  // Open the data file
  std::ifstream dataFile(_fileName.data());
  if (!dataFile.is_open())
    {
      std::cout <<  "ERROR::DATAFILE : Unable to open file " << _fileName << std::endl;
      exit(-1);
    }
  // Pour stocker chaque ligne
  std::string line;
  // Run through the dataFile to find the parameters
  while (getline(dataFile, line))
    {
      // Clean line
      std::string proper_line(cleanLine(line));
      if (proper_line.find("ResultsDir") != std::string::npos)
        {
          dataFile >> _resultsDir;
        }
      if (proper_line.find("Friction?") != std::string::npos)
      {
        dataFile >> _friction_yes;
      }
      if (proper_line.find("SaveFinalResultOnly") != std::string::npos)
        {
          dataFile >> _isSaveFinalTimeOnly;
        }
      if (proper_line.find("SaveFrequency") != std::string::npos)
        {
          dataFile >> _saveFrequency;
        }
      if (proper_line.find("IsTestCase") != std::string::npos)
        {
          dataFile >> _isTestCase;
        }
      if (proper_line.find("WhichTestCase") != std::string::npos)
        {
          dataFile >> _testCase;
        }
      if (proper_line.find("InitialCondition") != std::string::npos)
        {
          dataFile >> _initialCondition;
          std::cout << _initialCondition << '\n';
        }
      if (proper_line.find("InitialHeight") != std::string::npos)
        {
          dataFile >> _initialHeight;
        }
      if (proper_line.find("InitialDischarge") != std::string::npos)
        {
          dataFile >> _initialDischarge;
        }
      if (proper_line.find("xmin") != std::string::npos)
        {
          dataFile >> _xmin;
        }
      if (proper_line.find("xmax") != std::string::npos)
        {
          dataFile >> _xmax;
        }
      if (proper_line.find("dx") != std::string::npos)
        {
          dataFile >> _dx;
        }
      if (proper_line.find("Order") != std::string::npos)
        {
          dataFile >> _schemeOrder;
        }
      if (proper_line.find("InitialTime") != std::string::npos)
        {
          dataFile >> _initialTime;
        }
      if (proper_line.find("FinalTime") != std::string::npos)
        {
          dataFile >> _finalTime;
        }
      if (proper_line.find("TimeStep") != std::string::npos)
        {
          dataFile >> _timeStep;
        }
      if (proper_line.find("CFL") != std::string::npos)
        {
          dataFile >> _CFL;
        }
      if (proper_line.find("GravityAcceleration") != std::string::npos)
        {
          dataFile >> _g;
        }
      if (proper_line.find("eta") != std::string::npos)
        {
          dataFile >> _eta;
        }
      if (proper_line.find("kappa") != std::string::npos)
        {
          dataFile >> _kappa;
        }
      if (proper_line.find("LeftBoundaryCondition") != std::string::npos)
        {
          dataFile >> _leftBC;
        }
      if (proper_line.find("RightBoundaryCondition") != std::string::npos)
        {
          dataFile >> _rightBC;
        }
      if (proper_line.find("LeftBoundaryImposedHeight") != std::string::npos)
        {
          dataFile >> _leftBCImposedHeight;
        }
      if (proper_line.find("LeftBoundaryImposedDischarge") != std::string::npos)
        {
          dataFile >> _leftBCImposedDischarge;
        }
      if (proper_line.find("RightBoundaryImposedHeight") != std::string::npos)
        {
          dataFile >> _rightBCImposedHeight;
        }
      if (proper_line.find("RightBoundaryImposedDischarge") != std::string::npos)
        {
          dataFile >> _rightBCImposedDischarge;
        }
      if (proper_line.find("IsTopography") != std::string::npos)
        {
          dataFile >> _isTopography;
        }
      if (proper_line.find("TopographyType") != std::string::npos)
        {
          dataFile >> _topographyType;
        }
    }

  // Making a teporary directory in which to copy the initFile
  system("mkdir -p ./temp");
  system(("cp -T " + _initFile + " ./temp/initial_condition.txt 2> /dev/null").c_str());
  _initFile = "temp/initial_condition.txt";

  // Création et nettoyage du dossier de résultats


  system(("mkdir -p ./" +_resultsDir).c_str());
  system(("rm -f ./" +_resultsDir + "/solution*").c_str());
  system(("cp -r ./" + _fileName + " ./" + _resultsDir + "/parameters.txt").c_str());


  // Si pas de topo --> impose un fond plat
  if (!_isTopography)
    {
      _topographyType = "FlatBottom";
    }
  // Si pas de cas test
  if (!_isTestCase)
    {
      _testCase = "None";
    }


  _Nx = int(ceil((_xmax - _xmin)/_dx));
  _dx = (_xmax - _xmin)/_Nx;

}
