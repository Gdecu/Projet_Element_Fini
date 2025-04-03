#ifndef UTILS_H
#define UTILS_H

#include "fem.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Interpolation
double hermiteInterpolation(double d, double d_star, double h0, double h_star);
double geoSize(double x, double y);

// Maillage
void geoMeshGenerate(void);

// Assemblage du système
void femElasticityAssembleElements(femProblem *theProblem);
void femElasticityAssembleNeumann(femProblem *theProblem);

// Résolution FEM
double *femElasticitySolve(femProblem *theProblem);
double *femElasticityForces(femProblem *theProblem);

// Fonction utilitaire d’intégration
double fun(double x, double y);

#endif
