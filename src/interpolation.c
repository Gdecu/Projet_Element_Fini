/*
 *  main.c
 *  Library for EPL1110 : Finite Elements for dummies
 *  Utilisation de l'API de GMSH pour cr�er un maillage
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */
 
#include "glfem.h"
#include "interpolation.h"
#include "utils.h"
#include "fem.h"
#define M_PI 3.14159265358979323846


void run_interpolation()
{  
    /*
    FILE *outputFile = fopen("../data/problem.txt", "w");
    if (outputFile == NULL) {
        perror("Erreur lors de l'ouverture du fichier problem.txt");
        exit(EXIT_FAILURE);
    }
    freopen("../data/problem.txt", "w", stdout);

    fprintf(outputFile, "Analyser la déformation d'une section horizontale d'un barrage hydraulique sous l'effet de la pression hydrostatique de l'eau.\n");
    fprintf(outputFile, "Dans le cadre du projet du cours LEPL1100 éléments finis.\n");
    fprintf(outputFile, "Le maillage est généré à l'aide de GMSH et le problème est résolu avec la bibliothèque FEM.\n");
    */
      
    int ierr;
    
    geoInitialize();
    femGeo* theGeometry = geoGetGeometry();
    theGeometry->r_in = 50;
    theGeometry->r_out = 62.5;
    theGeometry->h = 2.5;
    theGeometry->h_in = 12.5;
    theGeometry->h_out = 1.25;
    theGeometry->height = 50;
    theGeometry->angle = M_PI/2;
   
    geoMeshGenerate();
    geoMeshImport();
    printf("Number of detected domains: %d\n", theGeometry->nDomains);

    
    geoSetDomainName(0,"Outer Disk");
    geoSetDomainName(1,"Right");
    geoSetDomainName(2,"Inner Disk");
    geoSetDomainName(3,"Left");

    char filename[] = "../data/mesh.txt";
    geoMeshWrite(filename);

    
//
// 1.1 création du probleme
//

double E   = 35e9;
double nu  = 0.2;
double rho = 2.3e3; 
double g   = 9.81;
double Force = rho * g * theGeometry->height + 101325; // 101325 Pa = 1 atm
printf(" ==== Force applied on the disk : %14.7e [N/m] \n",Force);
femProblem* theProblem = femElasticityCreate(theGeometry,E,nu,rho,g,PLANAR_STRAIN);
femElasticityAddBoundaryCondition(theProblem,"Left",DIRICHLET_X,0.0);
femElasticityAddBoundaryCondition(theProblem,"Right",DIRICHLET_X,0.0);
femElasticityAddBoundaryCondition(theProblem,"Left",DIRICHLET_Y,0.0);
femElasticityAddBoundaryCondition(theProblem,"Right",DIRICHLET_Y,0.0);
//femElasticityAddBoundaryCondition(theProblem,"Inner Disk",DIRICHLET_Y,0.0);
femElasticityAddBoundaryCondition(theProblem,"Outer Disk",NEUMANN_Y,-Force);
femElasticityPrint(theProblem);


//
//  -3- Resolution du probleme et calcul des forces
//

double *theSoluce = femElasticitySolve(theProblem);
double *theForces = femElasticityForces(theProblem);
double area = femElasticityIntegrate(theProblem, fun);   


//
//  -4- Deformation du maillage pour le plot final
//      Creation du champ de la norme du deplacement
//
    
femNodes *theNodes = theGeometry->theNodes;
double deformationFactor = 1.5e3;
double *normDisplacement = malloc(theNodes->nNodes * sizeof(double));
double *forcesX = malloc(theNodes->nNodes * sizeof(double));
double *forcesY = malloc(theNodes->nNodes * sizeof(double));

for (int i=0; i<theNodes->nNodes; i++){
    theNodes->X[i] += theSoluce[2*i+0]*deformationFactor;
    theNodes->Y[i] += theSoluce[2*i+1]*deformationFactor;
    normDisplacement[i] = sqrt(theSoluce[2*i+0]*theSoluce[2*i+0] + 
                               theSoluce[2*i+1]*theSoluce[2*i+1]);
    forcesX[i] = theForces[2*i+0];
    forcesY[i] = theForces[2*i+1]; }

double hMin = femMin(normDisplacement,theNodes->nNodes);  
double hMax = femMax(normDisplacement,theNodes->nNodes);  
printf(" ==== Minimum displacement          : %14.7e [m] \n", hMin);
printf(" ==== Maximum displacement          : %14.7e [m] \n", hMax);

//
//  -5- Calcul de la force globaleresultante
//

double theGlobalForce[2] = {0, 0};
for (int i=0; i<theProblem->geometry->theNodes->nNodes; i++) {
    theGlobalForce[0] += theForces[2*i+0];
    theGlobalForce[1] += theForces[2*i+1]; }
printf(" ==== Global horizontal force       : %14.7e [N] \n", theGlobalForce[0]);
printf(" ==== Global vertical force         : %14.7e [N] \n", theGlobalForce[1]);
printf(" ==== Weight                        : %14.7e [N] \n", area * rho * g);

//
//  -6- Visualisation du maillage
//  

int mode = 1; 
int domain = 0;
int freezingButton = FALSE;
double t, told = 0;
char theMessage[MAXNAME];


GLFWwindow* window = glfemInit("EPL1110 : Recovering forces on constrained nodes");
glfwMakeContextCurrent(window);

do {
    int w,h;
    glfwGetFramebufferSize(window,&w,&h);
    glfemReshapeWindows(theGeometry->theNodes,w,h);

    t = glfwGetTime();  
    if (glfwGetKey(window,'D') == GLFW_PRESS) { mode = 0;}
    if (glfwGetKey(window,'V') == GLFW_PRESS) { mode = 1;}
    if (glfwGetKey(window,'X') == GLFW_PRESS) { mode = 2;}
    if (glfwGetKey(window,'Y') == GLFW_PRESS) { mode = 3;}
    if (glfwGetKey(window,'N') == GLFW_PRESS && freezingButton == FALSE) { domain++; freezingButton = TRUE; told = t;}
    if (t-told > 0.5) {freezingButton = FALSE; }
    
    if (mode == 0) {
        domain = domain % theGeometry->nDomains;
        glfemPlotDomain( theGeometry->theDomains[domain]); 
        sprintf(theMessage, "%s : %d ",theGeometry->theDomains[domain]->name,domain);
        glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
    if (mode == 1) {
        glfemPlotField(theGeometry->theElements,normDisplacement);
        glfemPlotMesh(theGeometry->theElements); 
        sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
        glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
    if (mode == 2) {
        glfemPlotField(theGeometry->theElements,forcesX);
        glfemPlotMesh(theGeometry->theElements); 
        sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
        glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
    if (mode == 3) {
        glfemPlotField(theGeometry->theElements,forcesY);
        glfemPlotMesh(theGeometry->theElements); 
        sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
        glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
     glfwSwapBuffers(window);
     glfwPollEvents();
} while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
         glfwWindowShouldClose(window) != 1 );
        
// Check if the ESC key was pressed or the window was closed

/*
fclose(outputFile);
freopen("/dev/tty", "w", stdout); // Restaure stdout vers le terminal
*/

free(normDisplacement);
free(forcesX);
free(forcesY);
femElasticityFree(theProblem) ; 
geoFinalize();
glfwTerminate(); 

exit(EXIT_SUCCESS);
return 0; 
}