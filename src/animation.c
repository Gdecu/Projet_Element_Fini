/*
 *  main.c
 *  Library for EPL1110 : Finite Elements for dummies
 *  Utilisation de l'API de GMSH pour créer un maillage
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 */

 #include "glfem.h"
 #include "animation.h"
 #include "fem.h"
 #include "utils.h"
 #define M_PI 3.14159265358979323846

 
 void run_animation()
 {
     printf("\n\n    V : Mesh and size mesh field \n");
     printf("    D : Domains \n");
     printf("    N : Next domain highlighted\n");
 
     // Tableau des hauteurs pour l'animation
     double *heigts = malloc(6 * sizeof(double));
     heigts[0] = 1;
     heigts[1] = 10;
     heigts[2] = 20;
     heigts[3] = 30;
     heigts[4] = 40;
     heigts[5] = 50;
 
     // Allocation initiale des pointeurs (déclarés ici pour usage global dans main)
     double *normDisplacement = NULL;
     double *forcesX = NULL;
     double *forcesY = NULL;
 
     // Initialisation géométrie
     
     femGeo *theGeometry = geoGetGeometry();
     theGeometry->r_in = 2;
     theGeometry->r_out = 2.5;
     theGeometry->h = 0.1;
     theGeometry->height = 1;
     theGeometry->angle = M_PI / 2;
 
     // Paramètres physiques
     double E = 35e9;
     double nu = 0.2;
     double rho = 2.3e3;
     double g = 9.81;
     double deformationFactor = 1.5e2;
 
     
 
     // Variables d'animation
     int currentFrame = 0;
     int numFrames = 6;
     double frameTime = 1.0;
     double tOld = glfwGetTime() - frameTime;
 
     // Variables affichage
     int mode = 1;
     int domain = 0;
     int freezingButton = FALSE;
     double t, told = 0;
     char theMessage[MAXNAME];
 
    
     // Initialisation OpenGL
     GLFWwindow *window = glfemInit("EPL1110 : Recovering forces on constrained nodes");
     glfwMakeContextCurrent(window);
 
    for(int i = 0; i < numFrames; i++) {
        geoInitialize();
        femProblem *theProblem = NULL;
        theGeometry->height = heigts[i];
        printf("Frame %d\n", i);
        // Génération du maillage
        glClearColor(1.0, 1.0, 1.0, 1.0); // Fond blanc
        glClear(GL_COLOR_BUFFER_BIT);


        geoMeshGenerate();
        geoMeshImport();
        printf("Number of detected domains: %d\n", theGeometry->nDomains);
 
        // Définition des noms de domaine
        geoSetDomainName(0, "Outer Disk");
        geoSetDomainName(1, "Right");
        geoSetDomainName(2, "Inner Disk");
        geoSetDomainName(3, "Left");
 
        // Écriture du maillage dans un fichier
        char filename[100];
        sprintf(filename, "../data/mesh_%d.txt", i);
        geoMeshWrite(filename);
 
        // Création du problème FEM
        theProblem = femElasticityCreate(theGeometry, E, nu, rho, g, PLANAR_STRAIN);
 
        // Application des conditions aux limites
        double Force = rho * g * theGeometry->height * theGeometry->r_out * theGeometry->angle;
        femElasticityAddBoundaryCondition(theProblem,"Left",DIRICHLET_X,0.0);
        femElasticityAddBoundaryCondition(theProblem,"Right",DIRICHLET_X,0.0);
        femElasticityAddBoundaryCondition(theProblem,"Left",DIRICHLET_Y,0.0);
        femElasticityAddBoundaryCondition(theProblem,"Right",DIRICHLET_Y,0.0);
        femElasticityAddBoundaryCondition(theProblem,"Outer Disk",NEUMANN_Y,-Force);
 
        double *theSoluce = femElasticitySolve(theProblem);
        double *theForces = femElasticityForces(theProblem);
        double area = femElasticityIntegrate(theProblem, fun);  
        
        femNodes *theNodes = theGeometry->theNodes;
        normDisplacement = malloc(theNodes->nNodes * sizeof(double));
        forcesX = malloc(theNodes->nNodes * sizeof(double));
        forcesY = malloc(theNodes->nNodes * sizeof(double));

 
        for (int i=0; i<theNodes->nNodes; i++){
            theNodes->X[i] += theSoluce[2*i+0]*deformationFactor;
            theNodes->Y[i] += theSoluce[2*i+1]*deformationFactor;
            normDisplacement[i] = sqrt(theSoluce[2*i+0]*theSoluce[2*i+0] + 
                                       theSoluce[2*i+1]*theSoluce[2*i+1]);
            forcesX[i] = theForces[2*i+0];
            forcesY[i] = theForces[2*i+1]; }
        

        int w,h;
        glfwGetFramebufferSize(window,&w,&h);
        glfemReshapeWindows(theGeometry->theNodes,w,h);
        glfemPlotField(theGeometry->theElements,normDisplacement);
        glfemPlotMesh(theGeometry->theElements); 
        sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
        glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); 
        sprintf(theMessage, "Height : %.2f", theGeometry->height);
        glColor3f(0.0, 0.0, 1.0); // Bleu
        glfemDrawMessage(20, 40, theMessage);


        glfwSwapBuffers(window);
        glfwPollEvents();

        double tStart = glfwGetTime();
        while (glfwGetTime() - tStart < 1.0) {
            glfwPollEvents();
        }

        
        
        free(normDisplacement); normDisplacement = NULL;
        free(forcesX);          forcesX = NULL;
        free(forcesY);          forcesY = NULL;

        //free(theSoluce);
        //free(theForces);

        femElasticityFree(theProblem);
        geoFinalize();

        theProblem = NULL;
        currentFrame++;
        

            }

        while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS &&
            glfwWindowShouldClose(window) != 1) {
            glfwPollEvents();  // ← essentiel pour détecter la fermeture
     }
     
    free(heigts);
    glfwTerminate();
 }
 