#include "fem.h"
#include <math.h>
#include "utils.h"
#include "flag.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double hermiteInterpolation(double d, double d_star, double h0, double h_star) {
    if (d >= d_star) return h_star;
    double t = d / d_star;
    return h0 + (h_star - h0) * (3 * t * t - 2 * t * t * t);
}

double geoSize(double x, double y) {
    //printf("geoSize called with x=%f, y=%f\n", x, y);
    femGeo* theGeometry = geoGetGeometry();
    if (theGeometry == NULL) {
        printf("geoSize: ERROR – geometry not initialized!\n");
        return 1.0;  // Valeur par défaut sécurisée
    }

    double h = theGeometry->h;

    if (flag == 0) {
        // Interpolation pour le maillage
        return h;
    }


    double r_in = theGeometry->r_in;
    double r_out = theGeometry->r_out;
    double angle_max = theGeometry->angle / 2;

    double damWidth = fabs(r_out - r_in);
    double h0 = theGeometry->h_in;
    double h1 = theGeometry->h_out;


    double r = sqrt( x * x + y * y);
    double angle = atan2(y,x);
    
    double dR_in = fabs(r - r_in);
    double dR_out = fabs(r_out - r);
    //double dRight = fabs(r * sin(angle_max - angle));
    //double dLeft = fabs(r * sin ( - angle_max + angle));


    // Taille interpolée autour de l'encoche et du trou
    double hIn = (dR_in < damWidth) ? hermiteInterpolation(dR_in, damWidth, h0, h) : h;
    double hOut = (dR_out < damWidth) ? hermiteInterpolation(dR_out, damWidth, h1, h) : h;

    double finalSize = fmin(h, fmin(hIn, hOut));
    
    return finalSize;
}

void geoMeshGenerate() {
    femGeo* theGeometry = geoGetGeometry();

    double ri = theGeometry->r_in;
    double ro = theGeometry->r_out;
    double h0 = theGeometry->h;
    double angle = theGeometry->angle;

    
    int ierr;

    // 1. Création du grand disque
    int idBigDisk = gmshModelOccAddDisk(0, 0, 0, ro, ro, -1, NULL, 0, NULL, 0, &ierr);
    ErrorGmsh(ierr);

    //2.1 Création du disque à retirer
    int idSmallDisk = gmshModelOccAddDisk(0, 0, 0, ri, ri, -1, NULL, 0, NULL, 0, &ierr);
    ErrorGmsh(ierr);

    int bigDisk[] = {2, idBigDisk};
    int smallDisk[] = {2, idSmallDisk};

    //2.2 Soustraction du disque à retirer
    gmshModelOccCut(bigDisk,2,smallDisk,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr); 
    ErrorGmsh(ierr);

    // 2.3 Définition des points de découpe sur le contour du disque extérieur
    int idPoint1, idPoint2;
    idPoint1 = gmshModelOccAddPoint(ro * cos((M_PI - angle)*0.5), ro * sin((M_PI - angle)*0.5), 0, h0, -1, &ierr);
    ErrorGmsh(ierr);
    idPoint2 = gmshModelOccAddPoint(ro * cos((M_PI + angle)*0.5), ro * sin((M_PI + angle)*0.5), 0, h0, -1, &ierr);
    ErrorGmsh(ierr);

    // 2.4 Création des lignes de découpe partant du centre vers les points définis
    int idOrigin = gmshModelOccAddPoint(0, 0, 0, h0, -1, &ierr);
    ErrorGmsh(ierr);

    int idLine1 = gmshModelOccAddLine(idOrigin, idPoint1, -1, &ierr);
    int idLine2 = gmshModelOccAddLine(idOrigin, idPoint2, -1, &ierr);
    ErrorGmsh(ierr);

    // Ajout d'un arc reliant idPoint1 à idPoint2
    int idArc = gmshModelOccAddCircleArc(idPoint1, idOrigin, idPoint2, -1, idOrigin, &ierr);
    ErrorGmsh(ierr);

    // 2.5 Création d'une surface triangulaire reliant les points
    int curveLoop[] = {idLine1, idArc, idLine2};
    int idSurface = gmshModelOccAddCurveLoop(curveLoop, 3, -1, &ierr);
    ErrorGmsh(ierr);

    int idPlaneSurface = gmshModelOccAddPlaneSurface(&idSurface, 1, -1, &ierr);
    ErrorGmsh(ierr);

    // 2.6 Découper le disque pour ne garder que le quart supérieur
    int cuttingSurface[] = {2, idPlaneSurface};

    gmshModelOccIntersect(bigDisk, 2, cuttingSurface, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    ErrorGmsh(ierr);

    

    // 4. Synchronisation et génération du maillage
    geoSetSizeCallback(geoSize);
    gmshModelOccSynchronize(&ierr);
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshModelMeshGenerate(2, &ierr);

}

void femElasticityAssembleElements(femProblem *theProblem){
    femFullSystem  *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete    *theSpace = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theMesh = theGeometry->theElements;
    femMesh        *theEdges = theGeometry->theEdges;
    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
    int iElem,iInteg,iEdge,i,j,d,map[4],mapX[4],mapY[4];
    int nLocal = theMesh->nLocalNode;
    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;      
    double rho = theProblem->rho;
    double g   = theProblem->g;
    double **A = theSystem->A;
    double *B  = theSystem->B;
    
    
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j=0; j < nLocal; j++) {
            map[j]  = theMesh->elem[iElem*nLocal+j];
            mapX[j] = 2*map[j];
            mapY[j] = 2*map[j] + 1;
            x[j]    = theNodes->X[map[j]];
            y[j]    = theNodes->Y[map[j]];} 
        
        for (iInteg=0; iInteg < theRule->n; iInteg++) {    
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];  
            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
            
            double dxdxsi = 0.0;
            double dxdeta = 0.0;
            double dydxsi = 0.0; 
            double dydeta = 0.0;
            for (i = 0; i < theSpace->n; i++) {  
                dxdxsi += x[i]*dphidxsi[i];       
                dxdeta += x[i]*dphideta[i];   
                dydxsi += y[i]*dphidxsi[i];   
                dydeta += y[i]*dphideta[i]; }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            
            for (i = 0; i < theSpace->n; i++) {    
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; }            
            for (i = 0; i < theSpace->n; i++) { 
                for(j = 0; j < theSpace->n; j++) {
                    A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + 
                                            dphidy[i] * c * dphidy[j]) * jac * weight;                                                                                            
                    A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + 
                                            dphidy[i] * c * dphidx[j]) * jac * weight;                                                                                           
                    A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + 
                                            dphidx[i] * c * dphidy[j]) * jac * weight;                                                                                            
                    A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + 
                                            dphidx[i] * c * dphidx[j]) * jac * weight; }}
             for (i = 0; i < theSpace->n; i++) {
                B[mapY[i]] -= phi[i] * g * rho * jac * weight; }}} 
}


void femElasticityAssembleNeumann(femProblem *theProblem){
    
    femFullSystem  *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->ruleEdge;
    femDiscrete    *theSpace = theProblem->spaceEdge;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theEdges = theGeometry->theEdges;

    int iBnd, iElem, iInteg, i, j, map[2], mapU[2];
    double x[2], y[2], phi[2];
    double *B  = theSystem->B;

    // Parcours des conditions aux limites
    for(iBnd = 0; iBnd < theProblem->nBoundaryConditions; iBnd++){
        
        femBoundaryCondition *theCondition = theProblem->conditions[iBnd];

        if (theCondition->type != NEUMANN_X && theCondition->type != NEUMANN_Y) continue;

        double value = theCondition->value;
        femDomain *domain = theCondition->domain;


        for(iElem = 0; iElem < domain->nElem; iElem++){
            
            int edge = domain->elem[iElem];
            for (j = 0; j < 2; j++){
                map[j] = theEdges->elem[edge*2+j];
                x[j]   = theNodes->X[map[j]];
                y[j]   = theNodes->Y[map[j]];
                mapU[j] = (theCondition->type == NEUMANN_X) ? 2*map[j] : 2*map[j] + 1;
            }

            double length = sqrt((x[1]-x[0])*(x[1]-x[0]) + (y[1]-y[0])*(y[1]-y[0]));
            double jacobian = length / 2.0;

            for(iInteg = 0; iInteg < theRule->n; iInteg++){

                double xsi    = theRule->xsi[iInteg];
                double weight = theRule->weight[iInteg];

                femDiscretePhi(theSpace, xsi, phi);

                for(i = 0; i < 2; i++){
                    B[mapU[i]] += phi[i] * value * jacobian * weight;
                }
            }
        }
    }
}




double *femElasticitySolve(femProblem *theProblem){

    // Récupérer le système à résoudre
    femFullSystem *theSystem = theProblem->system;

    femFullSystemInit(theSystem);


    femElasticityAssembleElements(theProblem);

    femElasticityAssembleNeumann(theProblem);

    int iBnd, iElem, j, node;
    femBoundaryCondition *theCondition;
    femDomain *domain;
    for (iBnd = 0; iBnd < theProblem->nBoundaryConditions; iBnd++) {

        theCondition = theProblem->conditions[iBnd];
        domain = theCondition->domain;
        double value = theCondition->value;

        if (theCondition->type == DIRICHLET_X || theCondition->type == DIRICHLET_Y){
            for (iElem = 0; iElem < domain->nElem; iElem++){
                int elem = domain->elem[iElem];
                for (j = 0; j < domain->mesh->nLocalNode; j++){
                    node = domain->mesh->elem[elem*domain->mesh->nLocalNode + j];

                    if (theCondition->type == DIRICHLET_X){
                        femFullSystemConstrain(theSystem, 2*node, value);
                    }
                    if (theCondition->type == DIRICHLET_Y){
                        femFullSystemConstrain(theSystem, 2*node+1, value);
                    }
                }
            }
        }
    }

    double *solution = femFullSystemEliminate(theSystem);

    int size = theSystem->size;
    for (int i = 0; i < size; i++){
        theProblem->soluce[i] = solution[i];
    }

    return theProblem->soluce;
}


double *femElasticityForces(femProblem *theProblem){        

    femFullSystem *theSystem = theProblem->system;
    int size = theSystem->size;


    femFullSystemInit(theSystem);

    femElasticityAssembleElements(theProblem);
    femElasticityAssembleNeumann(theProblem);

    for (int i = 0; i < size; i++){
        theProblem->residuals[i] = 0.0;
        for (int j = 0; j < size; j++){
            theProblem->residuals[i] += theSystem->A[i][j] * theProblem->soluce[j];
        }
        theProblem->residuals[i] -= theSystem->B[i];
    }
    for (int i = 0; i < size; i++){
        theProblem->residuals[i] = - theProblem->residuals[i];
    }

    return theProblem->residuals;
}

double fun(double x, double y) 
{
    return 1;
}

