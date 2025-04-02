#include "fem.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double hermiteInterpolation(double d, double d_star, double h0, double h_star) {
    if (d >= d_star) return h_star;
    double t = d / d_star;
    return h0 + (h_star - h0) * (3 * t * t - 2 * t * t * t);
}

double geoSize(double x, double y) {
    femGeo* theGeometry = geoGetGeometry();

    double h = 0.1;
    double r = 1;
    // double h0 = theGeometry->hNotch;
    // double d0 = theGeometry->dNotch;

    // double dDisk = sqrt(x * x + y * y) - r;
    // double hDisk = (dDisk < d0) ? hermiteInterpolation(dDisk, d0, h0, h) : h;

    // return fmin(h, hDisk);
    return h;
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


