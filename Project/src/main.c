/*
 *  main.c
 *  Projet 2022-2023
 *  Elasticite lineaire plane
 *
 *  Code de calcul
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */
 
#include "fem.h"

int main(void)
{
//
    double E   = 211.e9;
    double nu  = 0.3;
    double rho = 7.85e3;
    double q = 100;
    double g   = 9.81;
//
    femGeo* theGeometry = geoGetGeometry();   
    geoMeshRead("../../data/mesh.txt");
    femProblem* theProblem = femElasticityRead(theGeometry,"../../data/problem.txt");
    femElasticityPrint(theProblem);
    double *theSoluce = femElasticitySolve(theProblem);
//
/*
    double *X = theProblem->geometry->theNodes->X;
    double *Y = theProblem->geometry->theNodes->Y;

    printf(" ==== Soluce : \n");
    for (int i=0; i<theProblem->geometry->theNodes->nNodes*2; i++) {
        if (i % 2 == 0) printf("U%d = %14.7e", i/2, -theSoluce[i]);
        else printf("V%d = %14.7e", i/2, theSoluce[i]);

        if (theProblem->constrainedNodes[i] == -1) {
            if (i%2==0) printf("\t expected : %14.7e\n", q*nu*(3 - X[i/2])/E);
            else        printf("\t expected : %14.7e\n", -q*Y[i/2]/E);
        }
        else printf("\t exepcted :  0.0\n");
    }
//*/

    femNodes *theNodes = theGeometry->theNodes;
    femFieldWrite(theNodes->nNodes,2,&theSoluce[0],"../../data/U.txt");
    femFieldWrite(theNodes->nNodes,2,&theSoluce[1],"../../data/V.txt");
    femElasticityFree(theProblem); 
    geoFree();
    return 0;  
}

 
