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
#include <time.h>

int main(int argc, char* argv[])
{
    femGeo* theGeometry = geoGetGeometry();   
    geoMeshRead("../data/mesh.txt");
    femProblem* theProblem = femElasticityRead(theGeometry,"../data/problem.txt");
    if(argc >= 2 && strcmp(argv[1],"basic") == 0)
        theProblem->constrainedNodes[0] = 0;
    //femElasticityPrint(theProblem);
    double *theSoluce = femElasticitySolve(theProblem);
    //femFullSystemPrint(theProblem->system);

    femNodes *theNodes = theGeometry->theNodes;
    femFieldWrite(theNodes->nNodes,2,&theSoluce[0],"../data/U.txt");
    femFieldWrite(theNodes->nNodes,2,&theSoluce[1],"../data/V.txt");
    femElasticityFree(theProblem); 
    geoFree();
    return 0;  
}

 
