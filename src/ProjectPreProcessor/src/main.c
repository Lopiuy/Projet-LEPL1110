/*
 *  main.c
 *  Projet 2022-2023
 *  Elasticite lineaire plane
 *
 *  Preprocesseur
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */
 
#include "glfem.h"

int main(int argc, char* argv[])
{  

//
//  -1- Construction de la geometrie 
//

    geoInitialize();
    femGeo* theGeometry = geoGetGeometry();
    theGeometry->elementType = FEM_TRIANGLE;

    double meshSize = 0.11;

    
//    geoMeshGenerateGeo();   // Utilisation de outils de GMSH
                            // Attention : les entit�s sont diff�rentes !
                            // On a aussi invers� la g�omtrie pour rire !
                            
//    geoMeshGenerateGeoFile("../data/mesh.geo");   // Lecture fichier geo

//
//  -2- Definition des problemes
//

    femProblem* theProblem;

// -2.0- Definition du problème luge

    if(argc >= 2 && strcmp(argv[1], "luge") == 0){
        designLuge(10, 40, 1.6, 7, 1.5, 1, meshSize);
        geoMeshImport();
        geoSetDomainName(2,"Attache");
        geoSetDomainName(13,"Assise");
        geoSetDomainName(20,"Patin");
        geoSetDomainName(16,"Rocher");


        geoMeshWrite("../../../data/mesh.txt");
        double E   = 211.e9;
        double nu  = 0.3;
        double rho = 7.85e3;
        double g   = 9.81;
        theProblem = femElasticityCreate(theGeometry,E,nu,rho,g,PLANAR_STRESS,FEM_BAND, FEM_XNUM);
        //femElasticityAddBoundaryCondition(theProblem,"Patin",DIRICHLET_X,0.0);
        femElasticityAddBoundaryCondition(theProblem,"Patin",DIRICHLET_Y,0.0);
        femElasticityAddBoundaryCondition(theProblem,"Assise",NEUMANN_Y,-2.61e3);
        femElasticityAddBoundaryCondition(theProblem,"Attache",NEUMANN_X,-0.07);
        femElasticityAddBoundaryCondition(theProblem,"Rocher",DIRICHLET_X,0.0);
        femElasticityWrite(theProblem,"../../../data/problem.txt");
    }

// -2.1- Definition du probleme section de poutre

    if(argc >= 2 && strcmp(argv[1],"poutre") == 0){
        double Lx = 1.0;
        double Ly = 1.0;
        theGeometry->LxPlate     =  Lx;
        theGeometry->LyPlate     =  Ly;
        theGeometry->h           =  Lx * meshSize;

        geoMeshGenerate();      // Utilisation de OpenCascade
        geoMeshImport();
        //geoMeshPrint();
        geoSetDomainName(0,"Symetry");
        geoSetDomainName(7,"Bottom");
        geoMeshWrite("../../../data/mesh.txt");

        double E   = 211.e9;
        double nu  = 0.3;
        double rho = 7.85e3;
        double g   = 9.81;
        theProblem = femElasticityCreate(theGeometry,E,nu,rho,g,PLANAR_STRAIN, FEM_BAND, FEM_YNUM);
        femElasticityAddBoundaryCondition(theProblem,"Symetry",DIRICHLET_X,0.0);
        femElasticityAddBoundaryCondition(theProblem,"Bottom",DIRICHLET_Y,0.0);
    }

//  -2.2- Definition du probleme basic elasticity

    if(argc >= 2 && strcmp(argv[1],"basic") == 0){
        geoBasicElasticityProblem();
        geoMeshImport();
        geoSetDomainName(1, "Top");
        geoSetDomainName(3, "Bottom");
        geoMeshWrite("../../../data/mesh.txt");
        double E   = 211.e9;
        double nu  = 0.3;
        double rho = 7.85e3;
        double g   = 0;
        theProblem = femElasticityCreate(theGeometry,E,nu,rho,g,PLANAR_STRESS, FEM_BAND, FEM_YNUM);
        double q = 100;
        femElasticityAddBoundaryCondition(theProblem, "Bottom", DIRICHLET_Y, 0.0);
        femElasticityAddBoundaryCondition(theProblem, "Top", NEUMANN_Y, -q);
    }


// -2.3- Definition du probleme maison

    if (argc >= 2 && strcmp(argv[1],"maison") == 0){
        designHouse(20.0, 20.0, 4.0, 6.0, 4.0, 5.0, 7.2, meshSize);
        geoMeshImport();
        geoSetDomainName(2,"BottomL");
        geoSetDomainName(4,"LeftDoor");
        geoSetDomainName(5,"RoofLeft");
        geoSetDomainName(6,"BottomDoor");
        geoSetDomainName(8,"RightDoor");
        geoSetDomainName(9,"RoofRight");
        geoSetDomainName(10,"BottomR");

        geoMeshWrite("../../../data/mesh.txt");
        double E   = 14e11;
        double nu  = 0.22;
        double rho = 1.8e3;
        double g   = 9.81;
        theProblem = femElasticityCreate(theGeometry,E,nu,rho,g,PLANAR_STRAIN,FEM_BAND, FEM_YNUM);
        femElasticityAddBoundaryCondition(theProblem,"BottomL",DIRICHLET_Y,0.0);
        femElasticityAddBoundaryCondition(theProblem,"BottomR",DIRICHLET_Y,0.0);
        femElasticityAddBoundaryCondition(theProblem,"BottomL",DIRICHLET_X,0.0);
        femElasticityAddBoundaryCondition(theProblem,"BottomR",DIRICHLET_X,0.0);
        femElasticityAddBoundaryCondition(theProblem,"RoofLeft",NEUMANN_N,-100000.0);
        femElasticityWrite(theProblem,"../../../data/problem.txt");
    }

    //femElasticityPrint(theProblem);
    femElasticityWrite(theProblem,"../../../data/problem.txt");


//
//  -3- Champ de la taille de r�f�rence du maillage
//

    double *meshSizeField = malloc(theGeometry->theNodes->nNodes*sizeof(double));
    femNodes *theNodes = theGeometry->theNodes;
    for(int i=0; i < theNodes->nNodes; ++i)
        meshSizeField[i] = theGeometry->geoSize(theNodes->X[i], theNodes->Y[i]);
    double hMin = femMin(meshSizeField,theNodes->nNodes);  
    double hMax = femMax(meshSizeField,theNodes->nNodes);  
    printf(" ==== Global requested h : %14.7e \n",theGeometry->h);
    printf(" ==== Minimum h          : %14.7e \n",hMin);
    printf(" ==== Maximum h          : %14.7e \n",hMax);
    
//
//  -4- Visualisation 
//

    if(argc >= 3 && strcmp(argv[2], "view") == 0){
        int mode = 1;
        int domain = 0;
        int freezingButton = FALSE;
        double t, told = 0;
        char theMessage[MAXNAME];


        GLFWwindow* window = glfemInit("EPL1110 : Project 2022-23 ");
        glfwMakeContextCurrent(window);

        do {
            int w,h;
            glfwGetFramebufferSize(window,&w,&h);
            glfemReshapeWindows(theGeometry->theNodes,w,h);

            t = glfwGetTime();
            if (glfwGetKey(window,'D') == GLFW_PRESS) { mode = 0;}
            if (glfwGetKey(window,'V') == GLFW_PRESS) { mode = 1;}
            if (glfwGetKey(window,'N') == GLFW_PRESS && freezingButton == FALSE) { domain++; freezingButton = TRUE; told = t;}

            if (t-told > 0.5) {freezingButton = FALSE; }
            if (mode == 1) {
                glfemPlotField(theGeometry->theElements,meshSizeField);
                glfemPlotMesh(theGeometry->theElements);
                sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
                glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
            if (mode == 0) {
                domain = domain % theGeometry->nDomains;
                glfemPlotDomain( theGeometry->theDomains[domain]);
                sprintf(theMessage, "%s : %d ",theGeometry->theDomains[domain]->name,domain);
                glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);  }

            glfwSwapBuffers(window);
            glfwPollEvents();
        } while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
                 glfwWindowShouldClose(window) != 1 );
    }

            
    // Check if the ESC key was pressed or the window was closed

    free(meshSizeField);
    femElasticityFree(theProblem) ; 
    geoFree();
    glfwTerminate(); 
    
    exit(EXIT_SUCCESS);
    return 0;  
}


 
