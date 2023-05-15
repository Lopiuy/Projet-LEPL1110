#include "fem.h"

// Il faut un fifrelin generaliser ce code.....
//  (1) Ajouter l'axisymétrique !    (mandatory)
//  (2) Ajouter les conditions de Neumann !   (mandatory)  
//  (3) Ajouter les conditions en normal et tangentiel !   (strongly advised)
//  (4) Et remplacer le solveur plein par un truc un fifrelin plus subtil  (mandatory)



double *femElasticitySolve(femProblem *theProblem)
{

    femFullSystem  *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete    *theSpace = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theMesh = theGeometry->theElements;
    femMesh        *theEdges = theGeometry->theEdges;
    
    
    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
    int iElem,iInteg,iEdge,iCondition,i,j,d,map[4],mapX[4],mapY[4];
    
    int nLocal = theMesh->nLocalNode;   // = 3 if triangle, = 4 if quads

    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;      
    double rho = theProblem->rho;
    double g   = theProblem->g;
    double **A = theSystem->A;
    double *B  = theSystem->B;
    
    
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j=0; j < nLocal; j++) {
            map[j]  = theMesh->elem[iElem*nLocal+j];    // les indices des nlocalnodes(j) associé a iElem
            mapX[j] = 2*map[j];                         // ??
            mapY[j] = 2*map[j] + 1;                     // ??
            x[j]    = theNodes->X[map[j]];              // x du node
            y[j]    = theNodes->Y[map[j]];              // y du node
        }
        
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
            double xLoc = 0.0;
            for (i = 0; i < theSpace->n; i++) {     // n = 3 si triangle, 4 si quad
                dxdxsi += x[i]*dphidxsi[i];       
                dxdeta += x[i]*dphideta[i];   
                dydxsi += y[i]*dphidxsi[i];   
                dydeta += y[i]*dphideta[i];
                xLoc += x[i]*phi[i];
            }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            
            for (i = 0; i < theSpace->n; i++) {    
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
            }

            // Planar case
            if (theProblem->planarStrainStress != AXISYM) {
                for (i = 0; i < theSpace->n; i++) {
                    for(j = 0; j < theSpace->n; j++) {                      // Note elasticity, page 10
                        A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] +
                                                dphidy[i] * c * dphidy[j]) * jac * weight;
                        A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] +
                                                dphidy[i] * c * dphidx[j]) * jac * weight;
                        A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] +
                                                dphidx[i] * c * dphidy[j]) * jac * weight;
                        A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] +
                                                dphidx[i] * c * dphidx[j]) * jac * weight;
                    }
                }
                for (i = 0; i < theSpace->n; i++) {
                    B[mapY[i]] -= phi[i] * g * rho * jac * weight;
                }
            }

            // Axisymetric case
            else if(theProblem->planarStrainStress == AXISYM){

                for (i = 0; i < theSpace->n; i++) {
                    for(j = 0; j < theSpace->n; j++) {                      // Note elasticity, page 10
                        A[mapX[i]][mapX[j]] += (dphidx[i] * a * xLoc * dphidx[j] +
                                                dphidy[i] * c * xLoc * dphidy[j] +
                                                phi[i] * ((b * dphidx[j]) + (a * phi[j] / xLoc)) +
                                                dphidx[i] * b * phi[j]) * jac * weight;
                        A[mapX[i]][mapY[j]] += (dphidx[i] * b * xLoc * dphidy[j] +
                                                dphidy[i] * c * xLoc * dphidx[j] +
                                                phi[i] * b * dphidy[j]) * jac * weight;
                        A[mapY[i]][mapX[j]] += (dphidy[i] * b * xLoc * dphidx[j] +
                                                dphidx[i] * c * xLoc * dphidy[j] +
                                                dphidy[i] * b * phi[j]) * jac * weight;
                        A[mapY[i]][mapY[j]] += (dphidy[i] * a * xLoc * dphidy[j] +
                                                dphidx[i] * c * xLoc * dphidx[j]) * jac * weight;
                    }
                }

                for (i = 0; i < theSpace->n; i++) {
                    B[mapY[i]] -= phi[i] * xLoc * g * rho * jac * weight;
                }

            } else {
                Error("Wrong planarStrainStress value");
            }
        }
    }

    // Neumann conditions
    for (iCondition = 0; theProblem->nBoundaryConditions; iCondition++) {
        femBoundaryCondition *theCondition = theProblem->conditions[iCondition];
        femDomain *theDomain = theCondition->domain;
        if( theCondition->type == NEUMANN_X || theCondition->type == NEUMANN_Y) {
            for(int iBoundaryElem = 0; iBoundaryElem < theDomain->nElem; iBoundaryElem++){
                int BoundaryElemSegment = theDomain->elem[iBoundaryElem];
                int BoundaryElemNode[2] = {theEdges->elem[2*BoundaryElemSegment], theEdges->elem[2*BoundaryElemSegment+1]};

                double x1 = theNodes->X[BoundaryElemNode[0]];
                double y1 = theNodes->Y[BoundaryElemNode[0]];
                double x2 = theNodes->X[BoundaryElemNode[1]];
                double y2 = theNodes->Y[BoundaryElemNode[1]];

                double length = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1)); // length/2 = integral

                if(theProblem->planarStrainStress != AXISYM){
                    if(theCondition->type == NEUMANN_X){
                        B[BoundaryElemNode[0]] += theCondition->value * length / 2;
                        B[BoundaryElemNode[1]] += theCondition->value * length / 2;
                    } else if(theCondition->type == NEUMANN_Y){
                        B[BoundaryElemNode[0]+1] += theCondition->value * length / 2;
                        B[BoundaryElemNode[1]+1] += theCondition->value * length / 2;
                    } else {
                        Error("Wrong boundary condition type");
                    }
                }
                else if (theProblem->planarStrainStress == AXISYM){
                    double xsi[2] = {-1/sqrt(3), 1/sqrt(3)};
                    double weight[2] = {1, 1};
                    double phi[2][2] = {{0.5*(1-xsi[0]), 0.5*(1+xsi[0])},
                                        {0.5*(1-xsi[1]), 0.5*(1+xsi[1])}};
                    double xLoc;
                    for (int j = 0; j < 2; j++) {

                        for (int k = 0; k < 2; k++){
                            xLoc = phi[k][0] * x1 + phi[k][1] * x2;

                            double integ = length/2*xLoc*weight[k]*phi[k][j];

                            if(theCondition->type == NEUMANN_X){
                                B[BoundaryElemNode[j]] += theCondition->value * integ;
                            } else if(theCondition->type == NEUMANN_Y){
                                B[BoundaryElemNode[j]+1] += theCondition->value * integ;
                            } else {
                                Error("Wrong boundary condition type");
                            }
                        }
                    }

                } else {
                    Error("Wrong planarStrainStress value");
                }
            }
        }

    }
  
    int *theConstrainedNodes = theProblem->constrainedNodes;     
    for (int i=0; i < theSystem->size; i++) {
        if (theConstrainedNodes[i] != -1) {
            double value = theProblem->conditions[theConstrainedNodes[i]]->value;
            femFullSystemConstrain(theSystem,i,value); }}
                            
    return femFullSystemEliminate(theSystem);
}
