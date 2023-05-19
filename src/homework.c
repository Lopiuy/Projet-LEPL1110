#include "fem.h"

// Il faut un fifrelin generaliser ce code.....
//  (1) Ajouter l'axisymétrique !    (mandatory)
//  (2) Ajouter les conditions de Neumann !   (mandatory)  
//  (3) Ajouter les conditions en normal et tangentiel !   (strongly advised)
//  (4) Et remplacer le solveur plein par un truc un fifrelin plus subtil  (mandatory)



double *femElasticitySolve(femProblem *theProblem)
{

    femFullSystem  *theFullSystem;
    femBandSystem  *theBandSystem;
    double** A; double* B; int size;

    switch (theProblem->solverType) {
        case FEM_FULL:
            theFullSystem = (femFullSystem*) theProblem->system;
            A = theFullSystem->A;
            B = theFullSystem->B;
            size = theFullSystem->size;
            break;
        case FEM_BAND:
            theBandSystem = (femBandSystem*) theProblem->system;
            A = theBandSystem->A;
            B = theBandSystem->B;
            size = theBandSystem->size;
            break;
        default:
            Error("Unexpected solver type"); break;
    }

    femIntegration *theRule = theProblem->rule;
    femDiscrete    *theSpace = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theMesh = theGeometry->theElements;
    femMesh        *theEdges = theGeometry->theEdges;
    
    
    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
    int iElem,iInteg,iEdge,iCondition,i,j,d,map[4],mapX[4],mapY[4];
    
    int nLocal = theMesh->nLocalNode;

    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;      
    double rho = theProblem->rho;
    double g   = theProblem->g;
    
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j=0; j < nLocal; j++) {
            map[j]  = theNodes->number[theMesh->elem[iElem*nLocal+j]];
            mapX[j] = 2*map[j];
            mapY[j] = 2*map[j] + 1;
            x[j]    = theNodes->X[theMesh->elem[iElem*nLocal+j]];
            y[j]    = theNodes->Y[theMesh->elem[iElem*nLocal+j]];
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
            for (i = 0; i < theSpace->n; i++) {
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
                    for(j = 0; j < theSpace->n; j++) {
                        if(theProblem->solverType == FEM_BAND){
                            // Upper Triangle of the matrix
                            A[mapX[i]][mapX[j]] += (mapX[j] >= mapX[i]) ? (dphidx[i] * a * dphidx[j] +
                                                    dphidy[i] * c * dphidy[j]) * jac * weight : 0.0;
                            A[mapX[i]][mapY[j]] += (mapY[j] >= mapX[i]) ? (dphidx[i] * b * dphidy[j] +
                                                    dphidy[i] * c * dphidx[j]) * jac * weight : 0.0;
                            A[mapY[i]][mapX[j]] += (mapX[j] >= mapY[i]) ? (dphidy[i] * b * dphidx[j] +
                                                    dphidx[i] * c * dphidy[j]) * jac * weight : 0.0;
                            A[mapY[i]][mapY[j]] += (mapY[j] >= mapY[i]) ? (dphidy[i] * a * dphidy[j] +
                                                    dphidx[i] * c * dphidx[j]) * jac * weight : 0.0;
                        }else{
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
                }
                for (i = 0; i < theSpace->n; i++) {
                    B[mapY[i]] -= phi[i] * g * rho * jac * weight;
                }
            }

            // Axisymetric case
            else if(theProblem->planarStrainStress == AXISYM){
                for (i = 0; i < theSpace->n; i++) {
                    for(j = 0; j < theSpace->n; j++) {
                        if(theProblem->solverType == FEM_BAND){
                            // Upper Triangle of the matrix
                            A[mapX[i]][mapX[j]] += (mapX[j] >= mapX[i]) ? (dphidx[i] * a * xLoc * dphidx[j] +
                                                    dphidy[i] * c * xLoc * dphidy[j] +
                                                    phi[i] * ((b * dphidx[j]) + (a * phi[j] / xLoc)) +
                                                    dphidx[i] * b * phi[j]) * jac * weight : 0.0;
                            A[mapX[i]][mapY[j]] += (mapY[j] >= mapX[i]) ? (dphidx[i] * b * xLoc * dphidy[j] +
                                                    dphidy[i] * c * xLoc * dphidx[j] +
                                                    phi[i] * b * dphidy[j]) * jac * weight : 0.0;
                            A[mapY[i]][mapX[j]] += (mapX[j] >= mapY[i]) ? (dphidy[i] * b * xLoc * dphidx[j] +
                                                    dphidx[i] * c * xLoc * dphidy[j] +
                                                    dphidy[i] * b * phi[j]) * jac * weight : 0.0;
                            A[mapY[i]][mapY[j]] += (mapY[j] >= mapY[i]) ? (dphidy[i] * a * xLoc * dphidy[j] +
                                                    dphidx[i] * c * xLoc * dphidx[j]) * jac * weight : 0.0;
                        } else {
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
                }

                for (i = 0; i < theSpace->n; i++) {
                    B[mapY[i]] -= phi[i] * xLoc * g * rho * jac * weight;
                }

            } else {
                Error("Wrong planarStrainStress value");
            }
        }
    }


    // Neumann conditions in X and Y
    for (iCondition = 0; iCondition < theProblem->nBoundaryConditions; iCondition++) {
        femBoundaryCondition *theBoundary = theProblem->conditions[iCondition];
        if( theBoundary->type == NEUMANN_X || theBoundary->type == NEUMANN_Y) {
            femDomain *theDomain = theBoundary->domain;
            for(int iBoundaryElem = 0; iBoundaryElem < theDomain->nElem; iBoundaryElem++){
                int BoundaryElemSegment = theDomain->elem[iBoundaryElem];
                int iNode1 = theEdges->elem[2*BoundaryElemSegment];
                int iNode2 = theEdges->elem[2*BoundaryElemSegment+1];
                int BoundaryElemNode[2] = {theNodes->number[iNode1], theNodes->number[iNode2]};

                double x1 = theNodes->X[iNode1];
                double y1 = theNodes->Y[iNode1];
                double x2 = theNodes->X[iNode2];
                double y2 = theNodes->Y[iNode2];

                double length = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1)); // length/2 = integral

                if(theProblem->planarStrainStress != AXISYM){
                    if(theBoundary->type == NEUMANN_X){
                        B[2*BoundaryElemNode[0]] += theBoundary->value * length / 2;
                        B[2*BoundaryElemNode[1]] += theBoundary->value * length / 2;
                    } else if(theBoundary->type == NEUMANN_Y){
                        B[2*BoundaryElemNode[0]+1] += theBoundary->value * length / 2;
                        B[2*BoundaryElemNode[1]+1] += theBoundary->value * length / 2;
                    } else {
                        Error("Wrong boundary condition type");
                    }
                }
                else if (theProblem->planarStrainStress == AXISYM){
                    double xsi[2] = {-1/sqrt(3), 1/sqrt(3)};
                    double weightLine[2] = {1, 1};
                    double phiLine[2][2] = {{0.5 * (1 - xsi[0]), 0.5 * (1 + xsi[0])},
                                            {0.5*(1-xsi[1]), 0.5*(1+xsi[1])}};
                    double xLoc;
                    for (int l = 0; l < 2; l++) {

                        for (int k = 0; k < 2; k++){
                            xLoc = phiLine[k][0] * x1 + phiLine[k][1] * x2;

                            double integ = length / 2 * xLoc * weightLine[k] * phiLine[k][l];

                            if(theBoundary->type == NEUMANN_X){
                                B[2*BoundaryElemNode[l]] += theBoundary->value * integ;
                            } else if(theBoundary->type == NEUMANN_Y){
                                B[2*BoundaryElemNode[l] + 1] += theBoundary->value * integ;
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
        // Neumann conditions in N and T
        else if (theBoundary->type == NEUMANN_N || theBoundary->type == NEUMANN_T){
            femDomain *theDomain = theBoundary->domain;
            for(int iBoundaryElem = 0; iBoundaryElem < theDomain->nElem; iBoundaryElem++){
                int BoundaryElemSegment = theDomain->elem[iBoundaryElem];
                int iNode1 = theEdges->elem[2*BoundaryElemSegment];
                int iNode2 = theEdges->elem[2*BoundaryElemSegment+1];
                int BoundaryElemNode[2] = {theNodes->number[iNode1], theNodes->number[iNode2]};

                double x1 = theNodes->X[iNode1];
                double y1 = theNodes->Y[iNode1];
                double x2 = theNodes->X[iNode2];
                double y2 = theNodes->Y[iNode2];

                double length = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1)); // length/2 = integral

                double xComponent, yComponent, norm;
                switch (theBoundary->type) {
                    case NEUMANN_N :
                        xComponent = y1 - y2;
                        yComponent = x2 - x1;
                        break;
                    case NEUMANN_T :
                        xComponent = x2 - x1;
                        yComponent = y2 - y1;
                        break;
                    default:
                        break;
                }
                norm = sqrt(xComponent*xComponent + yComponent*yComponent);
                xComponent /= norm;
                yComponent /= norm;

                if(theProblem->planarStrainStress != AXISYM){
                    B[2*BoundaryElemNode[0]] += theBoundary->value * xComponent * length / 2;
                    B[2*BoundaryElemNode[1]] += theBoundary->value * xComponent * length / 2;
                    B[2*BoundaryElemNode[0]+1] += theBoundary->value * yComponent * length / 2;
                    B[2*BoundaryElemNode[1]+1] += theBoundary->value * yComponent * length / 2;
                }
                else if (theProblem->planarStrainStress == AXISYM){
                    double xsi[2] = {-1/sqrt(3), 1/sqrt(3)};
                    double weightLine[2] = {1, 1};
                    double phiLine[2][2] = {{0.5*(1-xsi[0]), 0.5*(1+xsi[0])},
                                        {0.5*(1-xsi[1]), 0.5*(1+xsi[1])}};
                    double xLoc;
                    for (int l = 0; l < 2; l++) {

                        for (int k = 0; k < 2; k++){
                            xLoc = phiLine[k][0] * x1 + phiLine[k][1] * x2;

                            double integ = length/2*xLoc*weightLine[k]*phiLine[k][l];

                            B[2*BoundaryElemNode[l]] += theBoundary->value * integ * xComponent;
                            B[2*BoundaryElemNode[l] + 1] += theBoundary->value * integ * yComponent;

                        }
                    }

                } else {
                    Error("Wrong planarStrainStress value");
                }
            }
        }
        else if (theBoundary->type == DIRICHLET_N || theBoundary->type == DIRICHLET_T){
            femDomain *theDomain = theBoundary->domain;

            int currBoundaryElem, prevBoundaryElem, iNode, iNodeNext, iNodePrev;
            double prevNormalX, prevNormalY, nextNormalX, nextNormalY,  normalX, normalY, norm;
            double prevLength, nextLength;
            double tangentX, tangentY;

            int nDomainElem = theDomain->nElem;
            for(int iBoundaryElem = 0; iBoundaryElem <= nDomainElem ; iBoundaryElem++){
                currBoundaryElem = iBoundaryElem == nDomainElem ?
                        femFindExBoundaryElem(theEdges, theDomain->elem[iBoundaryElem-1],
                                              theEdges->elem[2*theDomain->elem[nDomainElem-1]+1]) :
                        theDomain->elem[iBoundaryElem];
                prevBoundaryElem = iBoundaryElem == 0 ? femFindExBoundaryElem(theEdges, currBoundaryElem, iNode) :
                                       theEdges->elem[2*currBoundaryElem-1];
                iNode = theEdges->elem[2*currBoundaryElem];
                iNodeNext = theEdges->elem[2*currBoundaryElem+1];
                iNodePrev = theEdges->elem[2*prevBoundaryElem];

                prevNormalX = iBoundaryElem == 0 ? theNodes->Y[iNodePrev] - theNodes->Y[iNode] : nextNormalX;
                prevNormalY = iBoundaryElem == 0 ? theNodes->X[iNode] - theNodes->X[iNodePrev] : nextNormalY;
                nextNormalX = theNodes->Y[iNode] - theNodes->Y[iNodeNext];
                nextNormalY = theNodes->X[iNodeNext] - theNodes->X[iNode];

                double xNode[3] = {theNodes->X[iNodePrev], theNodes->X[iNode], theNodes->X[iNodeNext]};
                double yNode[3] = {theNodes->Y[iNodePrev], theNodes->Y[iNode], theNodes->Y[iNodeNext]};

                prevLength = iBoundaryElem == 0 ? sqrt((yNode[1] - yNode[0])*(yNode[1] - yNode[0]) + (xNode[1] - xNode[0])*(xNode[1] - xNode[0])) : nextLength;
                nextLength = sqrt((yNode[2] - yNode[1])*(yNode[2] - yNode[1]) + (xNode[2] - xNode[1])*(xNode[2] - xNode[1]));

                normalX = prevLength*prevNormalX + nextLength*nextNormalX;
                normalY = prevLength*prevNormalY + nextLength*nextNormalY;

                normalX /= prevLength + nextLength;
                normalY /= prevLength + nextLength;

                norm = sqrt(normalX*normalX + normalY*normalY);

                normalX /= norm;
                normalY /= norm;

                tangentX = normalY;
                tangentY = -normalX;

                double psi;
                psi = (tangentY < 0) ? asin(tangentX) : M_PI - asin(tangentX);

                double val = theBoundary->value;

                double dirichletX = (theBoundary->type == DIRICHLET_N) ? val * cos(psi) : val * sin(psi);
                double dirichletY = (theBoundary->type == DIRICHLET_N) ? val * sin(psi) : val * -cos(psi);

                if(theProblem->solverType == FEM_BAND){
                    if(!theProblem->ntConditions[iNode]){
                        femBandSystemConstrain(theBandSystem,2*theNodes->number[iNode],dirichletX);
                        femBandSystemConstrain(theBandSystem,2*theNodes->number[iNode]+1,dirichletY);
                        theProblem->ntConditions[iNode] = 1;
                    } else {
                        theBandSystem->B[2*theNodes->number[iNode]] += dirichletX;
                        theBandSystem->B[2*theNodes->number[iNode]+1] += dirichletY;
                    }
                }else {
                    if(!theProblem->ntConditions[iNode]){
                        femFullSystemConstrain(theFullSystem,2*theNodes->number[iNode],dirichletX);
                        femFullSystemConstrain(theFullSystem,2*theNodes->number[iNode]+1,dirichletY);
                        theProblem->ntConditions[iNode] = 1;
                    } else {
                        theFullSystem->B[2*theNodes->number[iNode]] += dirichletX;
                        theFullSystem->B[2*theNodes->number[iNode]+1] += dirichletY;
                    }
                }
            }
        }
    }
    double *solution;

    if(theProblem->solverType == FEM_BAND){

        int *theConstrainedNodes = theProblem->constrainedNodes;
        for (i = 0; i < theBandSystem->size; i++) {
            if (theConstrainedNodes[i] != -1) {
                double value = theProblem->conditions[theConstrainedNodes[i]]->value;
                femBandSystemConstrain(theBandSystem,2*theNodes->number[i/2] + i%2,value); }}
        solution = femBandSystemEliminate(theBandSystem);

    }else{

        int *theConstrainedNodes = theProblem->constrainedNodes;
        for (i = 0; i < theFullSystem->size; i++) {
            if (theConstrainedNodes[i] != -1) {
                double value = theProblem->conditions[theConstrainedNodes[i]]->value;
                femFullSystemConstrain(theFullSystem,2*theNodes->number[i/2] + i%2,value); }}
        solution = femFullSystemEliminate(theFullSystem);

    }
    if(theProblem->renumType == FEM_XNUM || theProblem->renumType == FEM_YNUM){
        femDenumber(theNodes,size,solution);
    }
    return solution;
}
