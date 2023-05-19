#include "fem.h"

//
// Ici, vous pouvez définir votre géométrie :-)
//  (1) Raffiner intelligemment.... (yes )
//  (2) Construire la geometrie avec OpenCascade 
//  (3) Construire la geometrie avec les outils de GMSH
//  (4) Obtenir la geometrie en lisant un fichier .geo de GMSH


double geoSize(double x, double y) {

    femGeo* theGeometry = geoGetGeometry();
    return theGeometry->h * (1.0 - 0.5 * x);
}


void geoMeshGenerate() {

    femGeo* theGeometry = geoGetGeometry();
    geoSetSizeCallback(geoSize);   


    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;
    
    int ierr;
    double r = w/4;
    int idRect = gmshModelOccAddRectangle(0.0,0.0,0.0,w,h,-1,0.0,&ierr);
    int idDisk = gmshModelOccAddDisk(w/2.0,h/2.0,0.0,r,r,-1,NULL,0,NULL,0,&ierr);
    int idSlit = gmshModelOccAddRectangle(w/2.0,h/2.0-r,0.0,w,2.0*r,-1,0.0,&ierr);
    int rect[] = {2,idRect};
    int disk[] = {2,idDisk};
    int slit[] = {2,idSlit};

    gmshModelOccCut(rect,2,disk,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr); 
    gmshModelOccCut(rect,2,slit,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr); 
    gmshModelOccSynchronize(&ierr); 

    if (theGeometry->elementType == FEM_QUAD) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshOptionSetNumber("Mesh.RecombineAll",1,&ierr);
        gmshOptionSetNumber("Mesh.Algorithm",8,&ierr);  
        gmshOptionSetNumber("Mesh.RecombinationAlgorithm",1.0,&ierr); 
        gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  
        gmshModelMeshGenerate(2,&ierr);  }
  
    if (theGeometry->elementType == FEM_TRIANGLE) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshModelMeshGenerate(2,&ierr);  }
 
    return;
}

void designLuge(double h, double w, double larg, double leg_h, double leg_w, double tail, double meshSizeFactor){
    /**
    * h = height of our sleigh
    * w = width of our sleigh
    * larg = extent of our sleigh
    * leg_h = height of our leg
    * leg_w = width of our leg
    * tail = height of our tail
    * meshSizeFactor = meshSize
    */

    femGeo* theGeometry = geoGetGeometry();
    geoSetSizeCallback(geoSize);

    int ierr;
    gmshClear(&ierr);


    // Add points
    double r = h/2.0;

    double hole_x = 0;
    double hole_y = 0;
    double picket_x1 = 0;
    double picket_y1 = 0;
    double picket_x2 = 0;
    double picket_y2 = 0;

    double x = -23.5;
    double y = 15;
    double z = 0;

    gmshModelOccAddPoint(x,y,z,0,1,&ierr);
    y -= larg;
    gmshModelOccAddPoint(x,y,z,0,2,&ierr);
    y += larg/2.0;
    gmshModelOccAddPoint(x,y,z,0,27,&ierr);
    x += larg/2.0;
    gmshModelOccAddPoint(x,y,z,0,26,&ierr);
    x -= larg/2.0;
    y -= larg/2.0 + r;
    gmshModelOccAddPoint(x,y,z,0,3,&ierr);
    x -= r;
    gmshModelOccAddPoint(x,y,z,0,4,&ierr);
    x -= larg;
    gmshModelOccAddPoint(x,y,z,0,5,&ierr);
    y -= h/10.0;
    gmshModelOccAddPoint(x,y,z,0,23,&ierr);
    x += larg;
    gmshModelOccAddPoint(x,y,z,0,24,&ierr);
    x += r;
    gmshModelOccAddPoint(x,y,z,0,25,&ierr);
    y -= r;
    gmshModelOccAddPoint(x,y,z,0,6,&ierr);
    y -= larg;
    gmshModelOccAddPoint(x,y,z,0,7,&ierr);
    x += w/4.0 - r;
    gmshModelOccAddPoint(x,y,z,0,8,&ierr);
    picket_x1 = x + leg_w/2.0;
    picket_y1 = y + larg;
    y -= leg_h;
    hole_x = x + leg_w;
    hole_y = y;
    gmshModelOccAddPoint(x,y,z,0,9,&ierr);
    x -= w/4.0;
    gmshModelOccAddPoint(x,y,z,0,10,&ierr);
    x -= 4*tail/5.0;
    y += tail;
    gmshModelOccAddPoint(x,y,z,0,11,&ierr);
    x -= tail/5.0;
    gmshModelOccAddPoint(x,y,z,0,12,&ierr);
    y -= tail;
    gmshModelOccAddPoint(x,y,z,0,13,&ierr);
    x += tail;
    y -= tail;
    gmshModelOccAddPoint(x,y,z,0,14,&ierr);
    x += w;
    gmshModelOccAddPoint(x,y,z,0,15,&ierr);
    y += tail;
    gmshModelOccAddPoint(x,y,z,0,16,&ierr);
    x -= w/4.0;
    gmshModelOccAddPoint(x,y,z,0,17,&ierr);
    y += leg_h;
    gmshModelOccAddPoint(x,y,z,0,18,&ierr);
    picket_x2 = x - leg_w/2.0;
    picket_y2 = y + larg;
    x += w/4.0;
    gmshModelOccAddPoint(x,y,z,0,19,&ierr);
    y += larg;
    gmshModelOccAddPoint(x,y,z,0,20,&ierr);
    gmshModelOccAddPoint(picket_x2,picket_y2,z,0,21,&ierr);
    gmshModelOccAddPoint(picket_x1,picket_y1,z,0,22,&ierr);



    // Add curves
    gmshModelOccAddCircleArc(1,3,5,1,&ierr);
    gmshModelOccAddLine(5,23,2,&ierr);
    gmshModelOccAddCircleArc(23,25,7,3,&ierr);
    gmshModelOccAddLine(7,8,4,&ierr);
    gmshModelOccAddLine(8,9,5,&ierr);
    gmshModelOccAddLine(9,10,6,&ierr);
    gmshModelOccAddLine(10,11,7,&ierr);
    gmshModelOccAddLine(11,12,8,&ierr);
    gmshModelOccAddLine(12,13,9,&ierr);
    gmshModelOccAddCircleArc(13,10,14,10,&ierr);
    gmshModelOccAddLine(14,15,11,&ierr);
    gmshModelOccAddLine(15,16,12,&ierr);
    gmshModelOccAddLine(16,17,13,&ierr);
    gmshModelOccAddLine(17,18,14,&ierr);
    gmshModelOccAddLine(18,19,15,&ierr);
    gmshModelOccAddLine(19,20,16,&ierr);
    gmshModelOccAddLine(20,21,17,&ierr);
    gmshModelOccAddLine(21,22,18,&ierr);
    gmshModelOccAddLine(22,6,19,&ierr);
    gmshModelOccAddCircleArc(6,25,24,20,&ierr);
    gmshModelOccAddLine(24,4,21,&ierr);
    gmshModelOccAddCircleArc(4,3,2,22,&ierr);
    gmshModelOccAddCircleArc(2,27,26,23,&ierr);
    gmshModelOccAddCircleArc(26,27,1,24,&ierr);


    // Add wire (closed curve)
    int curveTags[24] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24};
    gmshModelOccAddWire(curveTags, 24, 1, 1, &ierr);

    // Add surface
    int wireTags[1] = {1};
    gmshModelOccAddPlaneSurface(wireTags, 1, 100, &ierr);

    int tagHole = gmshModelOccAddRectangle(hole_x,hole_y,0,w/2.0-2*leg_w,leg_h,-1,0.0,&ierr);

    int sleighId[] = {2, 100};
    int holeId[] = {2, tagHole};

    gmshModelOccCut(sleighId,2,holeId,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);

    // Sync
    gmshModelOccSynchronize(&ierr);

    // // Create physical group for surface
    int surfaceTags[1] = {100};
    gmshModelAddPhysicalGroup(2, surfaceTags, 1, -1, "bulk", &ierr);

    // // Create physical group for clamped curves

    if (theGeometry->elementType == FEM_QUAD) {
        gmshOptionSetNumber("Mesh.MeshSizeFactor", meshSizeFactor, &ierr);
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshOptionSetNumber("Mesh.RecombineAll",1,&ierr);
        //gmshOptionSetNumber("Mesh.Algorithm",8,&ierr);
        //gmshOptionSetNumber("Mesh.RecombinationAlgorithm",1.0,&ierr);
        gmshModelGeoMeshSetRecombine(2,1,45,&ierr);
        gmshModelMeshGenerate(2,&ierr);  }

    if (theGeometry->elementType == FEM_TRIANGLE) {
        gmshOptionSetNumber("Mesh.MeshSizeFactor", meshSizeFactor, &ierr);
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshModelMeshGenerate(2,&ierr);  }
}

void designHouse(double w, double h, double r_w, double r_h, double w_w, double d_w, double d_h, double meshSizeFactor){
    /**
    * w = width of our house
    * h = height of our house
    * r_w = width of our roof
    * r_h = height of our roof
    * w_w = width of our window
    * d_w = width of our door
    * d_h = height of our door
    * meshSizeFactor = meshSize
    */

    femGeo* theGeometry = geoGetGeometry();
    geoSetSizeCallback(geoSize);

    int ierr;
    gmshClear(&ierr);

    // Add points

    double window_x = 0 + 1.0/5.0 * w;
    double window_y = 0 - 1.0/5.0 * h - w_w;

    double x = 0;
    double y = 0;
    double z = 0;

    gmshModelOccAddPoint(x,y,z,0,1,&ierr);
    y -= h;
    gmshModelOccAddPoint(x,y,z,0,2,&ierr);
    x += 3.0/5.0 * w;
    gmshModelOccAddPoint(x,y,z,0,3,&ierr);
    y += d_h;
    gmshModelOccAddPoint(x,y,z,0,4,&ierr);
    x += d_w;
    gmshModelOccAddPoint(x,y,z,0,5,&ierr);
    y -= d_h;
    gmshModelOccAddPoint(x,y,z,0,6,&ierr);
    x += 2.0/5.0 * w - d_w;
    gmshModelOccAddPoint(x,y,z,0,7,&ierr);
    y += h;
    gmshModelOccAddPoint(x,y,z,0,8,&ierr);
    x += r_w;
    gmshModelOccAddPoint(x,y,z,0,9,&ierr);
    y += h/20.0;
    gmshModelOccAddPoint(x,y,z,0,10,&ierr);
    x += -r_w - w/2.0 + w/20.0;
    y += r_h;
    gmshModelOccAddPoint(x,y,z,0,11,&ierr);
    x -= w/10.0;
    gmshModelOccAddPoint(x,y,z,0,12,&ierr);
    y -= r_h;
    x += -w/2.0 - r_w + w/20.0;
    gmshModelOccAddPoint(x,y,z,0,13,&ierr);
    y -= h/20.0;
    gmshModelOccAddPoint(x,y,z,0,14,&ierr);


    // Add curves
    gmshModelOccAddLine(1,2,1,&ierr);
    gmshModelOccAddLine(2,3,2,&ierr);
    gmshModelOccAddLine(3,4,3,&ierr);
    gmshModelOccAddLine(4,5,4,&ierr);
    gmshModelOccAddLine(5,6,5,&ierr);
    gmshModelOccAddLine(6,7,6,&ierr);
    gmshModelOccAddLine(7,8,7,&ierr);
    gmshModelOccAddLine(8,9,8,&ierr);
    gmshModelOccAddLine(9,10,9,&ierr);
    gmshModelOccAddLine(10,11,10,&ierr);
    gmshModelOccAddLine(11,12,11,&ierr);
    gmshModelOccAddLine(12,13,12,&ierr);
    gmshModelOccAddLine(13,14,13,&ierr);
    gmshModelOccAddLine(14,1,14,&ierr);


    // Add wire (closed curve)
    int curveTags[14] = {1,2,3,4,5,6,7,8,9,10,11, 12, 13, 14};
    gmshModelOccAddWire(curveTags, 14, 1, 1, &ierr);

    // Add surface
    int wireTags[1] = {1};
    gmshModelOccAddPlaneSurface(wireTags, 1, 100, &ierr);

    int tagWindow = gmshModelOccAddRectangle(window_x,window_y,0,w_w,w_w,-1,0.0,&ierr);

    int houseId[] = {2, 100};
    int windowId[] = {2, tagWindow};

    gmshModelOccCut(houseId,2,windowId,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);

    // Sync
    gmshModelOccSynchronize(&ierr);

    // // Create physical group for surface
    int surfaceTags[1] = {100};
    gmshModelAddPhysicalGroup(2, surfaceTags, 1, -1, "bulk", &ierr);

    // // Create physical group for clamped curves

    if (theGeometry->elementType == FEM_QUAD) {
        gmshOptionSetNumber("Mesh.MeshSizeFactor", meshSizeFactor, &ierr);
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshOptionSetNumber("Mesh.RecombineAll",1,&ierr);
        gmshOptionSetNumber("Mesh.Algorithm",8,&ierr);
        gmshOptionSetNumber("Mesh.RecombinationAlgorithm",1.0,&ierr);
        gmshModelGeoMeshSetRecombine(2,1,45,&ierr);
        gmshModelMeshGenerate(2,&ierr);  }

    if (theGeometry->elementType == FEM_TRIANGLE) {
        gmshOptionSetNumber("Mesh.MeshSizeFactor", meshSizeFactor, &ierr);
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshModelMeshGenerate(2,&ierr);  }
}

void geoBasicElasticityProblem() {

    int ierr;
    gmshClear(&ierr);

    int L = 3;
    int H = 1;

    gmshModelOccAddPoint(0,0,0,4,1,&ierr);
    gmshModelOccAddPoint(0,H,0,4,2,&ierr);
    gmshModelOccAddPoint(L,H,0,4,3,&ierr);
    gmshModelOccAddPoint(L,0,0,4,4,&ierr);
    gmshModelOccAddLine(1,2,1,&ierr);
    gmshModelOccAddLine(2,3,2,&ierr);
    gmshModelOccAddLine(3,4,3,&ierr);
    gmshModelOccAddLine(4,1,4,&ierr);
    int curves[4] = {1,2,3,4};
    gmshModelOccAddWire(curves, 4, 1, 1, &ierr);
    int wiretags[1] = {1};
    gmshModelOccAddPlaneSurface(wiretags, 1, 1, &ierr);

    gmshModelOccSynchronize(&ierr);
    gmshModelMeshGenerate(2, &ierr);

    return;
}


void geoMeshGenerateGeo() {

    femGeo* theGeometry = geoGetGeometry();
    geoSetSizeCallback(geoSize);   


    /*
    4 ------------------ 3
    |                    |
    |                    |
    5 ------- 6          |
               \         |
                )        |
               /         |
    8 ------- 7          |
    |                    |
    |                    |
    1 ------------------ 2
    */

    int ierr;
    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;
    double r = w/4;
    double lc = theGeometry->h;

    int p1 = gmshModelGeoAddPoint(-w/2, -h/2, 0., lc, 1, &ierr);
    int p2 = gmshModelGeoAddPoint( w/2, -h/2, 0., lc, 2, &ierr);
    int p3 = gmshModelGeoAddPoint( w/2,  h/2, 0., lc, 3, &ierr);
    int p4 = gmshModelGeoAddPoint(-w/2,  h/2, 0., lc, 4, &ierr);
    int p5 = gmshModelGeoAddPoint(-w/2,    r, 0., lc, 5, &ierr);
    int p6 = gmshModelGeoAddPoint(0.,      r, 0., lc, 6, &ierr);
    int p7 = gmshModelGeoAddPoint(0.,     -r, 0., lc, 7, &ierr);
    int p8 = gmshModelGeoAddPoint(-w/2,   -r, 0., lc, 8, &ierr);
    int p9 = gmshModelGeoAddPoint(0.,     0., 0., lc, 9, &ierr); // center of circle


    int l1 = gmshModelGeoAddLine(p1, p2, 1, &ierr);
    int l2 = gmshModelGeoAddLine(p2, p3, 2, &ierr);
    int l3 = gmshModelGeoAddLine(p3, p4, 3, &ierr);
    int l4 = gmshModelGeoAddLine(p4, p5, 4, &ierr);
    int l5 = gmshModelGeoAddLine(p5, p6, 5, &ierr);
    int l6 = gmshModelGeoAddCircleArc(p7, p9, p6, 6, 0., 0., 0., &ierr); // NB : the direction of the curve is reversed
    int l7 = gmshModelGeoAddLine(p7, p8, 7, &ierr);
    int l8 = gmshModelGeoAddLine(p8, p1, 8, &ierr);

    int lTags[] = {l1, l2, l3, l4, l5, -l6, l7, l8}; // NB : "-l6" because the curve is reversed 
    int c1[] = {1};
    c1[0] = gmshModelGeoAddCurveLoop(lTags, 8, 1, 0, &ierr);  
    int s1 = gmshModelGeoAddPlaneSurface(c1, 1, 1, &ierr);
    gmshModelGeoSynchronize(&ierr);


    if (theGeometry->elementType == FEM_QUAD) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshOptionSetNumber("Mesh.RecombineAll",1,&ierr);
        gmshOptionSetNumber("Mesh.Algorithm",8,&ierr);  
        gmshOptionSetNumber("Mesh.RecombinationAlgorithm",1.0,&ierr); 
        gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  
        gmshModelMeshGenerate(2,&ierr);  }
  
    if (theGeometry->elementType == FEM_TRIANGLE) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshModelMeshGenerate(2,&ierr);  }

 //   gmshFltkRun(&ierr);
}


void geoMeshGenerateGeoFile(const char *filename){
    femGeo* theGeometry = geoGetGeometry();
    int ierr;
    gmshOpen(filename, &ierr); ErrorGmsh(ierr);
    if (theGeometry->elementType == FEM_QUAD) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshOptionSetNumber("Mesh.RecombineAll",1,&ierr);
        gmshOptionSetNumber("Mesh.Algorithm",8,&ierr); 
        gmshOptionSetNumber("Mesh.RecombinationAlgorithm",1.0,&ierr); 
        gmshModelGeoMeshSetRecombine(2,1,45,&ierr); 
        gmshModelMeshGenerate(2,&ierr);  }
 
    if (theGeometry->elementType == FEM_TRIANGLE) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshModelMeshGenerate(2,&ierr);  }
    return;
}