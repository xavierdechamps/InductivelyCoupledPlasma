/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * mesh2nodal.c
 *
 * This file reads in the element node connectivity array of a mesh and writes
 * out its dual in the format suitable for Metis.
 *
 * Started 9/29/97
 * George
 *
 * $Id: mesh2nodal.c,v 1.2 2006/07/19 12:25:18 cvsusers Exp $
 *
 */

#include <metis.h>

#ifdef UPPERCASE_
#define MESH2NODAL MESH2NODAL_
#endif

#ifdef lowercase_
#define MESH2NODAL mesh2nodal_
#endif

#ifdef lowercase
#define MESH2NODAL mesh2nodal
#endif

/*************************************************************************
* Let the game begin
**************************************************************************/
void MESH2NODAL(int *iverb)
/*void MESH2NODAL(*argcp,argv)*/
{
  int i, j, ne, nn, etype, numflag=0,argc=2;
  idxtype *elmnts, *xadj, *adjncy;
  timer IOTmr, DUALTmr;
  char fileout[256], etypestr[4][5] = {"TRI", "TET", "HEX", "QUAD"};
  if (argc != 2) {
    printf("Usage: %s <meshfile>\n","mesh2nodal");
    exit(0);
  }
  cleartimer(IOTmr);
  cleartimer(DUALTmr);
  
  starttimer(IOTmr);
  elmnts = ReadMesh("grid.mts", &ne, &nn, &etype);
  stoptimer(IOTmr);
  if(*iverb >= 2)
  {
    printf("**********************************************************************\n");
    printf("%s", METISTITLE);
    printf("Mesh Information ----------------------------------------------------\n");
    printf("  Name: %s, #Elements: %d, #Nodes: %d, Etype: %s\n\n", "grid.mts", ne, nn, etypestr[etype-1]);
    printf("Forming Nodal Graph... ----------------------------------------------\n");
  }
  
  xadj = idxmalloc(nn+1, "main: xadj");
  adjncy = idxmalloc(20*nn, "main: adjncy");

  starttimer(DUALTmr);
  METIS_MeshToNodal(&ne, &nn, elmnts, &etype, &numflag, xadj, adjncy);
  stoptimer(DUALTmr);

  if(*iverb >= 2)
  {
    printf("  Nodal Information: #Vertices: %d, #Edges: %d\n", nn, xadj[nn]/2);
  }
  
  sprintf(fileout, "%s.ngraph", "grid.mts");
  starttimer(IOTmr);
  WriteGraph(fileout, nn, xadj, adjncy);
  stoptimer(IOTmr);

  if(*iverb >= 2)
  {
    printf("\nTiming Information --------------------------------------------------\n");
    printf("  I/O:          \t\t %7.3f\n", gettimer(IOTmr));
    printf("  Nodal Creation:\t\t %7.3f\n", gettimer(DUALTmr));
    printf("**********************************************************************\n");
  }
  
  GKfree(&elmnts, &xadj, &adjncy, LTERM);
}


