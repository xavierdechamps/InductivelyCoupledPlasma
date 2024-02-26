/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * partnmesh.c
 *
 * This file reads in the element node connectivity array of a mesh and 
 * partitions both the elements and the nodes using KMETIS on the dual graph.
 *
 * Started 9/29/97
 * George
 *
 * $Id: partnmesh.c,v 1.2 2006/07/19 12:25:18 cvsusers Exp $
 *
 */

#include <metis.h>

#ifdef UPPERCASE_ 
#define PARTNMESH PARTNMESH_ 
#endif

#ifdef lowercase_ 
#define PARTNMESH partnmesh_ 
#endif

#ifdef lowercase 
#define PARTNMESH partnmesh 
#endif

/*************************************************************************
* Let the game begin
**************************************************************************/
void PARTNMESH(int* npartsp, int* iverb)
{
  int i, j, ne, nn, etype, numflag=0, nparts, edgecut,argc=3;
  idxtype *elmnts, *epart, *npart;
  timer IOTmr, DUALTmr;
  char etypestr[4][5] = {"TRI", "TET", "HEX", "QUAD"};
  GraphType graph;

  if (argc != 3) {
    printf("Usage: %s <meshfile> <nparts>\n","partnmesh");
    exit(0);
  }

/*  nparts = atoi(argv[2]);*/
  nparts=*npartsp;
  if (nparts < 2) {
    printf("nparts must be greater than one.\n");
    exit(0);
  }
   
  cleartimer(IOTmr);
  cleartimer(DUALTmr);

  starttimer(IOTmr);
  elmnts = ReadMesh("grid.mts", &ne, &nn, &etype);
  stoptimer(IOTmr);

  epart = idxmalloc(ne, "main: epart");
  npart = idxmalloc(nn, "main: npart");

  if(*iverb >= 2)
  {
    printf("**********************************************************************\n");
    printf("%s", METISTITLE);
    printf("Mesh Information ----------------------------------------------------\n");
    printf("  Name: %s, #Elements: %d, #Nodes: %d, Etype: %s\n\n","grid.mts", ne, nn, etypestr[etype-1]);
    printf("Partitioning Nodal Graph... -----------------------------------------\n");
  }

  starttimer(DUALTmr);
  METIS_PartMeshNodal(&ne, &nn, elmnts, &etype, &numflag, &nparts, &edgecut, epart, npart);
  stoptimer(DUALTmr);

  if(*iverb >= 2)
  {
    printf("  %d-way Edge-Cut: %7d, Balance: %5.2f\n", nparts, edgecut, ComputeElementBalance(ne, nparts, epart));
  }
  
  starttimer(IOTmr);
  WriteMeshPartition("grid.mts", nparts, ne, epart, nn, npart);
  stoptimer(IOTmr);

  if(*iverb >= 2)
  {
    printf("\nTiming Information --------------------------------------------------\n");
    printf("  I/O:          \t\t %7.3f\n", gettimer(IOTmr));
    printf("  Partitioning: \t\t %7.3f\n", gettimer(DUALTmr));
    printf("**********************************************************************\n");
  }
/*
  graph.nvtxs = ne;
  graph.xadj = idxmalloc(ne+1, "xadj");
  graph.vwgt = idxsmalloc(ne, 1, "vwgt");
  graph.adjncy = idxmalloc(10*ne, "adjncy");
  graph.adjwgt = idxsmalloc(10*ne, 1, "adjncy");

  METIS_MeshToDual(&ne, &nn, elmnts, &etype, &numflag, graph.xadj, graph.adjncy);

  ComputePartitionInfo(&graph, nparts, epart);

  GKfree(&graph.xadj, &graph.adjncy, &graph.vwgt, &graph.adjwgt, LTERM);
*/

  GKfree(&elmnts, &epart, &npart, LTERM);
}


