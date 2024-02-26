/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * onmetis.c
 *
 * This file contains the driving routine for multilevel method
 *
 * Started 8/28/94
 * George
 *
 * $Id: onmetis.c,v 1.2 2006/07/19 12:25:18 cvsusers Exp $
 *
 */

#include <metis.h>

#ifdef UPPERCASE_ 
#define ONMETIS ONMETIS_ 
#endif

#ifdef lowercase_ 
#define ONMETIS onmetis_ 
#endif

#ifdef lowercase 
#define ONMETIS onmetis 
#endif

/*************************************************************************
* Let the game begin
**************************************************************************/
/*void ONMETIS(int argc, char *argv[])*/
void ONMETIS(int *iverb)
{
  int i, options[10],argc=2;
  idxtype *perm, *iperm;
  GraphType graph;
  char filename[256];
  int numflag = 0, wgtflag;
  timer TOTALTmr, METISTmr, IOTmr, SMBTmr;

  if (argc != 2) {
    printf("Usage: %s <GraphFile>\n","onmetis");
    exit(0);
  }
    
/*  strcpy(filename, argv[1]);*/
  

  cleartimer(TOTALTmr);
  cleartimer(METISTmr);
  cleartimer(IOTmr);
  cleartimer(SMBTmr);

  starttimer(TOTALTmr);
  starttimer(IOTmr);
/*  ReadGraph(&graph, filename, &wgtflag);*/
  ReadGraph(&graph, "grid.mts.ngraph", &wgtflag);
  if (graph.nvtxs <= 0) {
    printf("Empty graph. Nothing to do.\n");
    exit(0);
  }
  if (graph.ncon != 1) {
    printf("Ordering can only be applied to graphs with one constraint.\n");
    exit(0);
  }
  stoptimer(IOTmr);

  /* Ordering does not use weights! */
  GKfree(&graph.vwgt, &graph.adjwgt, LTERM);

  if(*iverb >= 2)
  {
    printf("**********************************************************************\n");
    printf("%s", METISTITLE);
    printf("Graph Information ---------------------------------------------------\n");
    printf("  Name: %s, #Vertices: %d, #Edges: %d\n\n", "grid.mts.ngraph", graph.nvtxs, graph.nedges/2);
    printf("Node-Based Ordering... ----------------------------------------------\n");
  }

  perm = idxmalloc(graph.nvtxs, "main: perm");
  iperm = idxmalloc(graph.nvtxs, "main: iperm");
  options[0] = 0;
  
  

  starttimer(METISTmr);
  METIS_NodeND(&graph.nvtxs, graph.xadj, graph.adjncy, &numflag, options, perm, iperm);
  stoptimer(METISTmr);

  starttimer(IOTmr);
  WritePermutation("grid.mts.ngraph", iperm, graph.nvtxs); 
  stoptimer(IOTmr);

  starttimer(SMBTmr);
  ComputeFillIn(&graph, iperm);
  stoptimer(SMBTmr);

  stoptimer(TOTALTmr);

  if(*iverb >= 2)
  {
    printf("\nTiming Information --------------------------------------------------\n");
    printf("  I/O:                     \t %7.3f\n", gettimer(IOTmr));
    printf("  Ordering:                \t %7.3f   (ONMETIS time)\n", gettimer(METISTmr));
    printf("  Symbolic Factorization:  \t %7.3f\n", gettimer(SMBTmr));
    printf("  Total:                   \t %7.3f\n", gettimer(TOTALTmr));
    printf("**********************************************************************\n");
  }

  GKfree(&graph.xadj, &graph.adjncy, &perm, &iperm, LTERM);
}  


