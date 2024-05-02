/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * rename.h
 *
 * This file contains header files
 *
 * Started 10/2/97
 * George
 *
 * $Id: rename.h 13933 2013-03-29 22:20:46Z karypis $
 *
 */


#ifndef _K_METIS_RENAME_H_
#define _K_METIS_RENAME_H_

/* balance.c */
#define Balance2Way			K_METIS__Balance2Way
#define Bnd2WayBalance			K_METIS__Bnd2WayBalance
#define General2WayBalance		K_METIS__General2WayBalance
#define McGeneral2WayBalance            K_METIS__McGeneral2WayBalance

/* bucketsort.c */
#define BucketSortKeysInc		K_METIS__BucketSortKeysInc

/* checkgraph.c */
#define CheckGraph                      K_METIS__CheckGraph
#define CheckInputGraphWeights          K_METIS__CheckInputGraphWeights
#define FixGraph                        K_METIS__FixGraph

/* coarsen.c */
#define CoarsenGraph			K_METIS__CoarsenGraph
#define Match_RM                        K_METIS__Match_RM
#define Match_SHEM                      K_METIS__Match_SHEM
#define Match_2Hop                      K_METIS__Match_2Hop
#define Match_2HopAny                   K_METIS__Match_2HopAny
#define Match_2HopAll                   K_METIS__Match_2HopAll
#define PrintCGraphStats                K_METIS__PrintCGraphStats
#define CreateCoarseGraph		K_METIS__CreateCoarseGraph
#define CreateCoarseGraphNoMask		K_METIS__CreateCoarseGraphNoMask
#define CreateCoarseGraphPerm		K_METIS__CreateCoarseGraphPerm
#define SetupCoarseGraph		K_METIS__SetupCoarseGraph
#define ReAdjustMemory			K_METIS__ReAdjustMemory

/* compress.c */
#define CompressGraph			K_METIS__CompressGraph
#define PruneGraph			K_METIS__PruneGraph

/* contig.c */
#define FindPartitionInducedComponents  K_METIS__FindPartitionInducedComponents   
#define IsConnected                     K_METIS__IsConnected
#define IsConnectedSubdomain            K_METIS__IsConnectedSubdomain
#define FindSepInducedComponents        K_METIS__FindSepInducedComponents
#define EliminateComponents             K_METIS__EliminateComponents
#define MoveGroupContigForCut           K_METIS__MoveGroupContigForCut
#define MoveGroupContigForVol           K_METIS__MoveGroupContigForVol

/* debug.c */
#define ComputeCut			K_METIS__ComputeCut
#define ComputeVolume			K_METIS__ComputeVolume
#define ComputeMaxCut			K_METIS__ComputeMaxCut
#define CheckBnd			K_METIS__CheckBnd
#define CheckBnd2			K_METIS__CheckBnd2
#define CheckNodeBnd			K_METIS__CheckNodeBnd
#define CheckRInfo			K_METIS__CheckRInfo
#define CheckNodePartitionParams	K_METIS__CheckNodePartitionParams
#define IsSeparable			K_METIS__IsSeparable
#define CheckKWayVolPartitionParams     K_METIS__CheckKWayVolPartitionParams

/* fm.c */
#define FM_2WayRefine                   K_METIS__FM_2WayRefine
#define FM_2WayCutRefine                K_METIS__FM_2WayCutRefine
#define FM_Mc2WayCutRefine              K_METIS__FM_Mc2WayCutRefine
#define SelectQueue                     K_METIS__SelectQueue
#define Print2WayRefineStats            K_METIS__Print2WayRefineStats

/* fortran.c */
#define Change2CNumbering		K_METIS__Change2CNumbering
#define Change2FNumbering		K_METIS__Change2FNumbering
#define Change2FNumbering2		K_METIS__Change2FNumbering2
#define Change2FNumberingOrder		K_METIS__Change2FNumberingOrder
#define ChangeMesh2CNumbering		K_METIS__ChangeMesh2CNumbering
#define ChangeMesh2FNumbering		K_METIS__ChangeMesh2FNumbering
#define ChangeMesh2FNumbering2		K_METIS__ChangeMesh2FNumbering2

/* graph.c */
#define SetupGraph			K_METIS__SetupGraph
#define SetupGraph_adjrsum              K_METIS__SetupGraph_adjrsum
#define SetupGraph_tvwgt                K_METIS__SetupGraph_tvwgt
#define SetupGraph_label                K_METIS__SetupGraph_label
#define SetupSplitGraph                 K_METIS__SetupSplitGraph
#define CreateGraph                     K_METIS__CreateGraph
#define InitGraph                       K_METIS__InitGraph
#define FreeRData                       K_METIS__FreeRData
#define FreeGraph                       K_METIS__FreeGraph

/* initpart.c */
#define Init2WayPartition		K_METIS__Init2WayPartition
#define InitSeparator			K_METIS__InitSeparator
#define RandomBisection			K_METIS__RandomBisection
#define GrowBisection			K_METIS__GrowBisection
#define McRandomBisection               K_METIS__McRandomBisection
#define McGrowBisection                 K_METIS__McGrowBisection
#define GrowBisectionNode		K_METIS__GrowBisectionNode

/* kmetis.c */
#define MlevelKWayPartitioning		K_METIS__MlevelKWayPartitioning
#define InitKWayPartitioning            K_METIS__InitKWayPartitioning

/* kwayfm.c */
#define Greedy_KWayOptimize		K_METIS__Greedy_KWayOptimize
#define Greedy_KWayCutOptimize		K_METIS__Greedy_KWayCutOptimize
#define Greedy_KWayVolOptimize          K_METIS__Greedy_KWayVolOptimize
#define Greedy_McKWayCutOptimize        K_METIS__Greedy_McKWayCutOptimize
#define Greedy_McKWayVolOptimize        K_METIS__Greedy_McKWayVolOptimize
#define IsArticulationNode              K_METIS__IsArticulationNode
#define KWayVolUpdate                   K_METIS__KWayVolUpdate

/* kwayrefine.c */
#define RefineKWay			K_METIS__RefineKWay
#define AllocateKWayPartitionMemory	K_METIS__AllocateKWayPartitionMemory
#define ComputeKWayPartitionParams	K_METIS__ComputeKWayPartitionParams
#define ProjectKWayPartition		K_METIS__ProjectKWayPartition
#define ComputeKWayBoundary		K_METIS__ComputeKWayBoundary
#define ComputeKWayVolGains             K_METIS__ComputeKWayVolGains
#define IsBalanced			K_METIS__IsBalanced

/* mcutil */
#define rvecle                          K_METIS__rvecle
#define rvecge                          K_METIS__rvecge
#define rvecsumle                       K_METIS__rvecsumle
#define rvecmaxdiff                     K_METIS__rvecmaxdiff
#define ivecle                          K_METIS__ivecle
#define ivecge                          K_METIS__ivecge
#define ivecaxpylez                     K_METIS__ivecaxpylez
#define ivecaxpygez                     K_METIS__ivecaxpygez
#define BetterVBalance                  K_METIS__BetterVBalance
#define BetterBalance2Way               K_METIS__BetterBalance2Way
#define BetterBalanceKWay               K_METIS__BetterBalanceKWay
#define ComputeLoadImbalance            K_METIS__ComputeLoadImbalance
#define ComputeLoadImbalanceDiff        K_METIS__ComputeLoadImbalanceDiff
#define ComputeLoadImbalanceDiffVec     K_METIS__ComputeLoadImbalanceDiffVec
#define ComputeLoadImbalanceVec         K_METIS__ComputeLoadImbalanceVec

/* mesh.c */
#define CreateGraphDual                 K_METIS__CreateGraphDual
#define FindCommonElements              K_METIS__FindCommonElements
#define CreateGraphNodal                K_METIS__CreateGraphNodal
#define FindCommonNodes                 K_METIS__FindCommonNodes
#define CreateMesh                      K_METIS__CreateMesh
#define InitMesh                        K_METIS__InitMesh
#define FreeMesh                        K_METIS__FreeMesh

/* meshpart.c */
#define InduceRowPartFromColumnPart     K_METIS__InduceRowPartFromColumnPart

/* minconn.c */
#define ComputeSubDomainGraph           K_METIS__ComputeSubDomainGraph
#define UpdateEdgeSubDomainGraph        K_METIS__UpdateEdgeSubDomainGraph
#define PrintSubDomainGraph             K_METIS__PrintSubDomainGraph
#define EliminateSubDomainEdges         K_METIS__EliminateSubDomainEdges
#define MoveGroupMinConnForCut          K_METIS__MoveGroupMinConnForCut
#define MoveGroupMinConnForVol          K_METIS__MoveGroupMinConnForVol

/* mincover.c */
#define MinCover			K_METIS__MinCover
#define MinCover_Augment		K_METIS__MinCover_Augment
#define MinCover_Decompose		K_METIS__MinCover_Decompose
#define MinCover_ColDFS			K_METIS__MinCover_ColDFS
#define MinCover_RowDFS			K_METIS__MinCover_RowDFS

/* mmd.c */
#define genmmd				K_METIS__genmmd
#define mmdelm				K_METIS__mmdelm
#define mmdint				K_METIS__mmdint
#define mmdnum				K_METIS__mmdnum
#define mmdupd				K_METIS__mmdupd


/* ometis.c */
#define MlevelNestedDissection		K_METIS__MlevelNestedDissection
#define MlevelNestedDissectionCC	K_METIS__MlevelNestedDissectionCC
#define MlevelNodeBisectionMultiple	K_METIS__MlevelNodeBisectionMultiple
#define MlevelNodeBisectionL2		K_METIS__MlevelNodeBisectionL2
#define MlevelNodeBisectionL1		K_METIS__MlevelNodeBisectionL1
#define SplitGraphOrder			K_METIS__SplitGraphOrder
#define SplitGraphOrderCC		K_METIS__SplitGraphOrderCC
#define MMDOrder			K_METIS__MMDOrder

/* options.c */
#define SetupCtrl                       K_METIS__SetupCtrl
#define SetupKWayBalMultipliers         K_METIS__SetupKWayBalMultipliers
#define Setup2WayBalMultipliers         K_METIS__Setup2WayBalMultipliers
#define PrintCtrl                       K_METIS__PrintCtrl
#define FreeCtrl                        K_METIS__FreeCtrl
#define CheckParams                     K_METIS__CheckParams

/* parmetis.c */
#define MlevelNestedDissectionP		K_METIS__MlevelNestedDissectionP
#define FM_2WayNodeRefine1SidedP        K_METIS__FM_2WayNodeRefine1SidedP
#define FM_2WayNodeRefine2SidedP        K_METIS__FM_2WayNodeRefine2SidedP

/* pmetis.c */
#define MlevelRecursiveBisection	K_METIS__MlevelRecursiveBisection
#define MultilevelBisect		K_METIS__MultilevelBisect
#define SplitGraphPart			K_METIS__SplitGraphPart

/* refine.c */
#define Refine2Way			K_METIS__Refine2Way
#define Allocate2WayPartitionMemory	K_METIS__Allocate2WayPartitionMemory
#define Compute2WayPartitionParams	K_METIS__Compute2WayPartitionParams
#define Project2WayPartition		K_METIS__Project2WayPartition

/* separator.c */
#define ConstructSeparator		K_METIS__ConstructSeparator
#define ConstructMinCoverSeparator	K_METIS__ConstructMinCoverSeparator

/* sfm.c */
#define FM_2WayNodeRefine2Sided         K_METIS__FM_2WayNodeRefine2Sided 
#define FM_2WayNodeRefine1Sided         K_METIS__FM_2WayNodeRefine1Sided
#define FM_2WayNodeBalance              K_METIS__FM_2WayNodeBalance

/* srefine.c */
#define Refine2WayNode			K_METIS__Refine2WayNode
#define Allocate2WayNodePartitionMemory	K_METIS__Allocate2WayNodePartitionMemory
#define Compute2WayNodePartitionParams	K_METIS__Compute2WayNodePartitionParams
#define Project2WayNodePartition	K_METIS__Project2WayNodePartition

/* stat.c */
#define ComputePartitionInfoBipartite   K_METIS__ComputePartitionInfoBipartite
#define ComputePartitionBalance		K_METIS__ComputePartitionBalance
#define ComputeElementBalance		K_METIS__ComputeElementBalance

/* timing.c */
#define InitTimers			K_METIS__InitTimers
#define PrintTimers			K_METIS__PrintTimers

/* util.c */
#define iargmax_strd                    K_METIS__iargmax_strd 
#define iargmax_nrm                     K_METIS__iargmax_nrm
#define iargmax2_nrm                    K_METIS__iargmax2_nrm
#define rargmax2                        K_METIS__rargmax2
#define InitRandom                      K_METIS__InitRandom
#define metis_rcode                     K_METIS__metis_rcode

/* wspace.c */
#define AllocateWorkSpace               K_METIS__AllocateWorkSpace                  
#define AllocateRefinementWorkSpace     K_METIS__AllocateRefinementWorkSpace
#define FreeWorkSpace                   K_METIS__FreeWorkSpace
#define wspacemalloc                    K_METIS__wspacemalloc
#define wspacepush                      K_METIS__wspacepush
#define wspacepop                       K_METIS__wspacepop
#define iwspacemalloc                   K_METIS__iwspacemalloc
#define rwspacemalloc                   K_METIS__rwspacemalloc
#define ikvwspacemalloc                 K_METIS__ikvwspacemalloc
#define cnbrpoolReset                   K_METIS__cnbrpoolReset
#define cnbrpoolGetNext                 K_METIS__cnbrpoolGetNext
#define vnbrpoolReset                   K_METIS__vnbrpoolReset
#define vnbrpoolGetNext                 K_METIS__vnbrpoolGetNext

#define errexit                         K_METIS__errexit
#endif
