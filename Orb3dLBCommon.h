//Abhishek - Extract common code from Orb and mslb_notopo
#ifndef _ORB3DHELPER_H
#define _ORB3DHELPER_H

#include <charm++.h>
#include <charm++.h>
#include "cklists.h"
#include "ParallelGravity.h"
#include "TopoManager.h"

#include "Refiner.h"
#include "MapStructures.h"
#include "TaggedVector3D.h"
#include "Vector3D.h"
#include "CentralLB.h"
#define  ORB3DLB_NOTOPO_DEBUG(X)
// #define  ORB3DLB_NOTOPO_DEBUG(X) CkPrintf X

void Orb_PrintLBStats(BaseLB::LDStats *stats, int numobjs);
void write_LB_particles(BaseLB::LDStats* stats, const char *achFileName, bool bFrom);

/// @brief Hold information about Pe load and number of objects.
class PeInfo {
  public:
  int idx;
  double load;
  double items;
  PeInfo(int id, double ld, int it) : idx(id), load(ld), items(it) {}
};

/// @brief Utility class for sorting processor loads.
class ProcLdGreater {
  public:
  bool operator()(PeInfo& p1, PeInfo& p2) {
    // This can be done based on load or number of tps assigned to a PE
    return (p1.load > p2.load);
  }
};

/// @brief Common methods among Orb3d class load balancers.
class Orb3dCommon{
  // pointer to stats->to_proc
  protected:		
    decltype(BaseLB::LDStats::to_proc) *mapping;
    decltype(BaseLB::LDStats::from_proc) *from;

    CkVec<float> procload;

    /// Take into account memory constraints by limiting the number of pieces
    /// per processor.
    double maxPieceProc;

    /// index of first processor of the group we are considering
    int nextProc;

    // Greedy strategy to assign TreePieces to PEs on a node.
    void orbPePartition(vector<Event> *events, vector<OrbObject> &tp, int node,
        BaseLB::LDStats *stats) {

      std::vector<PeInfo> peinfo;
      float totalLoad = 0.0;
      int firstProc = CkNodeFirst(node);
      int lastProc = firstProc + CkNodeSize(node) - 1;
      for (int i = firstProc; i <= lastProc; i++) {
        peinfo.push_back(PeInfo(i, 0.0, 0));
      }
      // Make a heap of processors belonging to this node
      std::make_heap(peinfo.begin(), peinfo.end(), ProcLdGreater());

      int nextProc;
      for(int i = 0; i < events[XDIM].size(); i++){
        Event &ev = events[XDIM][i];
        OrbObject &orb = tp[ev.owner];

        // Pop the least loaded PE from the heap and assign TreePiece to it
        PeInfo p = peinfo.front();
        pop_heap(peinfo.begin(), peinfo.end(), ProcLdGreater());
        peinfo.pop_back();

        nextProc = p.idx;

        if(orb.numParticles > 0){
          (*mapping)[orb.lbindex] = nextProc;
          procload[nextProc] += ev.load;
          p.load += ev.load;
          p.items += 1;
          totalLoad += ev.load;
        } else{
          int fromPE = (*from)[orb.lbindex];
          procload[fromPE] += ev.load;
        }

        peinfo.push_back(p);
        push_heap(peinfo.begin(), peinfo.end(), ProcLdGreater());
      }
    }

/// @brief Recursively partition treepieces among processors by
/// bisecting the load in orthogonal directions.
/// @param events Array of three (1 per dimension) Event vectors.
/// These are separate in each dimension for easy sorting.
/// @param box Spatial bounding box
/// @param nprocs Number of processors over which to partition the
/// Events. N.B. if node_partition is true, then this is the number of nodes.
/// @param tp Vector of TreePiece data.
    void orbPartition(vector<Event> *events, OrientedBox<float> &box, int nprocs,
        vector<OrbObject> & tp, BaseLB::LDStats *stats,
        bool node_partition=false){

        ORB3DLB_NOTOPO_DEBUG(("partition events %d %d %d nprocs %d\n", 
          events[XDIM].size(),
          events[YDIM].size(),
          events[ZDIM].size(),
          nprocs
                              ));
      int numEvents = events[XDIM].size();
      CkAssert(numEvents == events[YDIM].size());
      CkAssert(numEvents == events[ZDIM].size());

      if(numEvents == 0)
	return;

      if(nprocs == 1){
        ORB3DLB_NOTOPO_DEBUG(("base: assign %d tps to proc %d\n", numEvents, nextProc));
        if (!stats->procs[nextProc].available) {
          nextProc++;
          return;
        }

        // If we are doing orb partition at the node level, then call
        // orbPePartition to assign the treepieces to the PEs belonging to the node.
        if (node_partition) {
          orbPePartition(events, tp, nextProc, stats);
        } else {
          // direct assignment of tree pieces to processors
          //if(numEvents > 0) CkAssert(nprocs != 0);
          float totalLoad = 0.0;
          for(int i = 0; i < events[XDIM].size(); i++){
            Event &ev = events[XDIM][i];
            OrbObject &orb = tp[ev.owner];
            if(orb.numParticles > 0){
              (*mapping)[orb.lbindex] = nextProc;
              totalLoad += ev.load;
            }
            else{
              int fromPE = (*from)[orb.lbindex];
              if (fromPE < 0 || fromPE >= procload.size()) {
                CkPrintf("[%d] trying to access fromPe %d nprocs %lu\n", CkMyPe(), fromPE, procload.size());
                CkAbort("Trying to access a PE which is outside the range\n");
              }
              procload[fromPE] += ev.load;
            }
          }
          procload[nextProc] += totalLoad;
        }

        if(numEvents > 0) nextProc++;
        return;
      }

      // find longest dimension

      int longestDim = XDIM;
      float longestDimLength = box.greater_corner[longestDim] - box.lesser_corner[longestDim];
      for(int i = YDIM; i <= ZDIM; i++){
        float thisDimLength = box.greater_corner[i]-box.lesser_corner[i];
        if(thisDimLength > longestDimLength){
          longestDimLength = thisDimLength;
          longestDim = i;
        }
      }

      ORB3DLB_NOTOPO_DEBUG(("dimensions %f %f %f longest %d\n", 
          box.greater_corner[XDIM]-box.lesser_corner[XDIM],
          box.greater_corner[YDIM]-box.lesser_corner[YDIM],
          box.greater_corner[ZDIM]-box.lesser_corner[ZDIM],
          longestDim
                            ));

      int nlprocs = nprocs/2;
      int nrprocs = nprocs-nlprocs;

      float ratio = (1.0*nlprocs)/(1.0*(nlprocs+nrprocs));

      // sum background load on each side of the processor split
      float bglprocs = 0.0;
      for(int np = nextProc; np < nextProc + nlprocs; np++)
        bglprocs += stats->procs[np].bg_walltime;
      float bgrprocs = 0.0;
      for(int np = nextProc + nlprocs; np < nextProc + nlprocs + nrprocs; np++)
        bgrprocs += stats->procs[np].bg_walltime;

      ORB3DLB_NOTOPO_DEBUG(("nlprocs %d nrprocs %d ratio %f\n", nlprocs, nrprocs, ratio));

      int splitIndex = partitionRatioLoad(events[longestDim],ratio,bglprocs,
                                          bgrprocs);
      if(splitIndex == numEvents) {
        ORB3DLB_NOTOPO_DEBUG(("evenly split 0 load\n"));
        splitIndex = splitIndex/2;
      }
      int nleft = splitIndex;
      int nright = numEvents-nleft;

#if 0
      if(nright < nrprocs) {  // at least one piece per processor
        nright = nrprocs;
        nleft = splitIndex = numEvents-nright;
        CkAssert(nleft >= nlprocs);
      }
      else if(nleft < nlprocs) {
        nleft = splitIndex = nlprocs;
        nright = numEvents-nleft;
        CkAssert(nright >= nrprocs);
      }
#endif

      if(nleft > nlprocs*maxPieceProc) {
	  nleft = splitIndex = (int) (nlprocs*maxPieceProc);
	  nright = numEvents-nleft;
	  }
      else if (nright > nrprocs*maxPieceProc) {
	  nright = (int) (nrprocs*maxPieceProc);
	  nleft = splitIndex = numEvents-nright;
	  }
      CkAssert(splitIndex >= 0);
      CkAssert(splitIndex < numEvents);

      OrientedBox<float> leftBox;
      OrientedBox<float> rightBox;

      leftBox = rightBox = box;
      float splitPosition = events[longestDim][splitIndex].position;
      leftBox.greater_corner[longestDim] = splitPosition;
      rightBox.lesser_corner[longestDim] = splitPosition;

      // classify events
      for(int i = 0; i < splitIndex; i++){
        Event &ev = events[longestDim][i];
        CkAssert(ev.owner >= 0);
        CkAssert(tp[ev.owner].partition == INVALID_PARTITION);
        tp[ev.owner].partition = LEFT_PARTITION;
      }
      for(int i = splitIndex; i < numEvents; i++){
        Event &ev = events[longestDim][i];
        CkAssert(ev.owner >= 0);
        CkAssert(tp[ev.owner].partition == INVALID_PARTITION);
        tp[ev.owner].partition = RIGHT_PARTITION;
      }

      vector<Event> leftEvents[NDIMS];
      vector<Event> rightEvents[NDIMS];

      for(int i = 0; i < NDIMS; i++){
        if(i == longestDim){ 
          leftEvents[i].resize(nleft);
          rightEvents[i].resize(nright);
        }
        else{
          leftEvents[i].reserve(nleft);
          rightEvents[i].reserve(nright);
        }
      }
      // copy events of split dimension
      memcpy(&leftEvents[longestDim][0],&events[longestDim][0],sizeof(Event)*nleft);
      memcpy(&rightEvents[longestDim][0],&events[longestDim][splitIndex],sizeof(Event)*nright);

      // copy events of other dimensions
      for(int i = XDIM; i <= ZDIM; i++){
        if(i == longestDim) continue;
        for(int j = 0; j < numEvents; j++){
          Event &ev = events[i][j];
          CkAssert(ev.owner >= 0);
          OrbObject &orb = tp[ev.owner];
          CkAssert(orb.partition != INVALID_PARTITION);
          if(orb.partition == LEFT_PARTITION) leftEvents[i].push_back(ev);
          else if(orb.partition == RIGHT_PARTITION) rightEvents[i].push_back(ev);
        }
      }

      // cleanup
      // next, reset the ownership information in the
      // OrbObjects, so that the next invocation may use
      // the same locations for its book-keeping
      vector<Event> &eraseVec = events[longestDim];
      for(int i = 0; i < numEvents; i++){
        Event &ev = eraseVec[i];
        CkAssert(ev.owner >= 0);
        OrbObject &orb = tp[ev.owner];
        CkAssert(orb.partition != INVALID_PARTITION);
        orb.partition = INVALID_PARTITION;
      }

      // free events from parent node,
      // since they are not needed anymore
      // (we have partition all events into the
      // left and right event subsets)
      for(int i = 0; i < NDIMS; i++){
        //events[i].free();
        vector<Event>().swap(events[i]);
      }
      orbPartition(leftEvents,leftBox,nlprocs,tp, stats, node_partition);
      orbPartition(rightEvents,rightBox,nrprocs,tp, stats, node_partition);
    }

/// @brief Prepare structures for the ORB partition.
/// @param tpEvents Array of 3 (1 per dimension) Event vectors.
/// @param box Reference to bounding box (set here).
/// @param numobjs Number of tree pieces to partition.
/// @param stats Data from the load balancing framework.
/// @param node_partition Are we partitioning on nodes.
    void orbPrepare(vector<Event> *tpEvents, OrientedBox<float> &box, int
    numobjs, BaseLB::LDStats * stats, bool node_partition=false){

      int nmig = stats->n_migrateobjs;
      if(dMaxBalance < 1.0)
        dMaxBalance = 1.0;

      // If using node based orb partition, then the maxPieceProc is total
      // migratable objs / total number of node.
      if (node_partition) {
        maxPieceProc = dMaxBalance * nmig / CkNumNodes();
      } else {
        maxPieceProc = dMaxBalance*nmig/stats->nprocs();
      }

      if(maxPieceProc < 1.0)
        maxPieceProc = 1.01;

      CkAssert(tpEvents[XDIM].size() == numobjs);
      CkAssert(tpEvents[YDIM].size() == numobjs);
      CkAssert(tpEvents[ZDIM].size() == numobjs);

      mapping = &stats->to_proc;
      from = &stats->from_proc;

      CkPrintf("[Orb3dLB_notopo] sorting\n");
      for(int i = 0; i < NDIMS; i++){
        sort(tpEvents[i].begin(),tpEvents[i].end());
      }

      box.lesser_corner.x = tpEvents[XDIM][0].position;
      box.lesser_corner.y = tpEvents[YDIM][0].position;
      box.lesser_corner.z = tpEvents[ZDIM][0].position;

      box.greater_corner.x = tpEvents[XDIM][numobjs-1].position;
      box.greater_corner.y = tpEvents[YDIM][numobjs-1].position;
      box.greater_corner.z = tpEvents[ZDIM][numobjs-1].position;

      nextProc = 0;

      procload.resize(stats->nprocs());
      for(int i = 0; i < stats->nprocs(); i++){
        procload[i] = stats->procs[i].bg_walltime;
      }

    }

    void refine(BaseLB::LDStats *stats, int numobjs){
#ifdef DO_REFINE
      int *from_procs = Refiner::AllocProcs(stats->nprocs(), stats);
      int *to_procs = Refiner::AllocProcs(stats->nprocs(), stats);
#endif

      for(int i = 0; i < numobjs; i++){
#ifdef DO_REFINE
        int pe = stats->to_proc[i];
        from_procs[i] = pe;
        to_procs[i] = pe;
#endif
      }

      int numRefineMigrated = 0;
#ifdef DO_REFINE
      CkPrintf("[orb3dlb_notopo] refine\n");
      Refiner refiner(1.050);
      refiner.Refine(stats->nprocs(),stats,from_procs,to_procs);

      for(int i = 0; i < numobjs; i++){
        if(to_procs[i] != from_procs[i]) numRefineMigrated++;
        stats->to_proc[i] = to_procs[i];
      }
#endif

      Orb_PrintLBStats(stats, numobjs);

#ifdef DO_REFINE
      // Free the refine buffers
      Refiner::FreeProcs(from_procs);
      Refiner::FreeProcs(to_procs);
#endif

    }

/// @brief Given a vector of Events, find a split that partitions them
/// into two partitions with a given ratio of loads.
/// @param events Vector of Events to split
/// @param ratio Target ratio of loads in left partition to total load.
/// @param bglp Background load on the left processors.
/// @param bgrp Background load on the right processors.
/// @return Starting index of right partition.
///
    int partitionRatioLoad(vector<Event> &events, float ratio, float bglp, float bgrp){

      float approxBgPerEvent = (bglp + bgrp) / events.size();
      float totalLoad = bglp + bgrp;
      for(int i = 0; i < events.size(); i++){
        totalLoad += events[i].load;
      }
      //CkPrintf("************************************************************\n");
      //CkPrintf("partitionEvenLoad start %d end %d total %f\n", tpstart, tpend, totalLoad);
      float perfectLoad = ratio * totalLoad;
      ORB3DLB_NOTOPO_DEBUG(("partitionRatioLoad bgl %f bgr %f\n",
                            bglp, bgrp));
      int splitIndex = 0;
      float prevLoad = 0.0;
      float leftLoadAtSplit = 0.0;
      for(splitIndex = 0; splitIndex < events.size(); splitIndex++){

        leftLoadAtSplit += events[splitIndex].load + approxBgPerEvent;

        if (leftLoadAtSplit > perfectLoad) {
          if ( fabs(leftLoadAtSplit - perfectLoad) < fabs(prevLoad - perfectLoad) ) {
            splitIndex++;
          }
          else {
            leftLoadAtSplit = prevLoad;
          }
          break;
        }
        prevLoad = leftLoadAtSplit;
      }

      ORB3DLB_NOTOPO_DEBUG(("partitionEvenLoad mid %d lload %f rload %f ratio %f\n", splitIndex, leftLoadAtSplit, totalLoad - leftLoadAtSplit, leftLoadAtSplit / totalLoad));
      return splitIndex;
    }

}; //end class

#endif
