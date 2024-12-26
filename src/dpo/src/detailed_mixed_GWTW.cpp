///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (c) 2024, Zhiang Wang (UCSD)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Currently do not enable multiple threads

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
#include "detailed_mixed.h"

#include <boost/tokenizer.hpp>

#include "detailed_global.h"
#include "detailed_vertical.h"
#include "detailed_hpwl.h"
#include "detailed_manager.h"
#include "detailed_orient.h"
#include "detailed_segment.h"
#include "rectangle.h"
#include "utility.h"
#include "utl/Logger.h"
#include <random>
#include <iostream>
#include <fstream>

using utl::DPO;

namespace dpo {

double DetailedMixedSwap::updateSolution(
  const std::vector<int>& bottom, 
  const std::vector<int>& left,
  const std::vector<int>& segId)
{
  mgr_->removeAllCellsFromSegments();
  for (int i = 0; i < network_->getNumNodes(); i++) {
    Node* nd = network_->getNode(i);
    nd->setBottom(bottom[nd->getId()]);
    nd->setLeft(left[nd->getId()]);
  }

  mgr_->assignCellsToSegments(mgr_->getSingleHeightCells());
  mgr_->resortSegments();

  double hpwl_x, hpwl_y;
  double best_hpwl = Utility::hpwl(network_, hpwl_x, hpwl_y);
  return best_hpwl;  
}


void DetailedMixedSwap::getSolution(
  std::vector<int>& bottom, 
  std::vector<int>& left,
  std::vector<int>& segId)
{
  for (int i = 0; i < network_->getNumNodes(); i++) {
    const Node* nd = network_->getNode(i);
    bottom[nd->getId()] = nd->getBottom();
    left[nd->getId()] = nd->getLeft();
  }
}


void DetailedMixedSwap::mixedSwapGWTW()
{
  // define the variable to store the best solution (use initial solution as the best solution)
  std::vector<int> bestBottom;
  std::vector<int> bestLeft;
  std::vector<int> bestSegId;
  bestBottom.resize(network_->getNumNodes());
  bestLeft.resize(network_->getNumNodes());
  bestSegId.resize(network_->getNumNodes());
  getSolution(bestBottom, bestLeft, bestSegId);

  double hpwl_x, hpwl_y;
  double best_hpwl = Utility::hpwl(network_, hpwl_x, hpwl_y);

  // Always use greedy approach as the baseline solution 
  // save the starting solution
  std::vector<int> startBottom = bestBottom;
  std::vector<int> startLeft = bestLeft; 
  std::vector<int> startSegId = bestSegId;
  std::vector<float> coolingRateList;
  coolingRateList.push_back(0.0);
  coolingRateList.push_back(mgr_->getCoolingRate());

  // 0.0 means vertical only, 1.0 means global only
  std::vector<float> globalVerticalRatioList = {0.0, 0.25, 0.5, 0.75, 1.0};
  int seed = mgr_->getRandomSeed();
  float initT = mgr_->getInitT();
  double normHPWL = best_hpwl;
  int numSteps = mgr_->getNumLoops();
  numSteps = 10000;

  // Multi start  
  for (auto coolingRate : coolingRateList) {
    for (auto globalVerticalRatio : globalVerticalRatioList) {
      mgr_->checkOverlapInSegments();
      updateSolution(startBottom, startLeft, startSegId);
      mgr_->checkOverlapInSegments();
      double hpwl = worker(coolingRate, initT, numSteps, globalVerticalRatio, seed, normHPWL);
      if (hpwl < best_hpwl) {
        best_hpwl = hpwl;
        getSolution(bestBottom, bestLeft, bestSegId);
      } 
    } 
  }

  // update the best solution
  updateSolution(bestBottom, bestLeft, bestSegId);

  bool LSMCFlag = false;
  int LSMCLoop = 5;
  if (LSMCFlag == true) {
    for (int i = 0; i <= LSMCLoop; i++) {
      std::cout << "[LSMC] : loop id : " << i << std::endl;
      updateSolution(bestBottom, bestLeft, bestSegId);
      // perform random performutation
      randomKickMove();
      getSolution(startBottom, startLeft, startSegId);
      std::cout << "After random kick move" << std::endl;
      mgr_->checkOverlapInSegments();
      // Using multistart SA as the local searech algorithm
      for (auto coolingRate : coolingRateList) {
        for (auto globalVerticalRatio : globalVerticalRatioList) {
          mgr_->checkOverlapInSegments();
          updateSolution(startBottom, startLeft, startSegId);
          mgr_->checkOverlapInSegments();
          double hpwl = worker(coolingRate, initT, numSteps, globalVerticalRatio, seed, normHPWL);
          if (hpwl < best_hpwl) {
            best_hpwl = hpwl;
            getSolution(bestBottom, bestLeft, bestSegId);
          }
        } 
      }
    } 
  }  

  // update the best solution
  updateSolution(bestBottom, bestLeft, bestSegId);
}


// report the HPWL
double DetailedMixedSwap::worker(
  float coolingRate,  // if coolingRate <= 0.0,  then use greedy
  float initT,
  int numSteps,
  float globalVerticalRatio,
  int seed, // random seed
  double normHPWL)
{
  std::mt19937 rand_gen(seed);
  std::uniform_real_distribution<float> distribution(0.0, 1.0);
  bool greedyFlag = coolingRate <= 0.0;
  bool verticalOnlyFlag = globalVerticalRatio == 0.0;
  bool globalOnlyFlag = globalVerticalRatio == 1.0;
  float T = initT;

  std::cout << "[Mixed Swap] : numSteps : " << numSteps << std::endl;
  std::cout << "[Mixed Swap] : Initial Temperature : " << initT << std::endl;
  std::cout << "[Mixed Swap] : Cooling Rate : " << coolingRate << std::endl;
  std::cout << "[Mixed Swap] : Global Vertical Ratio : " << globalVerticalRatio << std::endl;
  std::cout << "[Mixed Swap] : Random Seed : " << seed << std::endl;
  std::cout << "[Mixed Swap] : greedyFlag : " << greedyFlag << std::endl;
  std::cout << "[Mixed Swap] : verticalOnlyFlag : " << verticalOnlyFlag << std::endl;
  std::cout << "[Mixed Swap] : globalOnlyFlag : " << globalOnlyFlag << std::endl;
  std::cout << "[Mixed Swap] : normHPWL : " << normHPWL << std::endl;

  int numVerticalMoves = 0;
  int numGlobalMoves = 0;
  for (int loopId = 0; loopId < numSteps; loopId++) {
    traversal_ = 0;
    edgeMask_.resize(network_->getNumEdges());
    std::fill(edgeMask_.begin(), edgeMask_.end(), 0);
    mgr_->resortSegments();
    T *= coolingRate;

    // Get candidate cells.
    std::vector<Node*> candidates = mgr_->getSingleHeightCells();
    mgr_->shuffle(candidates);

    // Wirelength objective.
    DetailedHPWL hpwlObj(network_);
    hpwlObj.init(mgr_, nullptr);  // Ignore orientation.
 
    double currHpwl = hpwlObj.curr();
    bool verticalFlag = false;
    // check the condition for vertical or global move
    if (verticalOnlyFlag) {
      verticalFlag = true;
    } else if (globalOnlyFlag) {
      verticalFlag = false;
    } else { // avoid the overflow issue
      if (loopId - numSteps * globalVerticalRatio > 30) {
        verticalFlag = true;
      } else if (loopId - numSteps * globalVerticalRatio < -30){
        verticalFlag = false;
      } else {
        float opNum = distribution(rand_gen);
        const float vertical_ratio = 1.0 / (1.0 + std::exp( loopId - numSteps * globalVerticalRatio));
        verticalFlag = (opNum < vertical_ratio) ? false : true;
      }
    }

    // Consider each candidate cell once.
    for (auto ndi : candidates) {      
      if (verticalFlag) {
        if (!generateVerticalSwap(ndi)) {
          continue;
        }
        numVerticalMoves++;
      } else {
        if (!generateGlobalSwap(ndi)) {
          continue;
        }
        numGlobalMoves++;
      }

      int originId = ndi->getId();

      double delta = hpwlObj.delta(mgr_->getNMoved(),
                                  mgr_->getMovedNodes(),
                                  mgr_->getCurLeft(),
                                  mgr_->getCurBottom(),
                                  mgr_->getCurOri(),
                                  mgr_->getNewLeft(),
                                  mgr_->getNewBottom(),
                                  mgr_->getNewOri());
    
      bool acceptFlag = false;
      if (delta > 0) {
        acceptFlag = true;
      } else if (greedyFlag == false) {
        const float num = distribution(rand_gen);
        const float prob = std::exp(delta / normHPWL / T);
        acceptFlag = (num <= prob) ? true : false;
      }

      if (acceptFlag) {
        mgr_->acceptMove();
        currHpwl = currHpwl - delta;
      } else {
        mgr_->rejectMove();
      }

      //int newId = mgr_->getMovedNodes()[originId];
      int newId = ndi->getId();
      if (originId != newId) {
        std::cout << "[Mixed Swap] : originId : " << originId << " newId : " << newId << std::endl;
      }

    }
  }

  std::cout << "[Mixed Swap] number of vertical moves : " << numVerticalMoves << std::endl;
  std::cout << "[Mixed Swap] number of global moves : " << numGlobalMoves << std::endl;


  double hpwl_x, hpwl_y;
  double best_hpwl = Utility::hpwl(network_, hpwl_x, hpwl_y);
  std::cout << "[Mixed Swap] Final HPWL : " << best_hpwl << std::endl;
  return best_hpwl;
}


bool DetailedMixedSwap::generateRandomMove(std::vector<Node*>& candidates)
{
  const int ydim = mgr_->getNumSingleHeightRows();
  double xwid = arch_->getRow(0)->getSiteSpacing();
  const int xdim
      = std::max(0, (int) ((arch_->getMaxX() - arch_->getMinX()) / xwid));

  xwid = (arch_->getMaxX() - arch_->getMinX()) / (double) xdim;
  double ywid = (arch_->getMaxY() - arch_->getMinY()) / (double) ydim;

  Node* ndi = candidates[mgr_->getRandom(candidates.size())];
  const int spanned_i = arch_->getCellHeightInRows(ndi);
  if (spanned_i != 1) {
    return false;
  }
  
  // Segments for the source.
  const std::vector<DetailedSeg*>& segs_i
      = mgr_->getReverseCellToSegs(ndi->getId());
  if (segs_i.size() != 1) {
    mgr_->getLogger()->error(
        DPO, 386, "Only working with single height cells currently.");
  }

  // For the window size.  This should be parameterized.
  const int rly = 10;
  const int rlx = 10;

  const int tries = 5;
  for (int t = 1; t <= tries; t++) {
    // Position of the source.
    const double yi = ndi->getBottom() + 0.5 * ndi->getHeight();
    const double xi = ndi->getLeft() + 0.5 * ndi->getWidth();

    // Segment for the source.
    const int si = segs_i[0]->getSegId();

    // Random position within a box centered about (xi,yi).
    const int grid_xi = std::min(
        xdim - 1, std::max(0, (int) ((xi - arch_->getMinX()) / xwid)));
    const int grid_yi = std::min(
        ydim - 1, std::max(0, (int) ((yi - arch_->getMinY()) / ywid)));

    const int rel_x = mgr_->getRandom(2 * rlx + 1);
    const int rel_y = mgr_->getRandom(2 * rly + 1);

    const int grid_xj
        = std::min(xdim - 1, std::max(0, (grid_xi - rlx + rel_x)));
    const int grid_yj
        = std::min(ydim - 1, std::max(0, (grid_yi - rly + rel_y)));

    // Position of the destination.
    const double xj = arch_->getMinX() + grid_xj * xwid;
    double yj = arch_->getMinY() + grid_yj * ywid;

    // Row and segment for the destination.
    int rj = (int) ((yj - arch_->getMinY()) / mgr_->getSingleRowHeight());
    rj = std::min(mgr_->getNumSingleHeightRows() - 1, std::max(0, rj));
    yj = arch_->getRow(rj)->getBottom();
    int sj = -1;
    for (int s = 0; s < mgr_->getNumSegsInRow(rj); s++) {
      const DetailedSeg* segPtr = mgr_->getSegsInRow(rj)[s];
      if (xj >= segPtr->getMinX() && xj <= segPtr->getMaxX()) {
        sj = segPtr->getSegId();
        break;
      }
    }

    // Need to determine validity of things.
    if (sj == -1 || ndi->getRegionId() != mgr_->getSegment(sj)->getRegId()) {
      // The target segment cannot support the candidate cell.
      continue;
    }

    if (mgr_->tryMove(ndi,
                      ndi->getLeft(),
                      ndi->getBottom(),
                      si,
                      (int) std::round(xj),
                      (int) std::round(yj),
                      sj)) {
      ++moves_;
      return true;
    }
    if (mgr_->trySwap(ndi,
                      ndi->getLeft(),
                      ndi->getBottom(),
                      si,
                      (int) std::round(xj),
                      (int) std::round(yj),
                      sj)) {
      ++swaps_;
      return true;
    }
  }
  return false;
}


void DetailedMixedSwap::randomKickMove()
{
  std::vector<Node*> candidates = mgr_->getSingleHeightCells();
  mgr_->shuffle(candidates);
  int numPerturbs = candidates.size() * random_ratio_;
  std::cout << "[Mixed Swap] : number of perturbs : " << numPerturbs << std::endl;
  for (int i = 0; i < numPerturbs; i++) {
    if (generateRandomMove(candidates)) {
      mgr_->acceptMove();
    }
  }
}


} // namespace dpo

