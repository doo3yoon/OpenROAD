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

///////////////////////////////////////////////////////////////////////////////
//
// Description:
// Essentially a zero temperature annealer that can use a variety of
// move generators, different objectives and a cost function in order
// to improve a placement.


#include <boost/tokenizer.hpp>
#include <stack>

#include "utility.h"
#include "utl/Logger.h"
// For detailed improvement.
#include "detailed_manager.h"
#include "detailed_orient.h"
#include "detailed_random.h"
#include "detailed_segment.h"
// Detailed placement objectives.
#include "detailed_abu.h"
#include "detailed_displacement.h"
#include "detailed_global.h"
#include "detailed_hpwl.h"
#include "detailed_objective.h"
#include "detailed_vertical.h"
#include <iostream>
#include <random>


using utl::DPO;

namespace dpo {

double DetailedRandom::updateSolution(
  const std::vector<int>& bottom, 
  const std::vector<int>& left)
{
  mgrPtr_->removeAllCellsFromSegments();
  for (int i = 0; i < network_->getNumNodes(); i++) {
    Node* nd = network_->getNode(i);
    nd->setBottom(bottom[nd->getId()]);
    nd->setLeft(left[nd->getId()]);
  }

  mgrPtr_->assignCellsToSegments(mgrPtr_->getSingleHeightCells());
  mgrPtr_->resortSegments(); 

  double hpwl_x, hpwl_y;
  double best_hpwl = Utility::hpwl(network_, hpwl_x, hpwl_y);
  return best_hpwl;
}

void DetailedRandom::getSolution(
  std::vector<int>& bottom, 
  std::vector<int>& left)
{
  for (int i = 0; i < network_->getNumNodes(); i++) {
    const Node* nd = network_->getNode(i);
    bottom[nd->getId()] = nd->getBottom();
    left[nd->getId()] = nd->getLeft();
  }
}

double DetailedRandom::goGWTW()
{
  // define the variable to store the best solution (use initial solution as the best solution)
  std::vector<int> bestBottom;
  std::vector<int> bestLeft;
  bestBottom.resize(network_->getNumNodes());
  bestLeft.resize(network_->getNumNodes());
  getSolution(bestBottom, bestLeft);

  // Always use greedy approach as the baseline solution 
  // save the starting solution
  std::vector<int> startBottom = bestBottom;
  std::vector<int> startLeft = bestLeft; 
  std::vector<float> coolingRateList;
  coolingRateList.push_back(0.0);
  coolingRateList.push_back((1.0 + mgrPtr_->getCoolingRate()) / 2);
  coolingRateList.push_back(mgrPtr_->getCoolingRate());
  
  std::vector<float> initTList;
  initTList.push_back(mgrPtr_->getInitT() * 1000);
  for (int i = 0; i < 5; i++) {
    initTList.push_back(initTList.back() * 0.1);
  }
  
  int seed = mgrPtr_->getRandomSeed();
  int numSteps = mgrPtr_->getNumLoops();
  numSteps = 1000000;

  deltaCost_.resize(objectives_.size());
  initCost_.resize(objectives_.size());
  currCost_.resize(objectives_.size());
  nextCost_.resize(objectives_.size());
  for (size_t i = 0; i < objectives_.size(); i++) {
    deltaCost_[i] = 0.;
    initCost_[i] = objectives_[i]->curr();
    currCost_[i] = initCost_[i];
    nextCost_[i] = initCost_[i];

    if (objectives_[i]->getName() == "abu") {
      auto ptr = dynamic_cast<DetailedABU*>(objectives_[i]);
      if (ptr != nullptr) {
        ptr->measureABU(true);
      }
    }
  }

  double normTotalCost = eval(currCost_, expr_);
  double bestTotalCost = normTotalCost;

  // Multi start  
  for (auto coolingRate : coolingRateList) {    
    for (auto initT : initTList) {
      double oldTotalCost = updateSolution(startBottom, startLeft);
      for (size_t i = 0; i < objectives_.size(); i++) {
        deltaCost_[i] = 0.;
        initCost_[i] = objectives_[i]->curr();
        currCost_[i] = initCost_[i];
        nextCost_[i] = initCost_[i];

        if (objectives_[i]->getName() == "abu") {
          auto ptr = dynamic_cast<DetailedABU*>(objectives_[i]);
          if (ptr != nullptr) {
            ptr->measureABU(true);
          }
        }

        std::cout << "[Random Swap] : Objective : " << objectives_[i]->getName() << std::endl;
        std::cout << "[Random Swap] : Initial Cost : " << initCost_[i] << std::endl;
      }

      std::cout << "[Random Swap] : initTotalCost : " << eval(initCost_, expr_) << std::endl;

      //updateSolution(startBottom, startLeft);
      mgrPtr_->checkOverlapInSegments();
      std::cout << "[Random Swap] : oldTotalCost : " << oldTotalCost << std::endl;
      double currTotalCost = worker(coolingRate, initT, numSteps, seed++, normTotalCost);
      std::cout << "[Random Swap] : currTotalCost : " << currTotalCost << std::endl;
      if (currTotalCost < bestTotalCost) {
        bestTotalCost = currTotalCost;
        getSolution(bestBottom, bestLeft);
      }

      if (coolingRate <= 0.0) {
        break;
      }
    }
  }
  
  std::cout << "[Random Swap] : Best Total Cost : " << bestTotalCost << std::endl;
  std::cout << "[Random Swap] : currTotalCost : " << normTotalCost << std::endl;

  // update the best solution
  updateSolution(bestBottom, bestLeft);
  return ((normTotalCost - bestTotalCost) / normTotalCost);
}

// report the HPWL
double DetailedRandom::worker(
  float coolingRate,  // if coolingRate <= 0.0,  then use greedy
  float initT,
  int numSteps,
  int seed, // random seed
  double normTotalCost)
{
  std::mt19937 rand_gen(seed);
  std::uniform_real_distribution<float> distribution(0.0, 1.0);
  bool greedyFlag = coolingRate <= 0.0;
  float T = initT;

  std::cout << "[Random Swap] : numSteps : " << numSteps << std::endl;
  std::cout << "[Random Swap] : Initial Temperature : " << initT << std::endl;
  std::cout << "[Random Swap] : Cooling Rate : " << coolingRate << std::endl;
  std::cout << "[Random Swap] : Random Seed : " << seed << std::endl;
  std::cout << "[Random Swap] : greedyFlag : " << greedyFlag << std::endl;
  std::cout << "[Random Swap] : normTotalCost : " << normTotalCost << std::endl;

  // Collect candidate cells.
  collectCandidates();
  mgrPtr_->shuffle(candidates_);

  double currTotalCost = normTotalCost;

  std::vector<int> gen_count(generators_.size());
  std::fill(gen_count.begin(), gen_count.end(), 0);
  
  for (int attempt = 0; attempt < numSteps; attempt++) {
    // Pick a generator at random.
    int g = (int) mgrPtr_->getRandom(generators_.size());
    ++gen_count[g];
    // Generate a move list.
    if (generators_[g]->generate(mgrPtr_, candidates_) == false) {
      // Failed to generate anything so just move on to the next attempt.
      T *= coolingRate;
      continue;
    }

    // The generator has provided a successful move which is stored in the
    // manager.  We need to evaluate that move to see if we should accept
    // or reject it.  Scan over the objective functions and use the move
    // information to compute the weighted deltas; an overall weighted delta
    // better than zero implies improvement.
    for (size_t i = 0; i < objectives_.size(); i++) {
      // XXX: NEED TO WEIGHT EACH OBJECTIVE!
      double change = objectives_[i]->delta(mgrPtr_->getNMoved(),
                                            mgrPtr_->getMovedNodes(),
                                            mgrPtr_->getCurLeft(),
                                            mgrPtr_->getCurBottom(),
                                            mgrPtr_->getCurOri(),
                                            mgrPtr_->getNewLeft(),
                                            mgrPtr_->getNewBottom(),
                                            mgrPtr_->getNewOri());

      deltaCost_[i] = change;
      nextCost_[i] = currCost_[i] - deltaCost_[i];  // -delta is +ve is less.
    }
    
    const double nextTotalCost = eval(nextCost_, expr_);
    const double deltaTotalCost = nextTotalCost - currTotalCost;
    T *= coolingRate;
    bool acceptFlag = false;
    if (deltaTotalCost < 0) {
      acceptFlag = true;
    } else if (greedyFlag == false) {
      const float num = distribution(rand_gen);
      const float prob = std::exp(-deltaTotalCost / normTotalCost / T);  
      if (num <= prob) {
        acceptFlag = true;
      }
    }      
      
    if (acceptFlag) {
      mgrPtr_->acceptMove();
      for (auto objective : objectives_) {
        objective->accept();
      }
        
      for (size_t i = 0; i < objectives_.size(); i++) {
        currCost_[i] = nextCost_[i];
      }
      currTotalCost = nextTotalCost;
    } else {
      mgrPtr_->rejectMove();
      for (auto objective : objectives_) {
        objective->reject();
      }
    }
  }     
  
  for (auto generator : generators_) {
    generator->stats();
  }

  for (size_t i = 0; i < objectives_.size(); i++) {
    double scratch = objectives_[i]->curr();
    nextCost_[i] = scratch;  // Temporary.

    if (objectives_[i]->getName() == "abu") {
      auto ptr = dynamic_cast<DetailedABU*>(objectives_[i]);
      if (ptr != nullptr) {
        ptr->measureABU(true);
      }
    }
  }
  
  const double nextTotalCost = eval(nextCost_, expr_);
  std::cout << "[Random Swap] : End of pass, Total cost is " << nextTotalCost << std::endl;
  return currTotalCost;
}

} // namespace dpo






