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

#pragma once

#include <memory>
#include <vector>

#include "detailed_generator.h"
#include "rectangle.h"

namespace dpo {

class Architecture;
class DetailedMgr;
class Edge;
class Network;
class RoutingParams;

class DetailedVerticalSwap;
class DetailedGlobalSwap;


// CLASSES ===================================================================
class DetailedMixedSwap : public DetailedGenerator
{
 public:
  DetailedMixedSwap(Architecture* arch, Network* network, RoutingParams* rt, double global_ratio);
  DetailedMixedSwap();

  // Intefaces for scripting.
  void run(DetailedMgr* mgrPtr, const std::string& command);
  void run(DetailedMgr* mgrPtr, const std::vector<std::string>& args);

  // Interface for move generation.
  void init(DetailedMgr* mgr) override;
  void stats() override;
 
  void mixedSwap();

  void mixedSwapGWTW();

  bool getRange(Node* nd, Rectangle& nodeBbox);

  bool calculateEdgeBB(const Edge* ed,
                       const Node* nd,
                       Rectangle& bbox);

  double delta(Node* ndi, double new_x, double new_y);
  double delta(Node* ndi, Node* ndj);

  bool generateVerticalSwap(Node* ndi);
  bool generateGlobalSwap(Node* ndi);
  bool generateMove(DetailedMgr* mgr,
                    std::vector<Node*>& candidates,
                    bool isVertical);
  bool generate(DetailedMgr* mgr, std::vector<Node*>& candidates) override 
  { return true; };

  double updateSolution(const std::vector<int>& bottom, 
                        const std::vector<int>& left, 
                        const std::vector<int>& segId);

  void getSolution(std::vector<int>& bottom,  std::vector<int>& left, std::vector<int>& segId);
  
  // if coolingRate <= 0.0,  then use greedy
  double worker(float coolingRate,  float initT, int numSteps,
    float globalVerticalRatio, int seed,  double normHPWL);

  void randomKickMove();
 
  bool generateRandomMove(std::vector<Node*>& candidates);
  
 private:
  // Standard stuff.
  DetailedMgr* mgr_;
  Architecture* arch_;
  Network* network_;
  RoutingParams* rt_;

  double global_ratio_ = 1.0;
  double random_ratio_ = 0.10;

  // Other.
  int skipNetsLargerThanThis_;
  std::vector<int> edgeMask_;
  int traversal_;

  std::vector<double> xpts_;
  std::vector<double> ypts_;

  // For use as a move generator.
  int attempts_;
  int moves_;
  int swaps_;
};

}  // namespace dpo
