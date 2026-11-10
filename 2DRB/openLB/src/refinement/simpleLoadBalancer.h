/*  This file is part of the OpenLB library
*
 *  Copyright (C) 2026 Kilian Mahler
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/

#ifndef OLB_REFINEMENT_SIMPLE_LOAD_BALANCER_H
#define OLB_REFINEMENT_SIMPLE_LOAD_BALANCER_H
#include "communication/loadBalancer.h"
#include "geometry/cuboidDecomposition.h"

namespace olb {

template<typename T>
class PredefinedLoadBalancer: public LoadBalancer<T> {
#ifdef PARALLEL_MODE_MPI
  singleton::MpiNonBlockingHelper _mpiNbHelper;
#endif

public:

  PredefinedLoadBalancer(CuboidDecomposition3D<T>& decomposition)
  {
    this->_glob.clear();
    int rank = 0;
    int size = 1;
    int nC = decomposition.size();
  #ifdef PARALLEL_MODE_MPI
    rank = singleton::mpi().getRank();
    size = util::max<int>(singleton::mpi().getSize(), 1);
  #endif


    std::vector<int> cuboidToThread(nC);
    std::vector<Cuboid3D<T>> cuboids = decomposition.cuboids();

    if (size == 1) {
      for (int i = 0; i < nC; ++i) {
        this->_glob.push_back(i);
        this->_loc[i] = i;
        this->_rank[i] = 0;
      };
      this->_size = nC;
      return;
    }
    if (rank == 0) {
      int count = 0;
      for (int i = 0; i < nC; ++i) {
        if (cuboids[i].getPreCalculatedProcess() == std::nullopt) {
          std::cout << cuboids[i].getOrigin() << std::endl;
          throw std::logic_error("not all cuboids processes are set");
        }
        if (cuboids[i].getPreCalculatedProcess().value() == 0) {
          this->_glob.push_back(i);
          this->_loc[i] = count;
          count++;
        };
        this->_rank[i] = cuboids[i].getPreCalculatedProcess().value();
        cuboidToThread[i] = cuboids[i].getPreCalculatedProcess().value();
      }
      this->_size = count;
    }

    // Send all threads their number of cuboids

  #ifdef PARALLEL_MODE_MPI
    if (rank == 0) {
      // Send all threads their respective cuboids
      _mpiNbHelper.free();
      _mpiNbHelper.allocate(size-1);
      for (int i = 1; i < size; i++) {
        singleton::mpi().iSend(&cuboidToThread.front(),
                               nC, i, &_mpiNbHelper.get_mpiRequest()[i-1], 0);
      }
      singleton::mpi().waitAll(_mpiNbHelper);
    }
    else {
      int *tmpCuboids = new int[nC];
      singleton::mpi().receive(tmpCuboids, nC, 0, 0);
      int count = 0;
      for (int i = 0; i < nC; ++i) {
        if (tmpCuboids[i] == rank) {
          this->_glob.push_back(i);
          this->_loc[i] = count;
          count++;
        };
        this->_rank[i] = tmpCuboids[i];
      }
      delete[] tmpCuboids;
      this->_size = count;
    }
#endif
  }

  PredefinedLoadBalancer(CuboidDecomposition2D<T>& decomposition)
  {
    this->_glob.clear();
    int rank = 0;
    int size = 1;
    int nC = decomposition.size();
#ifdef PARALLEL_MODE_MPI
    rank = singleton::mpi().getRank();
    size = util::max<int>(singleton::mpi().getSize(), 1);
#endif


    std::vector<int> cuboidToThread(nC);
    std::vector<Cuboid2D<T>> cuboids = decomposition.cuboids();

    if (size == 1) {
      for (int i = 0; i < nC; ++i) {
        this->_glob.push_back(i);
        this->_loc[i] = i;
        this->_rank[i] = 0;
      };
      this->_size = nC;
      return;
    }
    if (rank == 0) {
      int count = 0;
      for (int i = 0; i < nC; ++i) {
        if (!cuboids[i].getPreCalculatedProcess().has_value()) {
          throw std::logic_error("not all cuboids processes are set");
        }
        if (cuboids[i].getPreCalculatedProcess().value() == 0) {
          this->_glob.push_back(i);
          this->_loc[i] = count;
          count++;
        };
        this->_rank[i] = cuboids[i].getPreCalculatedProcess().value();
        cuboidToThread[i] = cuboids[i].getPreCalculatedProcess().value();
      }
      this->_size = count;
    }

    // Send all threads their number of cuboids

#ifdef PARALLEL_MODE_MPI
    if (rank == 0) {
      // Send all threads their respective cuboids
      _mpiNbHelper.free();
      _mpiNbHelper.allocate(size-1);
      for (int i = 1; i < size; i++) {
        singleton::mpi().iSend(&cuboidToThread.front(),
                               nC, i, &_mpiNbHelper.get_mpiRequest()[i-1], 0);
      }
      singleton::mpi().waitAll(_mpiNbHelper);
    }
    else {
      int *tmpCuboids = new int[nC];
      singleton::mpi().receive(tmpCuboids, nC, 0, 0);
      int count = 0;
      for (int i = 0; i < nC; ++i) {
        if (tmpCuboids[i] == rank) {
          this->_glob.push_back(i);
          this->_loc[i] = count;
          count++;
        };
        this->_rank[i] = tmpCuboids[i];
      }
      delete[] tmpCuboids;
      this->_size = count;
    }
#endif
  }

};


}

#endif //OLB_REFINEMENT_SIMPLE_LOAD_BALANCER_H
