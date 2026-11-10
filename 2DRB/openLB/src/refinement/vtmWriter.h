/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2026 Adrian Kummerlaender
 *
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

#ifndef REFINEMENT_VTM_WRITER_H
#define REFINEMENT_VTM_WRITER_H

#include <string>
#include <functional>
#include <memory>

#include "refinement/mesh.h"
#include "io/superVtmWriter3D.h"
#include "communication/mpiManager.h"

namespace olb {

template <typename T, typename DESCRIPTOR>
class RefinedLattice;

template <typename T, unsigned D> requires (D == 3)
class RefinedVTMwriter3D {
private:
  RefinedMesh<T, D>& _mesh;
  std::string _name;
  int _overlap;

  // Private lightweight wrapper passed to user callback
  struct WriterWrapper {
    SuperVTMwriter3D<T>& vtmWriter;
    bool isMaster;
    int locC;
    int refC;
    std::size_t timestep;

    void addFunctor(SuperF3D<T>& f) {
      vtmWriter.addFunctor(f);
    }

    void write() {
      if (isMaster) {
        vtmWriter.writePVD(timestep);
      } else {
        vtmWriter.writeVTI(timestep, locC, refC);
      }
    }
  };

public:
  RefinedVTMwriter3D(RefinedMesh<T, D>& mesh, std::string name, int overlap = 0)
    : _mesh(mesh), _name(name), _overlap(overlap)
  {}

  void createMasterFile()
  {
    if (singleton::mpi().getRank() == 0) {
      SuperVTMwriter3D<T> masterWriter(_mesh.getDecompositionBase(), _name, _overlap);
      masterWriter.createMasterFile();
    }
  }

  template <typename DESCRIPTOR, typename F>
  requires (DESCRIPTOR::d == D)
  void write(std::size_t iT, RefinedLattice<T,DESCRIPTOR>& refinedLattice, F&& vtkContext)
  {
    // 1. Write the master PVD file synchronously on Rank 0
    if (singleton::mpi().getRank() == 0) {
      SuperVTMwriter3D<T> masterWriter(_mesh.getDecompositionBase(), _name, _overlap);
      WriterWrapper wrapper{masterWriter, true, 0, 0, iT};

      vtkContext(wrapper, 0, refinedLattice.getLevel(0));
    }

    std::string name = _name;
    int overlap = _overlap;
    auto& mesh = _mesh;

    // 2. Schedule parallel asynchronous background VTI writing
    for (std::size_t iLevel = 0; iLevel < refinedLattice.getLevelCount(); ++iLevel) {
      auto& sLattice = refinedLattice.getLevel(iLevel);
      sLattice.setProcessingContext(ProcessingContext::Evaluation);

      sLattice.scheduleBackgroundOutput([name, overlap, &mesh, iLevel, iT, &refinedLattice, vtkContext](int locC) {
        auto& sLatticeLevel = refinedLattice.getLevel(iLevel);
        if (auto refC = sLatticeLevel.getCuboidDecomposition()
                                     .get(sLatticeLevel.getLoadBalancer().glob(locC))
                                     .getReferenceCuboid()) {
          if (   mesh.getDecompositionBase().get(*refC).getDeltaR()
              == sLatticeLevel.getUnitConverter().getPhysDeltaX()) {
            SuperVTMwriter3D<T> vtmWriter(name, overlap);
            WriterWrapper wrapper{vtmWriter, false, locC, *refC, iT};
            // Executes the user context asynchronously in the background thread context
            vtkContext(wrapper, iLevel, sLatticeLevel);
          }
        }
      });
    }
  }

};

template <typename T, unsigned D>
RefinedVTMwriter3D(RefinedMesh<T, D>&, std::string, int = 0) -> RefinedVTMwriter3D<T, D>;

}

#endif // REFINEMENT_VTM_WRITER_H
