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

#ifndef REFINEMENT_AMR_WRITER_H
#define REFINEMENT_AMR_WRITER_H

#include <string>
#include <functional>
#include <memory>
#include <cmath>
#include <fstream>

#include "refinement/mesh.h"
#include "io/superVtmWriter3D.h"
#include "io/fileName.h"
#include "communication/mpiManager.h"
#include "core/singleton.h"

namespace olb {

template <typename T, typename DESCRIPTOR>
class RefinedLattice;

template <typename T, unsigned D> requires (D == 3)
class RefinedAMRwriter {
private:
  RefinedMesh<T, D>& _mesh;
  std::string _name;
  int _overlap;
  bool _writePerLevel;

  struct WriterWrapper {
    SuperVTMwriter3D<T>& vtmWriter;
    int locC;
    int refC;
    std::size_t timestep;

    void addFunctor(SuperF3D<T>& f) {
      vtmWriter.addFunctor(f);
    }

    void write() {
      vtmWriter.writeVTI(timestep, locC, refC);
    }
  };

  void updateSeries(std::size_t iT, const std::string& fullNameSeries, const std::string& namePiece);

public:
  RefinedAMRwriter(RefinedMesh<T, D>& mesh, std::string name, int overlap = 0, bool writePerLevel = true);

  void createMasterFile();

  template <typename DESCRIPTOR, typename F>
  requires (DESCRIPTOR::d == D)
  void write(std::size_t iT, RefinedLattice<T,DESCRIPTOR>& refinedLattice, F&& vtkContext);
};

template <typename T, unsigned D>
RefinedAMRwriter(RefinedMesh<T, D>&, std::string, int = 0, bool = true) -> RefinedAMRwriter<T, D>;

} // namespace olb

#endif // REFINEMENT_AMR_WRITER_H
