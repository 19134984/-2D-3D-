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

#ifndef REFINEMENT_AMR_WRITER_HH
#define REFINEMENT_AMR_WRITER_HH

#include "refinement/amrWriter.h"

namespace olb {

template <typename T, unsigned D> requires (D == 3)
RefinedAMRwriter<T, D>::RefinedAMRwriter(RefinedMesh<T, D>& mesh, std::string name, int overlap, bool writePerLevel)
  : _mesh(mesh), _name(name), _overlap(overlap), _writePerLevel(writePerLevel)
{}

template <typename T, unsigned D> requires (D == 3)
void RefinedAMRwriter<T, D>::updateSeries(std::size_t iT, const std::string& fullNameSeries, const std::string& namePiece)
{
  std::string content;
  {
    std::ifstream fin(fullNameSeries);
    if (fin) {
      content.assign((std::istreambuf_iterator<char>(fin)), std::istreambuf_iterator<char>());
    }
  }
  if (!content.empty()) {
    std::size_t pos = content.rfind("]");
    if (pos != std::string::npos) {
      std::string prefix = content.substr(0, pos);
      bool isEmpty = (prefix.find_last_not_of(" \n\r\t") == prefix.rfind("["));
      std::ofstream fout(fullNameSeries, std::ios::trunc);
      fout << prefix;
      if (!isEmpty) fout << ",\n";
      fout << "    { \"name\" : \"" << namePiece << "\", \"time\" : " << iT << " }\n";
      fout << "  ]\n}\n";
    }
  }
}

template <typename T, unsigned D> requires (D == 3)
void RefinedAMRwriter<T, D>::createMasterFile()
{
  if (singleton::mpi().getRank() == 0) {
    std::string pathSeries = singleton::directories().getVtkOutDir() + createFileName(_name) + ".vthb.series";
    std::ofstream fout(pathSeries, std::ios::trunc);
    if (fout) {
      fout << "{\n"
           << "  \"file-series-version\" : \"1.0\",\n"
           << "  \"files\" : [\n"
           << "  ]\n"
           << "}\n";
    }

    if (_writePerLevel) {
      for (std::size_t iLevel = 0; iLevel < _mesh.getLevelCount(); ++iLevel) {
        std::string levelName = _name + "_level_" + std::to_string(iLevel);
        std::string pathPVD = singleton::directories().getVtkOutDir() + createFileName(levelName) + ".pvd";
        std::ofstream foutPVD(pathPVD, std::ios::trunc);
        if (foutPVD) {
          foutPVD << "<?xml version=\"1.0\"?>\n"
                  << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
                  << "<Collection>\n"
                  << "</Collection>\n"
                  << "</VTKFile>\n";
        }
      }
    }
  }
}

template <typename T, unsigned D> requires (D == 3)
template <typename DESCRIPTOR, typename F>
requires (DESCRIPTOR::d == D)
void RefinedAMRwriter<T, D>::write(std::size_t iT, RefinedLattice<T,DESCRIPTOR>& refinedLattice, F&& vtkContext)
{
  if (singleton::mpi().getRank() == 0) {
    std::string vthbName = createFileName(_name, iT) + ".vthb";
    std::string pathVTHB = singleton::directories().getVtkOutDir() + "data/" + vthbName;
    std::ofstream fout(pathVTHB, std::ios::trunc);
    if (fout) {
      auto& decomp0 = _mesh.getDecompositionBase();
      const T delta0 = decomp0.getMotherCuboid().getDeltaR();
      const Vector<T,D> origin0 = decomp0.getMotherCuboid().getOrigin();

      Vector<T,D> vtkOrigin;
      for (unsigned d = 0; d < D; ++d) {
        vtkOrigin[d] = origin0[d] - 0.5 * delta0;
      }

      fout << "<?xml version=\"1.0\"?>\n";
      fout << "<VTKFile type=\"vtkOverlappingAMR\" version=\"1.1\" byte_order=\"LittleEndian\">\n";
      fout << "  <vtkOverlappingAMR origin=\""
           << vtkOrigin[0] << " " << vtkOrigin[1] << " " << vtkOrigin[2]
           << "\" grid_description=\"XYZ\">\n";

      for (std::size_t iLevel = 0; iLevel < _mesh.getLevelCount(); ++iLevel) {
        auto& decomp_l = _mesh.getLevel(iLevel).getCuboidDecomposition();
        T delta_l = decomp_l.getMotherCuboid().getDeltaR();

        int refRatio = (iLevel == 0) ? 1 : 2;
        fout << "    <Block level=\"" << iLevel << "\" refinement_ratio=\"" << refRatio << "\" spacing=\""
             << delta_l << " " << delta_l << " " << delta_l << "\">\n";

        int dsIndex = 0;
        for (int iC = 0; iC < decomp_l.size(); ++iC) {
          auto& cuboid = decomp_l.get(iC);
          if (!cuboid.isCoarseBorder()) {
            auto refC = cuboid.getReferenceCuboid();
            if (refC) {
              Vector<T,D> Pc = cuboid.getOrigin();
              Vector<int,D> Nc = cuboid.getExtent();

              int imin[3], imax[3];
              for (unsigned d = 0; d < D; ++d) {
                imin[d] = static_cast<int>(std::round((Pc[d] - origin0[d]) / delta_l + 0.5 * (delta0 / delta_l - 1.0)));
                imax[d] = imin[d] + Nc[d] - 1;
              }

              std::string vtiFile = createFileName(_name, iT, *refC) + ".vti";

              fout << "      <DataSet index=\"" << dsIndex++ << "\" amr_box=\""
                   << imin[0] << " " << imax[0] << " "
                   << imin[1] << " " << imax[1] << " "
                   << imin[2] << " " << imax[2] << "\" "
                   << "file=\"" << vtiFile << "\" />\n";
            }
          }
        }
        fout << "    </Block>\n";
      }
      fout << "  </vtkOverlappingAMR>\n";
      fout << "</VTKFile>\n";
    }

    std::string pathSeries = singleton::directories().getVtkOutDir() + createFileName(_name) + ".vthb.series";
    updateSeries(iT, pathSeries, "data/" + vthbName);

    if (_writePerLevel) {
      for (std::size_t iLevel = 0; iLevel < _mesh.getLevelCount(); ++iLevel) {
        std::string levelName = _name + "_level_" + std::to_string(iLevel);
        std::string levelVtmName = createFileName(levelName, iT) + ".vtm";
        std::string pathVTM = singleton::directories().getVtkOutDir() + "data/" + levelVtmName;

        std::ofstream foutVTM(pathVTM, std::ios::trunc);
        if (foutVTM) {
          foutVTM << "<?xml version=\"1.0\"?>\n"
                  << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\">\n"
                  << "  <vtkMultiBlockDataSet>\n";

          auto& decomp_l = _mesh.getLevel(iLevel).getCuboidDecomposition();
          int dsIndex = 0;
          for (int iC = 0; iC < decomp_l.size(); ++iC) {
            auto& cuboid = decomp_l.get(iC);
            if (!cuboid.isCoarseBorder()) {
              auto refC = cuboid.getReferenceCuboid();
              if (refC) {
                std::string vtiFile = createFileName(_name, iT, *refC) + ".vti";
                foutVTM << "    <Block index=\"" << dsIndex++ << "\">\n"
                        << "      <DataSet index=\"0\" file=\"" << vtiFile << "\" />\n"
                        << "    </Block>\n";
              }
            }
          }
          foutVTM << "  </vtkMultiBlockDataSet>\n"
                  << "</VTKFile>\n";
        }

        std::string pathPVD = singleton::directories().getVtkOutDir() + createFileName(levelName) + ".pvd";
        std::ofstream foutPVD(pathPVD, std::ios::in | std::ios::out | std::ios::ate);
        if (foutPVD) {
          // Overwrite closing tags: "</Collection>\n</VTKFile>\n" which is exactly 25 bytes
          foutPVD.seekp(-25, std::ios::end);
          foutPVD << "<DataSet timestep=\"" << iT << "\" group=\"\" part=\"\" file=\"data/" << levelVtmName << "\"/>\n"
                  << "</Collection>\n"
                  << "</VTKFile>\n";
        }
      }
    }
  }

  std::string name = _name;
  int overlap = _overlap;
  auto& mesh = _mesh;

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
          SuperVTMwriter3D<T> vtmWriter(name, overlap, true, true, VtkDataMode::CellData);
          WriterWrapper wrapper{vtmWriter, locC, *refC, iT};
          vtkContext(wrapper, iLevel, sLatticeLevel);
        }
      }
    });
  }
}

} // namespace olb

#endif // REFINEMENT_AMR_WRITER_HH
