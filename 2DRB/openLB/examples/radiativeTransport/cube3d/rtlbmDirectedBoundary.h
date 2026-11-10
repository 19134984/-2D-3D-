#pragma once
#include <olb.h>

using namespace olb;
using namespace olb::descriptors;

// PostProcessor to implement a directed boundary at the emittor
template<typename T, typename DESCRIPTOR, int discreteNormalX, int discreteNormalY, int discreteNormalZ>
struct RtlbmDirectedBoundaryPostProcessor3D
{
static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  T intensity;

  int getPriority() const {
    return 0;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell) any_platform{
    //T intensity = parameters.template get<INTENSITY>();
    for(int iPop = 0; iPop<DESCRIPTOR::q; iPop++){
      const auto c = descriptors::c<DESCRIPTOR>(iPop);
      T k = c[0]*discreteNormalX + c[1]*discreteNormalY + c[2]*discreteNormalZ;
      T norm_c = std::sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
      T norm_n = std::sqrt(discreteNormalX*discreteNormalX + discreteNormalY*discreteNormalY + discreteNormalZ*discreteNormalZ);
      T cos_theta = k / (norm_c * norm_n);
      if (util::nearZero(cos_theta + 1.)) {
        cell[iPop] = (1 - descriptors::t<T,DESCRIPTOR>(0))*intensity - descriptors::t<T,DESCRIPTOR>(iPop);
      } else {
        cell[iPop] = - descriptors::t<T,DESCRIPTOR>(iPop);
      }
    }
  }
};

// Set RtlbmDirectedBoundary for any indicated cells inside the block domain
template<typename T, typename DESCRIPTOR>
void setRtlbmDirectedBoundary(BlockLattice<T,DESCRIPTOR>& _block, BlockIndicatorF3D<T>& indicator, T intensity)
{
  using namespace boundaryhelper;
  const auto& blockGeometryStructure = indicator.getBlockGeometry();
  std::vector<int> discreteNormal(4,0);
  blockGeometryStructure.forSpatialLocations([&](auto iX, auto iY, auto iZ) {
    if (indicator(iX, iY, iZ)) {
      discreteNormal = blockGeometryStructure.getStatistics().getType(iX, iY, iZ);
      if (discreteNormal[1]!=0 || discreteNormal[2]!=0 || discreteNormal[3]!=0) {
        _block.addPostProcessor(
          typeid(stage::PostCollide), {iX,iY,iZ},
          promisePostProcessorForNormal<T, DESCRIPTOR, RtlbmDirectedBoundaryPostProcessor3D>(
            Vector<int,3>(discreteNormal.data()+1)
          )
        );
        _block.template defineDynamics<NoCollideDynamics>({iX, iY, iZ});
      } else {
        _block.template defineDynamics<BounceBack>({iX, iY, iZ});
      }
    }
  });
}

// Initialising the setRtlbmDirectedBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setRtlbmDirectedBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF3D<T>>&& indicator, T intensity)
{
  OstreamManager clout(std::cout, "setRtlbmDirectedBoundary");

  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); iCloc++) {
    setRtlbmDirectedBoundary<T,DESCRIPTOR>(sLattice.getBlock(iCloc), indicator->getBlockIndicatorF(iCloc), intensity);
  }
}