/* ----------------------------------------------------------------------

    Contributing author and copyright for this file:

    Yu Chen; School of Civil Engineering, The University of Sydney; July 2025

------------------------------------------------------------------------- */

#ifdef COHESION_MODEL
COHESION_MODEL(COHESION_LUB,lub,3)
#else

#ifndef COHESION_MODEL_LUB_H_
#define COHESION_MODEL_LUB_H_

#include "contact_models.h"
#include "cohesion_model_base.h"
#include <cmath>
#include "neighbor.h"
#include "global_properties.h"
#include "fix_property_atom.h"

namespace MODEL_PARAMS
{
  inline static ScalarProperty* createFluidViscosityLub(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {
    ScalarProperty* fluidViscosityScalar = MODEL_PARAMS::createScalarProperty(registry, "fluidViscosity", caller); 
    return fluidViscosityScalar; 
  }
  inline static ScalarProperty* createRoughnessLub(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {
    ScalarProperty* RoughnessScalar = MODEL_PARAMS::createScalarProperty(registry, "Roughness", caller);
    return RoughnessScalar;
  }
  inline static ScalarProperty* createMinSeparationDistanceRatioLub(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {
    ScalarProperty* minSeparationDistanceRatioScalar = MODEL_PARAMS::createScalarProperty(registry, "minSeparationDistanceRatio", caller); 
    return minSeparationDistanceRatioScalar;
  }
}

namespace LIGGGHTS {
namespace ContactModels {
  using namespace std;
  using namespace LAMMPS_NS;

  template<>
  class CohesionModel<COHESION_LUB> : public CohesionModelBase {
  public:
    CohesionModel(LAMMPS * lmp, IContactHistorySetup * hsetup, class ContactModelBase * c) :
      CohesionModelBase(lmp, hsetup, c),
      fluidViscosity(0.0),
      Roughness(0.0),
      minSeparationDistanceRatio(0.0),
      tangentialReduce_(false)  // Initialize the member variable
    {

    }

    void registerSettings(Settings & settings)
    {
      settings.registerOnOff("tangential_reduce",tangentialReduce_,false);
    }

    inline void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb) {}

    void connectToProperties(PropertyRegistry & registry)
    {
      registry.registerProperty("fluidViscosity", &MODEL_PARAMS::createFluidViscosityLub);
      registry.connect("fluidViscosity", fluidViscosity,"cohesion_model lub");
      registry.registerProperty("Roughness", &MODEL_PARAMS::createRoughnessLub);
      registry.connect("Roughness", Roughness,"cohesion_model lub");
      registry.registerProperty("minSeparationDistanceRatio", &MODEL_PARAMS::createMinSeparationDistanceRatioLub);
      registry.connect("minSeparationDistanceRatio", minSeparationDistanceRatio,"cohesion_model lub"); 
      
      if(force->cg_active())
        error->cg(FLERR,"cohesion model lub");
      neighbor->register_contact_dist_factor(minSeparationDistanceRatio);
    }

    void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {
      const double ri = sidata.radi;
      const double radj = sidata.is_wall ? ri : sidata.radj;  // Fixed: using sidata instead of scdata
      const double rEff = ri*radj/(ri+radj);  // Fixed: using ri instead of radi
      const double mult = -6*M_PI*fluidViscosity*rEff*rEff/Roughness;
      const double Fn_lub = mult*sidata.vn;

      if(tangentialReduce_) sidata.Fn += Fn_lub;
      if(sidata.contact_flags) *sidata.contact_flags |= CONTACT_COHESION_MODEL;

      // Apply normal force
      const double fx = Fn_lub * sidata.en[0];
      const double fy = Fn_lub * sidata.en[1];
      const double fz = Fn_lub * sidata.en[2];

      i_forces.delta_F[0] += fx;
      i_forces.delta_F[1] += fy;
      i_forces.delta_F[2] += fz;

      j_forces.delta_F[0] -= fx;
      j_forces.delta_F[1] -= fy;
      j_forces.delta_F[2] -= fz;
    }

    inline void endSurfacesIntersect(SurfacesIntersectData & sidata, ForceData&, ForceData&)
    {}
    
    void beginPass(SurfacesIntersectData&, ForceData&, ForceData&)
    {}
    
    void endPass(SurfacesIntersectData&, ForceData&, ForceData&)
    {}
    
    void surfacesClose(SurfacesCloseData & scdata, ForceData & i_forces, ForceData & j_forces)
    {
      if(scdata.contact_flags) *scdata.contact_flags |= CONTACT_COHESION_MODEL;
      
      const double radi = scdata.radi; 
      const double radj = scdata.is_wall ? radi : scdata.radj;
      const double r = sqrt(scdata.rsq);
      const double dist = scdata.is_wall ? r - radi : r - (radi + radj);
      const double rEff = radi*radj / (radi+radj);

      const double rinv = 1.0 / r;
      const double dx = scdata.delta[0]; 
      const double dy = scdata.delta[1]; 
      const double dz = scdata.delta[2]; 
      const double enx = dx * rinv;
      const double eny = dy * rinv;
      const double enz = dz * rinv;

      const double vr1 = scdata.v_i[0] - scdata.v_j[0];
      const double vr2 = scdata.v_i[1] - scdata.v_j[1];
      const double vr3 = scdata.v_i[2] - scdata.v_j[2];
      const double vn = vr1 * enx + vr2 * eny + vr3 * enz;
      
      if (dist < Roughness)
      {
        if(scdata.contact_flags) *scdata.contact_flags &= ~CONTACT_COHESION_MODEL; 
        return;
      }
      
      const double mult = -6.*M_PI*fluidViscosity*rEff*rEff;
      const double Fn_lub = mult*vn/std::max(Roughness,dist); 
      const double fx = Fn_lub * enx;
      const double fy = Fn_lub * eny;
      const double fz = Fn_lub * enz;
      scdata.has_force_update = true; 
      
      i_forces.delta_F[0] += fx;
      i_forces.delta_F[1] += fy;
      i_forces.delta_F[2] += fz;

      j_forces.delta_F[0] -= fx;
      j_forces.delta_F[1] -= fy;
      j_forces.delta_F[2] -= fz;
    }

  private:
    double fluidViscosity;
    double Roughness;
    double minSeparationDistanceRatio;
    bool tangentialReduce_;  // Added missing member variable
  };

}  // namespace ContactModels
}  // namespace LIGGGHTS

#endif // COHESION_MODEL_LUB_H_
#endif // COHESION_MODEL