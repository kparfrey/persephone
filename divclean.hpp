#include "common.hpp"
#include "tensor_field.hpp"
#include "face_communicator.hpp"

class Physics;
class NumericalFlux;

namespace divClean
{
    void soln_to_flux(const real_t* const __restrict__ matrix, 
                      const real_t* const __restrict__ U, 
                            real_t* const __restrict__ Uf, 
                      const LengthBucket lb, const int dir);

    void fluxDeriv_to_soln(const real_t* const __restrict__ matrix, 
                           const real_t* const __restrict__ F, 
                                 real_t* const __restrict__ dF, 
                           const LengthBucket lb, const int dir);
    
    void bulk_fluxes(const real_t* const __restrict__ Uf,
                           real_t* const __restrict__ F ,
                     const VectorField                S ,
                     const Physics* const __restrict__ physics,
                     const LengthBucket lb, const int dir);

    void flux_divergence(const VectorField                dF,
                         const real_t* const __restrict__ Jrdetg,
                               real_t* const __restrict__ divF,
                         const LengthBucket lb);
    
    void add_geometric_sources(      real_t* const __restrict__  divF, 
                               const real_t* const __restrict__  U,
                               const Physics* const __restrict__ physics,                               
                               const int Ns);

    void fill_face_data(const real_t* const __restrict__ Uf,
                              FaceCommunicator           face,
                        const LengthBucket               lb);

    void external_numerical_flux(const FaceCommunicator           face,
                                       real_t* const __restrict__ F,
                                 const NumericalFlux*             F_numerical,
                                 const VectorField                S,
                                 const LengthBucket               lb);

    void internal_numerical_flux(const real_t* const __restrict__ Uf,
                                       real_t* const __restrict__ F,
                                 const NumericalFlux*             F_numerical,
                                 const VectorField                S,
                                 const VectorField                normal,
                                 const LengthBucket               lb,
                                 const int                        dir);    
}
