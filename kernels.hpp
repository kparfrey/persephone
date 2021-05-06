#include "common.hpp"
#include "tensor_field.hpp"
#include "face_communicator.hpp"

class ConservedToPrimitive;
class FluxesFromPrimitive;
class NumericalFlux;

class WaveSpeedsFromPrimitive;


namespace kernels
{
    /* Memory allocation on device --- don't use for host-side only memory */
    real_t* alloc(int n);     // Zero-initialize

    real_t* alloc_raw(int n); // Don't initialize



    void free(real_t* a);

    void add_2_vectors(real_t* v1,     real_t* __restrict__ v2, 
                       real_t  a1,     real_t               a2, 
                       real_t* result, const int N);

    void multiply_by_scalar(      real_t* __restrict__ v, 
                            const real_t               scalar,
                            const int                  N);

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
                     const ConservedToPrimitive*  U_to_P,
                     const FluxesFromPrimitive* F_from_P,
                     const LengthBucket lb, const int dir);

    void flux_divergence(const VectorField                dF,
                         const real_t* const __restrict__ Jrdetg,
                               real_t* const __restrict__ divF,
                         const LengthBucket lb);

    void scalar_field_source(      real_t* const __restrict__ divF,
                             const real_t* const __restrict__ U, 
                             const LengthBucket lb,
                             const real_t c_h, const real_t damping_rate);

    void store_divB(const real_t* const __restrict__ divF,
                          real_t* const __restrict__ divB, 
                    const LengthBucket lb);
    
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

    void external_interface_average(const FaceCommunicator           face,
                                          real_t* const __restrict__ Uf,
                                    const LengthBucket               lb);

    void internal_interface_average(      real_t* const __restrict__ Uf,
                                    const LengthBucket               lb,
                                    const int                        dir);    
     
    void viscous_flux(const VectorField                dU,
                            real_t* const __restrict__ Fvisc, 
                      const real_t viscosity,
                      const LengthBucket lb, const int dir);

    void fill_velocity_vector(      VectorField Vf,
                              const VectorField Uf,
                              const LengthBucket lb);    
    
    /* Not a final kernel, obviously */
    real_t local_timestep(const real_t* const __restrict__ Uf,
                      const VectorField timestep_transform,
                      const ConservedToPrimitive*  U_to_P,
                      const WaveSpeedsFromPrimitive* c_from_P,
                      const LengthBucket lb, const int dir);
}
