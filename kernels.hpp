#ifndef KERNELS_HPP
#define KERNELS_HPP

#include "common.hpp"
#include "tensor_field.hpp"
#include "face_communicator.hpp"

#include "physics_includes.hpp"
#include "numerical_flux.hpp"
#include "boundary_conditions.hpp"

#include <iostream>



namespace kernels
{
    /* Memory allocation on device --- don't use for host-side only memory */
    inline
    real_t* alloc(int n);     // Zero-initialize

    inline
    real_t* alloc_raw(int n); // Don't initialize



    inline
    void free(real_t* a);

    inline
    void add_2_vectors(real_t* v1,     real_t* __restrict__ v2, 
                       real_t  a1,     real_t               a2, 
                       real_t* result, const int N);

    inline
    void add_3_vectors(const real_t* const v1, const real_t* const v2, 
                       const real_t* const __restrict__ v3,
                       const real_t a1, const real_t a2, const real_t a3, 
                       real_t* const result, const int N);

    inline
    void add_vectors_inPlace(      real_t* const __restrict__ v,
                             const real_t* const __restrict__ v_add,
                             const int N);

    inline
    void add_scaled_vectors_inPlace(      real_t* const __restrict__ v,
                                    const real_t* const __restrict__ v_add,
                                    const real_t                     scalar,
                                    const int N);
 
#if 0
    void product_2_vectors(const real_t* const __restrict__ v1,
                           const real_t* const __restrict__ v2,
                                 real_t* const __restrict__ result,
                           const int N);
#endif

    inline
    void multiply_by_scalar(const real_t* const __restrict__ v, 
                            const real_t                     scalar,
                                  real_t* const __restrict__ result,
                            const int                        N);

    inline
    void multiply_by_scalar_inPlace(      real_t* const __restrict__ v, 
                                    const real_t                     scalar,
                                    const int                        N);

    inline
    void soln_to_flux(const real_t* const __restrict__ matrix, 
                      const real_t* const __restrict__ U, 
                            real_t* const __restrict__ Uf, 
                      const LengthBucket lb, const int dir);

    inline
    void fluxDeriv_to_soln(const real_t* const __restrict__ matrix, 
                           const real_t* const __restrict__ F, 
                                 real_t* const __restrict__ dF, 
                           const LengthBucket lb, const int dir);
    
    inline
    void bulk_fluxes(const real_t* const __restrict__ Uf,
                           real_t* const __restrict__ F ,
                     const VectorField                S ,
                     const Physics<PhysicsType>       physics,
                     //const ConservedToPrimitive*  U_to_P,
                     //const FluxesFromPrimitive* F_from_P,
                     const LengthBucket lb, const int dir);

    inline
    void flux_divergence(const VectorField                dF,
                         const real_t* const __restrict__ Jrdetg,
                               real_t* const __restrict__ divF,
                         const LengthBucket lb);
    
    inline
    void filter_field(      real_t* const __restrict__ U,
                      const VectorField filter_matrices,
                      const LengthBucket lb);

    inline
    void add_geometric_sources(      real_t* const __restrict__  divF, 
                               const real_t* const __restrict__  U,
                               const VectorField                 dP,
                               const Physics<PhysicsType>        physics,
                               const int Nfield, const int Ns);

#if 0
    void scalar_field_source(      real_t* const __restrict__ divF,
                             const real_t* const __restrict__ U, 
                             const LengthBucket lb,
                             const real_t c_h, const real_t damping_rate);

    void store_divB(const real_t* const __restrict__ divF,
                          real_t* const __restrict__ divB, 
                    const real_t c_h,
                    const LengthBucket lb);
#endif
    
    inline
    void fill_face_data(const real_t* const __restrict__ Uf,
                              FaceCommunicator           face,
                        const LengthBucket               lb);

    inline
    void dirichlet_boundary_conditions(const FaceCommunicator& face);
    
    inline
    void neumann_boundary_conditions(const FaceCommunicator& face);
    
    inline
    void external_numerical_flux(const FaceCommunicator           face,
                                       real_t* const __restrict__ F,
                                 const NumericalFlux<NumFluxType> F_numerical,
                                 const VectorField                S,
                                 const LengthBucket               lb);

    inline
    void internal_numerical_flux(const real_t* const __restrict__ Uf,
                                       real_t* const __restrict__ F,
                                 const NumericalFlux<NumFluxType> F_numerical,
                                 const VectorField                S,
                                 const VectorField                normal,
                                 const LengthBucket               lb,
                                 const int                        dir);    

    inline
    void external_interface_average(const FaceCommunicator           face,
                                          real_t* const __restrict__ Pf,
                                    const LengthBucket               lb);

    inline
    void internal_interface_average(      real_t* const __restrict__ Pf,
                                    const LengthBucket               lb,
                                    const int                        dir);    
     
    inline
    void gradient_ref_to_phys(const VectorField dP_ref,
                                    VectorField dP_phys,
                              const VectorField dxdr[3],
                              const LengthBucket lb);

    inline
    void phys_vector_to_ref_density(const real_t* const __restrict__ V_phys,
                                          real_t* const __restrict__ V_ref,
                                    const VectorField                S,
                                    const LengthBucket               lb,
                                    const int                        dir);

    inline
    void diffusive_flux(const real_t* const __restrict__ Pf,
                        const VectorField                dPf,
                              real_t* const __restrict__ F,
                        const Physics<PhysicsType>       physics,
                        const VectorField                S,
                        const LengthBucket               lb,
                        const int                        dir);

    inline
    void floors(      real_t* const __restrict__  U,
                const Physics<PhysicsType> physics,
                const LengthBucket               lb);

    inline
    void conserved_to_primitive_fluxpoints(      real_t* const __restrict__  UPf,
                                           const Physics<PhysicsType>        physics,
                                           const LengthBucket                lb,
                                           const int                         dir);

    inline
    void conserved_to_primitive_faces(      FaceCommunicator     face,
                                      const Physics<PhysicsType> physics,
                                      const LengthBucket         lb);

    inline
    void wall_BC_derivatives(const FaceCommunicator face,
                                   VectorField      dPf,
                             const Physics<PhysicsType> physics,                                   
                             const LengthBucket     lb);

    /* Not a final "kernel", obviously */
    inline
    real_t local_timestep(const real_t* const __restrict__ Uf,
                                real_t& vmax,
                          const VectorField timestep_transform,
                          const Physics<PhysicsType> physics,
                          //const ConservedToPrimitive*  U_to_P,
                          //const WaveSpeedsFromPrimitive* c_from_P,
                          const LengthBucket lb, const int dir);
}


#endif
