#ifndef NUMERICAL_FLUX_HPP
#define NUMERICAL_FLUX_HPP

#include "common.hpp"
#include "physics_includes.hpp"

class NumericalFlux
{
    public:
    ConservedToPrimitive*    U_to_P;
    WaveSpeedsFromPrimitive* c_from_P;
    FluxesFromPrimitive*     F_from_P;

    inline virtual void operator()(const real_t* const UL, 
                                   const real_t* const UR,
                                         real_t* const PL,
                                         real_t* const PR,
                                         real_t* const cL,
                                         real_t* const cR,
                                         real_t* const FL,
                                         real_t* const FR,
                                         real_t (*Fnum)[3],
                                   const int        Nfield) const = 0;
};


class LaxFriedrichs : public NumericalFlux
{
    public:
    ACCEL_DECORATOR
    inline virtual void operator()(const real_t* const UL, 
                                   const real_t* const UR,
                                         real_t* const PL,
                                         real_t* const PR,
                                         real_t* const cL,
                                         real_t* const cR,
                                         real_t (*Fnum)[3],
                                   const int        Nfield) const
    {
        (*U_to_P)(UL, PL);
        (*U_to_P)(UR, PR);

        (*c_from_P)(PL, cL);
        (*c_from_P)(PR, cR);

        (*F_from_P)(
        return;
    }
};
#endif
