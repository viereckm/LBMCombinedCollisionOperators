/**
 * @file CustomCollisionModels.h
 * @author Markus Viereck
 * @brief Collision models containing parameters of two LBM models (e.g. SRT and TRT) to be used in simulation
 * @version 0.1
 * @date 2024-07-01
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#pragma once 

#include "lbm/lattice_model/CollisionModel.h"

namespace walberla {
namespace lbm {
namespace collision_model {

struct SRT_TRT_tag {};

class SRT_TRT {
    public:
        using tag = SRT_TRT_tag;
        static const bool constant = true;
        static const real_t threeSixteenth;

        SRT_TRT( const real_t _omega, const real_t _lambda_d, uint_t _level = uint_t( 0 ) ) :
            omega_( _omega ), lambda_d_( _lambda_d ), level_( _level )
        {}

        /**
         * @brief constructs TRT specific parameters from SRT relaxation parameter according to refinement level
         * 
         * @param _omega SRT relaxation parameter, corresponds to lambda_e in TRT
         * @param _magicNumber parameter to determine lambda_d, can be used to tune staility / accuracy of TRT
         * @param _level refinement level of the processed block
         * @return SRT_TRT SRT_TRT instance with simulation parameters set w.r.t. specified magic number
         */
        static SRT_TRT constructWithMagicNumber( const real_t _omega, const real_t _magicNumber = threeSixteenth, const uint_t _level = uint_t(0) );

        void pack( mpi::SendBuffer & buffer ) const { buffer << omega_ << lambda_d_ << magicNumber_ << viscosity_ << level_; }
        void unpack( mpi::RecvBuffer & buffer ) { buffer >> omega_ >> lambda_d_ >> magicNumber_ >> viscosity_ >> level_; }

        /**
         * @brief adjusts simulation parameters to refinement level
         * 
         * @param block currently processed part of domain, determines refinement level
         * @param sbs 
         */
        void configure( IBlock & block, StructuredBlockStorage & sbs );
        void reset( const real_t _omega, const real_t _lambda_d, const uint_t _omega_level = uint_t(0) );
        void resetWithMagicNumber( const real_t _omega, const real_t _magicNumber = threeSixteenth, const uint_t _omega_level = uint_t(0) );

        real_t omega() const { return omega_; }
        real_t lambda_e() const { return omega_; }
        real_t lambda_d() const { return lambda_d_; }
        real_t viscosity() const { return viscosity_; }

        inline real_t omega_bulk() const { return omega(); }
        inline real_t omega_odd() const { return lambda_d(); }
        
        static real_t lambda_e( const real_t _omega ) { return _omega; }

        /**
         * @brief constructs TRT relaxation parameter for asymmetric populations according to bulk relaxation parameter and magic number
         * 
         * @param _omega bulk relaxation parameter
         * @param _magicNumber tuned parameter
         * @return real_t asymmetric relaxation parameter
         */
        static real_t lambda_d( const real_t _omega, const real_t _magicNumber = threeSixteenth ) { return ( real_t(4) - real_t(2) * _omega ) / ( real_t(4) * _magicNumber * _omega + real_t(2) - _omega ); }

        static real_t magicNumber( const real_t _lambda_e, const real_t _lambda_d ) { return ( ( real_t(2) - _lambda_e ) * ( real_t(2) - _lambda_d ) ) / ( real_t(4) * _lambda_e * _lambda_d ); }
        uint_t level() const { return level_; }

        real_t viscosity( const uint_t _level ) const
        {
            const real_t _lambda_e = levelDependentRelaxationParameter( _level, omega_, level_ );
            return viscosityFromOmega( _lambda_e );
        }

        real_t omega( const cell_idx_t /* x */, const cell_idx_t /* y */, const cell_idx_t /* z */, const Vector3< real_t > & /* velocity */ = Vector3< real_t >(), const real_t /* rho */ = real_t(1) ) const { return omega_; }
        real_t viscosity( const cell_idx_t /* x */, const cell_idx_t /* y */, const cell_idx_t /* z */ ) const { return viscosity_; }

    private:
        SRT_TRT() : omega_( real_t(0) ), lambda_d_( real_t(0) ), magicNumber_( real_t(0) ), viscosity_( real_t(0) ), level_( uint_t(0) ) {}

        void initWithMagicNumber( const real_t _omega, const real_t _magicNumber, const uint_t _level );
        
        // SRT relaxation parameter:
        real_t omega_;

        // TRT parameters:
        // real_t lambda_e_: this is omega_;
        real_t lambda_d_;
        real_t magicNumber_;

        // simulation parameters:
        real_t viscosity_;
        uint_t level_;
};

} // namespace collision_model
} // namespace lbm
} // namespace walberla