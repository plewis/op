#pragma once

#include <boost/random/gamma_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include "xtreedist.hpp"

using namespace std;

namespace treedist {

    class Lot {
        public:
                                            Lot();
                                            ~Lot();
            
            void                            setSeed(unsigned seed);
            double                          uniform();
            int                             randint(int low, int high);
            pair<unsigned, unsigned>        nchoose2(unsigned n);
            double                          normal();
            double                          gamma(double shape, double scale);
            double                          logUniform();
            
            typedef std::shared_ptr<Lot>    SharedPtr;

        private:
        
            typedef boost::variate_generator<boost::mt19937 &, boost::random::uniform_01<> >                uniform_variate_generator_t;
            typedef boost::variate_generator<boost::mt19937 &, boost::random::normal_distribution<> >       normal_variate_generator_t;
            typedef boost::variate_generator<boost::mt19937 &, boost::random::gamma_distribution<> >        gamma_variate_generator_t;
            typedef boost::variate_generator<boost::mt19937 &, boost::random::uniform_int_distribution<> >  uniform_int_generator_t;

            unsigned                                   _seed;
            boost::mt19937                             _generator;
            std::shared_ptr<uniform_variate_generator_t>    _uniform_variate_generator;
            std::shared_ptr<normal_variate_generator_t>     _normal_variate_generator;
            std::shared_ptr<gamma_variate_generator_t>      _gamma_variate_generator;
            std::shared_ptr<uniform_int_generator_t>        _uniform_int_generator;

            double                                          _gamma_shape;
            int                                             _low;
            int                                             _high;
    };
    
    // member function bodies go here
    
    inline Lot::Lot() : _seed(0), _gamma_shape(1.0), _low(0), _high(100) {
        _generator.seed(static_cast<unsigned int>(time(0)));
        _uniform_variate_generator = std::shared_ptr<uniform_variate_generator_t>(new uniform_variate_generator_t(_generator, boost::random::uniform_01<>()));
        _normal_variate_generator = std::shared_ptr<normal_variate_generator_t>(new normal_variate_generator_t(_generator, boost::random::normal_distribution<>()));
        _gamma_variate_generator = std::shared_ptr<gamma_variate_generator_t>(new gamma_variate_generator_t(_generator, boost::random::gamma_distribution<>(_gamma_shape)));
        _uniform_int_generator = std::shared_ptr<uniform_int_generator_t>(new uniform_int_generator_t(_generator, boost::random::uniform_int_distribution<>(_low, _high)));
    }
        
    inline Lot::~Lot() {
        _uniform_variate_generator.reset();
        _normal_variate_generator.reset();
        _gamma_variate_generator.reset();
        _uniform_int_generator.reset();
    }
        
    inline void Lot::setSeed(unsigned seed) {
        _seed = seed;
        _generator.seed(_seed > 0 ? _seed : static_cast<unsigned int>(time(0)));
    }
        
    inline double Lot::uniform() {
        double u = (*_uniform_variate_generator)();
        while (u <= 0.0)
            u = (*_uniform_variate_generator)();
        return u;
    }

    inline double Lot::logUniform() {
        double u = (*_uniform_variate_generator)();
        while (u <= 0.0)
            u = (*_uniform_variate_generator)();
        return log(u);
    }
    
    inline double Lot::normal() {
        return (*_normal_variate_generator)();
    }

    inline double Lot::gamma(double shape, double scale) {
        assert(shape > 0.0);
        assert(scale > 0.0);
        if (shape != _gamma_shape) {
            _gamma_shape = shape;
            _gamma_variate_generator.reset(new gamma_variate_generator_t(_generator, boost::random::gamma_distribution<>(_gamma_shape,1)));
        }
        double deviate = (*_gamma_variate_generator)();
        return scale*deviate;
    }

    inline int Lot::randint(int low, int high) {
        if (low != _low || high != _high) {
            _low  = low;
            _high = high;
            _uniform_int_generator.reset(new uniform_int_generator_t(_generator, boost::random::uniform_int_distribution<>(_low,_high)));
        }
        return (*_uniform_int_generator)();
    }
    
    inline pair<unsigned, unsigned> Lot::nchoose2(unsigned n) {
        if (n < 2)
            throw XTreeDist(format("nchoose2 called with n = %d") % n);
        int i = 0;
        int j = 1;
        if (n > 2) {
            i = randint(1, n);
            j = i + randint(1, n-1) - 1;
            i = i - 1;
            j = j % n;
        }
        
        return make_pair(i, j);
    }
}
