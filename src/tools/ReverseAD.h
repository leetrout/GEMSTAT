#ifndef GEMSTAT_REVERSE_AD_H
#define GEMSTAT_REVERSE_AD_H

/*
 * Minimal reverse-mode automatic differentiation.
 *
 * A Var records every arithmetic operation it takes part in on a tape: one
 * node per operation, holding for each operand its index and the partial
 * derivative of the result with respect to it.  After a computation,
 * Tape::backward() propagates adjoints from the outputs to the inputs in one
 * sweep, so the gradient of a scalar function of n inputs costs about as
 * much as a few evaluations of the function, independent of n.
 *
 * Nodes may have any number of operands.  The dynamic-programming
 * recurrences of GEMSTAT are dominated by sums of products
 * (sum += w_k * Z_k over a site's neighbours); Accumulator records such a
 * sum as a single node instead of two nodes per term, which makes the tape,
 * the forward pass and the backward sweep several times cheaper.
 *
 * Only what the recurrences need is provided: + - * /, unary minus, exp,
 * log, sqrt, max, comparisons (on values) and isnan/isinf.  Comparisons and
 * the ternary operator branch on values, so a piecewise function is
 * differentiated on the branch that was taken.
 *
 * The tape a Var records to is the *active* tape of the current thread
 * (Tape::activate()).  Each OpenMP thread uses its own tape.
 */

#include <vector>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>

namespace gemstat_ad {

template< class R >
class Tape
{
    public:
        Tape() { node_end.reserve( 1 << 16 ); arg_index.reserve( 1 << 17 ); arg_partial.reserve( 1 << 17 ); }

        void clear() { node_end.clear(); arg_index.clear(); arg_partial.clear(); }
        std::size_t size() const { return node_end.size(); }

        // a new independent variable
        int leaf() { node_end.push_back( (int)arg_index.size() ); return (int)node_end.size() - 1; }
        int unary( int a, R da )
        {
            arg_index.push_back( a ); arg_partial.push_back( da );
            node_end.push_back( (int)arg_index.size() ); return (int)node_end.size() - 1;
        }
        int binary( int a, R da, int b, R db )
        {
            arg_index.push_back( a ); arg_partial.push_back( da );
            arg_index.push_back( b ); arg_partial.push_back( db );
            node_end.push_back( (int)arg_index.size() ); return (int)node_end.size() - 1;
        }
        // a node with the given operands (index, partial); n may be 0
        int nary( const std::pair< int, R >* args, std::size_t n )
        {
            for ( std::size_t k = 0; k < n; k++ ) { arg_index.push_back( args[k].first ); arg_partial.push_back( args[k].second ); }
            node_end.push_back( (int)arg_index.size() ); return (int)node_end.size() - 1;
        }

        /*
         * Reverse sweep.  adjoint must have size() entries and already hold the
         * seeds (the derivative of the quantity of interest with respect to each
         * output node, zero elsewhere).  On return adjoint[i] is the derivative
         * with respect to node i, for every node, including the leaves.
         */
        void backward( std::vector< R >& adjoint ) const
        {
            if ( adjoint.size() != node_end.size() ) throw std::invalid_argument( "Tape::backward: adjoint vector has the wrong size" );
            const int* idx = arg_index.empty() ? 0 : &arg_index[0];
            const R* d = arg_partial.empty() ? 0 : &arg_partial[0];
            for ( int i = (int)node_end.size() - 1; i >= 0; i-- )
            {
                const R a_i = adjoint[i];
                if ( a_i == R(0) ) continue;
                const int begin = i > 0 ? node_end[i - 1] : 0;
                const int end = node_end[i];
                for ( int k = begin; k < end; k++ ) adjoint[ idx[k] ] += a_i * d[k];
            }
        }

        static Tape*& active() { static thread_local Tape* t = 0; return t; }
        void activate() { active() = this; }

    private:
        std::vector< int > node_end;      // node i's operands are arg_*[ node_end[i-1] .. node_end[i] )
        std::vector< int > arg_index;
        std::vector< R > arg_partial;
};

template< class R >
class Var
{
    public:
        typedef R value_type;

        Var() : v( R(0) ), i( -1 ) {}
        explicit Var( R value ) : v( value ), i( -1 ) {}          // a constant: not on the tape
        Var( R value, int index ) : v( value ), i( index ) {}

        // an independent variable, recorded on the active tape
        static Var input( R value ) { return Var( value, tape().leaf() ); }

        Var& operator=( R value ) { v = value; i = -1; return *this; }

        R value() const { return v; }
        int index() const { return i; }
        bool is_constant() const { return i < 0; }

        Var& operator+=( const Var& o ) { *this = *this + o; return *this; }
        Var& operator-=( const Var& o ) { *this = *this - o; return *this; }
        Var& operator*=( const Var& o ) { *this = *this * o; return *this; }
        Var& operator/=( const Var& o ) { *this = *this / o; return *this; }
        Var& operator+=( double o ) { *this = *this + o; return *this; }
        Var& operator-=( double o ) { *this = *this - o; return *this; }
        Var& operator*=( double o ) { *this = *this * o; return *this; }
        Var& operator/=( double o ) { *this = *this / o; return *this; }

        Var operator-() const { return is_constant() ? Var( -v ) : Var( -v, tape().unary( i, R(-1) ) ); }

        // ---- binary arithmetic ----
        friend Var operator+( const Var& x, const Var& y ) { return make2( x.v + y.v, x, R(1), y, R(1) ); }
        friend Var operator-( const Var& x, const Var& y ) { return make2( x.v - y.v, x, R(1), y, R(-1) ); }
        friend Var operator*( const Var& x, const Var& y ) { return make2( x.v * y.v, x, y.v, y, x.v ); }
        friend Var operator/( const Var& x, const Var& y ) { R q = x.v / y.v; return make2( q, x, R(1) / y.v, y, -q / y.v ); }

        friend Var operator+( const Var& x, double y ) { return make1( x.v + R(y), x, R(1) ); }
        friend Var operator-( const Var& x, double y ) { return make1( x.v - R(y), x, R(1) ); }
        friend Var operator*( const Var& x, double y ) { return make1( x.v * R(y), x, R(y) ); }
        friend Var operator/( const Var& x, double y ) { return make1( x.v / R(y), x, R(1) / R(y) ); }
        friend Var operator+( double x, const Var& y ) { return make1( R(x) + y.v, y, R(1) ); }
        friend Var operator-( double x, const Var& y ) { return make1( R(x) - y.v, y, R(-1) ); }
        friend Var operator*( double x, const Var& y ) { return make1( R(x) * y.v, y, R(x) ); }
        friend Var operator/( double x, const Var& y ) { R q = R(x) / y.v; return make1( q, y, -q / y.v ); }

        // ---- functions ----
        friend Var exp( const Var& x ) { R e = std::exp( x.v ); return make1( e, x, e ); }
        friend Var log( const Var& x ) { return make1( std::log( x.v ), x, R(1) / x.v ); }
        friend Var sqrt( const Var& x ) { R s = std::sqrt( x.v ); return make1( s, x, R(0.5) / s ); }
        friend Var max( const Var& x, const Var& y ) { return x.v >= y.v ? x : y; }
        friend Var max( double x, const Var& y ) { return R(x) >= y.v ? Var( R(x) ) : y; }
        friend Var max( const Var& x, double y ) { return x.v >= R(y) ? x : Var( R(y) ); }
        friend bool isnan( const Var& x ) { return std::isnan( x.v ); }
        friend bool isinf( const Var& x ) { return std::isinf( x.v ); }

        // ---- comparisons (on values) ----
        friend bool operator<( const Var& x, const Var& y ) { return x.v < y.v; }
        friend bool operator>( const Var& x, const Var& y ) { return x.v > y.v; }
        friend bool operator<=( const Var& x, const Var& y ) { return x.v <= y.v; }
        friend bool operator>=( const Var& x, const Var& y ) { return x.v >= y.v; }
        friend bool operator==( const Var& x, const Var& y ) { return x.v == y.v; }
        friend bool operator!=( const Var& x, const Var& y ) { return x.v != y.v; }
        friend bool operator<( const Var& x, double y ) { return x.v < R(y); }
        friend bool operator>( const Var& x, double y ) { return x.v > R(y); }
        friend bool operator<=( const Var& x, double y ) { return x.v <= R(y); }
        friend bool operator>=( const Var& x, double y ) { return x.v >= R(y); }
        friend bool operator==( const Var& x, double y ) { return x.v == R(y); }
        friend bool operator!=( const Var& x, double y ) { return x.v != R(y); }
        friend bool operator<( double x, const Var& y ) { return R(x) < y.v; }
        friend bool operator>( double x, const Var& y ) { return R(x) > y.v; }
        friend bool operator<=( double x, const Var& y ) { return R(x) <= y.v; }
        friend bool operator>=( double x, const Var& y ) { return R(x) >= y.v; }

        static Tape< R >& tape()
        {
            Tape< R >* t = Tape< R >::active();
            if ( !t ) throw std::logic_error( "gemstat_ad::Var used with no active Tape on this thread" );
            return *t;
        }

    private:
        R v;
        int i;

        static Var make1( R value, const Var& x, R dx )
        {
            if ( x.is_constant() ) return Var( value );
            return Var( value, tape().unary( x.i, dx ) );
        }
        static Var make2( R value, const Var& x, R dx, const Var& y, R dy )
        {
            if ( x.is_constant() ) return make1( value, y, dy );
            if ( y.is_constant() ) return make1( value, x, dx );
            return Var( value, tape().binary( x.i, dx, y.i, dy ) );
        }
};

/*
 * A running sum  s = init + sum_k a_k * b_k (+ sum_m c_m),  recorded as one
 * tape node.  For a plain floating point type it is exactly the sequential
 * loop  s += a_k * b_k,  so the values are those of the original code.
 */
template< class T >
class Accumulator
{
    public:
        explicit Accumulator( const T& init ) : s( init ) {}
        void add( const T& a, const T& b ) { s += a * b; }
        void add( const T& a ) { s += a; }
        T result() const { return s; }
    private:
        T s;
};

template< class R >
class Accumulator< Var< R > >
{
    public:
        typedef Var< R > V;
        explicit Accumulator( const V& init ) : s( init.value() )
        {
            args.reserve( 32 );
            if ( !init.is_constant() ) args.push_back( std::make_pair( init.index(), R(1) ) );
        }
        void add( const V& a, const V& b )
        {
            s += a.value() * b.value();
            if ( !a.is_constant() ) args.push_back( std::make_pair( a.index(), b.value() ) );
            if ( !b.is_constant() ) args.push_back( std::make_pair( b.index(), a.value() ) );
        }
        void add( const V& a )
        {
            s += a.value();
            if ( !a.is_constant() ) args.push_back( std::make_pair( a.index(), R(1) ) );
        }
        V result() const
        {
            if ( args.empty() ) return V( s );
            return V( s, V::tape().nary( &args[0], args.size() ) );
        }
    private:
        R s;
        std::vector< std::pair< int, R > > args;   // the operands recorded so far (several accumulators may be live at once)
};

// value() that also works for plain floating point types, so templated code
// can print or check a value whatever the scalar type
template< class R > inline R value_of( const Var< R >& x ) { return x.value(); }
inline double value_of( double x ) { return x; }
inline long double value_of( long double x ) { return x; }

} // namespace gemstat_ad

#endif
