

#include <cassert>
#include <cmath>

#include <array>
#include <random>
#include <algorithm>
#include <iostream>

#include "CLI11.hpp"
#include "IrreducibleGF3.hpp"


#define MATRIX_SZ	20

extern int sobol_dj[N_SOBOL_INIT_TAB_ENTRIES];
extern int sobol_sj[N_SOBOL_INIT_TAB_ENTRIES];
extern int sobol_aj[N_SOBOL_INIT_TAB_ENTRIES];
extern int sobol_mk[N_SOBOL_INIT_TAB_ENTRIES][32];


template<typename T>
void showmx( const vector<vector<T>>& mx, const unsigned size ) 
{
    unsigned m= std::min<unsigned>(size, mx.size());
    for(unsigned r=  0; r < m; r++)
    {
        for(unsigned c= 0; c < m; c++) 
            std::cout << mx[r][c] << " ";
        std::cout << endl;
    }
    std::cout << endl;
}


void fill_mx(int sobol_mk_index, vector<vector<int>>& mx)
{
    for (int i = 0; i < MATRIX_SZ; i++) 
    {
        int val = sobol_mk[sobol_mk_index][i];
        int len = i+1;
        
        vector<int> digits = IntegerDigits(val,3,len);
        for (int j = 0; j < len; j++)
            mx[len-j-1][i] = digits[j];
    }
}


//
typedef std::vector< vector<int> > sobol3_matrix;

// integer digits / base 3
struct integer3
{
    std::array<int8_t, 10> digits;
    static constexpr unsigned pow3_tab[]= { 1, 3, 9, 27, 81, 243, 729, 2187, 6561, 19683, 59049 };
    
    integer3( ) : digits() {}
    
    integer3( unsigned x ) : digits()
    {
        for(unsigned i= 0; i < digits.size(); i++)
            digits[i]= (x / pow3_tab[i]) % 3;
    }
    
    unsigned value( const unsigned m= 10 ) const
    {
        unsigned x= 0;
        for(unsigned i= 0; i < m; i++)
            x+= pow3_tab[i] * digits[i];
        
        return x;
    }
    
    double value_double( const unsigned m ) const
    {
        return double(value(m)) / double(pow3_tab[m]);
    }
    
    operator unsigned( ) const
    {
        return value();
    }
    
    // x % 3
    static int8_t mod( const int x )
    {
        static constexpr int8_t tab_mod3[]= { 0, 1, 2, 0, 1, 2 };
        //~ assert(x >= 0);
        //~ assert(x < 6 );
        return tab_mod3[x];
    }

    // ( a + (b*c)%3 )%3
    static int8_t fma( const int a, const int b, const int c ) 
    {
        //~ assert(a >= 0);
        //~ assert(a < 3);
        //~ assert(b >= 0);
        //~ assert(b < 3);
        //~ assert(c >= 0);
        //~ assert(c < 3);
        static constexpr int8_t tab_fma4[]= { 0, 0, 0, 0, 0, 1, 2, 0, 0, 2, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 2, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 2, 2, 2, 0, 2, 0, 1, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        return tab_fma4[a*4*4 + b*4 + c];
    }
};
    
void print( const integer3& i3, const int m= 1 )
{
    int k= int(i3.digits.size()) -1;
    for(; k >= m; k--)
        if(i3.digits[k])
            break;
    
    for(; k >=0; k--)
        printf("%d", int(i3.digits[k]));
}

std::ostream& operator<< ( std::ostream& out, const integer3& i3 )
{
    int k= int(i3.digits.size()) -1;
    for(; k >= 1; k--)
        if(i3.digits[k])
            break;
    
    for(; k >=0; k--)
        out << int(i3.digits[k]);
    
    return out;
}



// build the Gray code of an integer
// cf https://en.wikipedia.org/wiki/Gray_code#n-ary_Gray_code
integer3 graycode( const integer3 n )
{
    integer3 gray;
    
    int shift= 0;
    for(int i= n.digits.size() -1; i >= 0; i--)
    {
        // The Gray digit gets shifted down by the sum of the higher digits.
        gray.digits[i]= integer3::mod(n.digits[i] + shift);
        shift= integer3::mod(shift + 3 - gray.digits[i]);
    }
    
    return gray;
}


// generate a sobol point in Gray code order.
integer3 point3_graycode( const sobol3_matrix& matrix, const unsigned m, const integer3& ig3, integer3& pg3, integer3& x3 )
{
    // find modified digit : previous graycode(i-1) and graycode(i)
    unsigned g= 0;
    for(; g < m; g++)
        if(ig3.digits[g] != pg3.digits[g])
            break;
    
    int d= int(ig3.digits[g]) - int(pg3.digits[g]);
    d= integer3::mod(d + 3);
    assert(d == 1);
    
    // update previous graycode index
    pg3= ig3;
    
    // update previous point
    for(unsigned j= 0; j < m; j++)
        x3.digits[j]= integer3::mod(x3.digits[j] + matrix[m-1 - j][g]);    // reverse order !!
    
    return x3;
}


// generate a sobol point by incrementally modifying the previous point.
integer3 point3_digits( const sobol3_matrix& matrix, const unsigned m, const integer3& i3, integer3& p3, integer3& x3 )
{
    for(unsigned k= 0; k < m; k++)
    {
        // find modified digits : previous index (i-1) and index (i)
        if(p3.digits[k] != i3.digits[k])
        {
            int d= int(i3.digits[k]) - int(p3.digits[k]);
            d= integer3::mod(d + 3);
            
            // update previous point
            for(unsigned j= 0; j < m; j++)
                x3.digits[j]= integer3::fma(x3.digits[j], d, matrix[m-1 - j][k]);
        }
    }
    
    // update previous index
    p3= i3;
    return x3;
}


// counter based rng
// all small crush tests ok
struct FCRNG
{
    unsigned n;
    unsigned key;
    
    FCRNG( ) : n(), key() { seed(123451); }
    FCRNG( const unsigned s ) : n(), key() { seed(s); }
    void seed( const unsigned s ) { n= 0; key= (s << 1) | 1u; }
    
    FCRNG& index( const unsigned i ) { n= i; return *this;}
    unsigned sample( ) { return hash(++n * key); }
    
    unsigned sample_range( const unsigned range )
    {
        // Efficiently Generating a Number in a Range
        // cf http://www.pcg-random.org/posts/bounded-rands.html
        unsigned divisor= ((-range) / range) + 1; // (2^32) / range
        if(divisor == 0) return 0;
        
        while(true)
        {
            unsigned x= sample() / divisor;
            if(x < range) return x;
        }
    }
    
    // c++ interface
    unsigned operator() ( ) { return sample(); }
    static constexpr unsigned min( ) { return 0; }
    static constexpr unsigned max( ) { return ~unsigned(0); }
    typedef unsigned result_type;
    
    inline unsigned hash( unsigned x )
    {
        x ^= x >> 16;
        x *= 0x21f0aaad;
        x ^= x >> 15;
        x *= 0xd35a2d97;
        x ^= x >> 15;
        return x;
    }
    // cf "hash prospector" https://github.com/skeeto/hash-prospector/blob/master/README.md
};


// nested uniform scrambling / owen srambling
integer3 scramble_base3( const integer3& a3, const unsigned seed, const unsigned ndigits )
{
    // all random permutations of base-3 digits
    static constexpr int8_t scramble[6][3]= 
    {
        {0, 1, 2},
        {0, 2, 1},
        {1, 0, 2},
        {1, 2, 0},
        {2, 0, 1},
        {2, 1, 0},
    };
    
    // counter-based random number generator
    FCRNG rng(seed);
    
    integer3 b3;
    unsigned node_index= 0;
    for(unsigned i= 0; i < ndigits; i++)
    {
        unsigned flip= rng.index(node_index).sample_range(6);
        unsigned digit= a3.digits[ndigits-1 - i];
        b3.digits[ndigits-1 - i]= scramble[flip][digit];
        
        node_index= 3*node_index +1 + digit;
        // heap layout, root i= 0, childs 3i+1, 3i+2, 3i+3
    }

    return b3;
}


int main( int argc, char **argv )
{
    int D= 8;
    int m= 5; // 3^5= 243 points
    int N= std::pow(3, m);
    
    string init_file = "initIrreducibleGF3.dat";

    CLI::App app("samplerGF3");
    app.add_option("-d", D, "number of dimensions, default: "+ std::to_string(D)  );
    app.add_option("-m", m, "size of the matrices, to produce up to 3^m points, default: " + std::to_string(m));
    app.add_option("-i,--init_file", init_file, "filename which contains initialization data, default: " + init_file);
    CLI11_PARSE(app, argc, argv)
    
    // read irreducible polynoms
    if(load_mk(init_file, false) < D)
    {
        std::cout << "[error] reading '" << init_file << "'..." << std::endl;
        return 1;
    }
    
    // sanity check
    assert(m > 0);
    assert(m < 10);
    assert(D < 48);
    
    // build matrices
    std::vector<sobol3_matrix> Cd;
    for(int d= 1; d <= D; d++)
    {
        generate_mkGF3(sobol_aj[d], sobol_sj[d], sobol_mk[d], 3);
        
        sobol3_matrix C(MATRIX_SZ, vector<int>(MATRIX_SZ));
        fill_mx(d, C);
        
        std::cout << "d= " << d-1 << std::endl;
        showmx(C, m);
        
        Cd.push_back(C);
    }
    
    
    // generate points
    std::vector<double> samples;
    
    // 0. scrambling seeds
    std::random_device hwseed;
    std::default_random_engine rng( hwseed() );
    
    std::vector<unsigned> seeds;
    for(int d= 0; d < D; d++)
        seeds.push_back( rng() );
    
    // 1. first sobol point, all dimensions
    std::vector<integer3> x3(D);
    std::vector<integer3> p3(D);
    
    // first owen point, all dimensions
    for(int d= 0; d < D; d++)
    {
        integer3 y= scramble_base3(x3[d], seeds[d], m);
        samples.push_back( y.value_double(m) );
    }
    
    
    // 2. all points
    for(int i= 1; i < N; i++)
    {
        integer3 i3= i;
        
        for(int d= 0; d < D; d++)
        {
            integer3 x= point3_digits(Cd[d], m, i3, p3[d], x3[d]);
            integer3 y= scramble_base3(x, seeds[d], m);
            samples.push_back( y.value_double(m) );
        }
    }
    
    
    // 3. export points
    std::cout << D << " dimensions.\n";
    std::cout << N << " points.\n";
    
    for(unsigned i= 0; i + D < samples.size(); i+= D)
    {
        for(int d= 0; d < D; d++)
            std::cout << samples[i+d] << " ";
    
        std::cout << std::endl;
    }
    
    return 0;
}
