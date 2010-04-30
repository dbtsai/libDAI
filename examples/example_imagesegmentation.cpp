#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <dai/alldai.h>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <CImg.h>

using namespace std;
using namespace cimg_library;
using namespace dai;

typedef boost::numeric::ublas::vector<double> ublasvector;
typedef boost::numeric::ublas::compressed_matrix<double> ublasmatrix;
typedef ublasmatrix::value_array_type::const_iterator matrix_vcit;
typedef ublasmatrix::index_array_type::const_iterator matrix_icit;


class BinaryPairwiseGM {
    public:
        size_t N;
        ublasmatrix w;
        ublasvector th;
        double logZ0;

        BinaryPairwiseGM() {}
        BinaryPairwiseGM( const FactorGraph &fg );
        BinaryPairwiseGM( size_t _N, const ublasmatrix &_w, const ublasvector &_th, double _logZ0 ) : N(_N), w(_w), th(_th), logZ0(_logZ0) {}
        BinaryPairwiseGM( const BinaryPairwiseGM &x ) : N(x.N), w(x.w), th(x.th), logZ0(x.logZ0) {};
        BinaryPairwiseGM & operator=( const BinaryPairwiseGM &x ) {
            if( this != &x ) {
                N  = x.N;
                w  = x.w;
                th = x.th;
                logZ0 = x.logZ0;
            }
            return *this;
        }
        double doBP( size_t maxiter, double tol, size_t verbose, ublasvector &m );
        FactorGraph toFactorGraph();
};


// w should be upper triangular or lower triangular
void WTh2FG( const ublasmatrix &w, const vector<double> &th, FactorGraph &fg ) {
    vector<Var>    vars;
    vector<Factor> factors;

    size_t N = th.size();
    assert( (w.size1() == N) && (w.size2() == N) );

    vars.reserve(N);
    for( size_t i = 0; i < N; i++ )
        vars.push_back(Var(i,2));

    factors.reserve( w.nnz() + N );
    // walk through the sparse array structure
    // this is similar to matlab sparse arrays
    // index2 gives the column index
    // index1 gives the starting indices for each row
    size_t i = 0;
//    cout << w << endl;
    for( size_t pos = 0; pos < w.nnz(); pos++ ) {
        while( pos == w.index1_data()[i+1] )
            i++;
        size_t j = w.index2_data()[pos];
        double w_ij = w.value_data()[pos];
//        cout << "(" << i << "," << j << "): " << w_ij << endl;
        factors.push_back( createFactorIsing( vars[i], vars[j], w_ij ) );
    }
    for( size_t i = 0; i < N; i++ )
        factors.push_back( createFactorIsing( vars[i], th[i] ) );

    fg = FactorGraph(factors);
}


template<class T>
void Image2net( const CImg<T> &img, double J, double th_min, double th_plus, double th_tol, double p_background, BinaryPairwiseGM &net ) {
    size_t dimx = img.dimx();
    size_t dimy = img.dimy();

    net.N = dimx * dimy;
    net.w = ublasmatrix(net.N,net.N,4*net.N);
    net.th = ublasvector(net.N);
    for( size_t i = 0; i < net.N; i++ )
        net.th[i] = 0.0;
    net.logZ0 = 0.0;

    CImg<float> hist = img.get_channel(0).get_histogram(256,0,255);
    size_t cum_hist = 0;
    size_t level = 0;
    for( level = 0; level < 256; level++ ) {
        cum_hist += (size_t)hist(level);
        if( cum_hist > p_background * dimx * dimy )
            break;
    }

    double th_avg = (th_min + th_plus) / 2.0;
    double th_width = (th_plus - th_min) / 2.0;
    for( size_t i = 0; i < dimx; i++ )
        for( size_t j = 0; j < dimy; j++ ) {
            if( i+1 < dimx )
                net.w(i*dimy+j, (i+1)*dimy+j) = J;
            if( i >= 1 )
                net.w(i*dimy+j, (i-1)*dimy+j) = J;
            if( j+1 < dimy )
                net.w(i*dimy+j, i*dimy+(j+1)) = J;
            if( j >= 1 )
                net.w(i*dimy+j, i*dimy+(j-1)) = J;
            double x = img(i,j);
            net.th[i*dimy+j] = th_avg + th_width * tanh((x - level)/th_tol);
/*            if( x < level )
                x = x / level * 0.5;
            else
                x = 0.5 + 0.5 * ((x - level) / (255 - level));*/
/*            if( x < level )
                x = 0.01;
            else
                x = 0.99;
            th[i*dimy+j] = 0.5 * (log(x) - log(1.0 - x));*/
        }
}


template<class T>
FactorGraph img2fg( const CImg<T> &img, double J, double th_min, double th_plus, double th_tol, double p_background ) {
    vector<Var> vars;
    vector<Factor> factors;

    size_t dimx = img.dimx();
    size_t dimy = img.dimy();
    size_t N = dimx * dimy;

    // create variables
    vars.reserve( N );
    for( size_t i = 0; i < N; i++ )
        vars.push_back( Var( i, 2 ) );

    // build histogram
    CImg<float> hist = img.get_channel(0).get_histogram(256,0,255);
    size_t cum_hist = 0;
    size_t level = 0;
    for( level = 0; level < 256; level++ ) {
        cum_hist += (size_t)hist(level);
        if( cum_hist > p_background * dimx * dimy )
            break;
    }

    // create factors
    factors.reserve( 3 * N - dimx - dimy );
    double th_avg = (th_min + th_plus) / 2.0;
    double th_width = (th_plus - th_min) / 2.0;
    for( size_t i = 0; i < dimx; i++ )
        for( size_t j = 0; j < dimy; j++ ) {
            if( i >= 1 )
                factors.push_back( createFactorIsing( vars[i*dimy+j], vars[(i-1)*dimy+j], J ) );
            if( j >= 1 )
                factors.push_back( createFactorIsing( vars[i*dimy+j], vars[i*dimy+(j-1)], J ) );
            double x = img(i,j);
            factors.push_back( createFactorIsing( vars[i*dimy+j], th_avg + th_width * tanh((x - level)/th_tol) ) );
        }
    }

    return FactorGraph( factors.begin(), factors.end(), vars.begin(), vars.end(), factors.size(), vars.size() );
}


double myBP( BinaryPairwiseGM &net, size_t maxiter, double tol, size_t verbose, ublasvector &m, size_t dimx, size_t dimy, CImgDisplay &disp );
double myMF( BinaryPairwiseGM &net, size_t maxiter, double tol, size_t verbose, ublasvector &m, size_t dimx, size_t dimy, CImgDisplay &disp );
double doInference( FactorGraph &fg, string AlgOpts, size_t maxiter, double tol, size_t verbose, ublasvector &m, size_t dimx, size_t dimy, CImgDisplay &disp );

int main(int argc,char **argv) {
    // Display program usage, when invoked from the command line with option '-h'.
    cimg_usage("Usage: example_imagesegmentation -i <inputimage1> -j <inputimage2> -o <outputimage1> -p <outputimage2> -J <J> -t <t> -s <s> -u <u> -x <x>");
    const char* file_i = cimg_option("-i","","Input image 1");
    const char* file_j = cimg_option("-j","","Input image 2");
    const char* file_o = cimg_option("-o","","Output image (with BP)");
    const char* file_p = cimg_option("-p","","Output image (without BP)");
    const double J = cimg_option("-J",0.0,"Coupling strength");
    const double th_min = cimg_option("-t",0.0,"Local evidence strength background");
    const double th_plus = cimg_option("-s",0.0,"Local evidence strength foreground");
    const double th_tol = cimg_option("-u",0.0,"Sensitivity for fore/background");
    const double p_background = cimg_option("-x",0.0,"Percentage of background in image");

    CImg<unsigned char> image1 = CImg<>(file_i);
    CImg<unsigned char> image2 = CImg<>(file_j);

    CImg<int> image3(image1);
    image3 -= image2;
    image3.abs();
    image3.norm_pointwise(1); // 1 = L1, 2 = L2, -1 = Linf
    // normalize
    for( size_t i = 0; i < image3.dimx(); i++ ) {
        for( size_t j = 0; j < image3.dimy(); j++ ) {
            int avg = 0;
            for( size_t c = 0; c < image1.dimv(); c++ )
                avg += image1(i,j,c);
            avg /= image1.dimv();
            image3(i,j,0) /= (1.0 + avg / 255.0);
        }
    }
    image3.normalize(0,255);

    CImgDisplay disp1(image1,"Input 1",0);
    CImgDisplay disp2(image2,"Input 2",0);
    CImgDisplay disp3(image3,"Absolute difference of both inputs",0);

    //BinaryPairwiseGM net;
    //Image2net( image3, J, th_min, th_plus, th_tol, p_background, net );
    FactorGraph fg = img2fg( image3, J, th_min, th_plis, th_tol, p_background );

    size_t dimx = image3.dimx();
    size_t dimy = image3.dimy();
    CImg<unsigned char> image4(dimx,dimy,1,3);

    ublasvector m;
    //net.doBP( 0, 1e-2, 3, m );
    BP bp( fg, PropertySet("[updates=SEQFIX,maxiter=0,tol=1e-9,verbose=0,logdomain=0]") );
    bp.init();
    for( size_t i = 0; i < dimx; i++ )
        for( size_t j = 0; j < dimy; j++ ) {
            unsigned char g = (unsigned char)(bp.belief(fg.var(i*dimy+j))[1] * 255.0);
 //           unsigned char g = (unsigned char)((m[i*dimy+j] + 1.0) / 2.0 * 255.0);
            if( g > 127 ) {
                image4(i,j,0) = 255;
                image4(i,j,1) = 2 * (g - 127);
                image4(i,j,2) = 2 * (g - 127);
            } else {
                image4(i,j,0) = 0;
                image4(i,j,1) = 0;
                image4(i,j,2) = 2*g;
            }
        }
    CImgDisplay disp4(image4,"Local evidence",0);
    image4.save_jpeg(file_p,100);

    // solve the problem and show intermediate steps
    CImgDisplay disp5(dimx,dimy,"Beliefs during inference",0);
    if( 1 ) {
        //FactorGraph fg = net.toFactorGraph();
        fg.WriteToFile( "joris.fg" );

        doInference( fg, "BP[updates=SEQMAX,maxiter=1,tol=1e-9,verbose=0,logdomain=0]", 1000, 1e-5, 3, m, dimx, dimy, disp5 );
        // doInference( fg, "HAK[doubleloop=0,clusters=LOOP,init=UNIFORM,loopdepth=4,tol=1e-9,maxiter=1,verbose=3]", 1000, 1e-5, 3, m, dimx, dimy, disp5 );
        // doInference( fg, "HAK[doubleloop=0,clusters=BETHE,init=UNIFORM,maxiter=1,tol=1e-9,verbose=3]", 1000, 1e-5, 3, m, dimx, dimy, disp5 );
        // doInference( fg, "MF[tol=1e-9,maxiter=1,damping=0.0,init=RANDOM,updates=NAIVE]", 1000, 1e-5, 3, m, dimx, dimy, disp5 );
    } else {
    //  myBP( net, 1000, 1e-5, 3, m, dimx, dimy, disp5 );
    //  myMF( net, 1000, 1e-5, 3, m, dimx, dimy, disp5 );
    }

    for( size_t i = 0; i < dimx; i++ )
        for( size_t j = 0; j < dimy; j++ ) {
//          unsigned char g = (unsigned char)(bp.belief(fg.var(i*dimy+j))[1] * 255.0);
            unsigned char g = (unsigned char)((m[i*dimy+j] + 1.0) / 2.0 * 255.0);
            if( g > 127 ) {
                image4(i,j,0) = image2(i,j,0);
                image4(i,j,1) = image2(i,j,1);
                image4(i,j,2) = image2(i,j,2);
            } else
                for( size_t c = 0; c < (size_t)image4.dimv(); c++ )
                    image4(i,j,c) = 255;
/*              if( g > 127 ) {
                    image4(i,j,0) = image4(i,j,1) = image4(i,j,2) = 0;
                } else {
                    image4(i,j,0) = image4(i,j,1) = image4(i,j,2) = 255;
                }*/
        }
    CImgDisplay main_disp(image4,"Segmentation result",0);
    image4.save_jpeg(file_o,100);

    while( !main_disp.is_closed )
        cimg::wait( 40 );

    return 0;
}


double myBP( BinaryPairwiseGM &net, size_t maxiter, double tol, size_t verbose, ublasvector &m, size_t dimx, size_t dimy, CImgDisplay &disp ) {
    clock_t tic = toc();

    if( verbose >= 1 )
        cout << "Starting myBP..." << endl;

    size_t nr_messages = net.w.nnz();
    ublasmatrix message( net.w );
    for( size_t ij = 0; ij < nr_messages; ij++ )
        message.value_data()[ij] = 0.0;
    // NOTE: message(i,j) is \mu_{j\to i}

    m = ublasvector(net.N);

    size_t _iterations = 0;
    double max_diff = 1.0;
    for( _iterations = 0; _iterations < maxiter && max_diff > tol; _iterations++ ) {
        // walk through the sparse array structure
        // this is similar to matlab sparse arrays
        // index2 gives the column index (ir in matlab)
        // index1 gives the starting indices for each row (jc in matlab)
//        for( size_t t = 0; t < 3; t++ ) {
            size_t i = 0;
            max_diff = 0.0;
            for( size_t pos = 0; pos < nr_messages; pos++ ) {
                while( pos == net.w.index1_data()[i+1] )
                    i++;
                size_t j = net.w.index2_data()[pos];
                double w_ij = net.w.value_data()[pos];
                // \mu_{j\to i} = \atanh \tanh w_{ij} \tanh (\theta_j + \sum_{k\in\nb{j}\setm i} \mu_{k\to j})
                double field = sum(row(message,j)) - message(j,i) + net.th[j];
                double new_message = atanh( tanh( w_ij ) * tanh( field ) );
                double diff = fabs(message(i,j) - new_message);
                if( diff > max_diff )
                    max_diff = diff;
//                if( (pos % 3) == t )
                    message(i,j) = new_message;
            }
//        }

        if( verbose >= 3 )
            cout << "myBP:  maxdiff " << max_diff << " after " << _iterations+1 << " passes" << endl;

        for( size_t j = 0; j < net.N; j++ ) {
            // m_j = \tanh (\theta_j + \sum_{k\in\nb{j}} \mu_{k\to j})
            double field = sum(row(message,j)) + net.th[j];
            m[j] = tanh( field );
        }
        CImg<unsigned char> image(dimx,dimy,1,3);
        for( size_t i = 0; i < dimx; i++ )
            for( size_t j = 0; j < dimy; j++ ) {
                // unsigned char g = (unsigned char)(bp.belief(fg.var(i*dimy+j))[1] * 255.0);
                unsigned char g = (unsigned char)((m[i*dimy+j] + 1.0) / 2.0 * 255.0);
                if( g > 127 ) {
                    image(i,j,0) = 255;
                    image(i,j,1) = 2 * (g - 127);
                    image(i,j,2) = 2 * (g - 127);
                } else {
                    image(i,j,0) = 0;
                    image(i,j,1) = 0;
                    image(i,j,2) = 2*g;
                }
            }
        disp << image;
        char filename[30] = "/tmp/movie000.jpg";
        sprintf( &filename[10], "%03ld", (long)_iterations );
        strcat( filename, ".jpg" );
        image.save_jpeg(filename,100);
    }

    if( verbose >= 1 ) {
        if( max_diff > tol ) {
            if( verbose == 1 )
                cout << endl;
                cout << "myBP:  WARNING: not converged within " << maxiter << " passes (" << toc() - tic << " clocks)...final maxdiff:" << max_diff << endl;
        } else {
            if( verbose >= 3 )
                cout << "myBP:  ";
                cout << "converged in " << _iterations << " passes (" << toc() - tic << " clocks)." << endl;
        }
    }

    return max_diff;
}


double doInference( FactorGraph& fg, string AlgOpts, size_t maxiter, double tol, size_t verbose, ublasvector &m, size_t dimx, size_t dimy, CImgDisplay &disp ) {
    InfAlg* ia = newInfAlgFromString( AlgOpts, fg );
    ia->init();

    m = ublasvector( fg.nrVars() );
    CImg<unsigned char> image(dimx,dimy,1,3);

    size_t _iterations = 0;
    double max_diff = 1.0;
    for( _iterations = 0; _iterations < maxiter && max_diff > tol; _iterations++ ) {
        max_diff = ia->run();
        for( size_t i = 0; i < fg.nrVars(); i++ )
            m[i] = ia->beliefV(i)[1] - ia->beliefV(i)[0];
        for( size_t i = 0; i < dimx; i++ )
            for( size_t j = 0; j < dimy; j++ ) {
                // unsigned char g = (unsigned char)(ia->beliefV(i*dimy+j)[1] * 255.0);
                unsigned char g = (unsigned char)((m[i*dimy+j] + 1.0) / 2.0 * 255.0);
                if( g > 127 ) {
                    image(i,j,0) = 255;
                    image(i,j,1) = 2 * (g - 127);
                    image(i,j,2) = 2 * (g - 127);
                } else {
                    image(i,j,0) = 0;
                    image(i,j,1) = 0;
                    image(i,j,2) = 2*g;
                }
            }
        disp << image;
        /*
        char filename[30] = "/tmp/movie000.jpg";
        sprintf( &filename[10], "%03ld", (long)_iterations );
        strcat( filename, ".jpg" );
        image.save_jpeg(filename,100);
        */
        cout << "_iterations = " << _iterations << ", max_diff = " << max_diff << endl;
    }

    delete ia;

    return max_diff;
}


double myMF( BinaryPairwiseGM &net, size_t maxiter, double tol, size_t verbose, ublasvector &m, size_t dimx, size_t dimy, CImgDisplay &disp ) {
    clock_t tic = toc();

    if( verbose >= 1 )
        cout << "Starting myMF..." << endl;

    m = ublasvector(net.N);
    for( size_t i = 0; i < net.N; i++ )
        m[i] = 0.0;

    size_t _iterations = 0;
    double max_diff = 1.0;
    for( _iterations = 0; _iterations < maxiter && max_diff > tol; _iterations++ ) {
        max_diff = 0.0;
        for( size_t t = 0; t < net.N; t++ ) {
            size_t i = (size_t)(rnd_uniform() * net.N);
            double new_m_i = tanh(net.th[i] + inner_prod(row(net.w,i), m));
            double diff = fabs( new_m_i - m[i] );
            if( diff > max_diff )
                max_diff = diff;
            m[i] = new_m_i;
        }

        if( verbose >= 3 )
            cout << "myMF:  maxdiff " << max_diff << " after " << _iterations+1 << " passes" << endl;

        CImg<unsigned char> image(dimx,dimy,1,3);
        for( size_t i = 0; i < dimx; i++ )
            for( size_t j = 0; j < dimy; j++ ) {
                // unsigned char g = (unsigned char)(bp.belief(fg.var(i*dimy+j))[1] * 255.0);
                unsigned char g = (unsigned char)((m[i*dimy+j] + 1.0) / 2.0 * 255.0);
                if( g > 127 ) {
                    image(i,j,0) = 255;
                    image(i,j,1) = 2 * (g - 127);
                    image(i,j,2) = 2 * (g - 127);
                } else {
                    image(i,j,0) = 0;
                    image(i,j,1) = 0;
                    image(i,j,2) = 2*g;
                }
            }
        disp << image;
        char filename[30] = "/tmp/movie000.jpg";
        sprintf( &filename[10], "%03ld", (long)_iterations );
        strcat( filename, ".jpg" );
        image.save_jpeg(filename,100);
    }

    if( verbose >= 1 ) {
        if( max_diff > tol ) {
            if( verbose == 1 )
                cout << endl;
                cout << "myMF:  WARNING: not converged within " << maxiter << " passes (" << toc() - tic << " clocks)...final maxdiff:" << max_diff << endl;
        } else {
            if( verbose >= 3 )
                cout << "myMF:  ";
                cout << "converged in " << _iterations << " passes (" << toc() - tic << " clocks)." << endl;
        }
    }

    return max_diff;
}


BinaryPairwiseGM::BinaryPairwiseGM( const FactorGraph &fg ) {
    assert( fg.isPairwise() );
    assert( fg.isBinary() );

    // create w and th
    N = fg.nrVars();

    // count non_zeros in w
    size_t non_zeros = 0;
    for( size_t I = 0; I < fg.nrFactors(); I++ )
        if( fg.factor(I).vars().size() == 2 )
            non_zeros++;
    w = ublasmatrix(N, N, non_zeros * 2);

    th = ublasvector(N);
    for( size_t i = 0; i < N; i++ )
        th[i] = 0.0;
    
    logZ0 = 0.0;

    for( size_t I = 0; I < fg.nrFactors(); I++ ) {
        const Factor &psi = fg.factor(I);
        if( psi.vars().size() == 0 )
            logZ0 += dai::log( psi[0] );
        else if( psi.vars().size() == 1 ) {
            size_t i = fg.findVar( *(psi.vars().begin()) );
            th[i] += 0.5 * (dai::log(psi[1]) - dai::log(psi[0]));
            logZ0 += 0.5 * (dai::log(psi[0]) + dai::log(psi[1]));
        } else if( psi.vars().size() == 2 ) {
            size_t i = fg.findVar( *(psi.vars().begin()) );
            VarSet::const_iterator jit = psi.vars().begin();
            size_t j = fg.findVar( *(++jit) );

            double w_ij = 0.25 * (dai::log(psi[3]) + dai::log(psi[0]) - dai::log(psi[2]) - dai::log(psi[1]));
            w(i,j) += w_ij; 
            w(j,i) += w_ij; 

            th[i] += 0.25 * (dai::log(psi[3]) - dai::log(psi[2]) + dai::log(psi[1]) - dai::log(psi[0]));
            th[j] += 0.25 * (dai::log(psi[3]) - dai::log(psi[1]) + dai::log(psi[2]) - dai::log(psi[0]));

            logZ0 += 0.25 * (dai::log(psi[0]) + dai::log(psi[1]) + dai::log(psi[2]) + dai::log(psi[3]));
        }
    }
}


double BinaryPairwiseGM::doBP( size_t maxiter, double tol, size_t verbose, ublasvector &m ) {
    double tic = toc();

    if( verbose >= 1 )
        cout << "Starting BinaryPairwiseGM::doBP..." << endl;

    size_t nr_messages = w.nnz();
    ublasmatrix message( w );
    for( size_t ij = 0; ij < nr_messages; ij++ )
        message.value_data()[ij] = 0.0;
    // NOTE: message(i,j) is \mu_{j\to i}
    Real maxDiff = INFINITY;

    size_t _iterations = 0;
    for( _iterations = 0; _iterations < maxiter && maxDiff > tol; _iterations++ ) {
        // walk through the sparse array structure
        // this is similar to matlab sparse arrays
        // index2 gives the column index (ir in matlab)
        // index1 gives the starting indices for each row (jc in matlab)
        size_t i = 0;
        maxDiff = -INFINITY;
        for( size_t pos = 0; pos < nr_messages; pos++ ) {
            while( pos == w.index1_data()[i+1] )
                i++;
            size_t j = w.index2_data()[pos];
            double w_ij = w.value_data()[pos];
            // \mu_{j\to i} = \atanh \tanh w_{ij} \tanh (\theta_j + \sum_{k\in\nb{j}\setm i} \mu_{k\to j})
            double field = sum(row(message,j)) - message(j,i) + th[j];
            double new_message = atanh( tanh( w_ij ) * tanh( field ) );
            maxDiff = std::max( maxDiff, fabs(message(i,j) - new_message) );
            message(i,j) = new_message;
        }

        if( verbose >= 3 )
            cout << "BinaryPairwiseGM::doBP:  maxdiff " << maxDiff << " after " << _iterations+1 << " passes" << endl;
    }

    m = ublasvector(N);
    for( size_t j = 0; j < N; j++ ) {
        // m_j = \tanh (\theta_j + \sum_{k\in\nb{j}} \mu_{k\to j})
        double field = sum(row(message,j)) + th[j];
        m[j] = tanh( field );
    }

    if( verbose >= 1 ) {
        if( maxDiff > tol ) {
            if( verbose == 1 )
                cout << endl;
                cout << "BinaryPairwiseGM::doBP:  WARNING: not converged within " << maxiter << " passes (" << toc() - tic << " clocks)...final maxdiff:" << maxDiff << endl;
        } else {
            if( verbose >= 3 )
                cout << "BinaryPairwiseGM::doBP:  ";
                cout << "converged in " << _iterations << " passes (" << toc() - tic << " clocks)." << endl;
        }
    }

    return maxDiff;
}


FactorGraph BinaryPairwiseGM::toFactorGraph() {
    vector<Var> vars;
    vector<Factor> factors;

    // create variables
    vars.reserve( N );
    for( size_t i = 0; i < N; i++ )
        vars.push_back( Var( i, 2 ) );

    // create single-variable factors
    size_t nrE = w.nnz();
    factors.reserve( N + nrE / 2 );
    for( size_t i = 0; i < N; i++ )
        factors.push_back( createFactorIsing( vars[i], th[i] ) );

    // create pairwise factors
    // walk through the sparse array structure
    // this is similar to matlab sparse arrays
    size_t i = 0;
    for( size_t pos = 0; pos < nrE; pos++ ) {
        while( pos == w.index1_data()[i+1] )
            i++;
        size_t j = w.index2_data()[pos];
        double w_ij = w.value_data()[pos];
        if( i < j )
            factors.push_back( createFactorIsing( vars[i], vars[j], w_ij ) );
    }

    factors.front() *= dai::exp( logZ0 );

    return FactorGraph( factors.begin(), factors.end(), vars.begin(), vars.end(), factors.size(), vars.size() );
}
