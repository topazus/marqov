#ifndef MARQOV_H
#define MARQOV_H

#include <array>
#include <vector>
#include <iostream>
#include <string>
#include <functional>
#include <type_traits>
#include <utility>
#include <tuple>
#include <random>
#include <unistd.h> // provides usleep
#include <stdexcept>
#include "cachecontainer.h"

namespace MARQOV
{
template<typename Function, typename Object, typename Tuple, size_t ... I>
auto _call(Function f, Object& obj, Tuple t, std::index_sequence<I ...>) {
	return (obj.*f)(std::get<I>(t) ...);
}

template<typename Function, typename Object, typename Tuple>
auto _call(Function f, Object& obj, Tuple t) {
	static constexpr auto size = std::tuple_size<Tuple>::value;
	return _call(f, obj, t, std::make_index_sequence<size>{});
}

// ------- elementary state vector calculus

template <class StateVector>
StateVector operator + (StateVector lhs,  StateVector rhs)
{
    StateVector res(lhs);
    for(int i = 0; i < std::tuple_size<StateVector>::value; ++i)
    res[i] += rhs[i];
    return res;
}

template <class StateVector>
StateVector operator - (StateVector lhs,  StateVector rhs)
{
    StateVector res(lhs);
    for(int i = 0; i < std::tuple_size<StateVector>::value; ++i)
    res[i] -= rhs[i];
    return res;
}

inline double dot(const double& a, const double& b)
{
    return a*b;
}

template<class VecType>
inline typename VecType::value_type dot(const VecType& a, const VecType& b)
{
    typedef typename VecType::value_type FPType;
    return std::inner_product(begin(a), end(a), begin(b), 0.0);
}


template <class StateVector>
inline void reflect(StateVector& vec, const StateVector mirror)
{
	const int SymD = std::tuple_size<StateVector>::value;
	
	const double dotp = dot(vec,mirror);

	for (int i=0; i<SymD; i++) vec[i] -= 2*dotp*mirror[i];
}	

template <class Container>
inline void normalize(Container& a)
{
	typename Container::value_type tmp_abs=std::sqrt(dot(a, a));

	for (int i = 0; i < a.size(); ++i) a[i] /= tmp_abs;
}

// --------------------------- MARQOV CLASS -------------------------------

template <class Grid, class Hamiltonian>
class Marqov 
{
	public:
		typedef typename Hamiltonian::StateVector StateVector;
		typedef typename Hamiltonian::redStateVector redStateVector; // reduced StateVector 
		typedef StateVector* StateSpace;

//Local classes. We gain access to all Types of Marqov        
        
        template <typename T>
        struct ObsRetType
        {
            typedef decltype(// this line and the next two determine the value that we will pass to the vector from the return type of the measure function of the observable...
            std::declval<T>().template measure<StateSpace, Grid>(
            std::declval<StateSpace>(), std::declval<Grid>() )) RetType;
        };

template <int N, typename Tup>
struct ObsCacheTupleIter
{
    typedef decltype(std::tuple_cat(
                std::declval<typename ObsCacheTupleIter<N-1, Tup>::RetType>(),
                std::make_tuple(std::declval<CacheContainer<typename ObsRetType<typename std::tuple_element<N, Tup>::type>::RetType>>())
    )) RetType;
    
    static auto getargtuple(H5::H5File& h5file, Tup& t){
        return std::tuple_cat( 
        ObsCacheTupleIter<N-1, Tup>::getargtuple(h5file, t),
        std::make_tuple(CacheContainerArgs(h5file, std::get<N>(t).name)));}
};

template <typename Tup>
struct ObsCacheTupleIter<0, Tup>
{
    typedef std::tuple<CacheContainer<
    typename ObsRetType<typename std::tuple_element<0, Tup>::type>::RetType
    > > RetType;
    
    static auto getargtuple(H5::H5File& h5file, Tup& t){return std::make_tuple(CacheContainerArgs(h5file, std::get<0>(t).name));}
};

template <typename Tup>
struct ObsTupleToObsCacheTuple
{
    typedef typename ObsCacheTupleIter<std::tuple_size<Tup>::value-1, Tup>::RetType
    RetType;
    static auto getargtuple(H5::H5File& h5file, Tup&& t){
        return ObsCacheTupleIter<std::tuple_size<Tup>::value - 1, Tup>::getargtuple(h5file, t);
    }
};

		std::vector<std::vector<std::vector<double>>> check;
		std::vector<int> checkidxs;

		/** The initial constructor call.
         * First we have the parameters for the MARQOV class, then follows the arbitrary number of
         * arguments for a particular Hamiltonian.
         * @param lattice The instantiated lattice object
         * @param outfile Where to create the output file
         * @param mybeta the temperature that governs the Metropolis dynamics
         * @param args A template parameter pack for the Hamiltonian
         */

        template <class ...Ts>
		Marqov(Grid& lattice, std::string outfile, double mybeta, Ts&& ... args) : ham(std::forward<Ts>(args) ... ),
													grid(lattice), 
													rng(0, 1), 
													beta(mybeta),
													metro(rng), 
													dump(outfile, H5F_ACC_TRUNC ),
													obscache(ObsTupleToObsCacheTuple<ObsTs>::getargtuple(dump, ham.getobs()))
		{
//			rng.seed(15); cout << "seed is fixed!" << endl << endl;
			rng.seed(time(NULL)+std::random_device{}());
			rng.set_integer_range(lattice.size());
			statespace = new typename Hamiltonian::StateVector[lattice.size()];
		}

		// Destructor
		~Marqov() {
            delete [] statespace; dump.close();
        }
        //FIXME: Fix assignment and copying...
        Marqov(const Marqov& rhs) = delete;
        Marqov& operator=(const Marqov& rhs) = delete;
        Marqov(Marqov&& other) = default;
        Marqov& operator=(Marqov&& other) = default;

        //For reference: this template is just to cool to forget....
// 		template<size_t N = 0, typename... Ts, typename... Args>
// 		inline typename std::enable_if_t<N == sizeof...(Ts), void>
// 		marqov_measure(std::tuple<Ts...>& t, Args... args) {}
// 		
// 		template<size_t N = 0, typename... Ts, typename... Args>
// 		inline typename std::enable_if_t<N < sizeof...(Ts), void>
// 		marqov_measure(std::tuple<Ts...>& t, Args... args)
// 		{
// 		     auto retval = _call(&std::tuple_element<N, 
// 							 std::tuple<Ts...> >::type::template measure<Args...>,
// 							 std::get<N>(t), 
// 							 std::make_tuple(std::forward<Args>(args)...) );
// 			marqov_measure<N + 1, Ts...>(t, args...);
//              std::get<N>(obscache)<<retval;
// 		}
		template<size_t N = 0, typename... Ts, typename S, typename G>
		inline typename std::enable_if_t<N == sizeof...(Ts), void>
		marqov_measure(std::tuple<Ts...>& t, S& s, G&& grid) {}
		
		template<size_t N = 0, typename... Ts, typename S, typename G>
		inline typename std::enable_if_t<N < sizeof...(Ts), void>
		marqov_measure(std::tuple<Ts...>& t, S& s, G&& grid)
		{
		     auto retval = _call(&std::tuple_element<N, 
							 std::tuple<Ts...> >::type::template measure<StateSpace, G>,
							 std::get<N>(t), 
                            std::forward_as_tuple(s, grid) );
			marqov_measure<N + 1, Ts...>(t, s, grid);
             std::get<N>(obscache)<<retval;
		}


		// ----------------- consistency check ------------------

		// perform consistency check according to 
		// Hasenbusch, J. Phys. A: Math. Gen. 34 8221 (2001)
		
		void perform_consistency_check(std::vector<int>& checkidxs)
		{
			std::vector<std::vector<double>> subcheck;

			for (int k=0; k<checkidxs.size(); k++)
			{
				const int checkidx = checkidxs[k];
				std::vector<double> subsubcheck;

				const auto checksite = statespace[checkidx];
				const auto nbrs = grid.getnbrs(0, checkidx);

				for (int i = 0; i < nbrs.size(); ++i)
				{
					const auto currentnbr = statespace[nbrs[i]];
					subsubcheck.push_back(dot(checksite,currentnbr));
				}

				const double selfdot = dot(checksite,checksite);

				subsubcheck.push_back(selfdot);
				subsubcheck.push_back((selfdot-1)*selfdot);

				subcheck.push_back(subsubcheck);
			}

			check.push_back(subcheck);

			// monitor memory consumption
			const int nsites = checkidxs.size();
			const int nmeasure = check.size(); 
			const int ncol = check[0][0].size(); 

			constexpr static double check_GB_limit = 4.0;

			if (nsites*nmeasure*ncol > check_GB_limit*1024*1024*1024/8) 
			{
				throw std::overflow_error(
				"\n Running out of memory during consistency check! \n Decrease either lattice size or number of measurements!");
			}
		}


		// evaluate check and display results
		void finalize_consistency_check()
		{
			const int SymD = std::tuple_size<StateVector>::value;
			const int ncol = 8;
			const int nmeasure = check.size();
			const int nsites = check[0].size();

			std::vector<double> sum(ncol,0);

			// compute averages in each column
			for (int k=0; k<nmeasure; k++)
			{
				for (int i=0; i<nsites; i++)
				{
					for (int j=0; j<ncol; j++)
					{
						sum[j] += check[k][i][j];
					}
				}
			}
			
			for (int j=0; j<ncol; j++) 
			{
				sum[j] = sum[j] / double(nmeasure) / double(nsites);
				std::cout << sum[j] << " ";
			}
			std::cout << std::endl;

			// summation formula
			double retval = 0.5*ham.beta*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]) - sum[6] - 2*ham.lambda*sum[7] + 0.5*SymD;
			std::cout << retval << "\n\n";

		}
					

		// specificy site indices which will enter the check
		// default: all (good statistic, but requires a somewhat large amount of memory)
		void prepare_consistency_check(std::vector<int>& checkidxs)
		{
			for (int i=0; i<grid.size(); i++)
			{
				checkidxs.push_back(i);
			}
		}


		
		// ----------------- consistency check end ------------------




		double elementaryMCstep(const int ncluster, const int nsweeps);
	    
	    	void gameloop(const int nsteps, const int ncluster, const int nsweeps, int myid)
		{

//			prepare_consistency_check(checkidxs);

			double avgclustersize = 0;
			for (int k=0; k<10; k++)
			{

				if (myid == 0) std::cout << "." << std::flush;
				for (int i=0; i<nsteps/10; ++i)
				{
					avgclustersize += elementaryMCstep(ncluster, nsweeps);
					auto obs = ham.getobs();
					marqov_measure(obs, statespace, grid);
//					perform_consistency_check(checkidxs);
				}
			}

			if (myid == 0) std::cout << "|\n" << avgclustersize/nsteps << std::endl;
//			finalize_consistency_check();
		}
	
	    	void wrmploop(const int nsteps, const int ncluster, const int nsweeps, int myid)
		{
			if (myid == 0) std::cout << "|";
			for (int k=0; k<10; k++)
			{
				if (myid == 0) std::cout << "." << std::flush;
				for (int i=0; i<nsteps/10; ++i) elementaryMCstep(ncluster, nsweeps);
			}
			if (myid == 0) std::cout << "|";
		}
	



		// -------------- special purpose functions ----------------


	    	void debugloop(const int nsteps, const int ncluster, const int nsweeps)
		{
			double avgclustersize = 0;
			for (int i=0; i<nsteps; ++i)
			{
				avgclustersize += elementaryMCstep(ncluster, nsweeps);
			}
			std::cout << avgclustersize/nsteps << std::endl;
		}



		// use carefully, might generate a lot of terminal output
	    	void gameloop_liveview(int nframes = 100, int nsweepsbetweenframe = 5)
		{
			const int ncluster = 0;
			const int nsweeps = 1;
			for (int i = 0; i < nframes; ++i)
			{

				for (int j=0; j<nsweepsbetweenframe; j++) elementaryMCstep(ncluster, nsweeps);
					
				unsigned int microsec = 30000; 
				usleep(microsec);
				system("tput reset");
				
				visualize_state_2d();
			}
		}
	


		void visualize_state_2d(int dim=2, double threshold=0.3)
		{
			std::cout << "_";
			for(int i = 0; i < grid.length; ++i) std::cout << " _";
			std::cout <<"\n";
			for(int i = 0; i < grid.length; ++i)
			{
				std::cout << "|";
				for(int j = 0; j < grid.length; ++j)
				{
					int curridx = grid.length*i+j;
					double current = statespace[curridx][dim];

					if (current > threshold) std::cout << "O ";
					else if (current < -threshold) std::cout << "  ";
					else if (current > 0) std::cout << "o ";
					else if (current < 0) std::cout << ". ";
				}
				std::cout << "|\n";
			}
			std::cout << "‾";
			for(int i = 0; i < grid.length; ++i) std::cout << " ‾";
			std::cout <<"\n\n";
		}


	
		 void init_cold_Ising_like()
		 {
		 	const int SymD = std::tuple_size<StateVector>::value;
			for(int i = 0; i < grid.size(); ++i)
			{
				for(int j = 0; j < SymD; ++j)
				{
					statespace[i][j] = 1;
				}
			}
		 }

		 void init_cold_Heisenberg()
		 {
			for(int i = 0; i < grid.size(); ++i)
			{
				statespace[i][0] = -1;
				statespace[i][1] = 0;
				statespace[i][2] = 0;
			}
		 }

		 void init_hot()
		 {
		 	const int SymD = std::tuple_size<StateVector>::value;
//			for(decltype(grid.size()) i = 0; i < grid.size(); ++i)
			for(int i = 0; i < grid.size(); ++i)
			{
				statespace[i] = rnddir<RND, typename StateVector::value_type, SymD>(rng);
			}
		 }

	private:



	inline int metropolisstep(int rsite);

	template <typename callable1, typename callable2>
	inline int metropolisstep(int rsite, callable1 filter_ref, callable2 filter_copy, int comp);

	inline int wolffstep(int rsite, const StateVector& rdir);
	inline int wolffstep_Ising(int rsite);
	inline int wolffstep_Heisenberg(int rsite, const StateVector& rdir);
	template <typename DirType>
	inline int wolffstep_general(int rsite, const DirType& rdir);

	StateSpace statespace;
	Hamiltonian ham;
    typedef decltype(std::declval<Hamiltonian>().getobs()) ObsTs;

    H5::H5File dump;///< The handle for the HDF5 file. must be before the obscaches
    typename ObsTupleToObsCacheTuple<ObsTs>::RetType obscache;
	Grid& grid;
	RND rng;
    double beta;

	//Get the MetroInitializer from the user, It's required to have one template argument left, the RNG.
	typename Hamiltonian::template MetroInitializer<RND> metro;//C++11

	// number of EMCS
	static constexpr int nstep = 250;
};

#include "update.h"
#include "emcs.h"
}
#endif
