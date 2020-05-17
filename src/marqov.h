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
    struct MARQOVConfig
    {
        //
	   // The standard constructor. It requires an outpath, the rest of the positional parameters are optional.
        //

        MARQOVConfig(std::string op, 
	   			 int i = 0, 
				 int s = 0, 
				 int ugli = 10, 
				 int nst = 250, 
				 int ws = 100, 
				 int gls = 200, 
				 int nc = 20, 
				 int nsw = 10) : outpath(op), 
				 			  id(i), 
							  seed(s), 
							  gli(ugli), 
							  warmupsteps(ws), 
							  gameloopsteps(gls), 
							  ncluster(nc), 
							  nsweeps(nsw) {}

        MARQOVConfig(const MARQOVConfig& rhs) = default; //< FIXME: Think about wether we can get rid of it.
        MARQOVConfig& operator=(const MARQOVConfig& rhs) = delete;
        MARQOVConfig(MARQOVConfig&& other) = default;
        MARQOVConfig& operator=(MARQOVConfig&& other) = default;


	   // Output
	   std::string outname; // the output filename; is empty but will be specified by a filter!
        std::string outpath; // the outpath; full filename will be "outpath/outfile.h5"
        std::string logpath; // the logpath. For lack of a better place it is currently stored here.


	   // MC variables
        int id;
        int seed; ///< Doing this correctly opens a whole can of worms.... At one point we need to dump the state of the RNG for restart.
        int gli; ///< The unknown gameloop integer
        int nsteps;
        int warmupsteps;
        int gameloopsteps;
        int ncluster;
        int nsweeps;


        /** A chain of setters to emulate the named parameter idiom.*/
        MARQOVConfig& setid(int i) {id = i; return *this;}
        MARQOVConfig& setseed(int s) {seed = s; return *this;}
        MARQOVConfig& setgli(int c) {gli = c; return *this;}
        MARQOVConfig& setnsteps(int ns) {nsteps = ns; return *this;}
        MARQOVConfig& setwarmupsteps(int w) {warmupsteps = w; return *this;}
        MARQOVConfig& setgameloopsteps(int g) {gameloopsteps = g; return *this;}
        MARQOVConfig& setncluster(int nc) {ncluster = nc; return *this;}
        MARQOVConfig& setnsweeps(int ns) {nsweeps = ns; return *this;}
    };
    namespace detail
    {
        template<class> struct type_sink { typedef void type; }; // consumes a type, and makes it `void`
        template<class T> using type_sink_t = typename type_sink<T>::type;
        template<class L, class=void, class=void> struct is_Lattice : std::false_type {};

        template<class Lattice> struct is_Lattice<Lattice,
        type_sink_t< decltype( std::declval<Lattice>().getnbrs(std::declval<int>(), std::declval<int>()) ) >,
        type_sink_t< decltype( std::declval<Lattice>().size() ) >
        > : std::true_type {};

        template <class L>
        class Ref
        {
        public:
            template<class... Args>
            Ref(L&& l, Args&& ... args) : grid(l) {}
            L& grid;
        };

        template <class L>
        class NonRef
        {
        public:
            template <class ...Args>
            NonRef(std::tuple<Args...>&& args ) : NonRef(std::forward<std::tuple<Args...>>(args), 
                                                 std::make_index_sequence<std::tuple_size<typename std::remove_reference<std::tuple<Args...>>::type>::value>()) {}
            template <class ...Args, size_t... S>
            NonRef(std::tuple<Args...>&& args, std::index_sequence<S...>) : grid(std::get<S>(std::forward<std::tuple<Args...>>(args))... ) {}
//             template <class ...Args>
//             NonRef(Args&&... args) : grid(args...) {}
            L grid;
        };
    };
    
    
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

template <class Grid, class Hamiltonian, template<class> class RefType = detail::Ref >
class Marqov : public RefType<Grid>
{
	public:
		typedef typename Hamiltonian::StateVector StateVector;
		typedef int redStateVector; // reduced StateVector (so far needed only for AT model, improve me!!!)
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
	



	
	/** ----- The original constructor call -----
	* First we have the parameters for the MARQOV class, then follows the arbitrary number of
	* arguments for a particular Hamiltonian.
	* @param lattice The instantiated lattice object
	* @param outfile Where to create the output file
	* @param mybeta the temperature that governs the Metropolis dynamics
	* @param args A template parameter pack for the Hamiltonian
	*/
	
	template <class ...HArgs>
	Marqov(Grid& lattice, MARQOVConfig mc, double mybeta, HArgs&& ... args) : 
		RefType<Grid>(std::forward<Grid>(lattice)),
		ham(std::forward<HArgs>(args) ... ),
		mcfg(mc),
		rng(0, 1), 
		beta(mybeta), 
		metro(rng),  
		dump(mc.outpath+mc.outname+".h5", H5F_ACC_TRUNC ),
		obscache(ObsTupleToObsCacheTuple<ObsTs>::getargtuple(dump, ham.getobs()))
	{
		//rng.seed(15); cout << "seed is fixed!" << endl << endl;
		rng.seed(time(NULL)+std::random_device{}());
		rng.set_integer_range(lattice.size());
		statespace = new typename Hamiltonian::StateVector[lattice.size()];
	}
		
		
	/** ----- Alternate constructor -----
	* if you require Marqov to instantiate and embed the lattice for you
	* @param outfile Where to create the output file
	* @param mybeta the temperature that governs the Metropolis dynamics
	* @param p A pair containing in the second Argument the lattice parameters and in the first the Hamiltonian parameters
	*/

	template <class ...HArgs, class ... LArgs>
	Marqov(std::tuple<LArgs...>& largs, MARQOVConfig mc, double mybeta, HArgs&& ... hargs) : 
		RefType<Grid>(std::forward<std::tuple<LArgs...>>(largs)),
		ham(std::forward<HArgs>(hargs) ... ),
		mcfg(mc),
		rng(0, 1), 
		beta(mybeta), 
		metro(rng),  
		dump(mc.outpath+mc.outname+".h5", H5F_ACC_TRUNC ),
		obscache(ObsTupleToObsCacheTuple<ObsTs>::getargtuple(dump, ham.getobs()))
	{
		// rng.seed(15); cout << "seed is fixed!" << endl << endl;
		rng.seed(time(NULL)+std::random_device{}());
		rng.set_integer_range(this->grid.size());
		statespace = new typename Hamiltonian::StateVector[this->grid.size()];
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
				const auto nbrs = this->grid.getnbrs(0, checkidx);

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
			for (int i=0; i<this->grid.size(); i++)
			{
				checkidxs.push_back(i);
			}
		}


		
		// ----------------- consistency check end ------------------




		double elementaryMCstep();
	    
        void gameloop()
		{
//			prepare_consistency_check(checkidxs);

			double avgclustersize = 0;
			for (int k=0; k < this->mcfg.gli; k++)
			{

				if (this->mcfg.id == 0) std::cout << "." << std::flush;
				for (int i=0; i < this->mcfg.gameloopsteps/10; ++i)
				{
					avgclustersize += elementaryMCstep();
					auto obs = ham.getobs();
					marqov_measure(obs, statespace, this->grid);
//					perform_consistency_check(checkidxs);
				}
			}

			if (this->mcfg.id == 0) std::cout << "|\n" << avgclustersize/this->mcfg.nsteps << std::endl;
//			finalize_consistency_check();
		}
	
	    	void wrmploop()
		{
			if (this->mcfg.id == 0) std::cout << "|";
			for (int k=0; k < this->mcfg.gli; k++)
			{
				if (this->mcfg.id == 0) std::cout << "." << std::flush;
				for (int i=0; i < this->mcfg.warmupsteps/10; ++i) elementaryMCstep();
			}
			if (this->mcfg.id == 0) std::cout << "|";
		}
	

		// -------------- special purpose functions ----------------


          void full_output_2D(int dim=0)
          {
               const int LL = this->grid.length;

               ofstream os;
               os.open("../log/fullout-"+std::to_string(LL)+"-"+std::to_string(beta)+".dat");

               for (int i=0; i<LL; i++)
               {
                    for (int j=0; j<LL; j++)
                    {
                         int curridx = LL*i+j;
                         double current = statespace[curridx][dim];
                         os << current << "\t";
                    }
                    os << endl;
               }
          }




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

				for (int j=0; j<nsweepsbetweenframe; j++) elementaryMCstep();
					
				unsigned int microsec = 30000; 
				usleep(microsec);
				system("tput reset");
				
				visualize_state_2d();
			}
		}
	


		void visualize_state_2d(int dim=2, double threshold=0.3)
		{
			std::cout << "_";
			for(int i = 0; i < this->grid.length; ++i) std::cout << " _";
			std::cout <<"\n";
			for(int i = 0; i < this->grid.length; ++i)
			{
				std::cout << "|";
				for(int j = 0; j < this->grid.length; ++j)
				{
					int curridx = this->grid.length*i+j;
					double current = statespace[curridx][dim];

					if (current > threshold) std::cout << "O ";
					else if (current < -threshold) std::cout << "  ";
					else if (current > 0) std::cout << "o ";
					else if (current < 0) std::cout << ". ";
				}
				std::cout << "|\n";
			}
			std::cout << "‾";
			for(int i = 0; i < this->grid.length; ++i) std::cout << " ‾";
			std::cout <<"\n\n";
		}
	
		 void init_cold_Ising_like()
		 {
		 	const int SymD = std::tuple_size<StateVector>::value;
			for(int i = 0; i < this->grid.size(); ++i)
			{
				for(int j = 0; j < SymD; ++j)
				{
					statespace[i][j] = 1;
				}
			}
		 }

		 void init_cold_Heisenberg()
		 {
			for(int i = 0; i < this->grid.size(); ++i)
			{
				statespace[i][0] = -1;
				statespace[i][1] = 0;
				statespace[i][2] = 0;
			}
		 }

		 void init_hot()
		 {
		 	const int SymD = std::tuple_size<StateVector>::value;
			for(int i = 0; i < this->grid.size(); ++i)
			{
				statespace[i] = rnddir<RND, typename StateVector::value_type, SymD>(rng);
			}
		 }

	private:
	inline int metropolisstep(int rsite);

	template <typename callable1, typename callable2>
	inline int metropolisstep(int rsite, callable1 filter_ref, callable2 filter_copy, int comp);

	template <typename DirType>
	inline int wolffstep(int rsite, const DirType& rdir);

	StateSpace statespace;
	Hamiltonian ham;
    MARQOVConfig mcfg;
    typedef decltype(std::declval<Hamiltonian>().getobs()) ObsTs;

    H5::H5File dump;///< The handle for the HDF5 file. must be before the obscaches
    typename ObsTupleToObsCacheTuple<ObsTs>::RetType obscache;
	RND rng;
    double beta;

	//Get the MetroInitializer from the user, It's required to have one template argument left, the RNG.
	typename Hamiltonian::template MetroInitializer<RND> metro;//C++11
};

template <class H, class L, class... LArgs, class... HArgs, size_t... S>
auto makeMarqov3(MARQOVConfig& mc, std::tuple<LArgs...>&& largs, std::tuple<HArgs...> hargs, std::index_sequence<S...> )
{
    return Marqov<L, H, detail::NonRef>(largs, mc, 
                                        std::get<S>(std::forward<std::tuple<HArgs...>>(hargs))...);
}

template <class H, class L, class... LArgs, class... HArgs>
auto makeMarqov(MARQOVConfig mc, std::pair<std::tuple<LArgs...>, std::tuple<HArgs...> >& p)
{
    return makeMarqov3<H, L>(mc, std::forward<decltype(p.first)>(p.first), p.second,
        std::make_index_sequence<std::tuple_size<typename std::remove_reference<std::tuple<HArgs...>>::type>::value>()
    );
}

template <class H, class L, class ...Args>
auto makeMarqov2(std::true_type, L&& latt, MARQOVConfig&& mc, Args&& ... args)
{
    //The first argument is a Lattice-like type -> from this we infer that 
    //We get a reference to sth. already allocated
    return Marqov<L, H, detail::Ref>(latt, mc, args...);
}

template <class H, class L, class ...Args>
auto makeMarqov(L&& latt, MARQOVConfig&& mc, Args&&... args)
{
    return makeMarqov2<H>(typename detail::is_Lattice<L>::type(), latt, std::forward<MARQOVConfig>(mc), args...);
}

#include "emcs.h"
}
#endif
