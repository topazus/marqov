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
#include <H5Cpp.h>
#include <H5File.h>
#include <unistd.h> // provides usleep
#include <stdexcept>

using std::cout;
using std::endl;
using std::flush;

template<typename Function, typename Object, typename Tuple, size_t ... I>
auto _call(Function f, Object& obj, Tuple t, std::index_sequence<I ...>) {
	return (obj.*f)(std::get<I>(t) ...);
}

template<typename Function, typename Object, typename Tuple>
auto _call(Function f, Object& obj, Tuple t) {
	static constexpr auto size = std::tuple_size<Tuple>::value;
	return _call(f, obj, t, std::make_index_sequence<size>{});
}

// --------------------------- MARQOV CLASS -------------------------------

template <class Grid, class Hamiltonian>
class Marqov 
{
	public:
		typedef typename Hamiltonian::StateVector StateVector;
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
struct TupleIter
{
    typedef decltype(std::tuple_cat(
                std::declval<typename TupleIter<N-1, Tup>::RetType>(),
                std::make_tuple(std::declval<std::vector<typename ObsRetType<typename std::tuple_element<N, Tup>::type>::RetType>>())
    )) RetType;
};

template <typename Tup>
struct TupleIter<0, Tup>
{
    typedef std::tuple<std::vector<
    typename ObsRetType<typename std::tuple_element<0, Tup>::type>::RetType
    > > RetType;
};

template <typename Tup>
struct TupleToTupleVector
{
    typedef typename TupleIter<std::tuple_size<Tup>::value-1, Tup>::RetType
    RetType;
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

		/**
         * Writes out the entire current cache of observable N
         */
		template <int N>
		void writecache ()
        {
            hsize_t num = cachepos[N];
            //figure out how to properly append in HDF5
            typedef typename std::tuple_element<N, TupleCacheType>::type::value_type OutType;
			constexpr int rank = H5Mapper<OutType>::rank;
			hsize_t dims[rank] = {num};
			H5::DataSpace mspace(rank, dims, NULL);
			hsize_t start[rank] = {dssize[N]};
			dssize[N] += num;
			dataset[N].extend(&dssize[N]);
			auto filespace = dataset[N].getSpace();
			hsize_t count[rank] = {num};
			filespace.selectHyperslab(H5S_SELECT_SET, count, start);
			dataset[N].write(
                std::get<N>(obscache).data(),
                H5Mapper<OutType>::H5Type(), mspace, filespace);
			// std::cout<<std::get<N>(t).name<<" "<<retval<<std::endl;
        }

template <int N, class M>
struct ObsCacheDestructor
{
    static void call(M& m)
    {
    m.template writecache<N>();
    ObsCacheDestructor<N-1, M>::call(m);
}
};

template <class M>
struct ObsCacheDestructor<0, M>
{
    static void call(M& m)
    {
    m.template writecache<0>();
}
};

		
		template<size_t N = 0, typename... Ts>
		inline typename std::enable_if_t<N == sizeof...(Ts), void>
		marqov_createds(std::tuple<Ts...>& t){}
		
		template<size_t N = 0, typename... Ts>
		inline typename std::enable_if_t<N < sizeof...(Ts), void>
		marqov_createds(std::tuple<Ts...>& t)
		{
			typedef decltype(std::get<N>(t).template measure<decltype(statespace), Grid>(statespace, grid)) OutType;
			constexpr int rank = H5Mapper<OutType>::rank;
			hsize_t fdims[rank] = {0}; // dim sizes of ds (on disk)
			hsize_t maxdims[rank] = {H5S_UNLIMITED};
			
			
			H5::DataSpace mspace1(rank, fdims, maxdims);
			H5::DSetCreatPropList cparms;
			auto fv = H5Mapper<OutType>::fillval;
			
			hsize_t chunk_dims[1] = {4096*1024/sizeof(OutType)};//4MB chunking
			cparms.setChunk( rank, chunk_dims );
			cparms.setDeflate(9);//Best (1-9) compression
			cparms.setFillValue(  H5Mapper<OutType>::H5Type(), &fv);
			dataset[N] = dump.createDataSet(std::get<N>(t).name, H5Mapper<OutType>::H5Type(), mspace1, cparms);
			dssize[N] = 0;
			marqov_createds<N + 1, Ts...>(t);
            std::get<N>(obscache).resize(maxcache);//allocate space for 1024 entries
            cachepos[N] = 0;
		}


		std::vector<std::vector<std::vector<double>>> check;
		std::vector<int> checkidxs;

		// Constructor
        template <class ...Ts>
		Marqov(Grid& lattice, std::string outfile, Ts&& ... args) : ham(std::forward<Ts>(args) ... ),
													grid(lattice), 
													rng(0, 1), 
													metro(rng), 
													dump(outfile, H5F_ACC_TRUNC ),
													ca(ObsTupleToObsCacheTuple<ObsTs>::getargtuple(dump, ham.getobs()))
		{
//			rng.seed(15); cout << "seed is fixed!" << endl << endl;
			rng.seed(time(NULL));
			rng.set_integer_range(lattice.size());
			statespace = new typename Hamiltonian::StateVector[lattice.size()];
			auto obs = ham.getobs();
			constexpr int nobs = std::tuple_size<ObsTs>::value;
			dataset = new H5::DataSet[nobs];
			dssize = new hsize_t[nobs];

			//Now we need to register the observables with HDF5...
			marqov_createds(obs);
		}

		// Destructor
		~Marqov() {
            ObsCacheDestructor<std::tuple_size<ObsTs>::value -1, decltype(*this)>::call(*this);
            delete [] statespace; delete [] dataset; dump.close();
        }

		template<size_t N = 0, typename... Ts, typename... Args>
		inline typename std::enable_if_t<N == sizeof...(Ts), void>
		marqov_measure(std::tuple<Ts...>& t, Args... args) {}
		
		template<size_t N = 0, typename... Ts, typename... Args>
		inline typename std::enable_if_t<N < sizeof...(Ts), void>
		marqov_measure(std::tuple<Ts...>& t, Args... args)
		{
		     auto retval = _call(&std::tuple_element<N, 
							 std::tuple<Ts...> >::type::template measure<Args...>,
							 std::get<N>(t), 
							 std::make_tuple(args...) );
			marqov_measure<N + 1, Ts...>(t, args...);
            std::get<N>(obscache)[cachepos[N]] = retval;
            cachepos[N] = cachepos[N] + 1;
            if (cachepos[N] >= maxcache)
            {
                writecache<N>();
                cachepos[N] = 0;
            }
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
				cout << sum[j] << " ";
			}
			cout << endl;

			// summation formula
			double retval = 0.5*ham.beta*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]) - sum[6] - 2*ham.lambda*sum[7] + 0.5*SymD;
			cout << retval << endl << endl;

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
	    
	    	void gameloop(const int nsteps, const int ncluster, const int nsweeps)
		{

//			prepare_consistency_check(checkidxs);

			double avgclustersize = 0;
			for (int k=0; k<10; k++)
			{
				cout << "." << flush;
				for (int i=0; i<nsteps/10; ++i)
				{
					avgclustersize += elementaryMCstep(ncluster, nsweeps);
					auto obs = ham.getobs();
					marqov_measure(obs, statespace, grid);
//					perform_consistency_check(checkidxs);
				}
			}

			cout << "|" << endl;
			cout << avgclustersize/nsteps << endl;
//			finalize_consistency_check();
		}
	
	    	void wrmploop(const int nsteps, const int ncluster, const int nsweeps)
		{
			cout << "|";
			for (int k=0; k<10; k++)
			{
				cout << "." << flush;
				for (int i=0; i<nsteps/10; ++i) elementaryMCstep(ncluster, nsweeps);
			}
			cout << "|";
		}
	



		// -------------- special purpose functions ----------------


	    	void debugloop(const int nsteps, const int ncluster, const int nsweeps)
		{
			double avgclustersize = 0;
			for (int i=0; i<nsteps; ++i)
			{
				avgclustersize += elementaryMCstep(ncluster, nsweeps);
			}
			cout << avgclustersize/nsteps << endl;
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
			cout << "_";
			for(int i = 0; i < grid.length; ++i) cout << " _";
			cout << endl;
			for(int i = 0; i < grid.length; ++i)
			{
				cout << "|";
				for(int j = 0; j < grid.length; ++j)
				{
					int curridx = grid.length*i+j;
					double current = statespace[curridx][dim];

					if (current > threshold) cout << "O ";
					else if (current < -threshold) cout << "  ";
					else if (current > 0) cout << "o ";
					else if (current < 0) cout << ". ";
				}
				cout << "|" << endl;
			}
			cout << "‾";
			for(int i = 0; i < grid.length; ++i) cout << " ‾";
			cout << endl << endl;
		}


	
		 // only for the Ising model so far!
		 void init_cold()
		 {
			for(int i = 0; i < grid.size(); ++i)
			{
				statespace[i][0] = -1;
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
			for(decltype(grid.size()) i = 0; i < grid.size(); ++i)
			{
				statespace[i] = rnddir<RND, typename StateVector::value_type, SymD>(rng);
			}
		 }

	private:



	inline int metropolisstep(int rsite);
	inline int wolffstep(int rsite, const StateVector& rdir);
	inline int wolffstep_Ising(int rsite);
	inline int wolffstep_Heisenberg(int rsite, const StateVector& rdir);
	template <typename DirType>
	inline int wolffstep_general(int rsite, const DirType& rdir);



	StateSpace statespace;
	Hamiltonian ham;
    typedef decltype(std::declval<Hamiltonian>().getobs()) ObsTs;
    typedef typename TupleToTupleVector<ObsTs>::RetType TupleCacheType;
    constexpr static int maxcache=1024;

    H5::H5File dump;///< The handle for the HDF5 file. must be before the obscaches
    typename ObsTupleToObsCacheTuple<ObsTs>::RetType ca;
	Grid& grid;
	RND rng;

	H5::DataSet* dataset;
    TupleCacheType obscache;
    int cachepos[std::tuple_size<ObsTs>::value];
	hsize_t* dssize;

	//Get the MetroInitializer from the user, It's required to have one template argument left, the RNG.
	typename Hamiltonian::template MetroInitializer<RND> metro;//C++11


	// obs now handled differently

	// number of EMCS
	static constexpr int nstep = 250;
};

#include "update.h"
#include "emcs.h"

#endif
