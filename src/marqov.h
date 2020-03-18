#ifndef MARQOV_H
#define MARQOV_H

#include <array>
#include <vector>
#include <iostream>
#include <string>
#include <functional>
#include <type_traits>
#include <H5Cpp.h>
#include <H5File.h>
#include <unistd.h> // provides usleep

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


// ---------------- HDF5 MAPPER -----------------

template <typename T>
class H5Mapper;

template <>
class H5Mapper<double>
{
	public:
		static constexpr double fillval = 0;
		static constexpr int rank = 1;
		static auto H5Type(){return H5::PredType::NATIVE_DOUBLE;}
};

template <>
class H5Mapper<int>
{
	public:
		static constexpr int fillval = 0;
		static constexpr int rank = 1;
		static auto H5Type(){return H5::PredType::NATIVE_INT;}
};

template <typename Tp>
class H5Mapper
{
	public:
		static constexpr int fillval = H5Mapper<typename Tp::value_type>::fillval;
		static constexpr int rank = std::tuple_size<Tp>::value;
		static auto H5Type(){return H5Mapper<typename Tp::value_type>::H5Type();}
};





// --------------------------- MARQOV CLASS -------------------------------

template <class Grid, class Hamiltonian>
class Marqov 
{
	public:
		typedef typename Hamiltonian::StateVector StateVector;
		typedef StateVector* StateSpace;
		
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
		}
        
        
		// Constructor
		Marqov(Grid& lattice, double mybeta)  : ham(mybeta),  
										grid(lattice), 
										rng(0, 1), 
										metro(rng), 
										dump("dump.h5", H5F_ACC_TRUNC )
		{
			rng.seed(42);
			rng.seed(time(NULL));
			rng.set_integer_range(lattice.size());
			statespace = new typename Hamiltonian::StateVector[lattice.size()];
			auto obs = ham.getobs();
			constexpr int nobs = std::tuple_size<decltype(obs)>::value;
			dataset = new H5::DataSet[nobs];
			dssize = new hsize_t[nobs];

			//Now we need to register the observables with HDF5...
			marqov_createds(obs);
		}

		// Destructor
		~Marqov() {delete [] statespace; delete [] dataset; dump.close();}

		// Definition of an EMCS
		void elementaryMCstep()
		{
			const int nsweeps = 3;

			for (int j=0; j<nsweeps; j++)
			{
				for(int i = 0; i < grid.size(); ++i)
				{
					const int rsite = rng.i(); // choose random site -> Moved one level above
					metropolisstep(rsite);
				}
			}
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
			typedef decltype(retval) OutType;
			//figure out how to properly append in HDF5
			constexpr int rank = H5Mapper<OutType>::rank;
			hsize_t dims[rank] = {1};
			H5::DataSpace mspace(rank, dims, NULL);
			++dssize[N];
			dataset[N].extend(&dssize[N]);
			auto filespace = dataset[N].getSpace();
			hsize_t start[rank] = {dssize[N]-1};
			hsize_t count[rank] = {1};
			filespace.selectHyperslab(H5S_SELECT_SET, count, start);
			dataset[N].write(&retval, H5Mapper<OutType>::H5Type(), mspace, filespace);
			// std::cout<<std::get<N>(t).name<<" "<<retval<<std::endl;
		}
	    
	    	void gameloop(int nsteps)
		{
			for (int i = 0; i < nsteps; ++i)
			{
				elementaryMCstep();
				auto obs = ham.getobs();
				marqov_measure(obs, statespace, grid);
			}
		}
	

		// use carefully, might generate a lot of terminal output
	    	void gameloop_liveview(int nframes = 100, int nsweepsbetweenframe = 5)
		{
			for (int i = 0; i < nframes; ++i)
			{

				for (int j=0; j<nsweepsbetweenframe; j++) elementaryMCstep();
					
				unsigned int microsec = 30000; 
				usleep(microsec);
				system("tput reset");
				
				visualize_state_2d();
			}
		}
	
		void visualize_state_2d(int dim=0, double threshold=0.3)
		{
			cout << "_";
			for(int i = 0; i < grid.length; ++i) cout << " _";
			cout << endl;
			for(int i = 0; i < grid.length; ++i)
			{
				cout << "|";
				for(int j = 0; j < grid.length; ++j)
				{
					double current = statespace[grid.length*i+j][dim];

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


	
		 void init_cold()
		 {
			for(int i = 0; i < grid.size(); ++i)
			{
				statespace[i][0] = -1;
			}
		 }

		 void init_hot()
		 {
			for(int i = 0; i < grid.size(); ++i)
			{
				statespace[i] = rnddir<RND, double, 3>(rng);
			}
		 }
	
	private:


	// Single Metropolis update step statevectors on a lattice
	// returns an integer which encodes whether the flip attempt was successful (1) or not (0)
	inline int metropolisstep(int rsite)
	{
	    StateVector& svold = statespace[rsite];
	    StateVector svnew = metro.newsv(svold);

	    // interaction part
	    double interactionenergydiff = 0;
	    for(int a = 0; a < ham.Nalpha; ++a)
	    {
	        auto nbrs = grid.getnbrs(a, rsite);
	        StateVector averagevector = {0};
	
	        for (int i = 0; i < nbrs.size(); ++i)
	        {
	            auto mynbr = nbrs[i];
	            auto myvec = ham.interactions[a]->operator()(statespace[mynbr]);
	            averagevector = averagevector + myvec;
	        }
	        interactionenergydiff += ham.interactions[a]->J * (dot(svnew - svold, averagevector));
	    }
	
	    // onsite energy part
	    double onsiteenergydiff = 0;
	    for (int b = 0; b < ham.Nbeta; ++b)
	    {
            // compute the difference
            auto diff = ham.onsite[b]->operator()(svnew) - ham.onsite[b]->operator()(svold);
            // multiply the constant
            onsiteenergydiff += dot(ham.onsite[b]->h, diff);
	    }
        
	    
	    // multi-site energy
	    double multisiteenergyold = 0;
	    double multisiteenergynew = 0;
	    for (int g = 0; g < ham.Ngamma; ++g)
	    {
	        multisiteenergynew += ham.multisite[g]->operator()(svnew, rsite, statespace);//FIXME: think about this...
	        multisiteenergyold += ham.multisite[g]->operator()(svold, rsite, statespace);//FIXME: think about this...
	        //forgot k_gamma
	    }
	    
	    double dE 	= interactionenergydiff + onsiteenergydiff + (multisiteenergynew - multisiteenergyold);
	
	
	    int retval = 0;
	    if ( dE <= 0 )
	    {
	        svold = svnew;
	        retval = 1;
	    }
	    else if (rng.d() < exp(-ham.beta*dE))
	    {
	        svold = svnew;
	        retval = 1;
	    }
	    
	    return retval;
	}


	void wolff()
	{
	}


	StateSpace statespace;
	Hamiltonian ham;
	Grid& grid;
	RND rng;

	H5::H5File dump;
	H5::DataSet* dataset;
	hsize_t* dssize;

	//Get the MetroInitializer from the user, It's required to have one template argument left, the RNG.
	typename Hamiltonian::template MetroInitializer<RND> metro;//C++11

	// number of EMCS
	static constexpr int nstep = 250;
};

#endif
