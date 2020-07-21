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
#include <chrono>
#include <unistd.h> // provides usleep
#include <stdexcept>
#include <random>
#include <ctime>
#include "cachecontainer.h"
#include "svmath.h"
#include "rngcache.h"

namespace MARQOV
{
    template <typename T>
    void dumpscalartoH5(H5::Group& h5loc, const T& s, std::string name)
    {

        H5::DataSpace dspace(H5S_SCALAR); // create a scalar data space
        H5::DataSet dset(h5loc.createDataSet(name.c_str(), H5Mapper<T>::H5Type(), dspace));
        dset.write(&s, H5Mapper<T>::H5Type());
    }

	struct MARQOVConfig
	{
		/**
		* The standard constructor. It requires an outpath, the rest of the positional parameters are optional.
		*/
		
		MARQOVConfig(	std::string op, 
					int i = 0, 
					int ri = 0, 
		  			int s = 0, 
		  			int ugli = 10, 
		  			int nst = 250, 
		  			int ws = 100, 
		  			int gls = 200, 
		  			int nc = 20, 
		  			int nsw = 10) : outpath(op), 
		  		 			  	 id(i), 
		  					  	 repid(ri),
		  					  	 seed(s), 
		  					  	 gli(ugli), 
		  					  	 warmupsteps(ws), 
		  					  	 gameloopsteps(gls), 
		  					  	 ncluster(nc), 
		  					  	 nsweeps(nsw) {}
		
		MARQOVConfig(const MARQOVConfig& rhs) = default; // < FIXME: Think about wether we can get rid of it.
		MARQOVConfig& operator=(const MARQOVConfig& rhs) = delete;
		MARQOVConfig(MARQOVConfig&& other) = default;
		MARQOVConfig& operator=(MARQOVConfig&& other) = default;
		
		
		// Output
		std::string outname; // the output filename; is empty but will be specified by a filter!
		std::string outpath; // the outpath; full filename will be "outpath/outfile.h5"
		std::string logpath; // the logpath. For lack of a better place it is currently stored here.
		
		
		// MC variables
		int id;
		int repid;
		int seed; ///< Doing this correctly opens a whole can of worms.... At one point we need to dump the state of the RNG for restart.
		int gli; ///< The unknown gameloop integer
		int nsteps;
		int warmupsteps;
		int gameloopsteps;
		int ncluster;
		int nsweeps;
		
		
		/** A chain of setters to emulate the named parameter idiom.*/
		MARQOVConfig& setid(int i) {id = i; return *this;}
		MARQOVConfig& setrepid(int ri) {repid = ri; return *this;}
		MARQOVConfig& setseed(int s) {seed = s; return *this;}
		MARQOVConfig& setgli(int c) {gli = c; return *this;}
		MARQOVConfig& setnsteps(int ns) {nsteps = ns; return *this;}
		MARQOVConfig& setwarmupsteps(int w) {warmupsteps = w; return *this;}
		MARQOVConfig& setgameloopsteps(int g) {gameloopsteps = g; return *this;}
		MARQOVConfig& setncluster(int nc) {ncluster = nc; return *this;}
		MARQOVConfig& setnsweeps(int ns) {nsweeps = ns; return *this;}
		void dumpparamstoH5(H5::Group& mcg) const
        {
            mcg.setComment("Here we store all parameters that are in the MARQOVconfig object. They are mostly method related numbers and strings");
            dumpscalartoH5(mcg, id, "id");
            dumpscalartoH5(mcg, repid, "repid");
            dumpscalartoH5(mcg, seed, "seed");
            dumpscalartoH5(mcg, gli, "gli");
            dumpscalartoH5(mcg, nsteps, "nsteps");
            dumpscalartoH5(mcg, warmupsteps, "warmupsteps");
            dumpscalartoH5(mcg, gameloopsteps, "gameloopsteps");
            dumpscalartoH5(mcg, ncluster, "ncluster");
            dumpscalartoH5(mcg, nsweeps, "nsweeps");
        };
	};
    
	template<class> 
	struct type_sink { typedef void type; }; // consumes a type, and makes it `void`
	
	template<class T> using type_sink_t = typename type_sink<T>::type;
    
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
				
				L grid;
		};
		        
		//A helper to decide whether a Hamiltonian provides a init function
		template<class StateSpace, class H, class L, class R, class=void, class... Ts> struct has_init : std::false_type {};
		template<class StateSpace, class H, class L, class RNG, class... Ts>
		struct has_init<StateSpace, H, L, RNG,
		type_sink_t< decltype( std::declval<H>().template initstatespace<StateSpace, L, RNG, Ts...>(std::declval<StateSpace&>(), std::declval<L&>(), std::declval<RNG&>(), std::declval<Ts>()... ) ) >, Ts... > : std::true_type {};
        
        //A helper to decide whether a Hamiltonian provides the paramname function
        template<class H, class = void> struct has_paramname : std::false_type {};
        
        template<class H>
		struct has_paramname<H,
		type_sink_t<decltype( std::declval<H>().paramname(std::declval<int>()) )> > : std::true_type {};
        
	};
    
    
	template<typename Function, typename Object, typename Tuple, size_t ... I>
	auto _call(Function f, Object& obj, Tuple t, std::index_sequence<I ...>) 
	{
		return (obj.*f)(std::get<I>(t) ...);
	}
	
	template<typename Function, typename Object, typename Tuple>
	auto _call(Function f, Object& obj, Tuple t) 
	{
		static constexpr auto size = std::tuple_size<Tuple>::value;
		return _call(f, obj, t, std::make_index_sequence<size>{});
	}




// --------------------------- MARQOV CLASS -------------------------------

template <class Grid, class Hamiltonian, template<class> class RefType = detail::Ref >
class Marqov : public RefType<Grid>
{
	public:
		typedef typename Hamiltonian::StateVector StateVector;
		typedef int redStateVector; // reduced StateVector (so far needed only for AT model, improve me!!!)
		typedef StateVector* StateSpace;

		// Local classes. We gain access to all Types of Marqov        
        
		template <typename T>
		struct ObsRetType
		{
			// this determines the value that we will pass to the vector from the return type of the measure function of the observable...
			typedef decltype(
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
		    
		    static auto getargtuple(H5::Group& h5loc, Tup& t){
		        return std::tuple_cat( 
		        ObsCacheTupleIter<N-1, Tup>::getargtuple(h5loc, t),
		        std::make_tuple(CacheContainerArgs(h5loc, std::get<N>(t).name)));}
		};
		
		template <typename Tup>
		struct ObsCacheTupleIter<0, Tup>
		{
		    typedef std::tuple<CacheContainer<
		    typename ObsRetType<typename std::tuple_element<0, Tup>::type>::RetType
		    > > RetType;
		    
		    static auto getargtuple(H5::Group& h5loc, Tup& t){return std::make_tuple(CacheContainerArgs(h5loc, std::get<0>(t).name));}
		};
		
		template <typename Tup>
		struct ObsTupleToObsCacheTuple
		{
		    typedef typename ObsCacheTupleIter<std::tuple_size<Tup>::value-1, Tup>::RetType
		    RetType;
		    static auto getargtuple(H5::Group& h5loc, Tup&& t){
		        return ObsCacheTupleIter<std::tuple_size<Tup>::value - 1, Tup>::getargtuple(h5loc, t);
		    }
		};



	
	/** ----- The original constructor call -----
	* First we have the parameters for the MARQOV class, then follows the arbitrary number of
	* arguments for a particular Hamiltonian.
	* @param lattice The instantiated lattice object
	* @param outfile Where to create the output file
	* @param mybeta the temperature that governs the Metropolis dynamics
	* @param args A template parameter pack for the Hamiltonian
	*/
	
	template <class ...HArgs>
	Marqov(Grid& lattice, MARQOVConfig mc, double mybeta, HArgs&& ... hargs) : 
		RefType<Grid>(std::forward<Grid>(lattice)),
		ham(std::forward<HArgs>(hargs) ... ),
		mcfg(mc),
		step(-1),
		beta(mybeta),
		dump(setupHDF5Container(mc, std::forward<HArgs>(hargs)...)),
		stategroup(dump.openGroup("/step"+std::to_string(step)+"/state")),
		obsgroup(dump.openGroup("/step"+std::to_string(step)+"/observables")),
		obscache(ObsTupleToObsCacheTuple<ObsTs>::getargtuple(obsgroup, ham.getobs())),
        obs(ham.getobs()),
		rngcache(time(NULL)+std::random_device{}()),
		metro(rngcache),
		statespace(setupstatespace(lattice.size()))
	{}
		
		
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
		step(0),
		beta(mybeta),
		dump(setupHDF5Container(mc, std::forward<HArgs>(hargs)...)),
		stategroup(dump.openGroup("/step"+std::to_string(step)+"/state")),
		obsgroup(dump.openGroup("/step"+std::to_string(step)+"/observables")),
		obscache(ObsTupleToObsCacheTuple<ObsTs>::getargtuple(obsgroup, ham.getobs())),
		obs(ham.getobs()),
		rngcache(time(NULL)+std::random_device{}()), 
		metro(rngcache),
                statespace(setupstatespace(this->grid.size()))
	{}

auto setupstatespace(int size)
{
    auto retval = new typename Hamiltonian::StateVector[size];
    if (step > 0)
    {
        std::cout<<"previous data found! Continuing simulation at step "<<step<<std::endl;
        //read in the state space
        {
        auto stateds = stategroup.openDataSet("hamiltonianstatespace");
        auto dataspace = stateds.getSpace();
        //read the data... For now we just hope that everything matches...
        int rank = dataspace.getSimpleExtentNdims();
        hsize_t dims_out[rank];
        int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);

        hsize_t fdims[rank] = {static_cast<hsize_t>(size)};
        hsize_t maxdims[rank] = {H5S_UNLIMITED};

        H5::DataSpace mspace1(rank, fdims, maxdims);

        hsize_t start[rank] = {0};
        hsize_t count[rank] = {static_cast<hsize_t>(size)};
        dataspace.selectHyperslab(H5S_SELECT_SET, count, start);
        stateds.read(statespace, H5Mapper<StateVector>::H5Type(), mspace1, dataspace);
        }
        
        //FIXME read in the state of the RNG
        std::vector<int64_t> rngstate;
        {
        auto stateds = stategroup.openDataSet("rngstate");
        auto dataspace = stateds.getSpace();
        //read the data... For now we just hope that everything matches...
        int rank = dataspace.getSimpleExtentNdims();
        hsize_t dims_out[rank];
        int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);
        
        rngstate.resize(dims_out[0]);
        
        hsize_t maxdims[rank];
        for(int i = 0; i < rank; ++i)
        {
           maxdims[i] = H5S_UNLIMITED;
        }

        H5::DataSpace mspace1(rank, dims_out, maxdims);

        hsize_t start[rank] = {0};
//        hsize_t count[rank] = {static_cast<hsize_t>(size)};
        dataspace.selectHyperslab(H5S_SELECT_SET, dims_out, start);
        stateds.read(rngstate.data(), H5Mapper<StateVector>::H5Type(), mspace1, dataspace);
        }
        rngcache.setstate(rngstate);
    }
    else
    {
        //FIXME: initial seed!!
    }
    return retval;
}

/* Helper function for HDF5
 * Operator function to find the last step
 */
static herr_t
findstep(hid_t loc_id, const char *name, const H5L_info_t *linfo, void *step)
{
    std::string gname(name);
    if (gname.substr(0, 4) == "step")
    {
    std::string rem(gname.substr(4));
    int c = std::stoi(rem);
     *static_cast<int*>(step) = std::max(*static_cast<int*>(step), c);
    }
    return 0;
}

	/** This function sets up the layout of the HDF5 Container.
     * @param mc. A MARQOVConfig object that we will dump to the respective path
     * @return An object for the HDF5 File
     */
    template <class ...HArgs>
    H5::H5File setupHDF5Container(const MARQOVConfig& mc, HArgs&& ...hargs)
    {
        std::string filepath = mc.outpath+mc.outname + ".h5";
        auto flag = H5F_ACC_TRUNC;
        if (std::ifstream(filepath).good() && H5::H5File::isHdf5(filepath)) flag = H5F_ACC_RDWR;
        H5::H5File retval(filepath, flag);
        if (flag == H5F_ACC_RDWR) // abuse flag
        {// We have to iterate through the root group and find the last step.
            H5Literate(retval.getId(), H5_INDEX_NAME, H5_ITER_NATIVE, NULL, findstep, &step);
        }
        step = step + 1; // If there is no file we start at 0, else we increment.
        createstep(retval, step, mc, std::forward<HArgs>(hargs)... );
        return retval;//hopefully the refcounting of HDF5 works....
    }
    
    std::string createparamname(std::false_type, int c)
    {
        return "param" + std::to_string(c);
    }
    
    std::string createparamname(std::true_type, int c)
    {
        return ham.paramname(c);
    }

    template <class ...HArgs>
    void dumphamparamstoH5(H5::Group& h5loc, HArgs&&... hargs)
    {
        H5::StrType strdatatype(H5::PredType::C_S1, ham.name.size());
         H5::DataSpace dspace(H5S_SCALAR); // create a scalar data space
         H5::DataSet dset(h5loc.createDataSet("Model", strdatatype, dspace));
         dset.setComment("This is the Model. This should correspond to a class name");
         dset.write(ham.name.c_str(), strdatatype);
         
        h5loc.setComment("These parameters are peculiar to the considered Hamiltonian.");
        dumpscalartoH5(h5loc, beta, "beta");
        //Let's dump the unknown number of unknown parameters of the Hamiltonian....
        int paramnr = 0;
        (void) std::initializer_list<int>{((void) dumpscalartoH5(h5loc, hargs,
            createparamname(typename MARQOV::detail::has_paramname<Hamiltonian>::type(), paramnr++)
        ), 0)... };
    }
    void dumplatparamstoH5(H5::Group& h5loc)
    {
        h5loc.setComment("These parameters are peculiar to the lattice at hand.");
        return;
    }
    void dumpenvtoH5(H5::Group& h5loc)
    {
     h5loc.setComment("Here we have various parameters of the host system.");
     //Acquire hostname
     char hostname[HOST_NAME_MAX];
     int err = gethostname(hostname, HOST_NAME_MAX);
     if (err == 0) // We were able to retrieve a hostname, hence we dump it
     {
         std::string hn(hostname);
         H5::StrType strdatatype(H5::PredType::C_S1, hn.size());
         H5::DataSpace dspace(H5S_SCALAR); // create a scalar data space
         H5::DataSet dset(h5loc.createDataSet("hostname", strdatatype, dspace));
         dset.write(hn.c_str(), strdatatype);
     }
         std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
         std::time_t start_time = std::chrono::system_clock::to_time_t(now);
         char timestamp[100];
         struct tm buf;
         buf = *(std::localtime(&start_time));
         std::strftime(timestamp, sizeof(timestamp), "%A %Y-%m-%d %H:%M:%S", &buf);
         std::string ts(timestamp);

         H5::StrType strdatatype(H5::PredType::C_S1, ts.size());
         H5::DataSpace dspace(H5S_SCALAR); // create a scalar data space
         H5::DataSet dset(h5loc.createDataSet("startingdate", strdatatype, dspace));
         dset.write(ts.c_str(), strdatatype);

    }
    /** This creates a single step in the HDF5 File.
     * @param file the HDF5 file where to create the step
     * @param s the number of the step
     * @param mc
     */
    template <class... HArgs>
    void createstep(H5::H5File& file, int s, const MARQOVConfig& mc, HArgs&& ...hargs)
    {
        std::string stepname = "step" + std::to_string(s);
        H5::Group step(file.createGroup(stepname));
        step.setComment("A single step encapsulates the initial config, the observable time series and the final state.");
          H5::Group environment(step.createGroup("environment"));
          dumpenvtoH5(environment);

          H5::Group config(step.createGroup("config"));
            config.setComment("Here we have all configuration related parameters that are required to start a simulation.");
            H5::Group marqovconfig(config.createGroup("marqovconfig"));
            mc.dumpparamstoH5(marqovconfig);

            H5::Group latticeconfig(config.createGroup("lattice"));
            dumplatparamstoH5(latticeconfig);
            H5::Group hamconfig(config.createGroup("hamiltonian"));
            dumphamparamstoH5(hamconfig, std::forward<HArgs>(hargs)...);


          H5::Group s1(step.createGroup("state"));
          H5::Group s2(step.createGroup("observables"));
    }

	// default state space initializer
	void init_hot()
	{
		const int SymD = std::tuple_size<StateVector>::value;
		for (decltype(this->grid.size()) i=0; i<this->grid.size(); ++i)
		{
			statespace[i] = rnddir<RNGCache<RNGType>, typename StateVector::value_type, SymD>(rngcache);
		}
	}
		
	template <typename StateSpace, class Lattice, class H, typename... Ts>
	auto haminit_helper(std::true_type, StateSpace& statespace, const Lattice& grid, H& ham, Ts&& ... ts)
	{
		return ham.initstatespace(statespace, grid, rngcache, std::forward<Ts>(ts) ...);
	}
        
	// If there's no user defined function we do a random initialization
	template <typename StateSpace, class Lattice, class H, typename... Ts>
	auto haminit_helper(std::false_type, StateSpace& statespace, const Lattice&, H& ham, Ts&& ... ts) -> void
	{
		this->init_hot();
	}
	
	template <typename... Ts>
	auto init(Ts&& ... ts)
	{
		return haminit_helper(typename detail::has_init<StateSpace, Hamiltonian, Grid, RNGCache<RNGType>, Ts... >::type(), this->statespace, this->grid, this->ham, std::forward<Ts>(ts)...);
	}
	
	void dumprng()
    {
        //FIXME: We currently lack a proper way of writing out the *type* of the RNG.
        auto rngstate = rngcache.dumpstate();
        //We interpret the rng state space as a time series of 64bit integers
        constexpr int rank = 1;
        auto len = rngstate.size();
        std::array<hsize_t, 1> fdims, maxdims;
        fdims.fill(static_cast<hsize_t>(len));
        maxdims.fill(len);

        H5::DataSpace mspace1(rank, fdims.data(), maxdims.data());
        H5::DSetCreatPropList cparms;
        auto fv = H5Mapper<int64_t>::fillval;
        
        //no compression
        
        cparms.setFillValue(H5Mapper<int64_t>::H5Type(), &fv);
        H5::DataSet dataset = stategroup.createDataSet("rngstate", H5Mapper<int64_t>::H5Type(), mspace1, cparms);
        
        auto filespace = dataset.getSpace();
        hsize_t start[rank] = {0};//This works for initialization
        std::array<hsize_t, rank> count;
        count.fill(static_cast<hsize_t>(len));
        filespace.selectHyperslab(H5S_SELECT_SET, count.data(), start);
        dataset.write(rngstate.data(), H5Mapper<int64_t>::H5Type(), mspace1, filespace);
    }
	/** A helper function to dump the entire statespace
     */
	void dumpstatespace()
    {        
        //We interpret the statespace as a time-series of lattice points points
        constexpr int rank = H5Mapper<StateVector>::rank;
        constexpr auto len = std::tuple_size<StateVector>::value;
        std::array<hsize_t, 1> fdims, maxdims, chunk_dims;
        fdims.fill(static_cast<hsize_t>(this->grid.size()));
        maxdims.fill(H5S_UNLIMITED);

        H5::DataSpace mspace1(rank, fdims.data(), maxdims.data());
        H5::DSetCreatPropList cparms;
        auto fv = H5Mapper<StateVector>::fillval;
        chunk_dims.fill(4096*1024/H5Mapper<StateVector>::bytecount);//4MB chunking
        
        cparms.setChunk( rank, chunk_dims.data() );
        cparms.setDeflate(9);//Best (1-9) compression
        cparms.setShuffle();
        cparms.setFillValue(H5Mapper<StateVector>::H5Type(), &fv);
        H5::DataSet dataset = stategroup.createDataSet("hamiltonianstatespace", H5Mapper<StateVector>::H5Type(), mspace1, cparms);
        
        
        auto filespace = dataset.getSpace();
        hsize_t start[rank] = {0};//This works for initialization
        std::array<hsize_t, rank> count;
        count.fill(static_cast<hsize_t>(this->grid.size()));
        filespace.selectHyperslab(H5S_SELECT_SET, count.data(), start);
        dataset.write(statespace, H5Mapper<StateVector>::H5Type(), mspace1, filespace);
    }
	// Destructor
	~Marqov() 
	{
        dumprng();
        dumpstatespace();
        delete [] statespace;
        dump.close();
	}

	//FIXME: Fix assignment and copying...
	Marqov(const Marqov& rhs) = delete;
	Marqov& operator=(const Marqov& rhs) = delete;
	Marqov(Marqov&& other) = default;
	Marqov& operator=(Marqov&& other) = default;

////		For reference: this template is just to cool to forget....
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



	// ------------------ update --------------------

	double elementaryMCstep();
	    
	void gameloop()
	{
		double avgclustersize = 0;
		for (int k=0; k < this->mcfg.gli; k++)
		{

			if (this->mcfg.id == 0) std::cout << "." << std::flush;
			for (int i=0; i < this->mcfg.gameloopsteps/10; ++i)
			{
				avgclustersize += elementaryMCstep();
				marqov_measure(obs, statespace, this->grid);
			}
		}

		if (this->mcfg.id == 0) std::cout << "|\n" << avgclustersize/this->mcfg.gameloopsteps << std::endl;
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
		this->mcfg.setnsweeps(nsweeps);
		this->mcfg.setncluster(ncluster);


		double avgclustersize = 0;
		for (int i=0; i<nsteps; ++i)
		{
			avgclustersize += elementaryMCstep();
		}
		std::cout << avgclustersize/nsteps << std::endl;
	}



	// use carefully, might generate a lot of terminal output
    	void gameloop_liveview(int nframes = 100, int nsweepsbetweenframe = 5)
	{
		visualize_state_2d();

		for (int i = 0; i < nframes; ++i)
		{

			for (int j=0; j<nsweepsbetweenframe; j++) elementaryMCstep();
				
			unsigned int microsec = 20000; 
			usleep(microsec);
			system("tput reset");
			
			visualize_state_2d();
		}
	}



	void visualize_state_2d(int dim=2, double threshold=0.3)
	{
		std::cout << "_";
		for(int i = 0; i < this->grid.len; ++i) std::cout << " _";
		std::cout <<"\n";
		for(int i = 0; i < this->grid.len; ++i)
		{
			std::cout << "|";
			for(int j = 0; j < this->grid.len; ++j)
			{
				int curridx = this->grid.len*i+j;
				double current = statespace[curridx][dim];

				if (current > threshold) std::cout << "O ";
				else if (current < -threshold) std::cout << "  ";
				else if (current > 0) std::cout << "o ";
				else if (current < 0) std::cout << ". ";
			}
			std::cout << "|\n";
		}
		std::cout << "‾";
		for(int i = 0; i < this->grid.len; ++i) std::cout << " ‾";
		std::cout <<"\n\n";
	}



	private:

		inline int metropolisstep(int rsite);
		int step; ///< the current step of the simulation. Used for HDF5 paths.

		template <typename callable1, typename callable2>
		inline int metropolisstep(int rsite, callable1 filter_ref, callable2 filter_copy, int comp);

		template <typename DirType>
		inline int wolffstep(int rsite, const DirType& rdir);

		StateSpace statespace;
		Hamiltonian ham;
		MARQOVConfig mcfg;
		typedef decltype(std::declval<Hamiltonian>().getobs()) ObsTs;
		double beta; ///< The inverse temperature
		H5::H5File dump; ///< The handle for the HDF5 file. must be before the obscaches
		H5::Group obsgroup; ///< The HDF5 Group of all our observables
		H5::Group stategroup; ///< The HDF5 Group where to dump the statespace
		typename ObsTupleToObsCacheTuple<ObsTs>::RetType obscache;//The HDF5 caches for each observable
        ObsTs obs; //the actual observables

		typedef std::ranlux48_base RNGType;
		RNGCache<RNGType> rngcache;///< The caching RNG


		//Get the MetroInitializer from the user, It's required to have one template argument left, the RNG.
		typename Hamiltonian::template MetroInitializer<RNGCache<RNGType> > metro;//C++11
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
