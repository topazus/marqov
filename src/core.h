/* This file is part of MARQOV:
 * A modern framework for classical spin models on general topologies
 * Copyright (C) 2020-2021, The MARQOV Project
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef CORE_H
#define CORE_H

#include <array>
#include <vector>
#include <iostream>
#include <string>
#include <functional>
#include <type_traits>
#include <utility>
#include <tuple>
#include <stdexcept>
#include <random>
#include <ctime>
#include <chrono>
#include <mutex>
#include "cachecontainer.h"
#include "svmath.h"
#include "rngcache.h"
#include "timetracker.h"

namespace MARQOV
{
    /** Marqov Config.
     * 
     * We have a global marqov object that collects parameters that are special
     * for MARQOV. Hamiltonian and lattice parameters are elsewhere.
     */
	struct Config
	{
        /** Constructs a MARQOV::Config object.
         * 
         * The standard constructor. It requires an outpath, the rest of the
         * positional parameters are optional.
         * 
         * @param i id.
         * @param ri replica id.
         * @param s random number seed. Will be ignored if restarted
         * @param ugli unknown game loop integer
         * @param nst number of steps
         * @param ws warmup steps
         * @param gls gameloop steps
         * @param nc number of cluster updates
         * @param nsw number of sweeps
		 */
		Config(	std::string op, 
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
		  					  	 nsteps(nst),
		  					  	 warmupsteps(ws), 
		  					  	 gameloopsteps(gls), 
		  					  	 ncluster(nc), 
		  					  	 nsweeps(nsw) {}
		
		Config(const Config& rhs) = default; // < FIXME: Think about wether we can get rid of it.
		Config& operator=(const Config& rhs) = delete;
		Config(Config&& other) = default;
		Config& operator=(Config&& other) = default;
		
		
		// Output
		std::string outname; ///< the output filename; is empty but will be specified by a filter!
		std::string outpath; ///< the outpath; full filename will be "outpath/outfile.h5"
		std::string logpath; ///< the logpath. For lack of a better place it is currently stored here.
		

		// MC variables
		int id; ///< id
		int repid; ///< replica id
		int seed; ///< Doing this correctly opens a whole can of worms... We now dump the RNG state.
		int gli; ///< The unknown gameloop integer.
		int nsteps; ///< The number of elementary Monte Carlo steps.
		int warmupsteps; ///< The number of steps to do for warmups.
		int gameloopsteps; ///< gameloop steps.
		int ncluster; ///< number of cluster updates.
		int nsweeps; ///< number of sweeps.

		Config& setid(int i) {id = i; return *this;}

		/** Set the replica id.
         */
		Config& setrepid(int ri) {repid = ri; return *this;}
		/** Set the seed.
         */
		Config& setseed(int s) {seed = s; return *this;}

		Config& setgli(int c) {gli = c; return *this;}
		/** Set the number of steps.
         */
		Config& setnsteps(int ns) {nsteps = ns; return *this;}

		/** Set the number of warmup steps.
         */
		Config& setwarmupsteps(int w) {warmupsteps = w; return *this;}

		/** Set the number of gameloop steps
         */
		Config& setgameloopsteps(int g) {gameloopsteps = g; return *this;}
		/** Set the number of cluster updates.
         */
		Config& setncluster(int nc) {ncluster = nc; return *this;}

		/**Set the number of sweeps.
         */
		Config& setnsweeps(int ns) {nsweeps = ns; return *this;}

		/** Dump parameters to HDF5 Group.
         * 
         * @param mcg the HDF5 group where to store the parameters.
         */
		void dumpparamstoH5(H5::Group& mcg) const
        {
            mcg.setComment("Here we store all parameters that are in the MARQOV::Config object. They are mostly method related numbers and strings");
            dumpscalartoH5(mcg, "id", id);
            dumpscalartoH5(mcg, "repid", repid);
            dumpscalartoH5(mcg, "seed", seed);
            dumpscalartoH5(mcg, "gli", gli);
            dumpscalartoH5(mcg, "nsteps", nsteps);
            dumpscalartoH5(mcg, "warmupsteps", warmupsteps);
            dumpscalartoH5(mcg, "gameloopsteps", gameloopsteps);
            dumpscalartoH5(mcg, "ncluster", ncluster);
            dumpscalartoH5(mcg, "nsweeps", nsweeps);
        };
	};
    
    /**  
     * This function gathers information about the environment and dumps it into 
     * the specified HDF5 Group.
     * @see marqov.cpp
     * @param h5loc the HDF5 group where we generate all the information.
     */
    void dumpEnvironmenttoHDF5Group(H5::Group& h5loc);
    
	namespace detail
	{
        /** C++17 type sink.
         * 
         * A generic type sink from C++17.
         * It consumes a type, and makes it `void`.
         */
		template<class> struct type_sink { typedef void type; };
		template<class T> using type_sink_t = typename type_sink<T>::type;
        
        /** Helper to figure out whether a type is a lattice.
         * 
         * We have defined a lattice as
         * something that has a .nbrs() and a .size() function.
         */
		template<class L, class=void, class=void> struct is_Lattice : std::false_type {};
		
		template<class Lattice> struct is_Lattice<Lattice,
		type_sink_t< decltype( std::declval<Lattice>().nbrs(std::declval<int>(), std::declval<int>()) ) >,
		type_sink_t< decltype( std::declval<Lattice>().size() ) >
		> : std::true_type {};
		
        /** Internal Base class if MARQOV gets the lattice by reference.
         * 
         * A base class that gets used if MARQOV::Core gets the lattice by reference.
         * @tparam L the lattice that will be used
         */
		template <class L>
		class Ref
		{
			public:
				template<class... Args>
				Ref(const L&& l, Args&& ... args) : grid(l) {}
				const L& grid; ///< A reference to the external lattice. Note that the name is the same as in NonRef.
		};
		
        /** Internal base class that is used if MARQOV creates a lattice.
         * 
         * A base class that gets used if MARQOV::Core gets the parameters
         * of the lattice and constructs the lattice by itself.
         * @tparam L the lattice that will be used
         */
		template <class L>
		class NonRef
		{
			public:
				template <class ...Args>
				NonRef(std::tuple<Args...>&& args ) : NonRef(std::forward<std::tuple<Args...>>(args), 
				                                                 std::make_index_sequence<std::tuple_size<typename std::remove_reference<std::tuple<Args...>>::type>::value>()) {}
				template <class ...Args, size_t... S>
				NonRef(std::tuple<Args...>&& args, std::index_sequence<S...>) : grid(std::get<S>(std::forward<std::tuple<Args...>>(args))... ) {}
				
				const L grid;///< The storage of the Lattice. Note that the name is the same as in Ref.
		};

        /** Helper to decide whether a Hamiltonian provides an init function.
         * 
         * This init function is ued for setting up the initial state space.
         * @tparam StateSpace the type of the state space
         * @tparam H the Hamiltonian type
         * @tparam L the Lattice type
         * @tparam R the RNG
         * @tparam Ts arguments for the Hamiltonian
         */
		template<class StateSpace, class H, class L, class R, class=void, class... Ts> struct has_init : std::false_type {};
		template<class StateSpace, class H, class L, class RNG, class... Ts>
		struct has_init<StateSpace, H, L, RNG,
		type_sink_t< decltype( std::declval<H>().template initstatespace<StateSpace, L, RNG, Ts...>(std::declval<StateSpace&>(), std::declval<L&>(), std::declval<RNG&>(), std::declval<Ts>()... ) ) >, Ts... > : std::true_type {};
        
        /** A helper to decide whether a Hamiltonian provides the paramname function.
         * 
         * This is used for naming the parameter names properly in the HDF5 file.
         * @tparam H The Hamiltonian
         */
        template<class H, class = void> struct has_paramname : std::false_type {};
        
        template<class H>
		struct has_paramname<H,
		type_sink_t<decltype( std::declval<H>().paramname(std::declval<int>()) )> > : std::true_type {};
        
        /** Helper to decide whether an observable provides a description.
         * 
         * This is used for providing observable descriptions in the HDF5 file.
         * @tparam O the observable that we check.
         */
        template<class O, class = void> struct obs_has_desc : std::false_type {};
        
        template<class O>
		struct obs_has_desc<O,
		type_sink_t<decltype( std::declval<O>().desc )> > : std::true_type {};
     
        /** Create proper constructor call of the cache container.
         * 
         * The user has provided description.
         * @see createCArgTuple
         * @tparam i the index of the current observable in Tup.
         * @tparam Tup a tuple of observables.
         * @returns A tuple with proper arguments for the Cache Container.
         */
        template <int i, class Tup>
        auto createCArgTuple_impl(H5::Group& h5loc, Tup& t, std::true_type) {return std::make_tuple(CacheContainerArgs(h5loc, std::get<i>(t).name, std::get<i>(t).desc));}

        /** Create proper constructor call of the cache container.
         * 
         * The user has not provided description.
         * @see createCArgTuple
         * @tparam i the index of the current observable in Tup.
         * @tparam Tup a tuple of observables.
         * @returns A tuple with proper arguments for the Cache Container.
         */
        template <int i, class Tup>
        auto createCArgTuple_impl(H5::Group& h5loc, Tup& t, std::false_type) {return std::make_tuple(CacheContainerArgs(h5loc, std::get<i>(t).name));}
        
        /** Create proper constructor call of the cache container.
         * 
         * This selects the right call, whether the user provides description or not.
         * @see createCArgTuple_impl
         * @tparam i the index of the current observable in Tup.
         * @tparam Tup a tuple of observables.
         * @returns A tuple with proper arguments for the Cache Container.
         */        
        template <int i, class Tup>
        auto createCArgTuple(H5::Group& h5loc, Tup& t) {return detail::createCArgTuple_impl<i>(h5loc, t, typename obs_has_desc<typename std::tuple_element<i, Tup>::type>::type() );}
	};
    
	template<typename Function, typename Object, typename Tuple, size_t ... I>
	auto _call(Function f, Object& obj, Tuple t, std::index_sequence<I ...>) 
	{
		return (obj.*f)(std::get<I>(t) ...);
	}
	
	/** A helper function to call an objects member function with a tuple of arguments.
     * 
     * @tparam Function Which member function to call.
     * @tparam Object The type of the object
     * @tparam Tuple The type of the Tuple of arguments.
     * 
     * @param f the name of the function.
     * @param obj the object.
     * @param t the arguments to f.
     * @returns the return value of the function call.
     */
	template<typename Function, typename Object, typename Tuple>
	auto _call(Function f, Object& obj, Tuple t) 
	{
		static constexpr auto size = std::tuple_size<Tuple>::value;
		return _call(f, obj, t, std::make_index_sequence<size>{});
	}

// --------------------------- MARQOV::Core class -------------------------------

/** The MARQOV Core class.
 *
 * This class fuses all the parts that the user has specified and calls them to
 * life if needed. The Hamiltonian and the lattice will be instantiated
 * and depending on these the respective Monte Carlo moves will be generated.
 * I/O for the Observables will be generated.
 * The RNGCache will be initialized.
 * @tparam Grid the Lattice that we should use
 * @tparam Hamiltonian The Hamiltonian that we should use.
 * @tparam RefType used internally to distinguish, whether MARQOV::Core should
 *                 create the lattice or whether it is user provided.
 */
template <class Grid, class Hamiltonian, template<class> class RefType = detail::Ref >
class Core : public RefType<Grid>
{
	public:
        typedef Hamiltonian HamiltonianType; ///< Via HamiltonianType the used Hamiltonian is accessible.
        typedef Grid Lattice; ///< The Type of the Lattice
		typedef typename Hamiltonian::StateVector StateVector; ///< The type of the StateVector as retrieved from the Hamiltonian.
		typedef int redStateVector; // reduced StateVector (so far needed only for AT model, improve me!!!)
		typedef StateVector* StateSpace;///< the type of the state space.

		timetracker marqovtime; ///< The TimeTracker for tracking times.

		// Local classes. We gain access to all Types of MARQOV::Core
		
		/** Get the return type of the measure function of an observable.
         * 
         * @tparam T The observable.
         */
		template <typename T>
		struct ObsRetType
		{
			typedef decltype(
				std::declval<T>().template measure<StateSpace, Grid>(
					std::declval<StateSpace>(), std::declval<Grid>() )) RetType;///< This determines the value that we will pass to the vector from the return type of the measure function of the observable...
		};

        //The following three class templates help to generate from a tuple of Observables the tuple of ObscacheContainers.
        /** Map Tuple of Observables to Tuple of ObservableCaches.
         * 
         * This peels of one level of the recursion.
         * @tparam N recursion level counter
         * @tparam Tup A Tuple of observables.
         */
		template <int N, typename Tup>
		struct ObsCacheTupleIter
		{
		    typedef decltype(std::tuple_cat(
		                std::declval<typename ObsCacheTupleIter<N-1, Tup>::RetType>(),
		                std::make_tuple(std::declval<CacheContainer<typename ObsRetType<typename std::tuple_element<N, Tup>::type>::RetType>>())
		    )) RetType;
		    
		    static auto getargtuple(H5::Group& h5loc, Tup& t){
		        return std::tuple_cat(
		        ObsCacheTupleIter<N-1, Tup>::getargtuple(h5loc, t)
                , detail::createCArgTuple<N>(h5loc, t));}
		};
        /** Map Tuple of Observables to Tuple of ObservableCaches.
         * 
         * Recursion end.
         * @tparam N recursion level counter
         * @tparam Tup A Tuple of observables.
         */
		template <typename Tup>
		struct ObsCacheTupleIter<0, Tup>
		{
		    typedef std::tuple<CacheContainer<
		    typename ObsRetType<typename std::tuple_element<0, Tup>::type>::RetType
		    > > RetType;///< At recursion end this is an ObsCache with the type of the observable.

		    static auto getargtuple(H5::Group& h5loc, Tup& t){return detail::createCArgTuple<0>(h5loc, t);}
		};

        /** Map tuple of observables to tuple of ObservableCaches.
         * 
         * This template is the entry point to map a tuple of observables
         * to a tuple of ObservableCaches at Compile Time.
         * This starts the recursion.
         * @tparam Tup A tuple of observables.
         */
		template <typename Tup>
		struct ObsTupleToObsCacheTuple
		{
		    typedef typename ObsCacheTupleIter<std::tuple_size<Tup>::value-1, Tup>::RetType
		    RetType;///< this holds the type of the tuple of ObservableCaches.
		    /** Construct the tuple of Observable Caches.
             * 
             * @param h5loc the HDF5 group where we store the observables.
             * @param t the tuple of observables.
             * @returns the final tuple of observable caches.
             */
		    static auto getargtuple(H5::Group& h5loc, Tup& t){
		        return ObsCacheTupleIter<std::tuple_size<Tup>::value - 1, Tup>::getargtuple(h5loc, t);
		    }
		};

	/** Construct MARQOV  with a predefined lattice
     * 
     * First we have the parameters for the MARQOV::Core class, then follows the arbitrary number of
     * arguments for a particular Hamiltonian.
     * @tparam HArgs The arguments to the Hamiltonian.
     * 
     * @param lattice A reference to the instantiated lattice object. You are responsible for managing its lifetime.
     * @param outfile Where to create the output file.
     * @param mybeta The temperature that governs the Metropolis dynamics.
     * @param hargs A template parameter pack for the Hamiltonian.
     */
	template <class ...HArgs>
	Core(Grid&& lattice, Config mc, std::mutex& mtx, double mybeta, HArgs&& ... hargs) : 
		RefType<Grid>(std::forward<Grid>(lattice)),
		beta(mybeta),
		ham(std::forward<HArgs>(hargs) ... ),
		mcfg(mc),
		step(-1),
		statespace(setupstatespace(lattice.size())),
		hdf5lock(mtx),
		dump(setupHDF5Container(mc, std::forward<HArgs>(hargs)...)),
		obsgroup(dump.openGroup("/step"+std::to_string(step)+"/observables")),
		stategroup(dump.openGroup("/step"+std::to_string(step)+"/state")),
		obscache(ObsTupleToObsCacheTuple<ObsTs>::getargtuple(obsgroup, ham.observables)),
		obs(ham.observables),
		rngcache(time(NULL)+std::random_device{}()),
		metro(rngcache)
		{
			hdf5lock.unlock();

			marqovtime.add_clock("cluster");
			marqovtime.add_clock("local");
			marqovtime.add_clock("measurements");
			marqovtime.add_clock("other");
			//marqovtime.status();
			marqovtime.run("other");

		}

	/** Construct MARQOV and let MARQOV create the lattice.
     * 
	* If you require MARQOV::Core to instantiate and embed the lattice for you.
    * @tparam HArgs the Arguments of the Hamiltonian.
    * @tparam LArgs The Arguments of the Lattice.
	* @param outfile Where to create the output file.
	* @param mybeta the temperature that governs the Metropolis dynamics.
	* @param p A pair containing in the second argument the lattice parameters and in the first the Hamiltonian parameters.
	*/
	template <class ...HArgs, class ... LArgs>
	Core(std::tuple<LArgs...>& largs, Config mc, std::mutex& mtx, double mybeta, HArgs&& ... hargs) : 
		RefType<Grid>(std::forward<std::tuple<LArgs...>>(largs)),
		beta(mybeta),
		ham(std::forward<HArgs>(hargs) ... ),
		mcfg(mc),
		step(0),
		statespace(setupstatespace(this->grid.size())),
		hdf5lock(mtx),
		dump(setupHDF5Container(mc, std::forward<HArgs>(hargs)...)),
		obsgroup(dump.openGroup("/step"+std::to_string(step)+"/observables")),
		stategroup(dump.openGroup("/step"+std::to_string(step)+"/state")),
		obscache(ObsTupleToObsCacheTuple<ObsTs>::getargtuple(obsgroup, ham.observables)),
		obs(ham.observables),
		rngcache(time(NULL)+std::random_device{}()), 
		metro(rngcache)
		{
			hdf5lock.unlock();

			marqovtime.add_clock("cluster");
			marqovtime.add_clock("local");
			marqovtime.add_clock("measurements");
			marqovtime.add_clock("other");
			//	marqovtime.status();
			marqovtime.run("other");

		}

		/** Set up and/or reinitialize state space.
         * 
         * This sets up the state space.
         * If a previous step is present we reread that from the files.
         * @param size the number of statevectors that we want to have.
         * @returns A pointer to the state space memory.
         */
		auto setupstatespace(int size)
        {
            auto retval = new typename Hamiltonian::StateVector[size];
            if (step > 0)
            {
                std::cout<<"Previous data found! Continuing simulation at step "<<step<<std::endl;
                //read in the state space
                {
                    auto stateds = stategroup.openDataSet("hamiltonianstatespace");
                    auto dataspace = stateds.getSpace();
                    //read the data... For now we just hope that everything matches...
                    int rank = dataspace.getSimpleExtentNdims();
                    hsize_t fdims[rank], maxdims[rank], start[rank];

                    for(int i = 0; i < rank; ++i)
                    {
                        fdims[i] = static_cast<hsize_t>(size);
                        maxdims[i] = H5S_UNLIMITED;
                        start[i] = 0;
                    }

                    H5::DataSpace mspace1(rank, fdims, maxdims);
                    dataspace.selectHyperslab(H5S_SELECT_SET, fdims, start);//We have no separate count array since fdims contains identical information
                    stateds.read(statespace, H5Mapper<StateVector>::H5Type(), mspace1, dataspace);
                }
        
                //compare the used RNGs
                {
                    auto stateds = stategroup.openDataSet("RNG");
                    auto dataspace = stateds.getSpace();
                    //FIXME: proper reading of strings...
                }
        
                std::vector<int64_t> rngstate;
                {
                    auto stateds = stategroup.openDataSet("rngstate");
                    auto dataspace = stateds.getSpace();
                    //read the data... For now we just hope that everything matches...
                    int rank = dataspace.getSimpleExtentNdims();
                    hsize_t dims_out[rank];
                    dataspace.getSimpleExtentDims( dims_out, NULL);

                    rngstate.resize(dims_out[0]);
                    hsize_t maxdims[rank], start[rank];
                    for(int i = 0; i < rank; ++i)
                    {
                        maxdims[i] = H5S_UNLIMITED;
                        start[i] = 0;
                    }

                    H5::DataSpace mspace1(rank, dims_out, maxdims);

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
        /** Helper function for HDF5.
         * 
         * Operator function to find the last step in a file.
         * 
         * @param loc_id an HDF5 id.
         * @param name The group name that we are looking at.
         * @param linfo Other info from HDF5.
         * @param step the last step that is in the file.
         * @returns 0.
         */
        static herr_t findstep(hid_t loc_id, const char *name, const H5L_info_t *linfo, void *step)
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
         * 
         * @tparam HArgs The template pack of the Hamiltonian parameters.
         * 
         * @param mc. A MARQOV::Config object that we will dump to the respective path.
         * @param hargs The arguments of the Hamiltonian.
         * @return An object for the HDF5 File.
         */
        template <class ...HArgs>
        H5::H5File setupHDF5Container(const Config& mc, HArgs&& ...hargs)
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
    
        /** Create default name for Hamiltonian parameter.
         * 
         * Creates strings of the form param0, param1,...
         * @param c the index of the parameter.
         * @returns param + c
         */
        std::string createparamname(std::false_type, int c)
        {
            return "param" + std::to_string(c);
        }

        /** Use the user supplied for Hamiltonian parameter.
         * 
         * asks the hamiltonian for the name of the parameter.
         * @param c the index of the parameter.
         * @returns ham.paramname(c).
         */
        std::string createparamname(std::true_type, int c)
        {
            return ham.paramname(c);
        }
        
        /** Write hamiltonian parameters to file.
         * 
         * @tparam HArgs a template parameter pack with the variables for the Hamiltonian.
         * 
         * @param h5loc the HDF5 group where we store these names.
         * @param hargs the arguments of the hamiltonian.
         */
        template <class ...HArgs>
        void dumphamparamstoH5(H5::Group& h5loc, HArgs&&... hargs)
        {
            H5::StrType strdatatype(H5::PredType::C_S1, ham.name.size());
            H5::DataSpace dspace(H5S_SCALAR); // create a scalar data space
            H5::DataSet dset(h5loc.createDataSet("Model", strdatatype, dspace));
            dset.setComment("This is the Model. This should correspond to a class name");
            dset.write(ham.name.c_str(), strdatatype);

            h5loc.setComment("These parameters are peculiar to the considered Hamiltonian.");
            //Let's dump the unknown number of unknown parameters of the Hamiltonian....
            int paramnr = 0;
            (void) std::initializer_list<int>{((void) dumpscalartoH5(h5loc,
                createparamname(typename MARQOV::detail::has_paramname<Hamiltonian>::type(), paramnr++)
                , hargs), 0)... };
        }
        /** Write out the parameters of the lattice.
         * 
         * UNIMPLEMENTED
         * @param h5loc The group where to dump the information.
         */
        void dumplatparamstoH5(H5::Group& h5loc)
        {
            h5loc.setComment("These parameters are peculiar to the lattice at hand.");
            dumpscalartoH5(h5loc, "size", this->grid.size());
            return;
        }
        /** This function tries to dump as much useful information about the environment into 
         *  the container as possible.
         * @param h5loc The group where to dump the information.
         */
        void dumpenvtoH5(H5::Group& h5loc)
        {
            dumpEnvironmenttoHDF5Group(h5loc);
        }

        /** This creates a single step in the HDF5 File.
         *
         * @tparam HArgs the argument tuple of the Hamiltonian.
         * 
         * @param file the HDF5 file where to create the step
         * @param s the number of the step
         * @param mc the parameters of MARQOV::Core
         * @param hargs the argument tuple of the Hamiltonian
         */
        template <class... HArgs>
        void createstep(H5::H5File& file, int s, const Config& mc, HArgs&& ...hargs)
        {
            file.setComment("A calculation is made up by a series of steps. Each step can use as input the previous step.");
            std::string stepname = "step" + std::to_string(s);
            H5::Group step(file.createGroup(stepname));
            step.setComment("A single step encapsulates the initial config, the observable time series and the final state.");
            H5::Group environment(step.createGroup("environment"));
            dumpenvtoH5(environment);
            
            H5::Group config(step.createGroup("config"));
            config.setComment("Here we have all configuration related parameters that are required to start a simulation.");
            H5::Group marqovconfig(config.createGroup("marqovconfig"));
            dumpscalartoH5(marqovconfig, "beta", beta);
            mc.dumpparamstoH5(marqovconfig);

            H5::Group latticeconfig(config.createGroup("lattice"));
            dumplatparamstoH5(latticeconfig);
            H5::Group hamconfig(config.createGroup("hamiltonian"));
            dumphamparamstoH5(hamconfig, std::forward<HArgs>(hargs)...);

            H5::Group s1(step.createGroup("state"));
            H5::Group s2(step.createGroup("observables"));
        }

        /** The default StateSpace initializer.
         * 
         * This initializes with random directions.
         */ 
        void init_hot()
        {
            constexpr int SymD = std::tuple_size<StateVector>::value;
            for (decltype(this->grid.size()) i = 0; i < this->grid.size(); ++i)
            {
                statespace[i] = rnddir<RNGCache<RNGType>, typename StateVector::value_type, SymD>(rngcache);
            }
        }

        /** Select the user supplied initializer
         * 
         * @tparam StateSpace The type of the statespace
         * @tparam Lattice The Lattice
         * @tparam H The Hamiltonian
         * @tparam Ts optional arguments
         * 
         * @param statespace The statespace of MARQOV.
         * @param ham The Hamiltonian.
         * @param grid the currently used lattice.
         * @param ts optional arguments.
         */
        template <typename StateSpace, class Lattice, class H, typename... Ts>
        auto haminit_helper(std::true_type, StateSpace& statespace, const Lattice& grid, H& ham, Ts&& ... ts)
        {
            return ham.initstatespace(statespace, grid, rngcache, std::forward<Ts>(ts) ...);
        }

        /** Initializer if no user supplied initializer.
         * 
         * If there's no user defined function we do a random initialization.
         * @tparam StateSpace The type of the statespace
         * @tparam Lattice The Lattice
         * @tparam H The Hamiltonian
         * @tparam Ts optional arguments
         * 
         * @param statespace The statespace of MARQOV
         * @param ham The Hamiltonian
         * @param ts optional arguments
         */
        template <typename StateSpace, class Lattice, class H, typename... Ts>
        auto haminit_helper(std::false_type, StateSpace& statespace, const Lattice&, H& ham, Ts&& ... ts) -> void
        {
            this->init_hot();
        }
        
        /** Initialize the state space.
         * 
         * This initializes the state space.
         * If the user has supplied his/her own function we use that,
         * else we use randomized directions.
         * @tparam Ts optional parameters
         * 
         * @param ts optional parameters
         */
        template <typename... Ts>
        auto init(Ts&& ... ts)
        {
            return haminit_helper(typename detail::has_init<StateSpace, Hamiltonian, Grid, RNGCache<RNGType>, Ts... >::type(), this->statespace, this->grid, this->ham, std::forward<Ts>(ts)...);
        }

        /** Dump RNG State to HDF5.
         * 
         * This writes out the state of the RNG to the HDF5 group.
         * We assume that the state is written out as a vector of 64bit integers.
         */
        void dumprng()
        {
            H5::StrType strdatatype(H5::PredType::C_S1, RNGName<RNGType>().name.size());
            H5::DataSpace dspace(H5S_SCALAR); // create a scalar data space
            H5::DataSet dset(stategroup.createDataSet("RNG", strdatatype, dspace));
            dset.write(RNGName<RNGType>().name.c_str(), strdatatype);
        
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
        /** A helper function to dump the entire statespace.
         * 
         * The internal data layout is the same as for a time series.
         */
        void dumpstatespace()
        {
            //We interpret the statespace as a time-series of lattice points
            constexpr int rank = H5Mapper<StateVector>::rank;
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

        /** Destructor.
         * 
         * Uses the HDF5 Mutex to serialize the access to the HDF5 library and hence the output to the library.
         */
        ~Core() 
        {
            //locking is necessary since the dump functions contain HDF5 calls.
            hdf5lock.lock();
            dumprng();
            dumpstatespace();
            delete [] statespace;
            dump.close();
        }

        //FIXME: Fix assignment and copying...
        Core(const Core& rhs) = delete;
        Core& operator=(const Core& rhs) = delete;
        Core(Core&& other) = default;
        Core& operator=(Core&& other) = default;
        
        template<size_t N = 0, typename... Ts, typename S, typename G>
        inline typename std::enable_if_t<N == sizeof...(Ts), void>
        marqov_measure(std::tuple<Ts...>& t, S& s, G&& grid) {}

        /**
         * Measure observables
         * 
         * This functions measures observables and recurses into the observable tuple.
         * The measured value gets written to HDF5
         * 
         * @tparam N the current index of the observable.
         * @tparam Ts the types of all observables.
         * @tparam S an index sequence.
         * @tparam G the type of the lattice.
         * 
         * @param t a tuple with all observables
         * @param s dummy argument for the index sequence
         * @param grid the actual grid
         */
        template<size_t N = 0, typename... Ts, typename S, typename G>
        inline typename std::enable_if_t<N < sizeof...(Ts), void>
        marqov_measure(std::tuple<Ts...>& t, S& s, G&& grid)
        {
            auto retval = _call(&std::tuple_element<N, 
                            std::tuple<Ts...> >::type::template measure<StateSpace, G>,
                            std::get<N>(t), 
                            std::forward_as_tuple(s, grid) );
            marqov_measure<N + 1, Ts...>(t, s, grid);
            hdf5lock.lock();
            std::get<N>(obscache)<<retval;
            hdf5lock.unlock();
        }

        // ------------------ update --------------------

        /** The elementary Monte Carlo step.
         * 
         * @returns 
         */
        double elementaryMCstep();
        
        /** The basic Monte Carlo loop
         * 
         * A couple of cycles of updating and measuring.
         */
        void gameloop()
        {

            double avgclustersize = 0;
            for (int k=0; k < this->mcfg.gli; k++)
            {

                if (this->mcfg.id == 0) std::cout << "." << std::flush;
                for (int i=0; i < this->mcfg.gameloopsteps/10; ++i)
                {
                    avgclustersize += elementaryMCstep();
                    marqovtime.switch_clock("measurements");
                    marqov_measure(obs, statespace, this->grid);
                }
            }
            marqovtime.stop();

            marqovtime.status();
            if (this->mcfg.id == 0) std::cout << "|\n" << avgclustersize/this->mcfg.gameloopsteps << std::endl;
        }

        /** Warm up loop
         * 
         * We warm up the state space with a couple of cycles.
         */
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
        
        /** Visual output of 2D statespace to a file.
         * 
         * @param dim
         */
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

        /** The Debug loop.
         * 
         * @param nsteps The number of steps.
         * @param ncluster The number of cluster updates.
         * @param nsweeps The number of sweeps.
         */
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
        
        /** Liveview helper function
         * 
         * use carefully, might generate a lot of terminal output
         * @param nframes
         * @param nsweepsbetweenframe
         */
    	void gameloop_liveview(int nframes = 100, int nsweepsbetweenframe = 5)
        {
            visualize_state_2d();

            for (int i = 0; i < nframes; ++i)
            {

                for (int j=0; j<nsweepsbetweenframe; j++) elementaryMCstep();

                unsigned int microsec = 20000;
                std::this_thread::sleep_for(std::chrono::microseconds(microsec));

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
        /** A Metropolis step.
         * 
         * This function dispatches the call for a metropolis step 
         * to the Metropolis class.
         * @param rsite the randomly chosen site for the update.
         */
		inline int metropolisstep(int rsite);

        /** A step of the Wolff Cluster Algorithm.
         * 
         * The plane of reflection for the cluster update does not be aligned
         * to a coordinate axis, but can be chosen randomly in space.
         * @param rsite The random site where to start the cluster.
         * @param rdir The randomized plane of reflection for the cluster update.
         */
		template <typename DirType>
		inline int wolffstep(int rsite, const DirType& rdir);

		double beta; ///< The inverse temperature.
		Hamiltonian ham; ///< An instance of the user-defined Hamiltonian.
		Config mcfg; ///< An instance of all our MARQOV related parameters.
		int step; ///< The current step of the simulation. Used for HDF5 paths.
		StateSpace statespace; ///< The statespace. It holds the current configuration space.
        std::unique_lock<std::mutex> hdf5lock; ///< The global lock to synchronize access to the HDF5 *library*.
		H5::H5File dump; ///< The handle for the HDF5 file. Must be before the obscaches.
		H5::Group obsgroup; ///< The HDF5 Group of all our observables.
		H5::Group stategroup; ///< The HDF5 Group where to dump the statespace.
		typedef decltype(std::declval<Hamiltonian>().observables) ObsTs; ///< This type is mostly a tuple of other observables.
		typename ObsTupleToObsCacheTuple<ObsTs>::RetType obscache; ///< The HDF5 caches for each observable.
        ObsTs& obs; ///<the actual observables defined in the hamiltonians.

//		typedef std::ranlux48_base RNGType;
		typedef std::mt19937_64 RNGType; ///< Here the type of the RNG is set.
		RNGCache<RNGType> rngcache;///< The caching RNG


		//Get the MetroInitializer from the user, It's required to have one template argument left, the RNG.
		typename Hamiltonian::template MetroInitializer<RNGCache<RNGType> > metro;
};

/** Instantiate Core with our own lattice.
 * 
 * Internal function to unpack the template parameter pack
 * 
 * @tparam H the type of Hamiltonian
 * @tparam L the type of the lattice
 * @tparam LArgs The arguments of the lattice.
 * @tparam HArgs the arguments of the hamiltonian.
 * @tparam S a parameter pack of integers for unpacking the hamiltonian parameters
 * 
 * @param mc the MARQOVConfig object
 * @param mtx The mutex that synchronizes access to the HDF5 files.
 * @param largs The lattice arguments
 * @param hargs The hamiltonian arguments.
 */
template <class H, class L, class... LArgs, class... HArgs, size_t... S>
auto makeCore3(Config& mc, std::mutex&& mtx, std::tuple<LArgs...>&& largs, std::tuple<HArgs...> hargs, std::index_sequence<S...> )
{
    return Core<L, H, detail::NonRef>(largs, mc, mtx,
                                        std::get<S>(std::forward<std::tuple<HArgs...>>(hargs))...);
}

/** Instantiate Core with our own lattice.
 * 
 * In this call we get the arguments for the lattice 
 * and construct it ourselves.
 * 
 * @tparam H the type of Hamiltonian
 * @tparam L the type of the lattice
 * @tparam LArgs The arguments of the lattice.
 * @tparam HArgs the arguments of the hamiltonian.
 * 
 * @param mc the MARQOVConfig object
 * @param mtx The mutex that synchronizes access to the HDF5 files.
 * @param p A pair of lattice arguments and Hamiltonian arguments.
 */
template <class H, class L, class... LArgs, class... HArgs>
auto makeCore(Config mc, std::mutex&& mtx, std::pair<std::tuple<LArgs...>, std::tuple<HArgs...> >& p)
{
    return makeCore3<H, L>(mc, std::forward<std::mutex>(mtx), std::forward<decltype(p.first)>(p.first), p.second,
        std::make_index_sequence<std::tuple_size<typename std::remove_reference<std::tuple<HArgs...>>::type>::value>()
    );
}

/** Instantiate Core with a pre-allocated lattice.
 * 
 * @tparam H the type of Hamiltonian
 * @tparam L the type of the lattice
 * @tparam Args the other arguments.
 * 
 * @param latt the lattice
 * @param mc the MARQOVConfig object
 * @param mtx The mutex that synchronizes access to the HDF5 files.
 * @param args the remaining arguments.
 */
template <class H, class L, class ...Args>
auto makeCore2(std::true_type, L&& latt, Config&& mc, std::mutex&& mtx, Args&& ... args)
{
    //The first argument is a Lattice-like type -> from this we infer that 
    //We get a reference to sth. already allocated
    return Core<L, H, detail::Ref>(latt, mc, mtx, args...);
}

/** Instantiate the right core class.
 * 
 * This function figures out from the function arguments which is the right .
 * version of MARQOV::Core to take.
 * @tparam H the type of Hamiltonian
 * @tparam L the type of the lattice
 * @tparam Args the other arguments.
 * 
 * @param latt the lattice
 * @param mc the MARQOVConfig object
 * @param mtx The mutex that synchronizes access to the HDF5 files.
 * @param args the remaining arguments.
 */
template <class H, class L, class ...Args>
auto makeCore(L&& latt, Config&& mc, std::mutex&& mtx, Args&&... args)
{
    return makeCore2<H>(typename detail::is_Lattice<L>::type(), latt, std::forward<Config>(mc), std::forward<std::mutex>(mtx), args...);
}

}

#include "emcs.h"
#endif
