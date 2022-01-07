/* This file is part of MARQOV:
 * A modern framework for classical spin models on general topologies
 * Copyright (C) 2021, The MARQOV Project
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
#ifndef MARQOV_DETAIL_H
#define MARQOV_DETAIL_H
#include <type_traits>
#include <tuple>
namespace MARQOV
{
    namespace detail
    {
        /** C++17 type sink.
         * 
         * A generic type sink from C++17.
         * It consumes a type, and makes it `void`.
         */
        template<class> struct type_sink {
            typedef void type;///< convert everything to void.
        };
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
        
        /** Internal base class if MARQOV gets the lattice by reference.
         * 
         * A base class that gets used if MARQOV::Core gets the lattice by reference.
         * @tparam L the lattice that will be used
         */
        template <class L>
        class Ref
        {
        public:
            /** The constructor that takes a reference.
             * 
             * @tparam Args the parameter pack with all parameters
             *
             * @param lattice A reference to a lattice.
             * @param args additional arguments that just get eaten.
             */
            template<class... Args>
            Ref(const L& lattice, Args&& ... args) : grid(lattice) {}
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
            /** Constructor if we construct the lattice ourselves.
             * 
             * @tparam Args The parameter pack of the lattice parameters.
             * 
             * @param args the actual lattice Arguments.
             */
            template <class ...Args>
            NonRef(std::tuple<Args...> args ) : NonRef(args, 
                                                         std::make_index_sequence<std::tuple_size<
                                                         std::tuple<Args...>
                                                         >::value>()) {}
                                                         
            /** Constructor to do the parameter unpacking.
            * 
            * @tparam Args The parameter pack of the lattice parameters.
            * @tparam S an integer sequence.
            * 
            * @param args the actual lattice Arguments.
            */
            template <class ...Args, size_t... S>
            NonRef(std::tuple<Args...> args, std::index_sequence<S...>) : grid(std::get<S>(args)... ) {}

            const L grid;///< The storage of the Lattice. Note that the name of the variable is the same as in detail::Ref.
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
        
        /** Implementation function to call an objects' member function with a tuple of arguments.
         * 
         * This call actually performs the function call.
         * @tparam Function Which member function to call.
         * @tparam Object The type of the object
         * @tparam Tuple The type of the Tuple of arguments.
         * @tparam I an index sequence
         * 
         * @param f the name of the function.
         * @param obj the object.
         * @param t the arguments to f.
         * @returns the return value of the function call.
         */
        template<typename Function, typename Object, typename Tuple, size_t ... I>
        inline auto _call(Function f, Object& obj, Tuple t, std::index_sequence<I ...>) 
        {
            return (obj.*f)(std::get<I>(t) ...);
        }
        
        /** A helper function to call an objects' member function with a tuple of arguments.
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
        inline auto _call(Function f, Object& obj, Tuple t) 
        {
            static constexpr auto size = std::tuple_size<Tuple>::value;
            return _call(f, obj, t, std::make_index_sequence<size>{});
        }
    };
};

#endif
