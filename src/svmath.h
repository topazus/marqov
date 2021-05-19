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

#ifndef SVMATH_H
#define SVMATH_H
#include <tuple>
#include <cmath>
#include <numeric>
#include <array>

// ------- elementary state vector calculus

/** Add two StateVectors elementwise: c=a+b.
 * 
 * @param arg1 a
 * @param arg2 b
 * @returns @f[ \vec{c} = \vec{a} + \vec{b} @f]
 */
template <class T, std::size_t N>
inline std::array<T, N> operator + (std::array<T, N> arg1, std::array<T, N> arg2)
{
    std::array<T, N> res(arg1);
    for(std::size_t i = 0; i < N; ++i)
    res[i] += arg2[i];
    return res;
}

/** Subtract two StateVectors elementwise: c=a-b.
 * 
 * @param arg1 a
 * @param arg2 b
* @returns @f[ \vec{c} = \vec{a} - \vec{b} @f]
 */
template <class T, std::size_t N>
inline std::array<T, N> operator - (std::array<T, N> arg1, std::array<T, N> arg2)
{
    std::array<T, N> res(arg1);
    for(std::size_t i = 0; i < N; ++i)
    res[i] -= arg2[i];
    return res;
}

/** mult is a component-wise multiplication.
 * 
 * This is used e.g. in the Metropolis algorithm.
 * @see Metropolis
 * @param a
 * @param b
 * @returns @f[ c = a * b @f]
 */
inline double mult(const double& a, const double& b)
{
    return a*b;
}

/** mult is a component-wise multiplication.
 * 
 * This is used e.g. in the Metropolis algorithm.
 * @see Metropolis
 * @param a
 * @param b
 * @returns @f[ c = a * b @f]
 */
inline double mult(const double& a, const int& b)
{
    return a*double(b);
}

/** mult is a component-wise multiplication.
 * 
 * This is used e.g. in the Metropolis algorithm.
 * @see Metropolis
 * @param a
 * @param b
 * @returns @f[ c_i = a_i * b_i @f]
 */
template <class VecType, class T, std::size_t N>
inline std::array<T, N> mult(const VecType& a, const std::array<T, N>& b)
{
    std::array<T, N> retval(b);
    for(int i = 0; i < N; ++i)
    retval[i] *= a[i];
    return retval;
}

/** mult is a component-wise multiplication.
 * 
 * This is used e.g. in the Metropolis algorithm.
 * @see Metropolis
 * @param a
 * @param b
 * @returns @f[ c_i = a * b_i @f]
 */
template <class T, std::size_t N>
inline std::array<T, N> mult(const int& a, const std::array<T, N>& b)
{
	std::array<T, N> retval(b);
	for(std::size_t i = 0; i < N; ++i) retval[i] *= a;
	return retval;
}

/** mult is a component-wise multiplication.
 * 
 * This is used e.g. in the Metropolis algorithm.
 * @see Metropolis
 * @param a
 * @param b
 * @returns @f[ c_i = a * b_i @f]
 */
template <class T, std::size_t N>
inline std::array<T, N> mult(const double& a, const std::array<T, N>& b)
{
	std::array<T, N> retval(b);
	for(std::size_t i = 0; i < N; ++i) retval[i] *= a;
	return retval;
}

/** The dot/inner product.
 * 
 * @param a
 * @param b
 * @returns @f[ c=a*b @f]
 */
inline double dot(const double& a, const double& b)
{
    return a*b;
}

/** The dot/inner product.
 * 
 * @param a
 * @param b
 * @returns @f[ c=\sum_i a_i b_i @f]
 */
template<class T, std::size_t N>
inline T dot(const std::array<T, N>& a, const std::array<T, N>& b)
{
    return std::inner_product(begin(a), end(a), begin(b), 0.0);
}

/** Reflect a StateVector along a plane.
 * 
 * This takes a State vector and a plane defined by its normal vector
 * and performs a reflection along it
 * 
 * @tparam StateVector the Container used for storing the vectors.
 * @param vec the state vector that should be reflected.
 * @param mirror The normal vector of the plane used for reflecting.
 * @returns the modified, reflected vector.
 */
template <class StateVector>
inline void reflect(StateVector& vec, const StateVector mirror)
{
	const int SymD = std::tuple_size<StateVector>::value;
	
	const double dotp = dot(vec, mirror);

	for (int i=0; i<SymD; i++) vec[i] -= 2*dotp*mirror[i];
}

/** Normalize vector
 * 
 * This normalizes a vector such that it has unit length.
 * 
 * @tparam Container The Container used for storing the vector.
 * @param a the vector to normalize
 */
template <class Container>
inline void normalize(Container& a)
{
	typename Container::value_type tmp_abs=std::sqrt(dot(a, a));

	for (decltype(a.size()) i = 0; i < a.size(); ++i) a[i] /= tmp_abs;
}

/** Dump array to stream.
 * 
 * @tparam T the type of the array elements.
 * @tparam N the length of the array.
 * 
 * @param os the ouput stream
 * @param vec the array that we want to write to the screen.
 * @return A stream that now contains a textual representation of the array.
 */
template <typename T, std::size_t N>
std::ostream& operator<<(std::ostream& os, const std::array<T, N>& vec)
{
    for (std::size_t i = 0; i < N; i++) os << vec[i] << "\t";
    return os;
}

#endif
