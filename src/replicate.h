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

#ifndef REPLICATE_H
#define REPLICATE_H
#include <vector>
#include <tuple>
#include "core.h"
#include "helpers.h"


template <class Params>
inline std::vector<Params> replicator(const std::vector<Params>& params, int nreplicas, int sortmode=0)
{
	std::vector<Params> newparams;

	int mcid = 0;
	for(std::size_t i=0; i<params.size(); ++i)
	{
		auto&& l = std::get<0>(params[i]);
		auto mp = std::get<1>(params[i]);
		auto&& hp = std::get<2>(params[i]);

		for (int j = 0; j < nreplicas; ++j)
		{
			auto mpr(mp);
			mpr.setrepid(j);
			mpr.setid(mcid++);
			newparams.emplace_back(l, mpr, hp);
		}
	}

	return newparams;
}

template <class L, class HArgs>
inline auto finalize_parameter(L lp, const MARQOV::Config& mp, const std::vector<HArgs>& hp)
{
	typedef std::tuple<L, MARQOV::Config, HArgs> RetType;
	std::vector<RetType> params;

	for(std::size_t i=0; i < hp.size(); ++i) 
	{
		params.emplace_back(lp, mp, hp[i]);
	}

	return params;
}
/*
template <class LArgs, class MArgs, class HArgs>
inline auto finalize_parameter_triple(const LArgs& lp, const MArgs& mp, const std::vector<HArgs>& hp)
{
	std::vector<Triple<LArgs, MArgs, HArgs>> params;

	for(std::size_t i=0; i<hp.size(); ++i) 
	{
		params.push_back(make_triple(lp, mp, hp[i]));
	}

	return params;
}


template <class MArgs, class HArgs>
inline auto finalize_parameter_pair(const MArgs& mp, const std::vector<HArgs>& hp)
{
	std::vector<std::pair<MArgs, HArgs>> params;

	for(std::size_t i=0; i<hp.size(); ++i) 
	{
		params.push_back(std::make_pair(mp, hp[i]));
	}

	return params;
}*/


#endif
