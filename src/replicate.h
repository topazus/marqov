#ifndef REPLICATE_H
#define REPLICATE_H

#include "helpers.h" 

template <class Params>
std::vector<Params> replicator(std::vector<Params>& params, int nreplicas, int sortmode=0)
{
	std::vector<Params> newparams;

	int mcid = 0;
	for(std::size_t i=0; i<params.size(); ++i)
	{
		auto lp = params[i].first;
		auto mp = params[i].second;
		auto hp = params[i].third;

		for (std::size_t j=0; j<nreplicas; ++j)
		{
			auto mpr(mp);
			mpr.setrepid(j);
			mpr.setid(mcid++);
			auto newparam = make_triple(lp, mpr, hp);
			newparams.push_back(newparam);
		}
	}

	return newparams;
}


template <class Params>
std::vector<Params> replicator_pair(std::vector<Params>& params, int nreplicas, int sortmode=0)
{
	std::vector<Params> newparams;

	int mcid = 0;
	for(std::size_t i=0; i<params.size(); ++i)
	{
		auto mp = params[i].first;
		auto hp = params[i].second;

		for (std::size_t j=0; j<nreplicas; ++j)
		{
			auto mpr(mp);
			mpr.setrepid(j);
			mpr.setid(mcid++);
			auto newparam = std::make_pair(mpr, hp);
			newparams.push_back(newparam);
		}
	}

	return newparams;
}



template <class LArgs, class MArgs, class HArgs>
auto finalize_parameter_triple(const LArgs& lp, const MArgs& mp, const std::vector<HArgs>& hp)
{
	std::vector<Triple<LArgs, MArgs, HArgs>> params;

	for(std::size_t i=0; i<hp.size(); ++i) 
	{
		params.push_back(make_triple(lp, mp, hp[i]));
	}

	return params;
}


template <class MArgs, class HArgs>
auto finalize_parameter_pair(const MArgs& mp, const std::vector<HArgs>& hp)
{
	std::vector<std::pair<MArgs, HArgs>> params;

	for(std::size_t i=0; i<hp.size(); ++i) 
	{
		params.push_back(std::make_pair(mp, hp[i]));
	}

	return params;
}


#endif
