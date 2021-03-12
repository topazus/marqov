#ifndef FILTERS_H
#define FILTERS_H

	// -------------------- filters --------------------

	// filter to determine output file path and name
	// the filter _must_ set the outname

	auto defaultfilter = [](auto& latt, auto p)
	{
		auto& mp = p.first;		// Monte Carlo params
		auto& hp = p.second;	// Hamiltonian params
		
		std::string str_repid = std::to_string(mp.repid);
		std::string str_beta  = "beta"+std::to_string(std::get<0>(hp));
		mp.outname = str_beta+"_"+str_repid;
		return std::tuple_cat(std::forward_as_tuple(latt), p);
	};

	auto defaultfilter_triple = [](auto p)
	{
       	auto& lp = p.first;
       	auto& mp = p.second;
       	auto& hp = p.third;

 		auto str_repid = std::to_string(mp.repid);
		auto str_beta  = "beta"+std::to_string(std::get<0>(hp));
		auto str_L     = std::to_string(std::get<0>(lp));

		mp.outname = str_beta+"_"+str_repid;

		return p;
	};

	auto xxzfilter = [](auto& latt, auto p)
	{	
		auto& mp = p.first;		// Monte Carlo params
		auto& hp = p.second;	// Hamiltonian params
	
		std::string str_repid = std::to_string(mp.repid);
		std::string str_beta  = "beta"+std::to_string(std::get<0>(hp));
		std::string str_extf  = "extf"+std::to_string(std::get<1>(hp));
		mp.outname = str_beta+"_"+str_extf+"_"+str_repid;

		return std::tuple_cat(std::forward_as_tuple(latt), p);
	};


	auto sshfilter = [](auto& latt, auto p)
	{
		auto& mp = p.first;		// Monte Carlo params
		auto& hp = p.second;	// Hamiltonian params

		std::string str_repid = std::to_string(mp.repid);
		std::string str_k     = "k"+std::to_string(std::get<2>(hp));
		std::string str_g     = "g"+std::to_string(std::get<3>(hp));
		std::string str_mu    = "mu"+std::to_string(std::get<4>(hp));

		double dtau = std::get<5>(hp);
		int Ltime = std::get<6>(hp); 

		std::string str_Ltime = "Ltime"+std::to_string(Ltime);
		std::string str_dtau  = "dtau"+std::to_string(dtau); 

		mp.outname = str_Ltime+"_"+str_k+"_"+str_g+"_"+str_mu+"_"+str_dtau+"_"+str_repid;
        auto p3 = std::make_pair(mp, std::tuple_cat(hp, std::forward_as_tuple(latt)));
		return std::make_tuple(latt, p3 );

	};

#endif
