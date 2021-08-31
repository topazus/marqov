
/** Find out if a string starts with something.
 * @param longword we search in this string
 * @param shortword we look for this
 * @return truen if longword strarts with shortword, else false
 */
bool startswith(const std::string& longword, const std::string& shortword) noexcept
{
    return longword.find(shortword) == 0;
}

/** Check the vailidity of the replica configuration.
 * throws if invalid.
 * 
 * @param nr number of replicas
 * @param nL amount of lattice simulations
 */
void checkreplicaconfig(int nr, int nL)
{
    if ((nr != nL) && (nr != 1)) throw std::invalid_argument("[MARQOV] Invalid replica configuration!");
}

/** Removes previous simulations.
 * 
 * @param outbasedir The folder that we remove entirely
 */
void delete_directory(const std::string& outbasedir)
{
	// delete previous output // fixme: don't do that by default!
	#ifdef MPIMARQOV
	    int myrank = 0;
	    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	    if(myrank == 0)
	    {
	#endif
        std::cout<<"[MARQOV::main] Erasing previous data!!!!!!!!!!!!!!!"<<std::endl;
        std::string command = "rm -r " + outbasedir;
        system(command.c_str());
        makeDir(outbasedir);
	#ifdef MPIMARQOV
	    }
	#endif
}

/** Print some nice information
 * 
 * @param registry Where to get the data from.
 * @param ham the hamiltonian
 */
void printInfoandcheckreplicaconfig(RegistryDB& registry, const std::string& ham)
{
    const auto dim 	      = registry.Get<int>(ham+".ini", ham, "dim" );
	const auto nreplicas  = registry.Get<std::vector<int>>(ham+".ini", ham, "rep" );
	const auto nreplicas_str = registry.Get<std::string>(ham+".ini", ham, "rep" );
	const auto nL  	      = registry.Get<std::vector<int>>(ham+".ini", ham, "L" );
	const auto nLs 	      = registry.Get<std::string>(ham+".ini", ham, "L" );

	cout << endl;
	cout << "Hamiltonian: \t" << ham << endl;
	cout << "Dimension: \t" << dim << endl;
	cout << "Lattice sizes:\t" << nLs << endl;
	cout << "Replicas:\t" << nreplicas_str << endl;
    
    checkreplicaconfig(nreplicas.size(), nL.size());
}

/** Utility function to write the common part of the config file
 * @param os the stream associated with the new config file.
 * @param nmetro 
 * @param nclusteramp
 * @param nclusterexp
 * @param warmupsteps
 * @param measuresteps
 */
void createcfgfooter(std::ostream& os, int nmetro, double nclusteramp, int nclusterexp, int warmupsteps, int measuresteps)
{
    os<<"[MC]\n"<<"nmetro = "<<nmetro<<"\nnclusteramp = "<<nclusteramp<<"\nnclusterexp = "<<nclusterexp<<"\nwarmupsteps = "<<warmupsteps<<"\nmeasuresteps = "<<measuresteps<<"\n";
    os<<"[IO]\n"<<"outdir = ./out\n"<<"[END]"<<std::endl;
}

