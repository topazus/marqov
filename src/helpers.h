std::vector<double> create_range(double rangestart, double rangefinal, int steps, std::string type="linear", int endpoint=0)
{
	std::vector<double> range(steps);

	if (type == "linear")
	{
		double rangestep = (rangefinal-rangestart)/double(steps-endpoint);

		for (int i=0; i<steps; i++) range[i] = rangestart + i*rangestep;
	}
	else
	{
		// implement me
		std::cout << "not implemented!" << std::endl;
	}

	return range;
}


