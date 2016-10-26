#include <dlib/matrix.h>
#include <fstream>
#include "simulator.cpp"


typedef struct
{
	bool        is_producer;
	std::size_t loc;
	bool control_type;
	double      WI;
	double      Pbh;
	double      rw;
	double      QINJ[3];
	double      Rso;
}                                            Well;
std::vector<Well>   mWells;
std::size_t Np, Nj, nw;

void
read_from_file(const char * filename, dlib::matrix<double> &v)
{
	double tmp;
	size_t n = 0;
	std::ifstream strm(filename);
	while(strm >> tmp)
	{
		n++;
	}
	v.set_size(n, 1);
	strm.close();
	strm.open(filename);
	for (unsigned i = 0; i < n; i++)
	{
		strm >> v(i, 0);
	}
}


std::size_t read_well()
{
	std::ifstream strm("input/well_inform.dat");

	Np = 0;
	Nj = 0;
	strm >> nw;
	mWells.resize(nw);
	for (unsigned i = 0; i < nw; i++)
	{
		strm >> mWells[i].is_producer;
		(mWells[i].is_producer) ? Np++ : Nj++;
		strm >> mWells[i].loc;
		strm >> mWells[i].rw;
		strm >> mWells[i].control_type;
		if (mWells[i].control_type)
		{
			strm >> mWells[i].Pbh;
		}
		else
		{
			strm >> mWells[i].QINJ[0];
			strm >> mWells[i].QINJ[1];
			strm >> mWells[i].QINJ[2];
		}
		mWells[i].Rso = 0;
	}
	return nw;
}

void F_D(simulator &model, double peturb_scale, dlib::matrix<double> &x)
{
	unsigned id = 15;
	double original, peturbed, peturbed_size;
	dlib::matrix<double> m, sens, x_p;
	read_from_file("input/m_pri.dat", m);
	dlib::matrix<double> grad(x.size(), 1);
	original = model.run(m, x, sens, 1);
	system("pause");
	std::cout << sens << std::endl;
	for (int i = 0; i < x.size(); i++)
	{
		x_p = x;
		if (i<x.size() / 2)
		{
			peturbed_size = peturb_scale*(x(i) - 3000);
			//peturbed_size = peturb_scale*x(i);
			x_p(i, 0) += peturbed_size;
		}
		else
		{
			std::cout << i << std::endl;
			peturbed_size = peturb_scale*(x(i) - 3000);
			//peturbed_size = peturb_scale*x(i);
			x_p(i, 0) += peturbed_size;
		}
		peturbed = model.run(m, x_p, sens, 0);
		grad(i, 0) = (peturbed - original) / peturbed_size;
		std::cout << "Gradient " << i << " is " << grad(i, 0) << std::endl;
	}
	dlib::matrix<double> error = (grad - sens);
	std::cout << "Adjoint" << std::endl;
	for (int i = 0; i < error.size(); i++)
	{
		std::cout << i << "  " << sens(i) << std::endl;
	}
	std::cout << "F_D " << std::endl;
	for (int i = 0; i < error.size(); i++)
	{
		std::cout << i << "  " << grad(i) << std::endl;
	}
	std::cout << "error " << std::endl;
	for (int i = 0; i < error.size(); i++)
	{
		error(i) = error(i) / (grad(i) + 10e-5);
		std::cout << i << "  " << error(i) * 100 << std::endl;
	}
	std::ofstream out("output/validation/FD.dat"), out1("output/validation/adjoint.dat"), out2("output/validation/err.dat");
	out << grad << std::endl;
	out1 << sens << std::endl;
	out2 << error << std::endl;
	std::cout << error << std::endl;
	system("pause");
}





int main()
{
	dlib::matrix<double> m, v,grad;
	read_from_file("input/m_pri.dat", m);
	read_well();
	simulator orig;
	read_from_file("input/wellschedule.dat",v);
	F_D(orig, 0.01, v);

}