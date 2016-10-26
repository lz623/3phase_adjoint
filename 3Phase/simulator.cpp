#include "DiscreteProblem.hpp"
#include "Intel_Pardiso.hpp"
#include "NewtonSolver.hpp"
#include "R_Precond_Solver.hpp"
#include "Intel_Prec_GMRES.hpp"
#include "Intel_ILU0.hpp"
#include "fastl/containers/pod_vector_unbounded.hpp"
#include <fstream>
#include <dlib/matrix.h>

//typedef GENSOL::R_Precond_Solver< GENSOL::Intel_Prec_GMRES< GENSOL::Intel_ILU0 > >  LINEARSOLVER;
typedef GENSOL::Intel_Pardiso                       LINEARSOLVER;
//typedef GENSOL::Intel_Prec_GMRES< GENSOL::Intel_ILU0 >  LINEARSOLVER;
typedef GENSOL::NewtonSolver< DiscreteProblem, LINEARSOLVER > STDN;

class simulator
{
	typedef DiscreteProblem::wellrate   w_rate;
	typedef DiscreteProblem::rate		 rate;
private:
	const std::size_t MAX_NLNITER = 15;
	const double      DT_INIT =5;
	const double      DT_CUT = 0.5;
	const double      DT_GROW = 2.0;
	std::size_t NX, NY, NZ;
	double      LX, LY, LZ;
	std::vector<double> times, pro_time, pro_dt;
	double Qo, Qw;
	dlib::matrix<double> m;
	std::vector<unsigned> ob_si;
	std::vector<double> vPORO;
	std::vector<double> vKX;
	std::vector<double> vKY;
	std::vector<double> vKZ;
	unsigned N_tc;
	unsigned Max_nz;
	unsigned Num_var;

	typedef struct
	{
		dlib::matrix<double> BHP;
		dlib::matrix<double> qt;
		dlib::matrix<double> qw;
		dlib::matrix<double> qo;
		dlib::matrix<double> wct;
		dlib::matrix<double> qg;
		dlib::matrix<double> gor;
	}CS_gradient;


public:

	void
	read_from_file(const char * filename, std::vector<double> &v)
	{
		double tmp;
		std::ifstream strm(filename);
		if (strm.good())
		{
			for (std::size_t i = 0; i < v.size(); ++i)
			if (strm >> tmp) v[i] = tmp;
		}
	}

	void
		dump_solution(const char * filename, const DiscreteProblem::StateVector & v)
	{
			std::ofstream strm(filename);
			for (std::size_t i = 0; i < v.size(); ++i)
				strm << v[i].Po.value() << "\t"
				<< v[i].Sw.value() << "\t"
				<< v[i].Sg.value() << "\t"
				<< 1.0 - v[i].Sw.value() - v[i].Sg.value() << "\t"
				<< v[i].Rso.value() << std::endl;
			strm.close();
		}

	void
		dump_field(const char * filename,
		const std::vector<double> &phi,
		const std::vector<double> &kx,
		const std::vector<double> &ky,
		const std::vector<double> &kz)
	{
			std::ofstream strm(filename);
			for (std::size_t i = 0; i < phi.size(); ++i)
				strm << phi[i] << "\t"
				<< kx[i] << "\t"
				<< ky[i] << "\t"
				<< kz[i] << std::endl;
			strm.close();
		}
	void dump_rate(std::string  filename,
		const std::vector<std::vector<rate>> well_rate)
	{
		std::string filenameo = filename + "_oil.dat";
		std::string filenameg = filename + "_gas.dat";
		std::string filenamew = filename + "_water.dat";
		std::string filenamep = filename + "_BHP.dat";
		std::ofstream strmo(filenameo);
		std::ofstream strmw(filenamew);
		std::ofstream strmg(filenameg);
		std::ofstream strmp(filenamep);
		for (unsigned i = 0; i < well_rate.size(); i++)
		{
			for (unsigned w = 0; w < well_rate[0].size(); w++)
			{

				strmo << well_rate[i][w].qo <<"  ";
				strmw << well_rate[i][w].qw << "  ";
				strmg << well_rate[i][w].qg/178.11025 << "  ";
				strmp << well_rate[i][w].BHP << "  ";
			}
			strmo << std::endl;
			strmw << std::endl;
			strmg << std::endl;
			strmp << std::endl;
		}
	}


	void read_reservoir()
	{
		std::ifstream strm("./input/reservoir.dat");
		strm >> NX;
		strm >> NY;
		strm >> NZ;
		strm >> LX;
		strm >> LY;
		strm >> LZ;
	}
	void read_schedule(std::vector<double> &times)
	{
		double tmp;
		std::ifstream strm("./input/schedule.dat");
		while (strm >> tmp)
		{
			times.push_back(tmp);
		}
	}

	void m_to_vector(dlib::matrix<double> &si,
		std::vector<double> &phi,
		std::vector<double> &kx,
		std::vector<double> &ky,
		std::vector<double> &kz)
	{
		size_t n = si.size();
		size_t nc = n / 4;
		phi.resize(NX*NY*NZ);
		kx.resize(NX*NY*NZ);
		ky.resize(NX*NY*NZ);
		kz.resize(NX*NY*NZ);
		for (unsigned i = 0; i < nc; i++)
		{
			kx[i] = exp(si(i, 0));
			ky[i] = exp(si(i, 0));
			kz[i] = exp(si(i, 0));
			phi[i] = si(i + 3 * nc, 0);
		}
	}
	double NPV(std::vector<bool> isproducer, std::vector<double> &pro_time, std::vector<double> &pro_dt,
		std::vector<std::vector<rate>> &well_prod, std::vector<std::vector<CS_gradient>> &all_grad, dlib::matrix<double> &NPV_grad)
	{
		double ro, rw, rinj, b,rg;
		std::ifstream input("./input/economic.dat");
		input >> ro;
		input >> rw;
		input >> rinj;
		input >> rg;
		input >> b;
		std::ofstream output("./output/NPV.out");
		double total = 0;
		NPV_grad.set_size(Num_var, 1);
		//NPV_grad = all_grad[29][0].qg;
		NPV_grad = 0;
		unsigned nw = well_prod[0].size();
		for (unsigned t = 0; t < pro_time.size(); t++)
		{
			double tmp = 0;
			dlib::matrix<double> t_grad(Num_var, 0);
			for (unsigned w = 0; w < nw; w++)
			{
				if (isproducer[w])
				{
					tmp += (ro*well_prod[t][w].qo - rw*well_prod[t][w].qw + rg*well_prod[t][w].qg / 178.11025);
					t_grad += (ro*all_grad[t][w].qo - rw*all_grad[t][w].qw + rg*all_grad[t][w].qg / 178.11025);
				}
				else
				{
					tmp -= (rinj*well_prod[t][w].qw);
					t_grad -= (rinj*all_grad[t][w].qw);
				}
			}
			total += tmp*pro_dt[t] / pow((1 + b), pro_time[t] / 365.0);
			NPV_grad += t_grad*pro_dt[t] / pow((1 + b), pro_time[t] / 365.0);
		}
		output << total;
		output.close();
		return total;
	}


	std::vector<double> transpose_prod(DiscreteProblem::Sparse_matrix &J_m, std::vector<double> &ajv, unsigned nwc, unsigned nrh)
	{
		std::vector<double> result(nwc*nrh, 0);
		unsigned ny = ajv.size() / nrh;
		for (unsigned j = 0; j <nrh; j++)
		{
			for (unsigned i = 0; i < J_m.val.size(); i++)
			{
				//std::cout << J_m.col[i] + j*nwc << "  " << J_m.row[i] + j*nwc <<"  "<< J_m.val[i] << std::endl;
				result[J_m.col[i] + j*nwc] += J_m.val[i] * ajv[J_m.row[i] + j*ny];
			}
		}

		return result;
	}

	void compute_well_cs(const std::vector<DiscreteProblem::CS> &cons,
		std::vector<unsigned> &w_cs, unsigned &nrh)
	{
		nrh = 0;
		w_cs.resize(cons.size());
		//count the number of constrain for each well
		for (unsigned w = 0; w < cons.size(); w++)
		{
			if (cons[w].BHP_i)
			{
				w_cs[w]++;
				nrh++;
			}
			if (cons[w].qt_i)
			{
				w_cs[w]++;
				nrh++;
			}
			if (cons[w].qo_i)
			{
				w_cs[w]++;
				nrh++;
			}
			if (cons[w].qw_i)
			{
				w_cs[w]++;
				nrh++;
			}
			if (cons[w].qg_i)
			{
				w_cs[w]++;
				nrh++;
			}
			if (cons[w].wct_i)
			{
				w_cs[w]++;
				nrh++;
			}
			if (cons[w].gor_i)
			{
				w_cs[w]++;
				nrh++;
			}
		}
	}
	void conglo_v(std::vector<double> &ob, const std::vector<DiscreteProblem::CS> &cons,
		const unsigned nv, const unsigned &nrh)
	{
		unsigned id = 0;
		double temp;
		//////////////////////////////////////////////////////////////////////////////
		ob.resize(nrh*nv);
		ob.assign(nrh*nv,0);
		for (unsigned w = 0; w < cons.size(); w++)
		{
			if (cons[w].BHP_i)
			{
				std::vector<double> x(nv, 0);
				cons[w].BHP.extract(temp, x.data(), nv);
				for (unsigned i = 0; i < nv; i++)	ob[i + id*nv] = x[i];
				id++;
			}
			if (cons[w].qt_i)
			{
				std::vector<double> x(nv, 0);
				cons[w].qt.extract(temp, x.data(), nv);
				for (unsigned i = 0; i < nv; i++)	ob[i + id*nv] = x[i];
				id++;
			}
			if (cons[w].qo_i)
			{
				std::vector<double> x(nv, 0);
				cons[w].qo.extract(temp, x.data(), nv);
				for (unsigned i = 0; i < nv; i++)	ob[i + id*nv] = x[i];
				id++;
			}
			if (cons[w].qw_i)
			{
				std::vector<double> x(nv, 0);
				cons[w].qw.extract(temp, x.data(), nv);
				for (unsigned i = 0; i < nv; i++)	ob[i + id*nv] = x[i];
				id++;
			}
			if (cons[w].qg_i)
			{
				std::vector<double> x(nv, 0);
				cons[w].qg.extract(temp, x.data(), nv);
				for (unsigned i = 0; i < nv; i++)	ob[i + id*nv] = x[i];
				id++;
			}
			if (cons[w].wct_i)
			{
				std::vector<double> x(nv, 0);
				cons[w].wct.extract(temp, x.data(), nv);
				for (unsigned i = 0; i < nv; i++)	ob[i + id*nv] = x[i];
				id++;
			}
			if (cons[w].wct_i)
			{
				std::vector<double> x(nv, 0);
				cons[w].qg.extract(temp, x.data(), nv);
				for (unsigned i = 0; i < nv; i++)	ob[i + id*nv] = x[i];
				id++;
			}
		}
	}
	void assign_gradient(std::vector<DiscreteProblem::CS> &cons,
		std::vector<CS_gradient> &grad, std::vector<double> &sum, unsigned nv)
	{
		unsigned id = 0;
		grad.resize(cons.size());
		for (unsigned w = 0; w < cons.size(); w++)
		{
			if (cons[w].BHP_i)
			{
				grad[w].BHP.set_size(nv, 1);
				grad[w].BHP = 0;
				for (unsigned i = 0; i < nv; i++) grad[w].BHP(i, 0) = sum[i + id*nv];
				id++;
			}
			if (cons[w].qt_i)
			{
				grad[w].qt.set_size(nv, 1);
				grad[w].qt = 0;
				for (unsigned i = 0; i < nv; i++) grad[w].qt(i, 0) = sum[i + id*nv];
				id++;
			}
			if (cons[w].qo_i)
			{
				grad[w].qo.set_size(nv, 1);
				grad[w].qo = 0;
				for (unsigned i = 0; i < nv; i++) grad[w].qo(i, 0) = sum[i + id*nv];
				id++;
			}
			if (cons[w].qw_i)
			{
				grad[w].qw.set_size(nv, 1);
				grad[w].qw = 0;
				for (unsigned i = 0; i < nv; i++) grad[w].qw(i, 0) = sum[i + id*nv];
				id++;
			}
			if (cons[w].qg_i)
			{
				grad[w].qg.set_size(nv, 1);
				grad[w].qg = 0;
				for (unsigned i = 0; i < nv; i++) grad[w].qg(i, 0) = sum[i + id*nv];
				id++;
			}
			if (cons[w].wct_i)
			{
				grad[w].wct.set_size(nv, 1);
				grad[w].wct = 0;
				for (unsigned i = 0; i < nv; i++) grad[w].wct(i, 0) = sum[i + id*nv];
				id++;
			}
			if (cons[w].gor_i)
			{
				grad[w].gor.set_size(nv, 1);
				grad[w].gor = 0;
				for (unsigned i = 0; i < nv; i++) grad[w].gor(i, 0) = sum[i + id*nv];
				id++;
			}
		}
	}
	void Backward_CS_gradient(std::vector<unsigned> & time,
		std::vector<CSR> & al_j,
		std::vector<DiscreteProblem::Sparse_matrix> & ol_j,
		std::vector<DiscreteProblem::Sparse_matrix> & al_jm,
		std::vector<DiscreteProblem::CS> &cons,
		std::vector<DiscreteProblem::CS> &cons_m,
		std::vector<CS_gradient> &grad)
	{
		unsigned nrh;
		unsigned nos = time.back() - 1;
		std::vector<unsigned> w_cs;
		compute_well_cs(cons, w_cs, nrh);
		LINEARSOLVER lnsolve(N_tc, Max_nz, nrh);
		std::size_t n = al_j[0].N(), nv = Num_var;
		std::vector<double> ajv(n*nrh, 0), tmp(nv*nrh), tmp1(n*nrh), sum(nv*nrh, 0);
		std::vector<double> oby, obm;
		// frist time to compute ajoint vector backward;

		conglo_v(oby, cons, n, nrh);
		conglo_v(obm, cons_m, nv, nrh);
		//frist run;
		lnsolve.solvet(al_j[nos], ajv, oby);
		tmp = transpose_prod(al_jm[nos], ajv, nv, nrh);
		////////////////////////// can simplify//////////////////////////////////////
		for (unsigned i = 0; i < nv*nrh; i++)
		{
			sum[i] += tmp[i];
		}
		//////////////////////////////////////////////////////////////////////////////
		//compute ajoint vector for each simulation time
		for (unsigned j = nos; j >0; j--)
		{
			tmp1 = transpose_prod(ol_j[j - 1], ajv, n, nrh);
			lnsolve.solvet(al_j[j - 1], ajv, tmp1);
			tmp = transpose_prod(al_jm[j - 1], ajv, nv, nrh);
			for (unsigned i = 0; i < nv*nrh; i++)
			{
				sum[i] += tmp[i];
			}
		}
		//summation of ajoint vecotor product of Fm match at time
		//////////////////////////////////////////////////////////////////////////////
		for (unsigned i = 0; i < nv*nrh; i++)
		{
			sum[i] -= obm[i];
		}
		//////////////////////////////////////////////////////////////////////////////
		assign_gradient(cons, grad, sum, nv);
	}
	simulator()
	{
		read_reservoir();
		read_schedule(times);
	}


	double run(dlib::matrix<double> &m, dlib::matrix<double> &v, dlib::matrix<double> &NPV_grad,bool output_i)
	{
		double NPV_v=0;
		std::vector<double> times;
		read_reservoir();
		read_schedule(times);
		m_to_vector(m, vPORO, vKX, vKY, vKZ);
		std::ofstream out("./output/timestep.dat");
		const double OPT_NLNITER = (3 < MAX_NLNITER ? MAX_NLNITER : 3);
		DiscreteProblem model(NX, NY, NZ, LX, LY, LZ, vKX, vKY, vKZ, vPORO);
		DiscreteProblem model_m(NX, NY, NZ, LX, LY, LZ, vKX, vKY, vKZ, vPORO);
		model.loadwell_data(times, v);
		model_m.loadwell_independent(times, v);
		LINEARSOLVER lnsolver(model.N_tc, model.max_num_nnz());
		STDN newton(model, MAX_NLNITER, 1);
		DiscreteProblem::StateVector uOld, uNew, aOld, aNew;
		model.initialize_state(uOld);
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		Num_var = v.nr();
		std::vector<std::vector<rate>> well_prod;
		CSR		Jaco(model.N_tc, model.max_num_nnz());
		DiscreteProblem::Sparse_matrix		o_Jaco;
		std::vector<DiscreteProblem::Sparse_matrix> F_v;
		std::vector<DiscreteProblem::Sparse_matrix> F_yold;
		std::vector<CSR> F_y;
		DiscreteProblem::Sparse_matrix	j_v;
		std::vector<std::vector<CS_gradient>> all_grad;
		std::vector<CS_gradient> grad;
		aOld.resize(uOld.size());
		aNew.resize(uOld.size());
		N_tc = model.N_tc;
		Max_nz = model.max_num_nnz();
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		double DT = DT_INIT;
		double time = 0.0;
		unsigned j = 0, ob_t = 0;
		for (std::size_t i = 0; i < times.size(); i++)
		{
			while (time < times[i])
			{
				model.sum_t = time;
				model_m.sum_t = time;
				uNew = uOld;
				STDN::report_t stdsmry = newton.solve_timestep_new(uNew, uOld, DT, model, lnsolver, Jaco, o_Jaco);
				if (stdsmry.is_converged)
				{
					model_m.initialize_state(uOld, uNew, aOld, aNew);// assign state to model_m;
					model_m.extract_F_der(aOld, aNew, DT, j_v);
					F_v.push_back(j_v);
					F_y.push_back(Jaco);
					F_yold.push_back(o_Jaco);
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					uOld = uNew;
					time += DT;
					if (stdsmry.niter < OPT_NLNITER) DT *= DT_GROW;
					DT = std::min(DT, (times[i] - time));
					DT == 0 ? DT = DT_INIT : DT = DT;
					out << DT << std::endl;
					j++;
				}
				else
				{
					DT *= DT_CUT;
					std::cout << "FAILED  @ time "<<model.sum_t<<" and current time step is: " << DT<<std::endl;
				}
			}
			model.extract_obj_der(i);	//extract derivetive for objective function
			model_m.extract_obj_der(i);
			pro_time.push_back(time);
			well_prod.push_back(model.well_rate);
			pro_dt.push_back(DT);
			ob_si.push_back(j);
			Backward_CS_gradient(ob_si, F_y, F_yold, F_v, model.H_constrain[ob_t], model_m.H_constrain[ob_t], grad);
			ob_t++;
			all_grad.push_back(grad);
		}
		NPV_v = NPV(model.is_producer, pro_time, pro_dt, well_prod, all_grad, NPV_grad);
		if (output_i)
		{
			dump_solution("./output/results.dat", uNew);
			dump_rate("./output/well_rate", well_prod);
		}
		//NPV_v = NPV(model.is_producer, pro_time, pro_dt, well_prod, all_grad, NPV_grad);
		pro_time.clear();
		pro_dt.clear();
		ob_si.clear();
		return NPV_v;
	}




};