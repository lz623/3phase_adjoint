#ifndef __EDFM_HPP_INCLUDED_
#define __EDFM_HPP_INCLUDED_
#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <dlib/matrix.h>


using namespace std;





//vector<fractureinf> n_fracture;
//vector<fractureinf> h_fracture;
//vector<grid> grids;
//vector<frac> n_frac;
//vector<frac> h_frac;
////vector<Hor_Well> h_well;
//intersection inter_nn;
//unsigned nw;
//unsigned nh, nn;
//unsigned Nx, Ny,Nz;
//double Lx, Ly,Lz,Dx,Dy,Dz;
//unsigned frac_index = 0;

class EDFM
{
public:
	typedef struct point
	{
		double x;
		double y;
	};

	static inline
		void eq(point &_v, point &_w)
	{
			_w.x = _v.x;
			_w.y = _v.y;
		};

	static inline
		void sw(point &A, point &B)
	{
			point tmp;
			eq(A, tmp);
			eq(B, A);
			eq(tmp, B);
		};

	struct grid
	{
		double perm;
		double poro;
	};

	typedef struct
	{

		double rw;
		point start, end;
		std::vector<double> schedule;
		std::vector<double> Vol_std;
		std::vector<size_t> fracindex;
		std::vector<size_t> gridindex;
		std::vector<double> WI;
		bool isproducer;  //Well type,(0=Pro,1=inj)
		bool isPcontrol;
	} Hor_Well;


	struct frac
	{
		point start, end, coor;
		double perm, aper, por, d, l, trans, A;
		unsigned index, f_index, g_index;
	};



	struct intersection
	{
		std::vector<point> coor;
		std::vector<unsigned> i_index, j_index;
		std::vector<double> perm_i, perm_j, A_i, A_j, d_i, d_j;
	};


	struct fractureinf
	{
		point start, end;
		double perm, aper, por;
		unsigned f_index, upper, lower;
	};


	//typedef struct
	//{
	//	vector<size_t> fracindex;
	//	vector<size_t> gridindex;
	//	vector<double> schedule;
	//	adetl::ADvector Pw;
	//	vector<double> WI;
	//	int wType; //Well type,(0=Pro,1=inj)
	//	double rw;
	//	double s_time, e_time;
	//} Hor_Well;
	//

public:
	EDFM();
	void initialize();
	void re_orient_fracture(point &A, point &B);
	void read_h_fracture();
	void locate_start_end(point start, point end, unsigned &index_s, unsigned &index_e);
	void read_n_fracture();
	void read_reservoir();
	void index_to_cor(unsigned index, unsigned &x, unsigned &y);
	void compute_gridshape(unsigned index, point &large, point &small);
	void vector_to_matrix(vector<double> & x, dlib::matrix<double> & xm);
	dlib::matrix<double> dot(dlib::matrix<double> & xm, dlib::matrix<double> & xn);
	bool point_on_edge(point &start, point &small, point &large);
	template <typename T> vector<size_t> sort_indexes(const vector<T> &v);
	double Norm_distance(point c, point start, point end);
	point compute_centroid(vector<double> & vx, vector<double> &vy, double &A);
	double AvgNorm(point &start, point &end, point &small, point &large);
	void find_next_index(unsigned index_start, unsigned &index_end, point &start, point &end, unsigned index_ref, unsigned &next_index, double kf, point &ref, point &next, double &avernorm);
	double harm_mean(double k1, double k2);
	void read_grid_properties();
	void assign_frac(unsigned gridindex, unsigned frac_index, double perm, double por, double aper, double d, double l, point start, point end, frac &seg);
	double dis(point start, point end);
	void compute_fracture(fractureinf & fracture);
	double mult(point a, point b, point c);
	bool on_frac(point start, point end, point intersect);
	bool intersect(point aa, point bb, point cc, point dd);
	point fractureintersection(point fracture1_start, point fracture1_end, point fracture2_start, point fracture2_end);
	void compute_intersection(fractureinf);
	unsigned locate_seg(fractureinf fracture, point p);
	bool onSegment(point p, point q, point r);
	int orientation(point p, point q, point r);
	bool doIntersect(point p1, point q1, point p2, point q2);
	void find_inter_nn(vector<fractureinf> n_fracture);
	void find_inter_hn(vector<fractureinf> n_fracture, vector<fractureinf> h_fracture);
	void comput_interfrac_nn(fractureinf fracture);
	double compute_refer_k(double k, point &ref, point &large, point &small);
	void transfer_data();
	void compute_hw_f();
	void find_inter_hh(vector<fractureinf> n_fracture, vector<fractureinf> h_fracture);
	void  set_Hwell();
public:

	std::vector<frac> n_frac;
	intersection inter_nn;

private:
	vector<fractureinf> n_fracture;
	vector<fractureinf> h_fracture;
	std::vector<Hor_Well> h_well;
	//std::vector<frac> h_frac;
	vector<grid> grids;
	//vector<Hor_Well> h_well;

	unsigned nw;
	unsigned nh, nn;
	unsigned Nx, Ny, Nz;
	double Lx, Ly, Lz, Dx, Dy, Dz;
	unsigned frac_index = 0;
};



void EDFM::re_orient_fracture(point &A, point &B)
{
	point tmp;
	if (A.y > B.y)
	{
		sw(A, B);
	}
	if (A.x > B.x)
	{
		sw(A, B);
	}
}

void	EDFM::set_Hwell()
{
	std::ifstream strm("Input/horizontal_well.txt");
	Hor_Well tmp;
	unsigned nche,n;
	double temp;
	strm >> n;
	for (unsigned i = 0; i < n; i++)
	{
		strm >> tmp.isproducer;
		strm >> tmp.isPcontrol;
		strm >> tmp.rw;
		strm >> tmp.start.x;
		strm >> tmp.start.y;
		strm >> tmp.end.x;
		strm >> tmp.end.y;
		if (tmp.end.x < tmp.start.x)
		{
			std::swap(tmp.start.x, tmp.end.x);
			std::swap(tmp.start.y, tmp.end.y);
		}
		strm >> nche;
		tmp.Pw.resize(nche);
		tmp.Vol_std.resize(nche);
		for (unsigned i = 0; i < nche; i++)
		{
			strm >> temp;
			tmp.schedule.push_back(temp);
			strm >> temp;
		}
	}
	h_well.push_back(tmp);
}

void EDFM::read_h_fracture()
{
	point tempoint1, tempoint2;
	double tmp;
	ifstream strm("Input/hydraulicfracture.txt");
	strm >> nh;
	h_fracture.resize(nh);
	for (unsigned i = 0; i < nh; i++)
	{
		if (strm >> tmp) tempoint1.x = tmp;
		if (strm >> tmp) tempoint1.y = tmp;
		if (strm >> tmp) tempoint2.x = tmp;
		if (strm >> tmp) tempoint2.y = tmp;
		re_orient_fracture(tempoint1, tempoint2);
		h_fracture[i].start = tempoint1;
		h_fracture[i].end = tempoint2;
	}
	for (unsigned i = 0; i < nh; i++)
	{
		if (strm >> tmp) h_fracture[i].perm = tmp;
		if (strm >> tmp) h_fracture[i].aper = tmp;
		if (strm >> tmp) h_fracture[i].por = tmp;
		h_fracture[i].f_index = i;
	}
	strm.close();
}

void EDFM::read_n_fracture()
{
	point tempoint1, tempoint2;
	double tmp;
	ifstream strm("Input/fracturesystem.txt");
	strm >> nn;
	n_fracture.resize(nn);
	for (unsigned i = 0; i < nn; i++)
	{
		if (strm >> tmp) tempoint1.x = tmp;
		if (strm >> tmp) tempoint1.y = tmp;
		if (strm >> tmp) tempoint2.x = tmp;
		if (strm >> tmp) tempoint2.y = tmp;
		re_orient_fracture(tempoint1, tempoint2);
		n_fracture[i].start = tempoint1;
		n_fracture[i].end = tempoint2;
	}
	for (unsigned i = 0; i < nn; i++)
	{
		if (strm >> tmp) n_fracture[i].perm = tmp;
		if (strm >> tmp) n_fracture[i].aper = tmp;
		if (strm >> tmp) n_fracture[i].por = tmp;
		n_fracture[i].f_index = i;
	}
	strm.close();
}

void EDFM::locate_start_end(point start, point end, unsigned &index_s, unsigned &index_e)
{
	unsigned i, j,l,k;
	j = static_cast<unsigned>(start.y / Dy);
	i = static_cast<unsigned>(start.x / Dx);
	l = static_cast<unsigned>(end.y / Dy);
	k = static_cast<unsigned>(end.x / Dx);
	double slope = end.y - start.y;
	double rs= start.y / Dy - (int)(start.y / Dy);
	double re = end.x / Dx - (int)(end.x / Dx);
	//When the point on boundary;
	if (start.x == Lx)
		i = i - 1;
	if (end.y == Ly)
		l = l - 1;
	if (start.y == Ly)
		j = j - 1;
	//When point on edge;
	if (rs == 0 && slope <0 && start.y!=Ly)
	{
		j = j - 1;
	}
	if (re==0 && end.x!=0)
	{
		k = k - 1;
	}
	index_s = j*Nx + i;
	index_e = l*Nx + k;	
}

void EDFM::read_reservoir()
{
	ifstream strm("Input/Reservoir.txt");
	strm >> Nx;
	strm >> Ny;
	strm >> Nz;
	strm >> Lx;
	strm >> Ly;
	strm >> Lz;
	Dx = static_cast<double>(Lx / Nx);
	Dy = static_cast<double>(Ly / Ny);
	Dz = static_cast<double>(Lz / Nz);
}
void EDFM::index_to_cor(unsigned index, unsigned &x, unsigned &y)
{
	x = index % Nx;
	y = index / Nx;
}

void EDFM::compute_gridshape(unsigned index, point &large, point &small)
{
	unsigned x, y;
	index_to_cor(index,x,y);
	small.x = x*Dx;
	small.y = y*Dy;
	large.x = small.x + Dx;
	large.y = small.y + Dy;
}



double EDFM::compute_refer_k(double k, point &ref, point &large, point &small)
{
	double k_ref;
	if (k < 0)
	{
		k_ref = (ref.y - small.y) / (large.x-ref.x);
	}
	else
	{
		k_ref = (large.y-ref.y) / (large.x - ref.x);
	}
	return k_ref;
}

//void set_Hwell()
//{
//	ifstream strm("Input/horizontal_well.txt");
//	size_t f;
//	point tmp;
//	double ro, WI, hf;
//	strm >> nw;
//	Builder.nw = nw;
//	Builder.Hw.resize(nw);
//
//	for (unsigned i = 0; i < nw; i++)
//	{
//		strm >> Builder.Hw[i].wType;
//		strm >> Builder.Hw[i].rw;
//		strm >> Builder.Hw[i].Pw.value();
//		strm >> Builder.Hw[i].s_time;
//		strm >> Builder.Hw[i].e_time;
//		strm >> startp[i].x;
//		strm >> startp[i].y;
//		strm >> endp[i].x;
//		strm >> endp[i].y;
//		if (endp[i].x < startp[i].x)
//		{
//			tmp = endp[i];
//			endp[i] = startp[i];
//			startp[i] = tmp;
//		}
//	}
//
//	for (unsigned i = 0; i < nw; i++)
//	{
//		for (unsigned j = 0; j < start_end.nh; j++)
//		{
//			tmp = fractureintersection(startp[i], endp[i], start_end.start[j], start_end.end[j]);
//			if (tmp.x>-0.1)
//			{
//				f = locate_seg(j, tmp);
//				Builder.Hw[i].fracindex.push_back(f);
//				Builder.Hw[i].gridindex.push_back(fracture.gridindex[f]);
//				hf = grid.Dz[fracture.gridindex[f]];
//				ro = 0.14*sqrt(pow(fracture.length[f], 2) + pow(hf, 2));
//				WI = 2 * pi*fracture.aper[f] * fracture.perm[f] / log(ro / Builder.Hw[i].rw);
//				Builder.Hw[i].WI.push_back(WI);
//			}
//		}
//
//	}
//}






//point compute_centroid(vector<double> & x, vector<double> &y,double &A)
//{
//	point c;
//	double xc = 0, yc = 0;
//	A = 0;
//	for (unsigned i = 0; i < (x.size() - 1); i++)
//	{
//		A += 0.5*(x[i]*y[i+1] - x[i + 1]*y[i]);
//	}
//	for (unsigned i = 0; i < (x.size() - 1); i++)
//	{
//		xc += (x[i] + x[i + 1])*(x[i] * y[i + 1] - x[i + 1] * y[i]);
//		yc += (y[i] + y[i + 1])*(x[i] * y[i + 1] - x[i + 1] * y[i]);
//	}
//	c.x = xc / (6 * A);
//	c.y = yc / (6 * A);
//	return c;
//}

void EDFM::vector_to_matrix(vector<double> & x, dlib::matrix<double> & xm)
{
	unsigned n = x.size();
	xm.set_size(n,1);
	for (unsigned i = 0; i <n; i++)
	{
		xm(i, 0) = x[i];
	}
}

dlib::matrix<double> EDFM::dot(dlib::matrix<double> & xm, dlib::matrix<double> & xn)
{

	unsigned n = xm.size();
	dlib::matrix<double> x(n,1);
	for (unsigned i = 0; i <n; i++)
	{
		x(i) = xm(i)*xn(i);
	}

	return x;
}

//void tri_penta(point &start, point &end, point &tri, point &penta)
//{
//	vector<double> x_1, y_1,x_2,y_2;
//	double A_1, A_2;
//	
//	x_1.push_back(tri.x);
//	x_1.push_back(start.x);
//	x_1.push_back(end.x);
//
//	y_1.push_back(tri.y);
//	y_1.push_back(start.y);
//	y_1.push_back(end.y);
//
//	x_2.push_back(tri.x);
//	x_2.push_back(start.x);
//	x_2.push_back(end.x);
//	x_2.push_back(start.x);
//	x_2.push_back(end.x);
//
//
//}

//void checkpoint(point &start, point &end, point &small, point &large)
//{
//	double k= (end.y - start.y) / (end.x - start.x);
//	double kref;
//	if (small.x<start.x && start.x<large.x && small.y<start.y && start.y<large.y)
//	{
//		if (k >= 0)
//		{
//			kref = (small.y-end.y) / (small.x-end.x);
//			if (k>kref)
//			{
//
//			}
//			else
//			{
//
//			}
//		}
//		else
//		{
//			kref = (large.y-end.y) / (small.x-end.x);
//			if (k>kref)
//			{
//
//			}
//			else
//			{
//
//			}
//		}
//
//	}
//
//	if (small.x<end.x && end.x<large.x && small.y<end.y && end.y<large.y)
//	{
//
//
//	}
//
//}




bool EDFM::point_on_edge(point &start, point &small, point &large)
{

	if (small.x<start.x && start.x<large.x && small.y<start.y && start.y<large.y)
	{
		return false;
	}
	return true;
}

template <typename T>
vector<size_t> EDFM::sort_indexes(const vector<T> &v) {

	// initialize original index locations
	vector<size_t> idx(v.size());
	for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;

	// sort indexes based on comparing values in v
	sort(idx.begin(), idx.end(),
		[&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });

	return idx;
}


double EDFM::Norm_distance(point c, point start, point end)
{
	double l =pow(pow(start.x-end.x,2)+pow(start.y-end.y,2),0.5);
	double n = abs((end.y-start.y)*c.x-(end.x-start.x)*c.y+end.x*start.y-end.y*start.x);
	double tmp = n / l;
	return n / l;
}

EDFM::point EDFM::compute_centroid(vector<double> & vx, vector<double> &vy, double &A)
{
	dlib::matrix<double> x, y;
	point c;
	unsigned n = vx.size();
	vector<unsigned> index;
	dlib::matrix<double> dx(n, 1);
	dlib::matrix<double> dy(n, 1), xy(n, 1), xx(n, 1), yy(n, 1), ydx(n, 1), dxdx(n, 1),xdy,dxdy;
	double xc = 0, yc = 0;
	vector_to_matrix(vx, x);
	vector_to_matrix(vy, y);
	double mx = mean(x);
	double my = mean(y);
	x = x - mx;
	y = y - my;
	for (unsigned i = 1; i <n; i++)
	{
		index.push_back(i);
	}
	index.push_back(0);

	for (unsigned i = 0; i <n; i++)
	{
		dx(i,0) = x(index[i],0)-x(i,0) ;
		dy(i,0) = y(index[i],0)-y(i,0);
	}

	A = sum(trans(y)*dx - trans(x)*dy) / 2.0;
	xy = dot(x, y);
	xx = dot(x,x);
	yy = dot(y,y);
	xdy = dot(x,dy);
	dxdy = dot(dx,dy);
	ydx = dot(y, dx);
	dxdx = dot(dx,dx);
	//cout << dot(xy, dx) << endl;
	//cout << dot(xy, dx) << endl;
	//cout << dot(xy, dx) << endl;
	//cout << dot(xy, dx) << endl;
	double Axc = sum(6 * dot(xy, dx) - 3 * dot(xx, dy) + 3 * dot(ydx, dx) + dot(dxdx,dy)) / 12;
	double Ayc = sum(3 * dot(yy, dx) - 6 * dot(xy, dy) - 3 * dot(xdy,dy) - dot(dxdy,dy)) / 12;
	if (A < 0)
	{
		A = -A;
		Axc = -Axc;
		Ayc = -Ayc;
	}
	xc = Axc / A;
	yc = Ayc / A;
	c.x = xc + mx;
	c.y = yc + my;
	return c;
}

double EDFM::AvgNorm(point &start, point &end, point &small, point &large)
{
	vector<double> x(6), y(6), a(6), tx(6), ty(6), poly_1x, poly_2x, poly_1y, poly_2y;
	unsigned index_a, index_b;
	double A1, A2;
	//initialize x,y;
	x[0] = small.x;
	x[1] = small.x;
	x[2] = large.x;
	x[3] = large.x;
	x[4] = start.x;
	x[5] = end.x;
	y[0] = small.y;
	y[1] = large.y;
	y[2] = small.y;
	y[3] = large.y;
	y[4] = start.y;
	y[5] = end.y;
	double cx = (small.x + large.x) / 2;
	double cy = (small.y + large.y) / 2;
	for (unsigned i = 0; i < 6; i++)
	{
		a[i] = -1 * atan2(y[i] - cy, x[i] - cx);
	}
	vector<size_t> index = sort_indexes(a);
	for (unsigned i = 0; i < 6; i++)
	{
		tx[i] = x[index[i]];
		ty[i] = y[index[i]];
		if (tx[i] == start.x && ty[i] == start.y)
		{
			index_a = i;
		}
		if (tx[i] == end.x && ty[i] == end.y)
		{
			index_b = i;
		}
	}
	for (unsigned i = min(index_a, index_b); i<=max(index_a, index_b); i++)
	{
		poly_1x.push_back(tx[i]);
		poly_1y.push_back(ty[i]);
	}
	for (unsigned i = max(index_a, index_b); i<6; i++)
	{
		poly_2x.push_back(tx[i]);
		poly_2y.push_back(ty[i]);
	}
	for (unsigned i = 0; i<=min(index_a, index_b); i++)
	{
		poly_2x.push_back(tx[i]);
		poly_2y.push_back(ty[i]);
	}
	point cp_1 = compute_centroid(poly_1x, poly_1y, A1);
	point cp_2 = compute_centroid(poly_2x, poly_2y, A2);
	return  (Norm_distance(cp_1, start, end)*A1 + Norm_distance(cp_2, start, end)*A2 ) / (A1 + A2);

}

void EDFM::find_next_index(unsigned index_start, unsigned &index_end, point &start, point &end, unsigned index_ref, unsigned &next_index, double kf, point &ref, point &next, double &avernorm)
{
	point large, small;
	compute_gridshape(index_ref, large, small);
	double k_ref = compute_refer_k(kf, ref, large, small);

	if (abs(kf)>abs(k_ref))
	{
		next_index = index_ref + kf / abs(kf)*Nx;
		if (kf<0)
		{
			next.y = small.y;
			next.x = ref.x +(next.y - ref.y) / kf;
		}
		else
		{
			next.y = large.y;
			next.x = ref.x + (next.y - ref.y) / kf;
		}
	}
	else if (abs(kf)<abs(k_ref))
	{
		next_index = index_ref + 1;
		next.x = large.x;
		next.y = (next.x - ref.x)*kf + ref.y;
	}
	else
	{
		next_index = index_ref + kf / abs(kf)*(Nx + 1);
		next.x = large.x;
		next.y = kf*(next.x - ref.x) + ref.x;
	}

	if (index_ref == index_start)
	{
		small.x = start.x;
	}
	else if (index_ref == index_end)
	{
		large.x= end.x;
		eq(end, next);
	}
	//cout << " start  " << ref.x << " " << ref.y ;
	//cout << " end  " << next.x << " " << next.y << endl;

	avernorm = AvgNorm(ref, next, small, large);
	//cout << avernorm << endl;
}
double EDFM::harm_mean(double k1, double k2)
{
	return 2 * k1*k2 / (k1+k2);
}

void EDFM::read_grid_properties()
{
	
	ifstream input("input/Poro_Perm.txt");
	unsigned ng = Nx*Ny*Nz;
	grids.resize(ng);
	for (unsigned i = 0; i < ng;i++)
	{
		input >> grids[i].perm;
		input >> grids[i].poro;
	}
	
}

void EDFM::assign_frac(unsigned gridindex, unsigned frac_index, double perm, double por, double aper, double d, double l, point start, point end, frac &seg)
{
	seg.l = l;
	seg.aper = aper;
	seg.perm = perm;
	seg.f_index = frac_index;
	seg.g_index = gridindex;
	eq(start,seg.start);
	eq(end, seg.end);
	seg.por = por;
	seg.A = l*Dz;
	//assigned grid perm;
	seg.trans=harm_mean(grids[gridindex].perm,perm)*seg.A/d;
}

double EDFM::dis(point start, point end)
{
	return pow(pow(start.x - end.x, 2) + pow(start.y - end.y, 2), 0.5);
}





//ofstream output("Input/mf.txt");
void EDFM::compute_fracture( fractureinf & fracture)
{
	frac seg;
	point start, end;
	double norm_dis, trans,l,A;
	fracture.lower = frac_index;
	eq(fracture.start,start);
	eq(fracture.end, end);
	unsigned index_start, index_end;
	point ref, next;
	locate_start_end(start, end, index_start, index_end);
	eq(start, ref);
	unsigned index_ref = index_start,next_index;

	if (start.x == end.x)
	{
		point large_1, small_1;
		while (index_ref != index_end)
		{
			next.x = start.x;
			compute_gridshape(index_ref, large_1, small_1);
			next.y = large_1.y;
			index_ref = index_ref+Nx;
			if (index_ref==index_start)
			{
				small_1.y = start.y;
			}
			else if (index_ref ==index_end)
			{
				large_1.y = end.y;
			}
			norm_dis=AvgNorm(ref,next,small_1,large_1);
			l = dis(ref,next);
			assign_frac(index_ref, frac_index, fracture.perm, fracture.por, fracture.aper, norm_dis, l, ref, next, seg);
			n_frac.push_back(seg);
			//output << "index  " << index_ref + 1 << "  ";
			//output << "current point  " << ref.x << " " << ref.y << "  ";
			//output << "next point  " << next.x << " " << next.y << "  " << endl;
			//output << index_ref + 1 << "  " << frac_index + 1 << endl;
			eq(next,ref);
			frac_index++;
		}
	}
	else
	{
		double kf = (end.y - start.y) / (end.x - start.x);
		while (index_ref != index_end)
		{
			find_next_index(index_start,index_end,start,end,index_ref, next_index, kf, ref, next,norm_dis);
			l = dis(ref, next);
			assign_frac(index_ref, frac_index, fracture.perm, fracture.por, fracture.aper, norm_dis, l, ref, next, seg);
			n_frac.push_back(seg);
			//cout << norm_dis << endl;
			//output << index_ref + 1 << "  " << frac_index+1 << endl;
			eq(next, ref);
			index_ref = next_index;
			frac_index++;
		}
	}
	point large_1, small_1;
	compute_gridshape(index_ref, large_1, small_1);
	norm_dis = AvgNorm(ref, end, small_1, large_1);
	l = dis(ref, end);
	assign_frac(index_ref, frac_index, fracture.perm, fracture.por, fracture.aper, norm_dis, l, ref, end, seg);
	n_frac.push_back(seg);
	fracture.upper=frac_index;
	//output << index_ref + 1 << "  " << frac_index + 1 << endl;
	frac_index++;
}

double EDFM::mult(point a, point b, point c)
{
	return (a.x - c.x)*(b.y - c.y) - (b.x - c.x)*(a.y - c.y);
}

bool EDFM::on_frac(point start, point end, point intersect)
{
	double largex, smallx;
	if (intersect.y >= start.y && intersect.y <= end.y)
	{
		largex = std::max(start.x, end.x);
		smallx = std::min(start.x, end.x);
		if (intersect.x >= smallx && intersect.x <= largex)
		{
			return true;
		}
	}
	return false;
}


bool EDFM::intersect(point aa, point bb, point cc, point dd)
{
	if (std::max(aa.x, bb.x)<std::min(cc.x, dd.x))
	{
		return false;
	}
	if (std::max(aa.y, bb.y)<std::min(cc.y, dd.y))
	{
		return false;
	}
	if (std::max(cc.x, dd.x)<std::min(aa.x, bb.x))
	{
		return false;
	}
	if (std::max(cc.y, dd.y)<std::min(aa.y, bb.y))
	{
		return false;
	}
	if (mult(cc, bb, aa)*mult(bb, dd, aa)<0)
	{
		return false;
	}
	if (mult(aa, dd, cc)*mult(dd, bb, cc)<0)
	{
		return false;
	}
	return true;
}


EDFM::point EDFM::fractureintersection(point fracture1_start, point fracture1_end, point fracture2_start, point fracture2_end)
{
	point intersectp;
	double a;
	double b;
	double c;
	double d;
	if (intersect(fracture1_start, fracture1_end, fracture2_start, fracture2_end) == 1)
	{
		if (fracture1_start.x == fracture1_end.x)
		{
			c = (fracture2_start.y - fracture2_end.y) / (fracture2_start.x - fracture2_end.x);
			d = (fracture2_end.y*fracture2_start.x - fracture2_start.y*fracture2_end.x) / (fracture2_start.x - fracture2_end.x);
			intersectp.x = fracture1_start.x;
			intersectp.y = c*intersectp.x + d;
		}
		else if (fracture2_start.x == fracture2_end.x)
		{
			a = (fracture1_start.y - fracture1_end.y) / (fracture1_start.x - fracture1_end.x);
			b = (fracture1_end.y*fracture1_start.x - fracture1_start.y*fracture1_end.x) / (fracture1_start.x - fracture1_end.x);
			intersectp.x = fracture2_start.x;
			intersectp.y = a*intersectp.x + b;
		}
		else
		{
			a = (fracture1_start.y - fracture1_end.y) / (fracture1_start.x - fracture1_end.x);
			b = (fracture1_end.y*fracture1_start.x - fracture1_start.y*fracture1_end.x) / (fracture1_start.x - fracture1_end.x);
			c = (fracture2_start.y - fracture2_end.y) / (fracture2_start.x - fracture2_end.x);
			d = (fracture2_end.y*fracture2_start.x - fracture2_start.y*fracture2_end.x) / (fracture2_start.x - fracture2_end.x);
			intersectp.x = (d - b) / (a - c);
			intersectp.y = a*intersectp.x + b;

		}

	}
	else
	{
		intersectp.x = -1;
		intersectp.y = -1;
	}
	return intersectp;
}

void EDFM::compute_intersection(fractureinf)
{
	point tmp;

}

unsigned EDFM::locate_seg(fractureinf fracture, point p)
{
	for (unsigned i = fracture.lower; i <= fracture.upper; i++)
	{
		if (p.x > n_frac[i].start.x && p.x < n_frac[i].end.x)
			return i;
		if (p.x - n_frac[i].end.x < 0.01)
			return i;
	}
}

bool EDFM::onSegment(point p, point q, point r)
{
	if (q.x <= max(p.x, r.x) && q.x >= min(p.x, r.x) &&
		q.y <= max(p.y, r.y) && q.y >= min(p.y, r.y))
		return true;

	return false;
}

int EDFM::orientation(point p, point q, point r)
{
	int val = (q.y - p.y) * (r.x - q.x) -
		(q.x - p.x) * (r.y - q.y);

	if (val == 0) return 0;  // colinear

	return (val > 0) ? 1 : 2; // clock or counterclock wise
}

bool EDFM::doIntersect(point p1, point q1, point p2, point q2)
{
	// Find the four orientations needed for general and
	// special cases
	int o1 = orientation(p1, q1, p2);
	int o2 = orientation(p1, q1, q2);
	int o3 = orientation(p2, q2, p1);
	int o4 = orientation(p2, q2, q1);

	// General case
	if (o1 != o2 && o3 != o4)
		return true;

	// Special Cases
	// p1, q1 and p2 are colinear and p2 lies on segment p1q1
	if (o1 == 0 && onSegment(p1, p2, q1)) return true;

	// p1, q1 and p2 are colinear and q2 lies on segment p1q1
	if (o2 == 0 && onSegment(p1, q2, q1)) return true;

	// p2, q2 and p1 are colinear and p1 lies on segment p2q2
	if (o3 == 0 && onSegment(p2, p1, q2)) return true;

	// p2, q2 and q1 are colinear and q1 lies on segment p2q2
	if (o4 == 0 && onSegment(p2, q1, q2)) return true;

	return false; // Doesn't fall in any of the above cases
}


//ofstream output1("Input/f_f.txt");
void EDFM::find_inter_nn(vector<fractureinf> n_fracture)
{

	point tmp;
	double l1, l2,d;

	for (unsigned i = 0; i < n_fracture.size(); i++)
	{
		for (unsigned j = i + 1; j < n_fracture.size(); j++)
		{

			if (doIntersect(n_fracture[i].start, n_fracture[i].end, n_fracture[j].start, n_fracture[j].end))
			{

				tmp = fractureintersection(n_fracture[i].start, n_fracture[i].end, n_fracture[j].start, n_fracture[j].end);
				inter_nn.coor.push_back(tmp);
				unsigned x = locate_seg(n_fracture[i], tmp);
				unsigned y = locate_seg(n_fracture[j], tmp);
				//output1 << x << "  " << y << endl;
				inter_nn.A_i.push_back(n_frac[x].aper*Dx);
				inter_nn.A_j.push_back(n_frac[y].aper*Dx);
				inter_nn.perm_i.push_back(n_frac[x].perm);
				inter_nn.perm_j.push_back(n_frac[y].perm);
				l1 = dis(n_frac[x].start, tmp);
				l2 = dis(n_frac[x].end, tmp);
				d = (0.5*pow(l1, 2) + 0.5*pow(l2, 2)) / dis(n_frac[x].start, n_frac[x].end);
				inter_nn.d_i.push_back(d);
				l1 = dis(n_frac[y].start, tmp);
				l2 = dis(n_frac[y].end, tmp);
				d = (0.5*pow(l1, 2) + 0.5*pow(l2, 2)) / dis(n_frac[y].start, n_frac[y].end);
				inter_nn.d_j.push_back(d);
				inter_nn.i_index.push_back(x);
				inter_nn.j_index.push_back(y);
			}

		}
	}
}


void EDFM::find_inter_hn(vector<fractureinf> n_fracture, vector<fractureinf> h_fracture)
{

	point tmp;
	double l1, l2, d;

	for (unsigned i = 0; i < n_fracture.size(); i++)
	{
		for (unsigned j = 0; j < h_fracture.size(); j++)
		{
			if (doIntersect(n_fracture[i].start, n_fracture[i].end, h_fracture[j].start, h_fracture[j].end))
			{
				tmp = fractureintersection(n_fracture[i].start, n_fracture[i].end, h_fracture[j].start, h_fracture[j].end);
				inter_nn.coor.push_back(tmp);
				unsigned x = locate_seg(n_fracture[i], tmp);
				unsigned y = locate_seg(h_fracture[j], tmp);
				inter_nn.A_i.push_back(n_frac[x].aper*Dx);
				inter_nn.A_j.push_back(n_frac[y].aper*Dx);
				inter_nn.perm_i.push_back(n_frac[x].perm);
				inter_nn.perm_j.push_back(n_frac[y].perm);
				l1 = dis(n_frac[x].start, tmp);
				l2 = dis(n_frac[x].end, tmp);
				d = (0.5*pow(l1, 2) + 0.5*pow(l2, 2)) / dis(n_frac[x].start, n_frac[x].end);
				inter_nn.d_i.push_back(d);
				l1 = dis(n_frac[y].start, tmp);
				l2 = dis(n_frac[y].end, tmp);
				d = (0.5*pow(l1, 2) + 0.5*pow(l2, 2)) / dis(n_frac[y].start, n_frac[y].end);
				inter_nn.d_j.push_back(d);
				inter_nn.i_index.push_back(x);
				inter_nn.j_index.push_back(y);
			}
		}
	}
}


void EDFM::find_inter_hh(vector<fractureinf> n_fracture, vector<fractureinf> h_fracture)
{

	point tmp;
	double l1, l2, d;

	for (unsigned i = 0; i < h_fracture.size(); i++)
	{
		for (unsigned j = i+1; j < h_fracture.size(); j++)
		{
			if (doIntersect(h_fracture[i].start, h_fracture[i].end, h_fracture[j].start, h_fracture[j].end))
			{
				tmp = fractureintersection(h_fracture[i].start, h_fracture[i].end, h_fracture[j].start, h_fracture[j].end);
				inter_nn.coor.push_back(tmp);
				unsigned x = locate_seg(h_fracture[i], tmp);
				unsigned y = locate_seg(h_fracture[j], tmp);
				inter_nn.A_i.push_back(n_frac[x].aper*Dx);
				inter_nn.A_j.push_back(n_frac[y].aper*Dx);
				inter_nn.perm_i.push_back(n_frac[x].perm);
				inter_nn.perm_j.push_back(n_frac[y].perm);
				l1 = dis(n_frac[x].start, tmp);
				l2 = dis(n_frac[x].end, tmp);
				d = (0.5*pow(l1, 2) + 0.5*pow(l2, 2)) / dis(n_frac[x].start, n_frac[x].end);
				inter_nn.d_i.push_back(d);
				l1 = dis(n_frac[y].start, tmp);
				l2 = dis(n_frac[y].end, tmp);
				d = (0.5*pow(l1, 2) + 0.5*pow(l2, 2)) / dis(n_frac[y].start, n_frac[y].end);
				inter_nn.d_j.push_back(d);
				inter_nn.i_index.push_back(x);
				inter_nn.j_index.push_back(y);
			}
		}
	}
}



void
EDFM::comput_interfrac_nn(fractureinf fracture)
{
	unsigned n;
	double T;

	n = fracture.lower;
	while (n != fracture.upper)
	{
		inter_nn.i_index.push_back(n);
		inter_nn.j_index.push_back(n+1);
		//output1 << n << "  " << n + 1 << endl;
		inter_nn.A_i.push_back(n_frac[n].aper*Dx);
		inter_nn.A_j.push_back(n_frac[n+1].aper*Dx);
		inter_nn.coor.push_back(n_frac[n].end);
		inter_nn.d_i.push_back(n_frac[n].l / 2.0);
		inter_nn.d_j.push_back(n_frac[n+1].l / 2.0);
		inter_nn.perm_i.push_back(n_frac[n].perm);
		inter_nn.perm_j.push_back(n_frac[n+1].perm);
		n++;
	}
}
//template< typename T1>
//void EDFM::transfer_data(T1& model)
//{
//	unsigned n = inter_nn.coor.size();
//	model.n_frac.resize(n_frac.size());
//	model.inter_nn.i_index.resize(n);
//	model.inter_nn.j_index.resize(n);
//	model.inter_nn.A_i.resize(n);
//	model.inter_nn.A_j.resize(n);
//	model.inter_nn.coor.resize(n);
//	model.inter_nn.d_i.resize(n);
//	model.inter_nn.d_j.resize(n);
//	model.inter_nn.perm_i.resize(n);
//	model.inter_nn.perm_j.resize(n);
//	for (unsigned i = 0; i < n_frac.size();i++)
//	{
//		model.n_frac[i].A = n_frac[i].A;
//		model.n_frac[i].aper = n_frac[i].aper;
//		model.n_frac[i].perm = n_frac[i].perm;
//		model.n_frac[i].d = n_frac[i].d;
//		model.n_frac[i].l = n_frac[i].l;
//		model.n_frac[i].coor.x = n_frac[i].coor.x;
//		model.n_frac[i].coor.y = n_frac[i].coor.y;
//		model.n_frac[i].por = n_frac[i].por;
//		model.n_frac[i].trans = n_frac[i].trans;
//		model.n_frac[i].f_index = n_frac[i].f_index;
//		model.n_frac[i].g_index = n_frac[i].g_index;
//	}
//	for (unsigned i = 0; i < inter_nn.coor.size(); i++)
//	{
//		model.inter_nn.i_index[i] = inter_nn.i_index[i];
//		model.inter_nn.j_index[i] = inter_nn.j_index[i];
//		model.inter_nn.A_i[i] = inter_nn.A_i[i];
//		model.inter_nn.A_j[i] = inter_nn.A_j[i];
//		model.inter_nn.coor[i].x = inter_nn.coor[i].x;
//		model.inter_nn.coor[i].y = inter_nn.coor[i].y;
//		model.inter_nn.d_i[i] = inter_nn.d_i[i];
//		model.inter_nn.d_j[i] = inter_nn.d_j[i];
//		model.inter_nn.perm_i[i] = inter_nn.perm_i[i];
//		model.inter_nn.perm_j[i] = inter_nn.perm_j[i];
//	}
//}

void EDFM::compute_hw_f()
{
	point tmp;
	unsigned x, y;
	double ro,WI;
	set_Hwell();
	for (unsigned i = 0; i <h_well.size() ;i++)
	{
		for (unsigned j = 0; j < h_fracture.size();j++)
		{
			if (doIntersect(h_well[i].start, h_well[i].end, h_fracture[j].start, h_fracture[j].end))
			{
				tmp = fractureintersection(h_well[i].start, h_well[i].end, h_fracture[j].start, h_fracture[j].end);
				unsigned f = locate_seg(h_fracture[j], tmp);
				h_well[i].fracindex.push_back(f);
				locate_start_end(tmp, tmp, x, y);
				h_well[i].gridindex.push_back(x);
				ro = 0.14*sqrt(pow(n_frac[f].l, 2) + pow(Dz, 2));
				WI = 2 * 3.1415926*n_frac[f].aper * n_frac[f].perm / log(ro / h_well[i].rw);
				h_well[i].WI.push_back(WI);
			}
		}
	}

}

EDFM::EDFM()
{
	//load information
	read_reservoir();
	read_h_fracture();
	read_n_fracture();
	read_grid_properties();
	//compute m_f connection & compute intersection point
	for (unsigned i = 0; i < nn; i++)
	{
		compute_fracture(n_fracture[i]);
		comput_interfrac_nn(n_fracture[i]);
	}
	for (unsigned i = 0; i < nh;i++)
	{
	compute_fracture(h_fracture[i]);
	comput_interfrac_nn(h_fracture[i]);
	}
	find_inter_nn(n_fracture); 
	find_inter_hn(n_fracture, h_fracture);
	compute_hw_f();
}

#endif