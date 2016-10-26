#pragma once

#include <vector>
#include "adetl/systems/ADvector.hpp"
#include "PropertyCalculator.hpp"
#include "DiscreteDomain.hpp"
#include "CSR_Matrix.hpp"
#include <dlib/matrix.h>
//#include "EDFM.hpp"
typedef GENSOL::CSR_Matrix<double, int> CSR;
class DiscreteProblem
{
	friend class EDFM;
public:

	typedef struct
	{

		double rw;
		std::vector<double> schedule;
		std::vector<double> Vol_std;
		std::vector<size_t> fracindex;
		std::vector<size_t> gridindex;
		std::vector<double> WI_f;
		std::vector<double> WI;
		bool isproducer;  //Well type,(0=Pro,1=inj)
		bool isPcontrol;
		bool isvertical;
		bool ifbreakthrough = false;
	} Hor_well;

	typedef struct
	{
		adetl::ADscalar<> qo = 0;
		adetl::ADscalar<> qw = 0;
		adetl::ADscalar<> qg = 0;
		adetl::ADscalar<> P = 0;
		adetl::ADscalar<> WCT = 0;
		adetl::ADscalar<> GOR = 0;
		adetl::ADscalar<> qt = 0;

	}											wellrate;
	typedef struct
	{
		double qo = 0;
		double qw = 0;
		double qg = 0;
		double BHP = 0;
	}											rate;
	typedef struct
	{
		std::vector<double> val;
		std::vector<int> col, row;

	} Sparse_matrix;
	typedef struct
	{
		bool BHP_i, qt_i, qo_i,qg_i, qw_i,gor_i, wct_i;
		adetl::ADscalar<> BHP, qt, qo, qw,qg,gor, wct;
	} CS;

private: 
   typedef PropertyCalculator::Phase     PhaseID;
   typedef PropertyCalculator::StatusID  StatusID;




public:


	typedef struct
	{

		double rw;
		std::vector<double> schedule;
		std::vector<double> Vol_std;
		adetl::ADvector Pw;
		std::vector<size_t> fracindex;
		std::vector<size_t> gridindex;
		std::vector<double> WI;
		bool isproducer;  //Well type,(0=Pro,1=inj)
		bool isPcontrol;
	} Hor_Well;


   typedef struct
   {
      adetl::ADscalar<> Po;
      adetl::ADscalar<> Sw;
      adetl::ADscalar<> Sg;  // active for OWG
      adetl::ADscalar<> Rso; // active for OW
      StatusID   status;
   }                                            StateElement;

   typedef struct
   {
      bool is_producer;
	  bool is_vertical;
	  double rw;
	  std::vector<double> WI;
	  std::vector<unsigned> grid_index;
   }     
   Well;



public:
   typedef std::vector< StateElement >          StateVector;
   typedef struct
   {
      double MatBal[3];
      double NormSat[3];
   }                                            ConvergenceInfo;
private:
   typedef PropertyCalculator::CellProps        CellProps;
   typedef std::vector< CellProps >             CellVector;

   typedef struct
   {
      double            T;
      adetl::ADscalar<> L[3];
      adetl::ADscalar<> Pot[3];
      adetl::ADscalar<> Rso;
   }                                           FaceProps;
   typedef std::vector< FaceProps >            FaceVector;

   typedef DiscreteDomain                      MeshType;
   typedef EDFM									Fracturesystem;
public:

	/////////////////////////////////////////////////
	//void DiscreteProblem::extract_F_der(DiscreteProblem::StateVector &oldstate, DiscreteProblem::StateVector &state, double DT,
	// std::vector<CSR> &F_y, std::vector<CSR> &F_yold, std::vector<CSR> &F_v);
	void DiscreteProblem::extract_F_der(DiscreteProblem::StateVector &oldstate,
		DiscreteProblem::StateVector &state, double DT, DiscreteProblem::Sparse_matrix &M);

	void initialize_state(DiscreteProblem::StateVector &uold, DiscreteProblem::StateVector &unew,
		DiscreteProblem::StateVector &aold, DiscreteProblem::StateVector &anew);

	void DiscreteProblem::extract_obj_der(unsigned nsche);

	void make_all_independent(dlib::matrix<double> &x, adetl::ADvector &ind_v);

	void loadwell_data(const std::vector<double> &time, const dlib::matrix<double> &well_sche);

	void loadwell_independent(const std::vector<double> &time, const dlib::matrix<double> &well_sche);
	template< typename V, typename M >
	void extract_acc(V &r, M &m, std::size_t offset);

	void extract_acc(DiscreteProblem::Sparse_matrix &M);
	double compute_total_WI(Hor_well &h_well);

	////////////////////////////////////////////////
	DiscreteProblem(std::size_t NX, std::size_t NY, std::size_t NZ,
		double LX, double LY, double LZ,
		const std::vector<double> &vKX, const std::vector<double> &vKY, const std::vector<double> &vKZ,
		const std::vector<double> &Phi_ref);

   std::size_t  max_num_eqns ( ) const { return 3 * mMesh.size_cells( )+nw; }
   std::size_t  max_num_nnz() const { return 63 * mMesh.size_cells() + nw * 21; }

   void set_Hwell();
   void initialize();
   void initialize_state(StateVector & state);

   void bind_to_old_state( const StateVector &_old_state );

   bool discretize( const StateVector & state, double DT );
   bool discretize_m(const StateVector & state, double DT);
   bool is_converged( ConvergenceInfo & nrm );

   template< typename V, typename M >
   void extract_R_J(V &residual, M &jacobian, std::size_t i );

   template< typename V >
   void extract_dR_dDT(V &dR_dDT );
   
   template< typename R >
   void update_state( StateVector &_state, const R & update, bool do_safeguard );

   template<typename x>
   void load_hschedule(std::vector<std::vector<x>> &hwells, dlib::matrix<x> &s);

   void dump_residual( std::ofstream &out )
   {
      fastl::spy( out, mResidual);
   }
   std::size_t read_well();
protected:
   void setup_wells( std::size_t NX, std::size_t NY, std::size_t NZ, 
		     double LX, double LY, double LZ,
		     const std::vector<double> &vKX, 
		     const std::vector<double> &vKY, 
		     const std::vector<double> &vKZ  );

   std::size_t eqnID( std::size_t c, PhaseID ph ) const
   { 
      std::size_t incr; 
      switch ( ph )
      {
	 case PhaseID::O :
	    incr = 0;
	    break;
	 case PhaseID::W :
	    incr = 1;
	    break;
	 case PhaseID::G :
	    incr = 2;
	    break;
      }
      return 3*c + incr; 
   }

   std::size_t eqnID( std::size_t c, std::size_t ph ) const
   { 
      PhaseID incr; 
      switch ( ph )
      {
	 case 0 :
	    incr = PhaseID::O;
	    break;
	 case 1 :
	    incr = PhaseID::W;
	    break;
	 case 2:
	    incr = PhaseID::G ;
	    break;
      }
      return eqnID( c, incr ); 
   }

   std::size_t wtoc(std::size_t wellID) const
   {
	   return 3 * N_c + wellID;
   }


   void initialize_transmissibility( const std::vector<double> & KX,
				     const std::vector<double> & KY,
				     const std::vector<double> & KZ );
   double safeguard_MAC( double upd );
   void compute_accumulation( );
   void compute_flow( );
   void compute_wells ( );
   void compute_face_properties ( );
   void compute_cell_properties( const StateVector & state );

   ////////////////////////////////////////////////////////////////////

   void compute_BHP();
   void extract_inform(const std::vector<std::vector<unsigned>> &x, std::vector<std::vector<unsigned>> &inf);
   void compute_total_mobility(unsigned w, PhaseID I, adetl::ADscalar<> &T, adetl::ADscalar<> &Tp);

   void extract_R_Jm(std::vector<double> &r, CSR &m, std::size_t offset);

   void CSR_to_Sparse(CSR & x, DiscreteProblem::Sparse_matrix &M);

   void read_from_file_unsigned(std::string filename, dlib::matrix<unsigned> &v);

   void DiscreteProblem::compute_oil_water_ratio(std::vector<adetl::ADscalar<>> &BHP, unsigned w);

   void load_constrain();

   void assign_wellrate(std::vector<wellrate> &temp_v);

   void check_schedule_size(unsigned nw, unsigned nct, unsigned nwc);

   std::size_t well_location(std::size_t x, std::size_t y, std::size_t z);
   /////////////////////////////////////////////////////////////////////////
public:
	double  sum_t;
   std::size_t			nw, nc, Nj, Np, N_hj, N_hp, N_tc, nct, nwc;
   std::vector<Well>	 mWells;
   std::vector<double>  schedule;

   std::vector<std::vector<adetl::ADscalar<>>>   H_chedule;
   std::vector<std::vector<CS>>	H_constrain;
   std::vector<std::vector<unsigned>>	H_control;
   std::vector<wellrate>   Hwell_q;
   std::vector<rate>   well_rate;
   std::vector<bool> is_producer;
private:
   MeshType            mMesh;
   //EDFM					mFrac;
   PropertyCalculator  mPropCalc;
   CellVector          mCells;
   CellVector          fCells;
   FaceVector          mFaces;
   FaceVector		   ffFaces;
   FaceVector		   mfFaces;
   std::vector<double> mPhi_ref;
   std::vector<double> mAccum_old;
   adetl::ADvector     mAccum;   
   adetl::ADvector     mFlow;
   adetl::ADvector     mResidual;
   double              mDT;
   adetl::ADscalar<>   mTmpVars[3];
   std::size_t		  nhw, Nhp, Nhj;
   std::size_t        N_f, N_m, N_mm, N_fm, N_ff,N_var,N_c,Nx,Ny;
   std::vector<Hor_Well> h_well;
   std::size_t			nsche=0;
};

std::ostream & 
operator << ( std::ostream & ostr, 
	      const DiscreteProblem::ConvergenceInfo & _out );

#include "DiscreteProblem_IMPL.hpp"
