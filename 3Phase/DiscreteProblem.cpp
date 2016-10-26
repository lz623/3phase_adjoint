#include "DiscreteProblem.hpp"
#include "MeshBuilder.hpp"
#include "Point.hpp"
#include <fstream>

DiscreteProblem::DiscreteProblem( std::size_t NX, std::size_t NY, std::size_t NZ,
				  double LX, double LY, double LZ,
				  const std::vector<double> &vKX, const std::vector<double> &vKY, const std::vector<double> &vKZ,
				  const std::vector<double> &Phi_ref) :
				  Nx(NX),
				  Ny(NY),
   mMesh( MeshBuilder<UniformCartesian> ( NX,NY,NZ,LX,LY,LZ ) ),
   mPropCalc(  ),
   mCells( mMesh.size_cells() ),
   mFaces( mMesh.size_faces() ),
   mWells( ),
   mPhi_ref( Phi_ref ),
   mAccum_old( 3 * mMesh.size_cells() ),
   mAccum( 3 * mMesh.size_cells(), 0.0 ),
   mFlow( 3 * mMesh.size_cells(), 0.0 )
{
	setup_wells(NX, NY, NZ, LX, LY, LZ, vKX, vKY, vKZ);
	well_rate.resize(nw);
	N_c = mMesh.size_cells();
	N_tc = 3 * mMesh.size_cells() + nw;
	initialize_transmissibility( vKX, vKY, vKZ );
	mResidual.resize(N_tc);
};

std::size_t DiscreteProblem::well_location(std::size_t x, std::size_t y, std::size_t z)
{
	unsigned grid_index;
	return grid_index=x+y*Nx+z*Nx*Ny;
}

std::size_t DiscreteProblem::read_well()
{
	std::ifstream strm;
	strm.open("input/well_inform.dat");
	bool control_type;
	int idx=0, idy=0,idzs=0,idze=0;
	Np = 0;
	Nj = 0;
	strm >> nw;
	mWells.resize(nw);
	is_producer.resize(nw);
	for (unsigned i = 0; i < nw; i++)
	{
		strm >> mWells[i].is_vertical;
		strm >> mWells[i].is_producer;
		(mWells[i].is_producer) ? Np++ : Nj++;
		(mWells[i].is_producer) ? is_producer[i] = 1 : is_producer[i] = 0;
		strm >> control_type;
		strm >> mWells[i].rw;
		if (mWells[i].is_vertical)
		{
			strm >> idx;
			strm >> idy;
			strm >> idzs;
			strm >> idze;
			for (int k = idzs; k <= idze;k++)
			{
				mWells[i].grid_index.push_back(well_location(idx, idy, k));
			}
		}
		else
		{
			//////////////////////////code for horizontal well////////////////////////////////////////////////////
		}
	}
	return nw;
}

void DiscreteProblem::compute_total_mobility(unsigned w, PhaseID I, adetl::ADscalar<> &T, adetl::ADscalar<> &Tp)
{
	T = 0;
	Tp = 0;
	for (std::size_t m = 0; m < mWells[w].grid_index.size(); m++)
	{
		unsigned c = mWells[w].grid_index[m];
		double WI = mWells[w].WI[m];
		T += WI * mCells[c].Kr[I] * mCells[c].bmu[I];
		Tp += WI * mCells[c].Kr[I] * mCells[c].bmu[I] * mCells[c].P[0];
	}
}

void DiscreteProblem::compute_BHP()
{
	adetl::ADscalar<>  To, Top, Tw, Twp, T, Tp;
	for (std::size_t w = 0; w < nw; ++w)
	{
		if (mWells[w].is_producer)
		{
			if (H_control[w][nsche] == 0)
			{
				Hwell_q[w].P = H_chedule[w][nsche];
			}
			else if (H_control[w][nsche] == 1)
			{
				compute_total_mobility(w, PhaseID::O, To, Top);
				compute_total_mobility(w, PhaseID::W, Tw, Twp);
				T = To + Tw;
				Tp = Top + Twp;
				Hwell_q[w].P = (Tp - H_chedule[w][nsche]) / T;
			}
			else if (H_control[w][nsche] == 2)
			{
				compute_total_mobility(w, PhaseID::O, To, Top);
				Hwell_q[w].P = (Top - H_chedule[w][nsche]) / To;
			}
			else if (H_control[w][nsche] == 3)
			{
				compute_total_mobility(w, PhaseID::W, Tw, Twp);
				Hwell_q[w].P = (Twp - H_chedule[w][nsche]) / Tw;
			}
			else if (H_control[w][nsche] == 4)
			{
				compute_total_mobility(w, PhaseID::G, Tw, Twp);
				Hwell_q[w].P = (Twp - H_chedule[w][nsche]) / Tw;
			}
			else
			{
				std::cout << "Invalid well control for producer. " << std::endl;
				system("pause");
			}
		}
		else
		{
			if (H_control[w][nsche] == 0 || H_control[w][nsche] == 7)
			{
				Hwell_q[w].P = H_chedule[w][nsche];
			}
			else if (H_control[w][nsche] == 3)
			{
				compute_total_mobility(w, PhaseID::W, Tw, Twp);
				Hwell_q[w].P = (Twp - H_chedule[w][nsche]) / Tw;
			}
			else if (H_control[w][nsche] == 4)
			{
				compute_total_mobility(w, PhaseID::G, Tw, Twp);
				Hwell_q[w].P = (Twp - H_chedule[w][nsche]) / Tw;
			}
			else if (H_control[w][nsche] == 8)
			{
				compute_total_mobility(w, PhaseID::G, Tw, Twp);
				Hwell_q[w].P = (Twp - H_chedule[w][nsche]) / Tw;
			}
			else
			{
				std::cout << "Invalid well control for injector. " << std::endl;
				system("pause");
			}
		}
	}
}



void
DiscreteProblem::setup_wells(std::size_t NX, std::size_t NY, std::size_t NZ, 
			     double LX, double LY, double LZ,
			     const std::vector<double> &vKX, 
			     const std::vector<double> &vKY, 
			     const std::vector<double> &vKZ  )
{
	std::size_t nw;
	nw = read_well();
	Hwell_q.resize(nw);
	const double DX = LX / static_cast<double>(NX);
	const double DY = LY / static_cast<double>(NY);
	const double DZ = LZ / static_cast<double>(NZ);
	for (std::size_t i = 0; i < nw; i++)
	{
		mWells[i].WI.resize(mWells[i].grid_index.size());
		for (std::size_t j = 0; j < mWells[i].grid_index.size();j++)
		{
			const double Kh = DZ * sqrt(vKY[mWells[i].grid_index[j]] * vKX[mWells[i].grid_index[j]]);
			const double r1 = vKY[mWells[i].grid_index[j]] / vKX[mWells[i].grid_index[j]];
			const double r2 = vKX[mWells[i].grid_index[j]] / vKY[mWells[i].grid_index[j]];
			const double ro = 0.28 * std::sqrt(std::sqrt(r1)*DX*DX + std::sqrt(r2)*DY*DY) / (std::pow(r1, 0.25) + std::pow(r2, 0.25));
			mWells[i].WI[j] = 0.00708*Kh / log(ro / mWells[i].rw);
		}
	}
}

void 
DiscreteProblem::initialize_state(DiscreteProblem::StateVector &state)
{
   state.resize(N_c);
   std::ifstream strm("./input/initial_state.dat");
   double Pi, Sw,Sg,Rs = 0;
   strm >> Pi >> Sw >> Sg>>Rs;
   for (std::size_t c = 0; c < N_c; ++c ) 
   {
      state[c].Po     = Pi;
      state[c].Po.make_independent( eqnID(c,PhaseID::O) ) ;
      state[c].Sw     = Sw;
      state[c].Sw.make_independent( eqnID(c,PhaseID::W) ) ;
      state[c].Sg     = Sg;
	  if (Sg==0)
	  {
		  state[c].Rso = Rs*178.1076;
		  state[c].Rso.make_independent(eqnID(c, PhaseID::G));
		  state[c].status = StatusID::OW;
	  }
	  else
	  {
		  state[c].Sg.make_independent(eqnID(c, PhaseID::G));
		  state[c].Rso.make_constant();
		  state[c].status = StatusID::OWG;
	  }
   }
   for (std::size_t w = 0; w < nw; w++)
   {
	   std::size_t c = wtoc(w);
	   Hwell_q[w].P = Pi;
	   Hwell_q[w].P.make_independent(c);
   }
}


void DiscreteProblem::loadwell_independent(const std::vector<double> &time, const dlib::matrix<double> &well_sche)
{
	dlib::matrix<unsigned> ct, cs;  // input well constrain and well control
	read_from_file_unsigned("input/wellcontrol.dat", ct);
	read_from_file_unsigned("input/wellconstrain.dat", cs);
	check_schedule_size(nw, time.size(), well_sche.size());
	schedule = time;
	nct = time.size();
	nwc = well_sche.size();
	H_chedule.resize(nw);
	for (unsigned i = 0; i < nw; i++)
	{
		H_chedule[i].resize(nct);
		for (unsigned j = 0; j < time.size(); j++)
		{
			H_chedule[i][j].value() = well_sche(i*time.size() + j);
			H_chedule[i][j].make_independent(i*time.size() + j);
		}
	}
	load_hschedule(H_control, ct);
	load_constrain();
}


void 
DiscreteProblem::bind_to_old_state( const DiscreteProblem::StateVector &old_state )
{
   compute_cell_properties( old_state );
   compute_accumulation( );
   for ( std::size_t c=0; c<mMesh.size_cells(); ++c )
   {
      std::size_t eqn1 = eqnID(c,PhaseID::O);
      std::size_t eqn2 = eqnID(c,PhaseID::W);
      std::size_t eqn3 = eqnID(c,PhaseID::G);
      mAccum_old[ eqn1 ] = mAccum[ eqn1 ].value();
      mAccum_old[ eqn2 ] = mAccum[ eqn2 ].value();
      mAccum_old[ eqn3 ] = mAccum[ eqn3 ].value();
   }
}

bool
DiscreteProblem::discretize( const DiscreteProblem::StateVector &state, double DT )
{
   bool is_badvalue = false;
   mDT = DT;
   compute_cell_properties( state );
   compute_accumulation( );
   compute_flow( );
   for ( std::size_t c=0; c<mMesh.size_cells(); ++c )
   {
      std::size_t eqn1 = eqnID(c,PhaseID::O);
      std::size_t eqn2 = eqnID(c,PhaseID::W);
      std::size_t eqn3 = eqnID(c,PhaseID::G);
      mResidual[eqn1] = (mAccum[eqn1] - mAccum_old[eqn1]) + mDT * mFlow[eqn1];
      mResidual[eqn2] = (mAccum[eqn2] - mAccum_old[eqn2]) + mDT * mFlow[eqn2];
      mResidual[eqn3] = (mAccum[eqn3] - mAccum_old[eqn3]) + mDT * mFlow[eqn3];
      if ( !std::isfinite( mResidual[eqn1].value() ) ) is_badvalue = true;
      if ( !std::isfinite( mResidual[eqn2].value() ) ) is_badvalue = true;
      if ( !std::isfinite( mResidual[eqn3].value() ) ) is_badvalue = true;
   }   
   for (std::size_t w = 0; w<nw; ++w)
   {
	   std::size_t eqn1 = wtoc(w);
	   if (H_control[w][nsche] == 0 || H_control[w][nsche] == 7)
	   {
		   mResidual[eqn1] = H_chedule[w][nsche] - Hwell_q[w].P;
	   }
	   else if (H_control[w][nsche] == 1)
	   {
		   mResidual[eqn1] = mDT*(H_chedule[w][nsche] - Hwell_q[w].qt);
	   }
	   else if (H_control[w][nsche] == 2)
	   {
		   mResidual[eqn1] = mDT*(H_chedule[w][nsche] - Hwell_q[w].qo);
	   }
	   else if (H_control[w][nsche] == 3)
	   {
		   mResidual[eqn1] = mDT*(H_chedule[w][nsche] - Hwell_q[w].qw);
	   }
	   else if (H_control[w][nsche] == 4 || H_control[w][nsche] == 8)
	   {
		   mResidual[eqn1] = mDT*(H_chedule[w][nsche] - Hwell_q[w].qg);
	   }
	   else if (H_control[w][nsche] == 5)
	   {
		   mResidual[eqn1] = mDT*(H_chedule[w][nsche] - Hwell_q[w].WCT);
	   }
	   else if (H_control[w][nsche] == 6)
	   {
		   mResidual[eqn1] = mDT*(H_chedule[w][nsche] - Hwell_q[w].GOR);
	   }
   }

   return is_badvalue;
}

bool
DiscreteProblem::is_converged ( ConvergenceInfo & nrm )
{
   bool is_MATBAL_converged = true;
   bool is_NRMSAT_converged = true;
   bool is_well_converged = true;
   double tot_PV = 0.0;
   for ( std::size_t c=0; c<mMesh.size_cells( ); ++c )
      tot_PV += mCells[c].phi.value() * mMesh.cell_measure( c )*0.1781076;
   
   double sum_well = 0;
   for (std::size_t w = 0; w < nw; ++w)
   {
	   std::size_t c = wtoc(w);
	   sum_well += abs(mResidual[c].value());
   }
   if (abs(sum_well)> 1.0e-5) is_well_converged = false;

   for ( std::size_t phs = 0; phs < 3; ++phs )
   {
      double sum_R = 0.0;
      double avg_B = 0.0;
      double max_R_PV = 0.0;
      for ( std::size_t c=0; c<mMesh.size_cells( ); ++c )
      {
		sum_R    += mResidual[ eqnID( c, phs ) ].value();
		avg_B    += 1.0 / mCells[c].b[ phs ].value();
		double R_PV = std::abs(mResidual[ eqnID( c, phs ) ].value() / (mCells[c].phi.value() * mMesh.cell_measure( c )*0.1781076 ));
		if (R_PV > max_R_PV) max_R_PV = R_PV;
      }
      avg_B /= mMesh.size_cells( );
      nrm.MatBal[ phs ]  = std::abs( avg_B * sum_R / tot_PV );
      nrm.NormSat[ phs ] = avg_B * max_R_PV;
      if ( nrm.MatBal[ phs ] > 1.0e-7 )  is_MATBAL_converged = false;
      if ( nrm.NormSat[ phs ] > 0.001 )  is_NRMSAT_converged = false;
   }
   return (is_MATBAL_converged && is_NRMSAT_converged &&is_well_converged);
}

void 
DiscreteProblem::initialize_transmissibility( const std::vector<double> & KX,
					     const std::vector<double> & KY,
					     const std::vector<double> & KZ )
{
   Point KL, KR;
   for ( std::size_t f=0; f < mMesh.size_faces( ) ; ++f )
   {
      const std::size_t c1 = mMesh.face_adj_cell_1( f );
      const std::size_t c2 = mMesh.face_adj_cell_2( f );
      KL.p[0] = KX[c1];
      KL.p[1] = KY[c1];
      KL.p[2] = KZ[c1];
      KR.p[0] = KX[c2];
      KR.p[1] = KY[c2];
      KR.p[2] = KZ[c2];
      const Point nrml = mMesh.unit_normal( f );
      const double kl = dot_product( nrml, KL );
      const double kr = dot_product( nrml, KR );
      const double DX = norm( mMesh.cell_coord(c1) - mMesh.cell_coord( c2 ) );
      const double A  = mMesh.face_measure( f );
      const double beta = DX*0.5/A*( 1.0/kl + 1.0/kr );
      mFaces[f].T = 0.00112712 / beta;
   }
};

void 
DiscreteProblem::compute_cell_properties( const DiscreteProblem::StateVector &state )
{
   for (std::size_t c = 0; c < mMesh.size_cells(); ++c ) 
   {
      mPropCalc.calculate( mPhi_ref[c], state[c], mCells[c] );
   }
}

void 
DiscreteProblem::compute_accumulation( )
{
   for (std::size_t c = 0; c < mMesh.size_cells(); ++c ) 
   {
      const double VOLUME = mMesh.cell_measure( c )*0.1781076;
      mAccum[ eqnID(c,PhaseID::O) ] = VOLUME * mCells[c].phi * mCells[c].S[PhaseID::O] * mCells[c].b[PhaseID::O];
      mAccum[ eqnID(c,PhaseID::W) ] = VOLUME * mCells[c].phi * mCells[c].S[PhaseID::W] * mCells[c].b[PhaseID::W];
      mAccum[ eqnID(c,PhaseID::G) ] = VOLUME * mCells[c].phi * ( mCells[c].S[PhaseID::G] * mCells[c].b[PhaseID::G] +
                                                               mCells[c].Rso * mCells[c].S[PhaseID::O] * mCells[c].b[PhaseID::O] );
   }
}

void 
DiscreteProblem::compute_flow( )
{
   compute_face_properties( );

   for (std::size_t c = 0; c < (3*mMesh.size_cells()); ++c ) 
      mFlow[c] = 0.0;

   for (std::size_t f = 0; f < mMesh.size_faces(); ++f ) 
   {
      std::size_t c1 = mMesh.face_adj_cell_1( f );
      std::size_t c2 = mMesh.face_adj_cell_2( f );

      mTmpVars[PhaseID::W] = mFaces[f].T * mFaces[f].L[PhaseID::W] * mFaces[f].Pot[PhaseID::W];
      mTmpVars[PhaseID::O] = mFaces[f].T * mFaces[f].L[PhaseID::O] * mFaces[f].Pot[PhaseID::O];
      mTmpVars[PhaseID::G] = mFaces[f].T * mFaces[f].L[PhaseID::G] * mFaces[f].Pot[PhaseID::G] + 
         mFaces[f].Rso * mTmpVars[PhaseID::O];
      
      mFlow[ eqnID(c1,PhaseID::O) ] -= mTmpVars[ PhaseID::O ];
      mFlow[ eqnID(c1,PhaseID::W) ] -= mTmpVars[ PhaseID::W ];
      mFlow[ eqnID(c1,PhaseID::G) ] -= mTmpVars[ PhaseID::G ];

      mFlow[ eqnID(c2,PhaseID::O) ] += mTmpVars[ PhaseID::O ];
      mFlow[ eqnID(c2,PhaseID::W) ] += mTmpVars[ PhaseID::W ];
      mFlow[ eqnID(c2,PhaseID::G) ] += mTmpVars[ PhaseID::G ];
   }
   compute_wells( );
}

void 
DiscreteProblem::compute_face_properties( )
{
   for (std::size_t f = 0; f < mMesh.size_faces(); ++f ) 
   {
      std::size_t c1 = mMesh.face_adj_cell_1( f );
      std::size_t c2 = mMesh.face_adj_cell_2( f );
      double      dz = ( mMesh.cell_coord(c2) - mMesh.cell_coord(c1) ).p[2];

      for ( std::size_t ph=0; ph < 3; ++ph )
      {
		mFaces[f].Pot[ph] = (mCells[c2].P[ph] - mCells[c1].P[ph]) + 
			0.00694 * 0.5 * (mCells[c2].Rho[ph] + mCells[c1].Rho[ph]) * dz;
		mFaces[f].L[ph] = 0.5 * (mCells[c2].bmu[ph] + mCells[c1].bmu[ph]);
		if ( mFaces[f].Pot[ph].value() > 0.0 )
			mFaces[f].L[ph] *= mCells[c2].Kr[ph];
		else
			mFaces[f].L[ph] *= mCells[c1].Kr[ph];
      }
      mFaces[f].Rso = 0.5 * ( mCells[c2].Rso + mCells[c1].Rso);
   }
}


void 
DiscreteProblem::compute_wells( )
{
	if (schedule[nsche] == sum_t)
	{
		nsche++;
		std::cout << "The " << nsche << "th schedule at " << sum_t << " day." << std::endl;
	}
   for ( std::size_t w = 0; w<mWells.size(); ++w )
   {
	   adetl::ADscalar<> qo = 0, qw = 0,qg=0;
	   if (mWells[w].is_producer)
	   {
		   //compute_oil_water_ratio(BHP,w);
		   for (std::size_t m = 0; m < mWells[w].grid_index.size(); m++)
		   {
			   unsigned c = mWells[w].grid_index[m];

			   if (mCells[c].P[PhaseID::O].value() >= Hwell_q[w].P.value())
				   mTmpVars[0] = mCells[c].P[PhaseID::O] - Hwell_q[w].P;
			   else  
				   mTmpVars[0].value() = 0;
			   double WI = mWells[w].WI[m];
			   qo += WI * mTmpVars[0] * mCells[c].Kr[PhaseID::O] * mCells[c].bmu[PhaseID::O];
			   qw += WI * mTmpVars[0] * mCells[c].Kr[PhaseID::W] * mCells[c].bmu[PhaseID::W];
			   qg += WI * mTmpVars[0] * (mCells[c].Kr[PhaseID::G] * mCells[c].bmu[PhaseID::G] +
				   mCells[c].Rso * mCells[c].Kr[PhaseID::O] * mCells[c].bmu[PhaseID::O] );
			   mFlow[eqnID(c, PhaseID::O)] += WI * mTmpVars[0] * mCells[c].Kr[PhaseID::O] * mCells[c].bmu[PhaseID::O];
			   mFlow[eqnID(c, PhaseID::W)] += WI * mTmpVars[0] * mCells[c].Kr[PhaseID::W] * mCells[c].bmu[PhaseID::W];
			   mFlow[eqnID(c, PhaseID::G)] += WI * mTmpVars[0] * (mCells[c].Kr[PhaseID::G] * mCells[c].bmu[PhaseID::G] +
				   mCells[c].Rso * mCells[c].Kr[PhaseID::O] * mCells[c].bmu[PhaseID::O] );

		   }
		   Hwell_q[w].qo = qo;
		   Hwell_q[w].qw = qw;
		   Hwell_q[w].qg = qg;
		   Hwell_q[w].qt = qo + qw;
		   Hwell_q[w].WCT = qw / Hwell_q[w].qt;
		   Hwell_q[w].GOR = qg / qo;
	   }
	   else
	   {
		   adetl::ADscalar<> total_mob = 0;
		   if (H_control[w][nsche] == 0)
		   {
			   for (std::size_t m = 0; m < mWells[w].grid_index.size(); m++)
			   {
				   unsigned c = mWells[w].grid_index[m];
				   if (mCells[c].P[PhaseID::O].value()<=Hwell_q[w].P.value())
				   {
					   mTmpVars[0] = mCells[c].P[PhaseID::O] - Hwell_q[w].P;
				   }
				   else
				   {
					   mTmpVars[0].value() = 0;
				   }
				   double WI = mWells[w].WI[m];
				   total_mob = mCells[c].Kr[PhaseID::W] * mCells[c].bmu[PhaseID::W] +
					   mCells[c].Kr[PhaseID::O] * mCells[c].bmu[PhaseID::O] + 
					   mCells[c].Kr[PhaseID::G] * mCells[c].bmu[PhaseID::G];
				   qw += WI * mTmpVars[0] * total_mob;
				   mFlow[eqnID(c, PhaseID::W)] += WI * mTmpVars[0] * total_mob;
			   }
			   Hwell_q[w].qw = qw;
		   }
		   else if (H_control[w][nsche] == 3)
		   {
			   for (std::size_t m = 0; m < mWells[w].grid_index.size(); m++)
			   {
				   unsigned c = mWells[w].grid_index[m];
				   double WI = mWells[w].WI[m];
				   mTmpVars[0] = mCells[c].P[PhaseID::O] - Hwell_q[w].P;
				   total_mob = mCells[c].Kr[PhaseID::W] * mCells[c].bmu[PhaseID::W] +
					   mCells[c].Kr[PhaseID::O] * mCells[c].bmu[PhaseID::O] +
					   mCells[c].Kr[PhaseID::G] * mCells[c].bmu[PhaseID::G];
				   qw += WI * mTmpVars[0] * total_mob;
				   mFlow[eqnID(c, PhaseID::W)] += WI * mTmpVars[0] * total_mob;
			   }
			   Hwell_q[w].qw = qw;
		   }
		   else if (H_control[w][nsche] == 4)
		   {
			   for (std::size_t m = 0; m < mWells[w].grid_index.size(); m++)
			   {
				   unsigned c = mWells[w].grid_index[m];
				   double WI = mWells[w].WI[m];
				   mTmpVars[0] = mCells[c].P[PhaseID::O] - Hwell_q[w].P;
				   total_mob = mCells[c].Kr[PhaseID::W] * mCells[c].bmu[PhaseID::W] +
					   mCells[c].Kr[PhaseID::O] * mCells[c].bmu[PhaseID::O] +
					   mCells[c].Kr[PhaseID::G] * mCells[c].bmu[PhaseID::G];
				   qg += WI * mTmpVars[0] * total_mob;
				   mFlow[eqnID(c, PhaseID::G)] += WI * mTmpVars[0] * total_mob;
			   }
			   Hwell_q[w].qg = qg;
		   }
		   else if (H_control[w][nsche] == 7 || H_control[w][nsche] == 8)
		   {
			   for (std::size_t m = 0; m < mWells[w].grid_index.size(); m++)
			   {
				   unsigned c = mWells[w].grid_index[m];
				   double WI = mWells[w].WI[m];
				   if (mCells[c].P[PhaseID::O].value() <= Hwell_q[w].P.value())
				   {
					   mTmpVars[0] = mCells[c].P[PhaseID::O] - Hwell_q[w].P;
				   }
				   else
				   {
					   mTmpVars[0].value() = 0;
				   }
				   total_mob = mCells[c].Kr[PhaseID::W] * mCells[c].bmu[PhaseID::W] +
					   mCells[c].Kr[PhaseID::O] * mCells[c].bmu[PhaseID::O] +
					   mCells[c].Kr[PhaseID::G] * mCells[c].bmu[PhaseID::G];
				   qg += WI * mTmpVars[0] * total_mob;
				   mFlow[eqnID(c, PhaseID::G)] += WI * mTmpVars[0] * total_mob;
			   }
			   Hwell_q[w].qg = qg;
		   }
		   else
		   {
			   std::cout << H_control[w][nsche] << "   Invalid well control. " << std::endl;
		   }
	   }
   }
}

void
DiscreteProblem::initialize_state(DiscreteProblem::StateVector &uold, DiscreteProblem::StateVector &unew,
DiscreteProblem::StateVector &aold, DiscreteProblem::StateVector &anew)
{
	for (unsigned c = 0; c< N_c; c++)
	{
		aold[c].status = uold[c].status;
		anew[c].status = unew[c].status;
		aold[c].Po = uold[c].Po.value();
		anew[c].Po = unew[c].Po.value();
		aold[c].Sw = uold[c].Sw.value();
		anew[c].Sw = unew[c].Sw.value();
		aold[c].Sg = uold[c].Sg.value();
		anew[c].Sg = unew[c].Sg.value();
		aold[c].Rso = uold[c].Rso.value();
		anew[c].Rso = unew[c].Rso.value();
	}
}

void DiscreteProblem::extract_obj_der(unsigned nsche)
{
	unsigned ny = max_num_eqns();
	std::vector<double> obj_y(ny);
	adetl::ADscalar<> x;
	//assign_wellrate(temp_y);
	for (unsigned i = 0; i < nw; i++)
	{
		if (H_constrain[nsche][i].BHP_i)
		{
			H_constrain[nsche][i].BHP = -Hwell_q[i].P;
			well_rate[i].BHP = Hwell_q[i].P.value();
		}
		if (H_constrain[nsche][i].qt_i)
		{
			H_constrain[nsche][i].qt = -Hwell_q[i].qt;
		}
		if (H_constrain[nsche][i].qo_i)
		{
			H_constrain[nsche][i].qo = -Hwell_q[i].qo;
			well_rate[i].qo = Hwell_q[i].qo.value();
		}
		if (H_constrain[nsche][i].qw_i)
		{
			H_constrain[nsche][i].qw = -Hwell_q[i].qw;
			well_rate[i].qw = Hwell_q[i].qw.value();
		}
		if (H_constrain[nsche][i].qg_i)
		{
			H_constrain[nsche][i].qg = -Hwell_q[i].qg;
			well_rate[i].qg = Hwell_q[i].qg.value();
		}
		if (H_constrain[nsche][i].gor_i)
		{
			H_constrain[nsche][i].gor = -Hwell_q[i].GOR;
		}
		if (H_constrain[nsche][i].wct_i)
		{
			H_constrain[nsche][i].wct = -Hwell_q[i].WCT;
		}
	}
}



void DiscreteProblem::loadwell_data(const std::vector<double> &time, const dlib::matrix<double> &well_sche)
{
	dlib::matrix<unsigned> ct, cs;  // input well constrain and well control
	read_from_file_unsigned("input/wellcontrol.dat", ct);
	check_schedule_size(nw, time.size(), well_sche.size());
	schedule = time;
	nwc = ct.size();
	nct = time.size();
	H_chedule.resize(nw);
	for (unsigned i = 0; i < Np; i++)
	{
		H_chedule[i].resize(time.size());
		for (unsigned j = 0; j < time.size(); j++)
		{
			H_chedule[i][j] = well_sche(i*time.size() + j);
		}
	}
	for (unsigned i = Np; i < nw; i++)
	{
		H_chedule[i].resize(time.size());
		for (unsigned j = 0; j < time.size(); j++)
		{
			H_chedule[i][j] = well_sche(i*time.size() + j);
		}
	}
	load_hschedule(H_control, ct);
	load_constrain();
}


void DiscreteProblem::load_constrain()
{
	std::ifstream input("input/wellconstrain.dat");
	H_constrain.resize(nct);
	for (unsigned w = 0; w < nw; w++)
	{
		for (unsigned i = 0; i < nct; i++)
		{
			H_constrain[i].resize(nw);
			input >> H_constrain[i][w].BHP_i;
			input >> H_constrain[i][w].qt_i;
			input >> H_constrain[i][w].qo_i;
			input >> H_constrain[i][w].qw_i;
			input >> H_constrain[i][w].qg_i;
			input >> H_constrain[i][w].wct_i;
			input >> H_constrain[i][w].gor_i;
		}
	}


}


void DiscreteProblem::check_schedule_size(unsigned nw, unsigned nct, unsigned nwc)
{
	if (nw*nct != nwc)
	{
		std::cout << "The schedule's size is not match" << std::endl;
		system("pause");
	}
}

void
DiscreteProblem::read_from_file_unsigned(std::string filename, dlib::matrix<unsigned> &v)
{
	double tmp;
	size_t n = 0;
	std::ifstream strm(filename);
	while (strm >> tmp)
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




bool
DiscreteProblem::discretize_m(const DiscreteProblem::StateVector &state, double DT)
{
	bool is_badvalue = false;
	mDT = DT;
	compute_cell_properties(state);
	compute_accumulation();
	compute_BHP();
	compute_flow();
	for (std::size_t c = 0; c<N_c; ++c)
	{
		std::size_t eqn1 = eqnID(c, PhaseID::O);
		std::size_t eqn2 = eqnID(c, PhaseID::W);
		std::size_t eqn3 = eqnID(c, PhaseID::G);
		mResidual[eqn1] = (mAccum[eqn1] - mAccum_old[eqn1]) + mDT * mFlow[eqn1];
		mResidual[eqn2] = (mAccum[eqn2] - mAccum_old[eqn2]) + mDT * mFlow[eqn2];
		mResidual[eqn3] = (mAccum[eqn3] - mAccum_old[eqn3]) + mDT * mFlow[eqn3];
		if (!std::isfinite(mResidual[eqn1].value())) is_badvalue = true;
		if (!std::isfinite(mResidual[eqn2].value())) is_badvalue = true;
		if (!std::isfinite(mResidual[eqn3].value())) is_badvalue = true;
	}
	for (std::size_t w = 0; w<nw; ++w)
	{
		for (std::size_t w = 0; w<nw; ++w)
		{
			std::size_t eqn1 = wtoc(w);
			if (H_control[w][nsche] == 0 || H_control[w][nsche] == 7)
			{
				mResidual[eqn1] = H_chedule[w][nsche] - Hwell_q[w].P;
			}
			else if (H_control[w][nsche] == 1)
			{
				mResidual[eqn1] = mDT*(H_chedule[w][nsche] - Hwell_q[w].qt);
			}
			else if (H_control[w][nsche] == 2)
			{
				mResidual[eqn1] = mDT*(H_chedule[w][nsche] - Hwell_q[w].qo);
			}
			else if (H_control[w][nsche] == 3)
			{
				mResidual[eqn1] = mDT*(H_chedule[w][nsche] - Hwell_q[w].qw);
			}
			else if (H_control[w][nsche] == 4 || H_control[w][nsche] == 8)
			{
				mResidual[eqn1] = mDT*(H_chedule[w][nsche] - Hwell_q[w].qg);
			}
			else if (H_control[w][nsche] == 5)
			{
				mResidual[eqn1] = mDT*(H_chedule[w][nsche] - Hwell_q[w].WCT);
			}
			else if (H_control[w][nsche] == 6)
			{
				mResidual[eqn1] = mDT*(H_chedule[w][nsche] - Hwell_q[w].GOR);
			}
		}
	}

	return is_badvalue;
}




void DiscreteProblem::extract_F_der(DiscreteProblem::StateVector &oldstate, DiscreteProblem::StateVector &state, double DT, DiscreteProblem::Sparse_matrix &M)
{
	std::vector<double> r(N_c);
	CSR f_y;
	bind_to_old_state(oldstate);
	discretize_m(state, DT);
	extract_R_Jm(r, f_y, 0);
	CSR_to_Sparse(f_y, M);

}

void
DiscreteProblem::extract_acc(DiscreteProblem::Sparse_matrix &M)
{
	std::size_t offset = 0;
	std::vector<double> r(N_tc, 0);
	CSR m;
	mAccum.extract_CSR_non_sysmetric(r, m.rowptr(), m.colind(), m.value(), N_tc);
	for (std::size_t i = 0; i< m.rowptr().size(); ++i) m.rowptr()[i] += offset;
	for (std::size_t i = 0; i< m.colind().size(); ++i) m.colind()[i] += offset;
	m.check_size();
	CSR_to_Sparse(m, M);
}


void
DiscreteProblem::extract_R_Jm(std::vector<double> &r, CSR &m, std::size_t offset)
{
	mResidual.extract_CSR_non_sysmetric(r, m.rowptr(), m.colind(), m.value(), nwc);
	for (std::size_t i = 0; i< m.rowptr().size(); ++i) m.rowptr()[i] += offset;
	for (std::size_t i = 0; i< m.colind().size(); ++i) m.colind()[i] += offset;
	m.check_size();
	//if (r.size() != mResidual.size()) std::cout << "BUG IN JACOBIAN\t zero row found" << r.size() << "!=" << mResidual.size() << std::endl;
}

void DiscreteProblem::CSR_to_Sparse(CSR & x, DiscreteProblem::Sparse_matrix &M)
{
	unsigned n = x.value().size();
	unsigned  old = 0;
	M.row.clear();
	M.col = x.colind();
	M.val = x.value();
	for (unsigned i = 0; i < x.rowptr().size(); i++)
	{
		if (x.rowptr()[i] != 0)
		{
			for (unsigned j = 0; j < (x.rowptr()[i] - old); j++)
			{
				M.row.push_back(i);
			}
			old = x.rowptr()[i];
		}
	}
}





double
DiscreteProblem::safeguard_MAC( double upd )
{
	if (std::abs(upd) > 0.2)
	{
		return upd / std::abs(upd) * 0.2;
	}
   else
      return upd;

}

std::ostream & operator << ( std::ostream & ostr, 
			     const DiscreteProblem::ConvergenceInfo & _out )
{
   ostr << _out.NormSat[0]  << "\t" 
	<< _out.NormSat[1]  << "\t" 
	<< _out.NormSat[2]  << "\t"
	<< _out.MatBal[0]  << "\t" 
	<< _out.MatBal[1]  << "\t" 
	<< _out.MatBal[2];

   return ostr;
}

