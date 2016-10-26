#pragma once

template< typename V, typename M >
void
DiscreteProblem::extract_R_J(V &r, M &m, std::size_t offset)
{
	mResidual.extract_CSR(r, m.rowptr(), m.colind(), m.value());
	for (std::size_t i = 0; i< m.rowptr().size(); ++i) m.rowptr()[i] += offset;
	for (std::size_t i = 0; i< m.colind().size(); ++i) m.colind()[i] += offset;
	m.check_size();
	if (r.size() != mResidual.size()) std::cout << "BUG IN JACOBIAN\t zero row found" << r.size() << "!=" << mResidual.size() << std::endl;
}

template< typename V >
void 
DiscreteProblem::extract_dR_dDT( V &r )
{
   r.resize( mFlow.size() );
   for ( std::size_t i=0 ; i< mFlow.size() ; ++ i ) r[i] = mFlow[i].value();
}


template<typename x>
void DiscreteProblem::load_hschedule(std::vector<std::vector<x>> &hwells, dlib::matrix<x> &s)
{
	hwells.resize(nw);
	for (unsigned i = 0; i < Np; i++)
	{
		hwells[i].resize(nct);
		for (unsigned j = 0; j <nct; j++)
		{
			hwells[i][j] = s(i*nct + j);
		}
	}
	for (unsigned i = Np; i < nw; i++)
	{
		hwells[i].resize(nct);
		for (unsigned j = 0; j < nct; j++)
		{
			hwells[i][j] = s(i*nct + j);
		}
	}
}



template< typename R >
void 
DiscreteProblem::update_state( StateVector &state, const R & update, bool do_safeguard )
{
   for ( std::size_t c = 0; c < mMesh.size_cells(); ++c )
   {
      if ( state[c].status == StatusID::OW )
      {
		 state[c].Po  += update[ eqnID(c,PhaseID::O) ];
		if (do_safeguard)
			state[c].Sw  += safeguard_MAC( update[ eqnID(c,PhaseID::W) ] );
		else
			state[c].Sw  += update[ eqnID(c,PhaseID::W) ];
		state[c].Rso += update[ eqnID(c,PhaseID::G) ];
		state[c].Sg   = 0.0;

		 double Pb =  mPropCalc.calc_bubble_point(  state[c].Rso );
		if ( state[c].Po <= Pb )
		 {
//	    std::cout << c << "\t" << "OW to OWG" << std::endl;
			state[c].Po.value()  = Pb;
			state[c].Sg.value()  = 0.05;
			state[c].Sg.make_independent( eqnID(c,PhaseID::G) );
			state[c].Rso.make_constant( );
			state[c].status = StatusID::OWG;
		}
      }
      else
      {
	 state[c].Po += update[eqnID(c,PhaseID::O)];
	 if (do_safeguard)
	 {
	    state[c].Sw += safeguard_MAC( update[ eqnID(c,PhaseID::W) ] );
	    state[c].Sg += safeguard_MAC( update[ eqnID(c,PhaseID::G) ] );
	 }
	 else
	 {
	    state[c].Sw += update[ eqnID(c,PhaseID::W) ];
	    state[c].Sg += update[ eqnID(c,PhaseID::G) ];
	 }
	 double RSO_sat = mPropCalc.calc_bubble_Rso( state[c].Po );
	 state[c].Rso = RSO_sat;

	 if ( state[c].Sg.value() <= 0.0 )
	 {
//	    std::cout << c << "\t" << "OWG to OW" << std::endl;
	    state[c].Sg = 0.0;
	    std::vector<std::size_t> neibs;
	    mMesh.cell_neighbors( c, neibs );
	    double RSO_AVG_NEIB = 0.0;
	    for ( std::size_t i=0; i<neibs.size(); ++i )
	       RSO_AVG_NEIB += state[neibs[i]].Rso.value();
	    RSO_AVG_NEIB /= neibs.size();
	    state[c].Rso.value() = ( RSO_AVG_NEIB < RSO_sat ? RSO_AVG_NEIB : RSO_sat );
	    state[c].Rso.make_independent( eqnID(c,PhaseID::G) );
	    state[c].status = StatusID::OW;
	 }

      } 
   }

   for (std::size_t w = 0; w < nw; ++w)
   {
	   std::size_t c = wtoc(w);
	   Hwell_q[w].P += update[c];
   }
}
   
