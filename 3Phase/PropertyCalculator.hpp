#pragma once
#include <fstream>
#include <vector>
#include "adetl/scalars/ADscalar.hpp"
#include "func_table.hpp"
#include <iostream>


class PropertyCalculator
{
private:
   typedef adetl::ADscalar< >         scalar;
   typedef fastl::sorted_key<double>  TableKey;
   typedef fastl::interpolatable_range<double,
				       fastl::forward_linear,
				       fastl::centered_linear,
				       fastl::backward_linear>  TableValue;

public:
   enum Phase { O=0, W=1, G=2 };
   enum StatusID { OWG, OW };
   typedef struct
   {
      // W , O, G
      StatusID   status;
      scalar S[3];
      scalar P[3];
      scalar Pb;
      scalar Kr[3];
      scalar Rso;
      scalar b[3];
      scalar bmu[3];
      scalar Rho[3];
      scalar phi;
   }                                            CellProps;

public:
   PropertyCalculator() : CR(3.4e-4), PREF_ROCK(4000.0), 
			  RHO_W_S(63.02), RHO_O_S(45.0), RHO_G_S(0.0702),
			  PREF_W(3600.0), BW_REF(1.00341), CW(4.0e-4), MUW_REF(0.52341), CVW(1.2e-6)
   {

      read_PVDG_table( );
      read_PVTO1_table( );
      read_PVTO2_table( );
      read_SWOF_table( );
      read_SGOF_table( );
	  std::size_t i = SW_key(Swr);
	  KROW_val(i, SW_key, Swr, Krocw);
	  Pb_ini.status = OWG;
	  calc_oil_props(Pb_ini);
	  calc_gas_props(Pb_ini);
   };

   double calc_bubble_point( const scalar & Rso )
   {
      const scalar RSO_SCALED = Rso/178.1076;
      std::size_t i = RS_key( RSO_SCALED );
      PO_val( i, RS_key, RSO_SCALED, mTmp[0] );
      return mTmp[0].value();
   }

   double calc_bubble_Rso( const scalar & Pb )
   {
      std::size_t i = PO_key( Pb );
      RS_val( i, PO_key, Pb, mTmp[0] );
      return mTmp[0].value()*178.1076;
   }

   template< typename T >
   void calculate( double PHI, const T & _state, CellProps & p )
   {
      p.status = _state.status;
      p.P[O]   = _state.Po;
      p.S[W]   = _state.Sw;
      if ( p.status == OWG )
      {
	 p.S[G] = _state.Sg;
      }
      else
      {
		p.S[G] = 0.0;
		p.Rso  = _state.Rso;
      }
      p.S[O] = 1.0 - p.S[W] - p.S[G];

      calc_rock_fluid( p );
      calc_water_props( p );
      calc_gas_props( p );
      calc_oil_props( p );
      calc_phi( PHI, p );
   }


   template< typename T >
   void calculat_initial(double PHI, const T & _state, CellProps & p)
   {
	   scalar x;
	   p.P[O] = _state.Po;
	   p.P[G] = _state.Po;
	   p.S[W] = _state.Sw;
	   if (p.P[0].value()<Pb_ini.P[O].value())
	   {
		   p.status = OWG;
		   calc_phi(PHI, p);
		   calc_phi(PHI, Pb_ini);
		   calc_gas_props(p);
		   calc_oil_props(p);
		   x= ((Pb_ini.phi*(1 - p.S[W])*Pb_ini.Rso*Pb_ini.b[O]) - (p.phi*(1 - p.S[W])*p.Rso*p.b[O])) / p.phi / (p.b[G] - p.Rso*p.b[O]);
		   p.S[G] = x.value();
	   }
	   else
	   {
		   p.status = OW;
		   p.S[G] = 0.0;
	   }
   }

private:

   void calc_phi( double PHI, CellProps &p ) const
   {
      const scalar X = CR * ( p.P[O] - PREF_ROCK );
      p.phi = PHI * ( 1.0 + X + 0.5 * X * X );
   }

   void calc_rock_fluid( CellProps &p )
   {
      // mTmp[0]  1   2    3
      // pcow, krow, pcog, krog;

      std::size_t i = SW_key( p.S[W] );
      KRW_val(  i, SW_key, p.S[W], p.Kr[W] );
      KROW_val( i, SW_key, p.S[W], mTmp[1] );
	  PCOW_val(i, SW_key, p.S[W], mTmp[0]);


      p.P[W] = p.P[O] - mTmp[0];

      i = SG_key( p.S[G] );
      KRG_val(  i, SG_key, p.S[G], p.Kr[G] );
      KROG_val( i, SG_key, p.S[G], mTmp[3] );
      PCOG_val(  i, SG_key, p.S[G], mTmp[2] );
      p.P[G] = p.P[O] + mTmp[2];
	  //p.Kr[O] = (mTmp[1] + p.Kr[W]) * (mTmp[3] + p.Kr[G]) - p.Kr[W] - p.Kr[G];
	  p.Kr[O] = Krocw.value()*(mTmp[1] / Krocw.value() + p.Kr[W]) * (mTmp[3] / Krocw.value() + p.Kr[G])
		  - p.Kr[W] - p.Kr[G];
   }

   void calc_water_props( CellProps &p )
   { 
      mTmp[0]    = CW*(p.P[W]-PREF_W);
      p.b[W]     = ( 1.0 + mTmp[0] + 0.5 * mTmp[0] * mTmp[0] ) / BW_REF;
      mTmp[0]    = ( CW - CVW ) * ( p.P[W] - PREF_W );
      p.bmu[W]   = ( 1.0 + mTmp[0] + 0.5 * mTmp[0] * mTmp[0] ) / ( BW_REF * MUW_REF );
      p.Rho[W]   = RHO_W_S * p.b[W];
   }  

   void calc_gas_props( CellProps &p ) const
   { 
      std::size_t i = PG_key( p.P[G] );
      bg_val(   i, PG_key, p.P[G], p.b[G] );
	  p.b[G] *= 178.1076;
      bmug_val( i, PG_key, p.P[G], p.bmu[G] );
      p.bmu[G] *= 178.1076;
      p.Rho[G] = RHO_G_S * p.b[G];
   }  

   void calc_oil_props( CellProps &p )
   { 
      if ( p.status == OWG )
      {
	 std::size_t i = PO_key( p.P[O] );
	 p.Pb = p.P[O].value();
	 RS_val(   i, PO_key, p.P[O], p.Rso );
	 bo_val(   i, PO_key, p.P[O], p.b[O] );
	 bmuo_val( i, PO_key, p.P[O], p.bmu[O] );
         p.Rso *= 178.1076;
      }
      else
      {
	 // mTmp[]    0       1        2   3     4
	 //         bo_bub, bmuo_bub, DPO, DBO, DBMUO
         const scalar RSO = p.Rso / 178.1076;
		std::size_t i = RS_key( RSO );
		PO_val(     i, RS_key, RSO, p.Pb );
		bo_val(     i, RS_key, RSO, mTmp[0] );
		bmuo_val(   i, RS_key, RSO, mTmp[1] );
		mTmp[2] = p.P[O] / p.Pb;
		i = DPO_key( mTmp[2] );
		DBO_val(   i, DPO_key, mTmp[2], mTmp[3] );
		DBMUO_val( i, DPO_key, mTmp[2], mTmp[4] );
		p.b[O]   = mTmp[3] * mTmp[0];
		p.bmu[O] = mTmp[4] * mTmp[1];
      }
      p.Rho[O] = (RHO_O_S + RHO_G_S * p.Rso ) * p.b[O];
  }  

   void
   read_PVDG_table( )
   {
      // read PVDG gas table
      double tmp1, tmp2, tmp3;
      std::ifstream strm ( "./input/PVDG.dat" );
      while ( strm.good( ) )
      {
	 if (strm >> tmp1) PG_key.push_back( tmp1 );
	 if (strm >> tmp2) bg_val.push_back( 1.0/tmp2 );
	 if (strm >> tmp3) bmug_val.push_back( 1.0/(tmp2*tmp3) );
      }
      strm.close(); 
   }

   void
   read_PVTO1_table( )
   {
      // read PVTO oil tables 
      double tmp1, tmp2, tmp3, tmp4;
      std::ifstream strm ( "./input/PVTO_1.dat" );
      while ( strm.good( ) )
      {
	 if ( strm >> tmp1 )
	 {
	    RS_key.push_back( tmp1 );
	    RS_val.push_back( tmp1 );
	 }

	 if ( strm >> tmp2 )
	 {
	    PO_key.push_back( tmp2 );
	    PO_val.push_back( tmp2 );
	 }

	 if ( strm >> tmp3 ) bo_val.push_back( 1.0/tmp3 );

	 if (strm >> tmp4)
	 {
		 bmuo_val.push_back(1.0 / (tmp3*tmp4));
	 }
      }
      strm.close();
   }

   void 
   read_PVTO2_table ( )
   {
      double PB1, PB2, B1, B2, BMU1, BMU2;
      std::ifstream strm ( "./input/PVTO_2.dat" );
      if ( strm.good( ) )
      {
		strm >> PB1;
		Pb_ini.P[O]= PB1;
		strm >> B1;
		B1 = 1.0/B1;
		strm >> BMU1;
		BMU1 = B1 / BMU1;
      }
      DPO_key.push_back( 1.0 );
      DBO_val.push_back( 1.0);
      DBMUO_val.push_back( 1.0 );
      while ( strm.good( ) )
      {
	 if ( strm >> PB2 ) DPO_key.push_back( PB2 / PB1 );

	 if (strm >> B2)
	 {
	    B2 = 1.0 / B2;
	    DBO_val.push_back( B2 / B1 );
	 }

	 if ( strm >> BMU2 )
	 {
	    BMU2 = B2/BMU2;
	    DBMUO_val.push_back( BMU2 / BMU1 );
	 }
      }
      strm.close();      
   }

   void
   read_SWOF_table( )
   {
      // read SWOF table
      double sw,krw,kro,pcow;
      std::ifstream strm ( "./input/SWOF.dat" );
      while ( strm.good( ) )
      {
		  if (strm >> sw)
		  {
			  SW_key.push_back(sw);

		  }
		if (strm >>krw) KRW_val.push_back( krw );
		if (strm >> kro) KROW_val.push_back( kro );
		if (strm >> pcow) PCOW_val.push_back(pcow);
		if (krw == 0) Swr.value() = sw;
		//std::cout << sw << "  " << krw << "  " << kro << "  " << pcow << "  " << std::endl;
      }
      strm.close(); 
   }

   void
   read_SGOF_table( )
   {
      // read SGOF table
      double tmp;
      std::ifstream strm ( "./input/SGOF.dat" );
      while ( strm.good( ) )
      {
	 if (strm >> tmp) SG_key.push_back( tmp );
	 if (strm >> tmp) KRG_val.push_back( tmp );
	 if (strm >> tmp) KROG_val.push_back( tmp );
	 if (strm >> tmp) PCOG_val.push_back( tmp );
      }
      strm.close(); 
   }

private:
	scalar  Swr;
	scalar Krocw;
   const double CR;
   const double PREF_ROCK;
   const double RHO_W_S;
   const double RHO_O_S;
   const double RHO_G_S;
   const double PREF_W;
   const double BW_REF;
   const double CW;
   const double MUW_REF;
   const double CVW;
   CellProps Pb_ini;
   TableKey     PG_key;
   TableValue   bg_val;
   TableValue   bmug_val;
   TableKey     RS_key;
   TableValue   RS_val;
   TableKey     PO_key;
   TableValue   PO_val;
   TableValue   bo_val;
   TableValue   bmuo_val;
   TableKey     DPO_key;
   TableValue   DBO_val;
   TableValue   DBMUO_val;  
   TableKey     SW_key;
   TableValue   KRW_val;
   TableValue   KROW_val;
   TableValue   PCOW_val;
   TableKey     SG_key;
   TableValue   KRG_val;
   TableValue   KROG_val;
   TableValue   PCOG_val;
   scalar       mTmp[5];
};
