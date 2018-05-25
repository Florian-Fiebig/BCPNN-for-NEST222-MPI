/*
 *  mynest::aeif_cond_exp_multisynapse.cpp
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "aeif_cond_exp_multisynapse.h"
#include "nest_names.h"

#ifdef HAVE_GSL_1_11

#include "universal_data_logger_impl.h"

#include "exceptions.h"
#include "network.h"
#include "dict.h"
#include "integerdatum.h"
#include "doubledatum.h"
#include "dictutils.h"
#include "numerics.h"
#include <limits>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <cstdio>

/* ---------------------------------------------------------------- 
 * Recordables map
 * ---------------------------------------------------------------- */

nest::RecordablesMap<mynest::aeif_cond_exp_multisynapse> mynest::aeif_cond_exp_multisynapse::recordablesMap_;

namespace nest
{
  /*
   * template specialization must be placed in namespace
   *
   * Override the create() method with one call to RecordablesMap::insert_() 
   * for each quantity to be recorded.
   */
  template <>
  void RecordablesMap<mynest::aeif_cond_exp_multisynapse>::create()
  {
    // use standard names whereever you can for consistency!
    insert_(nest::names::V_m, &mynest::aeif_cond_exp_multisynapse::get_y_elem_<mynest::aeif_cond_exp_multisynapse::State_::V_M>);
    insert_(Name("g_AMPA"), &mynest::aeif_cond_exp_multisynapse::get_y_elem_<mynest::aeif_cond_exp_multisynapse::State_::G_AMPA>);
    insert_(Name("g_NMDA"), &mynest::aeif_cond_exp_multisynapse::get_y_elem_<mynest::aeif_cond_exp_multisynapse::State_::G_NMDA>);
    insert_(Name("g_NMDA_NEG"), &mynest::aeif_cond_exp_multisynapse::get_y_elem_<mynest::aeif_cond_exp_multisynapse::State_::G_NMDA_NEG>);
    insert_(Name("g_AMPA_NEG"), &mynest::aeif_cond_exp_multisynapse::get_y_elem_<mynest::aeif_cond_exp_multisynapse::State_::G_AMPA_NEG>);
    insert_(Name("g_GABA"), &mynest::aeif_cond_exp_multisynapse::get_y_elem_<mynest::aeif_cond_exp_multisynapse::State_::G_GABA>);
    insert_(nest::names::w, &mynest::aeif_cond_exp_multisynapse::get_y_elem_<mynest::aeif_cond_exp_multisynapse::State_::W>);
    insert_(Name("z_j"),    &mynest::aeif_cond_exp_multisynapse::get_y_elem_<mynest::aeif_cond_exp_multisynapse::State_::Z_J>);
    insert_(Name("e_j"),    &mynest::aeif_cond_exp_multisynapse::get_y_elem_<mynest::aeif_cond_exp_multisynapse::State_::E_J>);
    insert_(Name("p_j"),    &mynest::aeif_cond_exp_multisynapse::get_y_elem_<mynest::aeif_cond_exp_multisynapse::State_::P_J>);
    insert_(Name("bias"),   &mynest::aeif_cond_exp_multisynapse::get_bias_);
    insert_(Name("epsilon"),&mynest::aeif_cond_exp_multisynapse::get_epsilon_);
    insert_(Name("kappa"),  &mynest::aeif_cond_exp_multisynapse::get_kappa_);
    insert_(Name("I_AMPA"), &mynest::aeif_cond_exp_multisynapse::get_I_AMPA_); 
    insert_(Name("I_NMDA"), &mynest::aeif_cond_exp_multisynapse::get_I_NMDA_); 
    insert_(Name("I_AMPA_NEG"), &mynest::aeif_cond_exp_multisynapse::get_I_AMPA_NEG_); 
    insert_(Name("I_NMDA_NEG"), &mynest::aeif_cond_exp_multisynapse::get_I_NMDA_NEG_); 
    insert_(Name("I_GABA"), &mynest::aeif_cond_exp_multisynapse::get_I_GABA_); 
  }
}

extern "C"
int mynest::aeif_cond_exp_multisynapse_dynamics (double, const double y[], double f[], void* pnode)
{
  // a shorthand
  typedef mynest::aeif_cond_exp_multisynapse::State_ S;

  // get access to node so we can almost work as in a member function
  assert(pnode);
  mynest::aeif_cond_exp_multisynapse& node =  *(reinterpret_cast<mynest::aeif_cond_exp_multisynapse*>(pnode));

  // y[] here is---and must be---the state vector supplied by the integrator,
  // not the state vector in the node, node.S_.y[]. 
  
  // The following code is verbose for the sake of clarity. We assume that a
  // good compiler will optimize the verbosity away ...

  // This constant is used below as the largest admissible value for the exponential spike upstroke
  static const nest::double_t largest_exp=std::exp(10.);

  // shorthand for state variables
  const nest::double_t& V     = y[S::V_M  ];
  const nest::double_t& w     = y[S::W    ];

  //const nest::double_t I_AMPA = -y[S::G_AMPA] * ( V - node.P_.AMPA_E_rev );
  //const nest::double_t I_NMDA = -y[S::G_NMDA] * ( V - node.P_.NMDA_E_rev );
  //const nest::double_t I_NMDA_NEG = -y[S::G_NMDA_NEG] * ( V - node.P_.NMDA_NEG_E_rev );
  //const nest::double_t I_AMPA_NEG = -y[S::G_AMPA_NEG] * ( V - node.P_.AMPA_NEG_E_rev );
  //const nest::double_t I_GABA = -y[S::G_GABA] * ( V - node.P_.GABA_E_rev );
  //const nest::double_t I_bias = node.P_.gain * std::log(y[S::P_J]);

  //node.S_.I_AMPA_    = I_AMPA;
  //node.S_.I_NMDA_    = I_NMDA;
  //node.S_.I_NMDA_NEG_    = I_NMDA_NEG;
  //node.S_.I_AMPA_NEG_    = I_AMPA_NEG;
  //node.S_.I_GABA_    = I_GABA;
 
  node.S_.I_AMPA_     = -y[S::G_AMPA] * ( V - node.P_.AMPA_E_rev );
  node.S_.I_NMDA_     = -y[S::G_NMDA] * ( V - node.P_.NMDA_E_rev );
  node.S_.I_NMDA_NEG_ = -y[S::G_NMDA_NEG] * ( V - node.P_.NMDA_NEG_E_rev );
  node.S_.I_AMPA_NEG_ = -y[S::G_AMPA_NEG] * ( V - node.P_.AMPA_NEG_E_rev );
  node.S_.I_GABA_     = -y[S::G_GABA] * ( V - node.P_.GABA_E_rev );
  const nest::double_t I_bias = node.P_.gain * std::log(y[S::P_J]);

  //const nest::double_t I_syn =  I_AMPA + I_NMDA + I_NMDA_NEG + I_AMPA_NEG + I_GABA + node.B_.I_stim_ ;
  const nest::double_t I_syn =  node.S_.I_AMPA_ + node.S_.I_NMDA_ + node.S_.I_NMDA_NEG_ + node.S_.I_AMPA_NEG_ + node.S_.I_GABA_ + node.B_.I_stim_;
  // We pre-compute the argument of the exponential
  const nest::double_t exp_arg=(V - node.P_.V_th) / node.P_.Delta_T;
  // If the argument is too large, we clip it.
  const nest::double_t I_spike = (exp_arg>10.)? largest_exp : node.P_.Delta_T * std::exp(exp_arg);

  // dv/dt
  f[S::V_M  ] = ( -node.P_.g_L *( (V-node.P_.E_L) - I_spike) - w + node.P_.I_e + I_bias + I_syn) / node.P_.C_m;

  // Adaptation current w.
  f[S::W    ] = ( node.P_.a * (V - node.P_.E_L) - w ) / node.P_.tau_w;

  // Synapse dynamics dg_AMPA/dt
  f[ S::G_AMPA ] = -y[ S::G_AMPA ] / node.P_.AMPA_Tau_decay;
  f[ S::G_AMPA_NEG ] = -y[ S::G_AMPA_NEG ] / node.P_.AMPA_Tau_decay;

  // dg_NMDA/dt
  f[ S::G_NMDA ] = -y[ S::G_NMDA ] / node.P_.NMDA_Tau_decay;
  f[ S::G_NMDA_NEG ] = -y[ S::G_NMDA_NEG ] / node.P_.NMDA_Tau_decay;

  // dg_GABA_A/dt
  f[ S::G_GABA ] = -y[ S::G_GABA ] / node.P_.GABA_Tau_decay;

  f[S::Z_J] = (- y[S::Z_J] + node.P_.epsilon) / node.P_.tau_j;
  f[S::E_J] = (y[S::Z_J] - y[S::E_J]) / node.P_.tau_e;
  f[S::P_J] = node.P_.kappa * (y[S::E_J] - y[S::P_J]) / node.P_.tau_p;

  return GSL_SUCCESS;
}

/* ---------------------------------------------------------------- 
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */
    
mynest::aeif_cond_exp_multisynapse::Parameters_::Parameters_()
  : V_peak_    (   0.0 ), // mV 
    V_reset_   ( -60.0 ), // mV
    t_ref_     (   0.0 ), // ms
    g_L        (  30.0 ), // nS
    C_m        ( 281.0 ), // pF
    E_L        ( -70.6 ), // mV
    Delta_T    (   2.0 ), // mV
    tau_w      ( 144.0 ), // ms
    a          (   4.0 ), // nS
    b          (  80.5 ), // pA
    V_th       ( -50.4 ), // mVs
    I_e        (   0.0 ), // pA
    gsl_error_tol( 1e-6),
    AMPA_E_rev       (  0.0   ),   // mV
    AMPA_Tau_decay   (  2.0   ),   // ms
    NMDA_E_rev       (  0.0   ),   // mV
    NMDA_Tau_decay   (  150.0 ),   // ms
    NMDA_NEG_E_rev       (  -70.0   ),   // mV
    AMPA_NEG_E_rev       (  -70.0   ),   // mV
    GABA_E_rev       (-70.0    ),  // mV
    GABA_Tau_decay   (  2.0    ),  // ms
    tau_j   ( 10.0    ),  // ms
    tau_e   (100.0    ),  // ms
    tau_p   (1000.0   ),  // ms
    kappa   (1.0      ), // dopamine
    fmax    (20.0     ), 
    gain    (1.0     ), 
    bias    (0.0      ),
    epsilon (0.000001     )
{
    recordablesMap_.create();
}

mynest::aeif_cond_exp_multisynapse::State_::State_(const Parameters_ &p)
  : I_AMPA_(0.0),
  I_NMDA_(0.0),
  I_NMDA_NEG_(0.0),
  I_AMPA_NEG_(0.0),
  I_GABA_(0.0),
  r_(0),
  bias(0)
{
  y_[0] = p.E_L;
  for ( size_t i = 1; i <STATE_VEC_SIZE; ++i )
    y_[i] = 0;
  y_[Z_J] = 0.01;
  y_[E_J] = 0.01;
  y_[P_J] = 0.01;
}

mynest::aeif_cond_exp_multisynapse::State_::State_(const State_ &s)
  : I_AMPA_(  s.I_AMPA_  ),
    I_NMDA_(  s.I_NMDA_  ),
    I_NMDA_NEG_(  s.I_NMDA_NEG_  ),
    I_AMPA_NEG_(  s.I_AMPA_NEG_  ),
    I_GABA_(s.I_GABA_),
    r_(s.r_),
    bias(s.bias)
{
  for ( size_t i = 0; i < STATE_VEC_SIZE; ++i )
    y_[i] = s.y_[i];
}

mynest::aeif_cond_exp_multisynapse::State_& mynest::aeif_cond_exp_multisynapse::State_::operator=(const State_ &s)
{
  assert(this != &s);  // would be bad logical error in program
  
  for ( size_t i = 0; i < STATE_VEC_SIZE; ++i )
    y_[i] = s.y_[i];
  I_AMPA_    = s.I_AMPA_;
  I_NMDA_    = s.I_NMDA_;
  I_NMDA_NEG_    = s.I_NMDA_NEG_;
  I_AMPA_NEG_    = s.I_AMPA_NEG_;
  I_GABA_ = s.I_GABA_;
  r_ = s.r_;  
  bias = s.bias;
  return *this;
}

/* ---------------------------------------------------------------- 
 * Paramater and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void mynest::aeif_cond_exp_multisynapse::Parameters_::get(DictionaryDatum &d) const
{
  def<double>(d,nest::names::C_m,        C_m);
  def<double>(d,nest::names::V_th,       V_th);
  def<double>(d,nest::names::t_ref,      t_ref_);
  def<double>(d,nest::names::g_L,        g_L);
  def<double>(d,nest::names::E_L,        E_L); 
  def<double>(d,nest::names::V_reset,    V_reset_);
  def<double>(d,nest::names::a,          a);
  def<double>(d,nest::names::b,          b);
  def<double>(d,nest::names::Delta_T,    Delta_T);
  def<double>(d,nest::names::tau_w,      tau_w);
  def<double>(d,nest::names::I_e,        I_e);
  def<double>(d,nest::names::V_peak,     V_peak_);
  def<double>(d,nest::names::gsl_error_tol, gsl_error_tol);
  def<nest::double_t>(d, "AMPA_E_rev",         AMPA_E_rev);
  def<nest::double_t>(d, "AMPA_Tau_decay",     AMPA_Tau_decay);
  def<nest::double_t>(d, "NMDA_E_rev",         NMDA_E_rev);
  def<nest::double_t>(d, "NMDA_NEG_E_rev",     NMDA_NEG_E_rev);
  def<nest::double_t>(d, "AMPA_NEG_E_rev",     AMPA_NEG_E_rev);
  def<nest::double_t>(d, "NMDA_Tau_decay",     NMDA_Tau_decay);
  def<nest::double_t>(d, "GABA_E_rev",         GABA_E_rev);
  def<nest::double_t>(d, "GABA_Tau_decay",     GABA_Tau_decay);
  def<nest::double_t>(d, "tau_j",      tau_j);
  def<nest::double_t>(d, "tau_e",      tau_e);
  def<nest::double_t>(d, "tau_p",      tau_p);
  def<nest::double_t>(d, "kappa",      kappa);
  def<nest::double_t>(d, "bias",      bias);
  def<nest::double_t>(d, "gain",      gain);
  def<nest::double_t>(d, "fmax",      fmax);
  def<nest::double_t>(d, "epsilon",      epsilon);
}

void mynest::aeif_cond_exp_multisynapse::Parameters_::set(const DictionaryDatum &d)
{
  updateValue<double>(d,nest::names::V_th,    V_th);
  updateValue<double>(d,nest::names::V_peak,  V_peak_);
  updateValue<double>(d,nest::names::t_ref,   t_ref_);
  updateValue<double>(d,nest::names::E_L,     E_L);
  updateValue<double>(d,nest::names::V_reset, V_reset_);    
  updateValue<double>(d,nest::names::C_m, C_m);
  updateValue<double>(d,nest::names::g_L, g_L);
  updateValue<double>(d,nest::names::a,       a);
  updateValue<double>(d,nest::names::b,       b);
  updateValue<double>(d,nest::names::Delta_T, Delta_T);
  updateValue<double>(d,nest::names::tau_w,   tau_w);
  updateValue<double>(d,nest::names::I_e, I_e);
  updateValue<double>(d,nest::names::gsl_error_tol, gsl_error_tol);
  updateValue<nest::double_t>(d, "AMPA_E_rev",        AMPA_E_rev);
  updateValue<nest::double_t>(d, "AMPA_Tau_decay",    AMPA_Tau_decay);
  updateValue<nest::double_t>(d, "NMDA_E_rev",        NMDA_E_rev);
  updateValue<nest::double_t>(d, "NMDA_NEG_E_rev",        NMDA_NEG_E_rev);
  updateValue<nest::double_t>(d, "AMPA_NEG_E_rev",        AMPA_NEG_E_rev);
  updateValue<nest::double_t>(d, "NMDA_Tau_decay",    NMDA_Tau_decay);
  updateValue<nest::double_t>(d, "GABA_E_rev",     GABA_E_rev);
  updateValue<nest::double_t>(d, "GABA_Tau_decay", GABA_Tau_decay);
  updateValue<nest::double_t>(d, "tau_j",      tau_j);
  updateValue<nest::double_t>(d, "tau_e",      tau_e);
  updateValue<nest::double_t>(d, "tau_p",      tau_p);
  updateValue<nest::double_t>(d, "kappa",      kappa);
  updateValue<nest::double_t>(d, "gain",      gain);
  updateValue<nest::double_t>(d, "bias",      bias);
  updateValue<nest::double_t>(d, "fmax",      fmax);
  updateValue<nest::double_t>(d, "epsilon",      epsilon);

  if ( V_peak_ <= V_th )
    throw nest::BadProperty("V_peak must be larger than threshold.");

  if ( V_reset_ >= V_peak_ )
    throw nest::BadProperty("Ensure that: V_reset < V_peak .");
    
  if ( C_m <= 0 )
    throw nest::BadProperty("Ensure that C_m >0");
    
  if ( t_ref_ < 0 )
    throw nest::BadProperty("Ensure that t_ref >= 0");
      
  if ( AMPA_Tau_decay    <= 0 ||
       NMDA_Tau_decay    <= 0 ||
       GABA_Tau_decay <= 0 )
    throw nest::BadProperty("All time constants must be strictly positive.");

  if ( gsl_error_tol <= 0. )
    throw nest::BadProperty("The gsl_error_tol must be strictly positive.");
}

void mynest::aeif_cond_exp_multisynapse::State_::get(DictionaryDatum &d) const
{
  def<double>(d,nest::names::V_m,  y_[V_M]);
  def<double>(d,nest::names::w,    y_[W]);
  def<double>(d,Name("p_j"), y_[P_J]);
}

void mynest::aeif_cond_exp_multisynapse::State_::set(const DictionaryDatum &d, const Parameters_ &)
{
  updateValue<double>(d,nest::names::V_m,  y_[V_M]);
  updateValue<double>(d,nest::names::w,    y_[W]);
  updateValue<double>(d,Name("p_j"), y_[P_J]);
}

mynest::aeif_cond_exp_multisynapse::Buffers_::Buffers_(mynest::aeif_cond_exp_multisynapse &n)
  : logger_(n),
    s_(0),
    c_(0),
    e_(0)
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

mynest::aeif_cond_exp_multisynapse::Buffers_::Buffers_(const Buffers_ &, mynest::aeif_cond_exp_multisynapse &n)
  : logger_(n),
    s_(0),
    c_(0),
    e_(0)
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

/* ---------------------------------------------------------------- 
 * Default and copy constructor for node, and destructor
 * ---------------------------------------------------------------- */

mynest::aeif_cond_exp_multisynapse::aeif_cond_exp_multisynapse()
  : Archiving_Node(), 
    P_(), 
    S_(P_),
    B_(*this)
{
  recordablesMap_.create();
}

mynest::aeif_cond_exp_multisynapse::aeif_cond_exp_multisynapse(const aeif_cond_exp_multisynapse &n)
  : Archiving_Node(n), 
    P_(n.P_), 
    S_(n.S_),
    B_(n.B_, *this)
{
}

mynest::aeif_cond_exp_multisynapse::~aeif_cond_exp_multisynapse()
{
  // GSL structs may not have been allocated, so we need to protect destruction
  if ( B_.s_ ) gsl_odeiv_step_free(B_.s_);
  if ( B_.c_ ) gsl_odeiv_control_free(B_.c_);
  if ( B_.e_ ) gsl_odeiv_evolve_free(B_.e_);
}

/* ---------------------------------------------------------------- 
 * Node initialization functions
 * ---------------------------------------------------------------- */

void mynest::aeif_cond_exp_multisynapse::init_state_(const Node &proto)
{
  const mynest::aeif_cond_exp_multisynapse &pr = downcast<mynest::aeif_cond_exp_multisynapse>(proto);
  S_ = pr.S_;
}

void mynest::aeif_cond_exp_multisynapse::init_buffers_()
{
  B_.spikes_AMPA_.clear();       // includes resize
  B_.spikes_NMDA_.clear();       // includes resize
  B_.spikes_NMDA_NEG_.clear();       // includes resize
  B_.spikes_AMPA_NEG_.clear();       // includes resize
  B_.spikes_GABA_.clear();    // includes resize
  B_.currents_.clear();           // includes resize
  Archiving_Node::clear_history();

  B_.logger_.reset();

  B_.step_ = nest::Time::get_resolution().get_ms();

  // We must integrate this model with high-precision to obtain decent results
  B_.IntegrationStep_ = std::min(0.01, B_.step_);

  static const gsl_odeiv_step_type* T1 = gsl_odeiv_step_rkf45;
  
  if ( B_.s_ == 0 )
    B_.s_ = gsl_odeiv_step_alloc (T1, State_::STATE_VEC_SIZE);
  else 
    gsl_odeiv_step_reset(B_.s_);
    
  if ( B_.c_ == 0 )  
    B_.c_ = gsl_odeiv_control_yp_new (P_.gsl_error_tol,P_.gsl_error_tol);
  else
    gsl_odeiv_control_init(B_.c_, P_.gsl_error_tol, P_.gsl_error_tol, 0.0, 1.0);
    
  if ( B_.e_ == 0 )  
    B_.e_ = gsl_odeiv_evolve_alloc(State_::STATE_VEC_SIZE);
  else 
    gsl_odeiv_evolve_reset(B_.e_);
  
  B_.sys_.function  = mynest::aeif_cond_exp_multisynapse_dynamics; 
  B_.sys_.jacobian  = NULL;
  B_.sys_.dimension = State_::STATE_VEC_SIZE;
  B_.sys_.params    = reinterpret_cast<void*>(this);
  B_.I_stim_ = 0.0;
}

void mynest::aeif_cond_exp_multisynapse::calibrate()
{
  B_.logger_.init();  // ensures initialization in case mm connected after Simulate
  V_.RefractoryCounts_ = nest::Time(nest::Time::ms(P_.t_ref_)).get_steps();
  assert(V_.RefractoryCounts_ >= 0);  // since t_ref_ >= 0, this can only fail in error
}

/* ---------------------------------------------------------------- 
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

void mynest::aeif_cond_exp_multisynapse::update(const nest::Time &origin, const nest::long_t from, const nest::long_t to)
{
  assert ( to >= 0 && (nest::delay) from < nest::Scheduler::get_min_delay() );
  assert ( from < to );
  assert ( State_::V_M == 0 );

  for ( nest::long_t lag = from; lag < to; ++lag )
  {
    double t = 0.0;

    if ( S_.r_ > 0 )
      --S_.r_;

    // numerical integration with adaptive step size control:
    // ------------------------------------------------------
    // gsl_odeiv_evolve_apply performs only a single numerical
    // integration step, starting from t and bounded by step;
    // the while-loop ensures integration over the whole simulation
    // step (0, step] if more than one integration step is needed due
    // to a small integration step size;
    // note that (t+IntegrationStep > step) leads to integration over
    // (t, step] and afterwards setting t to step, but it does not
    // enforce setting IntegrationStep to step-t
    while ( t < B_.step_ )
    {
      const int status = gsl_odeiv_evolve_apply(B_.e_, B_.c_, B_.s_, 
						&B_.sys_,             // system of ODE
						&t,                   // from t
						B_.step_,             // to t <= step
						&B_.IntegrationStep_, // integration step size
						S_.y_);               // neuronal state
      
      if ( status != GSL_SUCCESS )
        throw nest::GSLSolverFailure(get_name(), status);

      // check for unreasonable values; we allow V_M to explode
      if ( S_.y_[State_::V_M] < -1e3 ||
	   S_.y_[State_::W  ] <    -1e6 || S_.y_[State_::W] > 1e6    )
	throw nest::NumericalInstability(get_name());
      
      // spikes are handled inside the while-loop
      // due to spike-driven adaptation
      if ( S_.r_ > 0 )
	S_.y_[State_::V_M] = P_.V_reset_;
      else if ( S_.y_[State_::V_M] >= P_.V_peak_ )
	{
	  S_.y_[State_::V_M]  = P_.V_reset_;
	  S_.y_[State_::W]   += P_.b; // spike-driven adaptation
	  S_.r_               = V_.RefractoryCounts_;
          S_.y_[State_::Z_J] += (1000.0/(P_.fmax*B_.step_) - S_.y_[State_::Z_J] + P_.epsilon) * B_.step_ / P_.tau_j;
          S_.y_[State_::E_J] += (S_.y_[State_::Z_J] - S_.y_[State_::E_J]) * B_.step_ / P_.tau_e;
          S_.y_[State_::P_J] += P_.kappa * (S_.y_[State_::E_J] - S_.y_[State_::P_J]) * B_.step_ / P_.tau_p;
	  
	  set_spiketime(nest::Time::step(origin.get_steps() + lag + 1));
	  nest::SpikeEvent se;
	  network()->send(*this, se, lag);
	}
    }  
    S_.y_[State_::G_AMPA]    += B_.spikes_AMPA_.get_value(lag);
    S_.y_[State_::G_NMDA]    += B_.spikes_NMDA_.get_value(lag);
    S_.y_[State_::G_NMDA_NEG]    += B_.spikes_NMDA_NEG_.get_value(lag);
    S_.y_[State_::G_AMPA_NEG]    += B_.spikes_AMPA_NEG_.get_value(lag);
    S_.y_[State_::G_GABA]    += B_.spikes_GABA_.get_value(lag);

    S_.bias = P_.gain * std::log(S_.y_[State_::P_J]);
    
    // set new input current
    B_.I_stim_ = B_.currents_.get_value(lag);
    
    // log state data
    B_.logger_.record_data(origin.get_steps() + lag);
  }
}
  
void mynest::aeif_cond_exp_multisynapse::handle(nest::SpikeEvent &e)
{
  assert ( e.get_delay() > 0 );
  assert(0 <= e.get_rport() && e.get_rport() < SUP_SPIKE_RECEPTOR - MIN_SPIKE_RECEPTOR);

  // If AMPA
  if (e.get_rport() == AMPA - MIN_SPIKE_RECEPTOR && e.get_weight() >= 0)
      B_.spikes_AMPA_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()), e.get_weight() * e.get_multiplicity() );

  else if (e.get_rport() == AMPA - MIN_SPIKE_RECEPTOR && e.get_weight() < 0)
      B_.spikes_AMPA_NEG_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()), -e.get_weight() * e.get_multiplicity() );

  // If NMDA
  else if (e.get_rport() == NMDA - MIN_SPIKE_RECEPTOR && e.get_weight() >= 0)
           B_.spikes_NMDA_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()), e.get_weight() * e.get_multiplicity() );

  else if (e.get_rport() == NMDA - MIN_SPIKE_RECEPTOR && e.get_weight() < 0)
           B_.spikes_NMDA_NEG_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()), -e.get_weight() * e.get_multiplicity() );  

  // If GABA_A
  else if (e.get_rport() == GABA - MIN_SPIKE_RECEPTOR)
           B_.spikes_GABA_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()), -e.get_weight() * e.get_multiplicity() );
}

void mynest::aeif_cond_exp_multisynapse::handle(nest::CurrentEvent &e)
{
  assert ( e.get_delay() > 0 );

  const nest::double_t c=e.get_current();
  const nest::double_t w=e.get_weight();

  // add weighted current; HEP 2002-10-04
  B_.currents_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()), 
			 w*c);
  assert(0 <= e.get_rport() && e.get_rport() < SUP_CURR_RECEPTOR - MIN_CURR_RECEPTOR);
}

void mynest::aeif_cond_exp_multisynapse::handle(nest::DataLoggingRequest &e)
{
  B_.logger_.handle(e);
}

#endif // HAVE_GSL_1_11
