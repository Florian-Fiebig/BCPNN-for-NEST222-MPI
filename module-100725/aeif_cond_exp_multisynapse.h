/*
 *  aeif_cond_exp_multisynapse.h
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

#ifndef AEIF_COND_EXP_MULTISYNAPSE_H
#define AEIF_COND_EXP_MULTISYNAPSE_H

#include "config.h"

#ifdef HAVE_GSL_1_11

#include "nest.h"
#include "event.h"
#include "archiving_node.h"
#include "ring_buffer.h"
#include "connection.h"
#include "universal_data_logger.h"
#include "recordables_map.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

/* BeginDocumentation
Name: aeif_cond_exp_multisynapse - Conductance based exponential integrate-and-fire neuron model according to Brette and Gerstner (2005).

Description: 

aeif_cond_exp_multisynapse is the adaptive exponential integrate and fire neuron
according to Brette and Gerstner (2005), with post-synaptic
conductances in the form of truncated exponentials.

This implementation uses the embedded 4th order Runge-Kutta-Fehlberg
solver with adaptive stepsize to integrate the differential equation.

The membrane potential is given by the following differential equation:
C dV/dt= -g_L(V-E_L)+g_L*Delta_T*exp((V-V_T)/Delta_T)-g_e(t)(V-E_e) -g_i(t)(V-E_i)-w +I_e

and

tau_w * dw/dt= a(V-E_L) -W


Note that the spike detection threshold V_peak is automatically set to
V_th+10 mV to avoid numerical instabilites that may result from
setting V_peak too high.

Parameters: 
The following parameters can be set in the status dictionary.

Dynamic state variables:
  V_m        double - Membrane potential in mV
  g_ex       double - Excitatory synaptic conductance in nS.
  g_in       double - Inhibitory synaptic conductance in nS.
  w          double - Spike-adaptation current in pA.

Membrane Parameters:
  C_m        double - Capacity of the membrane in pF
  t_ref      double - Duration of refractory period in ms. 
  V_reset    double - Reset value for V_m after a spike. In mV.
  E_L        double - Leak reversal potential in mV. 
  g_L        double - Leak conductance in nS.
  I_e        double - Constant external input current in pA.

Spike adaptation parameters:
  a          double - Subthreshold adaptation in nS.
  b          double - Spike-triggered adaptation in pA.
  Delta_T    double - Slope factor in mV
  tau_w      double - Adaptation time constant in ms
  V_t        double - Spike initiation threshold in mV
  V_peak     double - Spike detection threshold in mV.

Synaptic parameters
  E_ex       double - Excitatory reversal potential in mV.
  tau_syn_ex double - Rise time of excitatory synaptic conductance in ms (exp function).
  E_in       double - Inhibitory reversal potential in mV.
  tau_syn_in double - Rise time of the inhibitory synaptic conductance in ms (exp function).

Integration parameters
  gsl_error_tol  double - This parameter controls the admissible error of the GSL integrator.
                          Reduce it if NEST complains about numerical instabilities.

Author: Adapted from aeif_cond_alpha by Lyle Muller

Sends: SpikeEvent

Receives: SpikeEvent, CurrentEvent, DataLoggingRequest

References: Brette R and Gerstner W (2005) Adaptive Exponential Integrate-and-Fire Model as 
            an Effective Description of Neuronal Activity. J Neurophysiol 94:3637-3642

SeeAlso: iaf_cond_exp, aeif_cond_alpha
*/

namespace mynest
{
  /**
   * Function computing right-hand side of ODE for GSL solver.
   * @note Must be declared here so we can befriend it in class.
   * @note Must have C-linkage for passing to GSL. Internally, it is
   *       a first-class C++ function, but cannot be a member function
   *       because of the C-linkage.
   * @note No point in declaring it inline, since it is called
   *       through a function pointer.
   * @param void* Pointer to model neuron instance.
   */
  extern "C"
  int aeif_cond_exp_multisynapse_dynamics (double, const double*, double*, void*);

  class aeif_cond_exp_multisynapse:
    public nest::Archiving_Node
  {
    
  public:        
    
    aeif_cond_exp_multisynapse();
    aeif_cond_exp_multisynapse(const aeif_cond_exp_multisynapse&);
    ~aeif_cond_exp_multisynapse();

    /**
     * Import sets of overloaded virtual functions.
     * We need to explicitly include sets of overloaded
     * virtual functions into the current scope.
     * According to the SUN C++ FAQ, this is the correct
     * way of doing things, although all other compilers
     * happily live without.
     */
    using nest::Node::connect_sender;
    using nest::Node::handle;

    nest::port check_connection(nest::Connection &, nest::port);
    
    void handle(nest::SpikeEvent &);
    void handle(nest::CurrentEvent &);
    void handle(nest::DataLoggingRequest &);
    
    nest::port connect_sender(nest::SpikeEvent &, nest::port);
    nest::port connect_sender(nest::CurrentEvent &, nest::port);
    nest::port connect_sender(nest::DataLoggingRequest &, nest::port);

    void get_status(DictionaryDatum &) const;
    void set_status(const DictionaryDatum &);
    
  private:
    
    void init_state_(const Node &proto);
    void init_buffers_();
    void calibrate();
    void update(const nest::Time &, const nest::long_t, const nest::long_t);

    static const nest::port MIN_SPIKE_RECEPTOR = 1;
    enum SpikeSynapseTypes { AMPA=MIN_SPIKE_RECEPTOR, NMDA, GABA, SUP_SPIKE_RECEPTOR };
    static const nest::size_t NUM_SPIKE_RECEPTORS = SUP_SPIKE_RECEPTOR - MIN_SPIKE_RECEPTOR;
    static const nest::port MIN_CURR_RECEPTOR = SUP_SPIKE_RECEPTOR;
    enum CurrentSynapseTypes { CURR = MIN_CURR_RECEPTOR, SUP_CURR_RECEPTOR };
    static const nest::size_t NUM_CURR_RECEPTORS = SUP_CURR_RECEPTOR - MIN_CURR_RECEPTOR;

    // END Boilerplate function declarations ----------------------------

    // Friends --------------------------------------------------------

    // make dynamics function quasi-member
    friend int mynest::aeif_cond_exp_multisynapse_dynamics(double, const double*, double*, void*);

    // The next two classes need to be friends to access the State_ class/member
    friend class nest::RecordablesMap<aeif_cond_exp_multisynapse>;
    friend class nest::UniversalDataLogger<aeif_cond_exp_multisynapse>;

  private:
    // ---------------------------------------------------------------- 

    //! Independent parameters
    struct Parameters_
    {
      nest::double_t V_peak_;     //!< Spike detection threshold in mV
      nest::double_t V_reset_;    //!< Reset Potential in mV
      nest::double_t t_ref_;      //!< Refractory period in ms

      nest::double_t g_L;         //!< Leak Conductance in nS
      nest::double_t C_m;         //!< Membrane Capacitance in pF
      nest::double_t E_L;         //!< Leak reversal Potential (aka resting potential) in mV
      nest::double_t Delta_T;     //!< Slope faktor in ms.
      nest::double_t tau_w;       //!< adaptation time-constant in ms.
      nest::double_t a;           //!< Subthreshold adaptation in nS.
      nest::double_t b;           //!< Spike-triggered adaptation in pA
      nest::double_t V_th;        //!< Spike threshold in mV.
      nest::double_t t_ref;       //!< Refractory period in ms.
      nest::double_t I_e;         //!< Intrinsic current in pA.

      nest::double_t AMPA_E_rev;        //!< AMPA reversal Potential in mV
      nest::double_t AMPA_Tau_decay;    //!< Synaptic Time Constant AMPA Synapse in ms
      nest::double_t NMDA_E_rev;        //!< NMDA reversal Potential in mV
      nest::double_t NMDA_NEG_E_rev;        //!< NMDA reversal Potential in mV
      nest::double_t AMPA_NEG_E_rev;        //!< NMDA reversal Potential in mV
      nest::double_t NMDA_Tau_decay;    //!< Synaptic Time Constant NMDA Synapse in ms
      nest::double_t GABA_E_rev;        //!< GABAA reversal Potential in mV
      nest::double_t GABA_Tau_decay;    //!< Rise Time Constant GABAA Synapse in ms

      nest::double_t tau_j;
      nest::double_t tau_e;
      nest::double_t tau_p;
      nest::double_t bias;
      nest::double_t fmax;
      nest::double_t gain;
      nest::double_t epsilon;
      nest::double_t kappa;
  
      nest::double_t gsl_error_tol;   //!< error bound for GSL integrator
  
      Parameters_();  //!< Sets default parameter values

      void get(DictionaryDatum &) const;  //!< Store current values in dictionary
      void set(const DictionaryDatum &);  //!< Set values from dicitonary
    };

  public:
    // ---------------------------------------------------------------- 

    /**
     * State variables of the model.
     * @note Copy constructor and assignment operator required because
     *       of C-style array.
     */
    struct State_
    {
      /**
       * Enumeration identifying elements in state array State_::y_.
       * The state vector must be passed to GSL as a C array. This enum
       * identifies the elements of the vector. It must be public to be
       * accessible from the iteration function.
       */  
      enum StateVecElems
      {
	V_M   = 0,
	W        ,  // 3
	G_AMPA   ,
	G_NMDA   ,
	G_NMDA_NEG   ,
	G_AMPA_NEG   ,
	G_GABA   ,  
        Z_J      ,    
        E_J      ,    
        P_J      ,      
	STATE_VEC_SIZE
      };

      nest::double_t y_[STATE_VEC_SIZE];  //!< neuron state, must be C-array for GSL solver
      nest::int_t    r_;           //!< number of refractory steps remaining

      nest::double_t I_AMPA_;    //!< AMPA current; member only to allow recording
      nest::double_t I_NMDA_;    //!< NMDA current; member only to allow recording
      nest::double_t I_NMDA_NEG_;    //!< NMDA current; member only to allow recording
      nest::double_t I_AMPA_NEG_;    //!< NMDA current; member only to allow recording
      nest::double_t I_GABA_; //!< GABAA current; member only to allow recording

      nest::double_t bias;
      nest::double_t epsilon;
      nest::double_t kappa;

      State_(const Parameters_ &);  //!< Default initialization
      State_(const State_ &);
      State_& operator = (const State_ &);

      void get(DictionaryDatum &) const;
      void set(const DictionaryDatum &, const Parameters_ &);
    };

    // ---------------------------------------------------------------- 

    /**
     * Buffers of the model.
     */
    struct Buffers_
    {
      Buffers_(aeif_cond_exp_multisynapse &);                    //!<Sets buffer pointers to 0
      Buffers_(const Buffers_ &, aeif_cond_exp_multisynapse &);  //!<Sets buffer pointers to 0

      //! Logger for all analog data
      nest::UniversalDataLogger<aeif_cond_exp_multisynapse> logger_;

      /** buffers and sums up incoming spikes/currents */
      nest::RingBuffer spikes_AMPA_;
      nest::RingBuffer spikes_NMDA_;
      nest::RingBuffer spikes_NMDA_NEG_;
      nest::RingBuffer spikes_AMPA_NEG_;
      nest::RingBuffer spikes_GABA_;
      nest::RingBuffer currents_;

      /** GSL ODE stuff */
      gsl_odeiv_step*    s_;    //!< stepping function
      gsl_odeiv_control* c_;    //!< adaptive stepsize control function
      gsl_odeiv_evolve*  e_;    //!< evolution function
      gsl_odeiv_system   sys_;  //!< struct describing system
      
      // IntergrationStep_ should be reset with the neuron on ResetNetwork,
      // but remain unchanged during calibration. Since it is initialized with
      // step_, and the resolution cannot change after nodes have been created,
      // it is safe to place both here.
      nest::double_t step_;             //!< step size in ms
      double   IntegrationStep_;  //!< current integration time step, updated by GSL

      /** 
       * Input current injected by CurrentEvent.
       * This variable is used to transport the current applied into the
       * _dynamics function computing the derivative of the state vector.
       * It must be a part of Buffers_, since it is initialized once before
       * the first simulation, but not modified before later Simulate calls.
       */
      nest::double_t I_stim_;
    };

    // ---------------------------------------------------------------- 

    /**
     * Internal variables of the model.
     */
    struct Variables_
    {
      nest::double_t PSConInit_AMPA;
      nest::double_t PSConInit_NMDA;
      nest::double_t PSConInit_NMDA_NEG;
      nest::double_t PSConInit_AMPA_NEG;
      nest::double_t PSConInit_GABA;
      nest::int_t RefractoryCounts_;
    };

    // Access functions for UniversalDataLogger -------------------------------
    
    //! Read out state vector elements, used by UniversalDataLogger
    template <State_::StateVecElems elem>
    nest::double_t get_y_elem_() const { return S_.y_[elem]; }
    nest::double_t get_I_AMPA_() const { return S_.I_AMPA_;  }
    nest::double_t get_I_NMDA_() const { return S_.I_NMDA_;  }
    nest::double_t get_I_NMDA_NEG_() const { return S_.I_NMDA_NEG_;  }
    nest::double_t get_I_AMPA_NEG_() const { return S_.I_AMPA_NEG_;  }
    nest::double_t get_I_GABA_() const { return S_.I_GABA_;  }
    nest::double_t get_bias_() const { return S_.bias; }
    nest::double_t get_epsilon_() const { return S_.epsilon; }
    nest::double_t get_kappa_() const { return S_.kappa; }

    // ---------------------------------------------------------------- 

    Parameters_ P_;
    State_      S_;
    Variables_  V_;
    Buffers_    B_;

    //! Mapping of recordables names to access functions
    static nest::RecordablesMap<aeif_cond_exp_multisynapse> recordablesMap_;
  };

  inline  
  nest::port aeif_cond_exp_multisynapse::check_connection(nest::Connection &c, nest::port receptor_type)
  {
    nest::SpikeEvent e;
    e.set_sender(*this);
    c.check_event(e);
    return c.get_target()->connect_sender(e, receptor_type);
  }

  inline
  nest::port aeif_cond_exp_multisynapse::connect_sender(nest::SpikeEvent &, nest::port receptor_type)
  {
	// If receptor type is less than 1 =(MIN_SPIKE_RECEPTOR) or greater or equal to 4
	// (=SUP_SPIKE_RECEPTOR) then provided receptor type is not a spike receptor.
	if ( receptor_type < MIN_SPIKE_RECEPTOR || receptor_type >= SUP_SPIKE_RECEPTOR )
		// Unknown receptor type is less than 0 or greater than 6
		// (SUP_CURR_RECEPTOR).
		if ( receptor_type < 0 || receptor_type >= SUP_CURR_RECEPTOR )
			throw nest::UnknownReceptorType(receptor_type, get_name());
	// Otherwise it is a current receptor or receptor 0 (data logging request
	// not used here and therefore incompatible.
		else
			throw nest::IncompatibleReceptorType(receptor_type, get_name(), "SpikeEvent");
	// If we arrive here the receptor type is a spike receptor and either 1, 2 or 3 e.i.
	// greater or equal to MIN_SPIKE_RECEPTOR = 1, and less than SUP_SPIKE_RECEPTOR
	// = 4. Then 0, 1, or 2 is returned.
	return receptor_type - MIN_SPIKE_RECEPTOR;
  }

  inline
  nest::port aeif_cond_exp_multisynapse::connect_sender(nest::DataLoggingRequest& dlr, 
				     nest::port receptor_type)
  {
	// If receptor type does not equal 0 then it is not a data logging request
	// receptor.
	if ( receptor_type != 0 )
		// If not a spike or current receptor that is less than 0 or greater or
		//  equal to 4 (SUP_CURR_RECEPTOR).
		if ( receptor_type < 0 || receptor_type >= SUP_CURR_RECEPTOR )
			throw nest::UnknownReceptorType(receptor_type, get_name());
	// Otherwise it is a spike or current receptor type.
		else
			throw nest::IncompatibleReceptorType(receptor_type, get_name(), "DataLoggingRequest");
	// CHANGED
	//B_.logger_.connect_logging_device(dlr, recordablesMap_);
	//return 0;

	// TO
	return B_.logger_.connect_logging_device(dlr, recordablesMap_);
  }

  inline
  nest::port aeif_cond_exp_multisynapse::connect_sender(nest::CurrentEvent &, nest::port receptor_type)
  {
	// If receptor type is less than 4 (MIN_CURR_RECEPTOR) or greater or equal
	// to 5 (SUP_CURR_RECEPTOR) the provided receptor type is not current
	// receptor.
	if ( receptor_type < MIN_CURR_RECEPTOR || receptor_type >= SUP_CURR_RECEPTOR )
		// If receptor is not a current receptor but still a receptor type that is
		// the receptor type is greater or equal to 0 or less than 3
		// (MIN_CURR_RECEPTOR).
		if ( receptor_type >= 0 && receptor_type < MIN_CURR_RECEPTOR )
			throw nest::IncompatibleReceptorType(receptor_type, get_name(), "CurrentEvent");
	// Otherwise unknown receptor type.
		else
			throw nest::UnknownReceptorType(receptor_type, get_name());
	//MIN_CURR_RECEPTOR =4, If here receptor type equals 4  and 0 is returned.
	return receptor_type - MIN_CURR_RECEPTOR;
  }
 
  inline
  void aeif_cond_exp_multisynapse::get_status(DictionaryDatum &d) const
  {
    P_.get(d);
    S_.get(d);
    nest::Archiving_Node::get_status(d);

    (*d)[nest::names::recordables] = recordablesMap_.get_list();
  }

  inline
  void aeif_cond_exp_multisynapse::set_status(const DictionaryDatum &d)
  {
    Parameters_ ptmp = P_;  // temporary copy in case of errors
    ptmp.set(d);            // throws if BadProperty
    State_      stmp = S_;  // temporary copy in case of errors
    stmp.set(d, ptmp);      // throws if BadProperty

    // We now know that (ptmp, stmp) are consistent. We do not 
    // write them back to (P_, S_) before we are also sure that 
    // the properties to be set in the parent class are internally 
    // consistent.
    nest::Archiving_Node::set_status(d);

    // if we get here, temporaries contain consistent set of properties
    P_ = ptmp;
    S_ = stmp;
  }
  
} // namespace

#endif // HAVE_GSL_1_11
#endif // AEIF_COND_EXP_MULTISYNAPSE_H
