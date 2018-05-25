/*
 *  iaf_cond_exp_bias.h
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

#ifndef IAF_COND_EXP_BIAS_H
#define IAF_COND_EXP_BIAS_H

#include "config.h"

#ifdef HAVE_GSL

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
Name: iaf_cond_exp_bias - Simple conductance based leaky integrate-and-fire neuron model.

Description:
iaf_cond_exp_bias is an implementation of a spiking neuron using IAF dynamics with
conductance-based synapses. Incoming spike events induce a post-synaptic change 
of conductance modelled by an exponential function. The exponential function 
is normalised such that an event of weight 1.0 results in a peak conductance of 1 nS. 

Parameters: 
The following parameters can be set in the status dictionary.

V_m        double - Membrane potential in mV 
E_L        double - Leak reversal potential in mV.
C_m        double - Capacity of the membrane in pF
t_ref      double - Duration of refractory period in ms. 
V_th       double - Spike threshold in mV.
V_reset    double - Reset potential of the membrane in mV.
E_ex       double - Excitatory reversal potential in mV.
E_in       double - Inhibitory reversal potential in mV.
g_L        double - Leak conductance in nS;
tau_syn_ex double - Time constant of the excitatory synaptic exponential function in ms.
tau_syn_in double - Time constant of the inhibitory synaptic exponential function in ms.
I_e        double - Constant external input current in pA.

Sends: SpikeEvent

Receives: SpikeEvent, CurrentEvent, DataLoggingRequest

References: 

Meffin, H., Burkitt, A. N., & Grayden, D. B. (2004). An analytical
model for the large, fluctuating synaptic conductance state typical of
neocortical neurons in vivo. J.  Comput. Neurosci., 16, 159–175.

Author: Sven Schrader

SeeAlso: iaf_psc_delta, iaf_psc_exp, iaf_cond_exp_bias
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
  int iaf_cond_exp_bias_dynamics (double, const double*, double*, void*);
  
  class iaf_cond_exp_bias : public nest::Archiving_Node
  {
    
  public:        
    
    typedef Node base;
    
    iaf_cond_exp_bias();
    iaf_cond_exp_bias(const iaf_cond_exp_bias&);
    ~iaf_cond_exp_bias();

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

    nest::port check_connection(nest::Connection&, nest::port);
    
    void handle(nest::SpikeEvent &);
    void handle(nest::CurrentEvent &);
    void handle(nest::DataLoggingRequest &); 
    
    nest::port connect_sender(nest::SpikeEvent &, nest::port);
    nest::port connect_sender(nest::CurrentEvent &, nest::port);
    nest::port connect_sender(nest::DataLoggingRequest &, nest::port);
    
    void get_status(DictionaryDatum &) const;
    void set_status(const DictionaryDatum &);
    
  private:
    void init_state_(const Node& proto);
    void init_buffers_();
    void calibrate();
    void update(nest::Time const &, const nest::long_t, const nest::long_t);

    // END Boilerplate function declarations ----------------------------

    // Friends --------------------------------------------------------

    // make dynamics function quasi-member
    friend int iaf_cond_exp_bias_dynamics(double, const double*, double*, void*);

    // The next two classes need to be friends to access the State_ class/member
    friend class nest::RecordablesMap<iaf_cond_exp_bias>;
    friend class nest::UniversalDataLogger<iaf_cond_exp_bias>;

  private:

    // ---------------------------------------------------------------- 

    //! Model parameters
    struct Parameters_ {
      double_t V_th_;       //!< Threshold Potential in mV
      double_t V_reset_;    //!< Reset Potential in mV
      double_t t_ref_;      //!< Refractory period in ms
      double_t g_L;         //!< Leak Conductance in nS
      double_t C_m;         //!< Membrane Capacitance in pF
      double_t E_ex;        //!< Excitatory reversal Potential in mV
      double_t E_in;        //!< Inhibitory reversal Potential in mV
      double_t E_L;         //!< Leak reversal Potential (aka resting potential) in mV
      double_t tau_synE;    //!< Synaptic Time Constant Excitatory Synapse in ms
      double_t tau_synI;    //!< Synaptic Time Constant for Inhibitory Synapse in ms
      double_t I_e;         //!< Constant Current in pA
      nest::double_t tau_j;
      nest::double_t tau_e;
      nest::double_t tau_p;
      nest::double_t fmax;
      nest::double_t gain;
      nest::double_t epsilon;
      nest::double_t K;
    
      Parameters_();  //!< Sets default parameter values

      void get(DictionaryDatum&) const;  //!< Store current values in dictionary
      void set(const DictionaryDatum&);  //!< Set values from dicitonary
    };

  public:
    // ---------------------------------------------------------------- 

    /**
     * State variables of the model.
     * @note Copy constructor and assignment operator required because
     *       of C-style array.
     */
    struct State_ {

      //! Symbolic indices to the elements of the state vector y
      enum StateVecElems { V_M = 0,           
			   G_EXC,     
			   G_INH,
                           Z_J,    
                           E_J,    
                           P_J,    
                           I_BIAS,  
			   STATE_VEC_SIZE };

      nest::double_t y_[STATE_VEC_SIZE];  //!< neuron state, must be C-array for GSL solver
      nest::int_t    r_;                  //!< number of refractory steps remaining
      nest::double_t epsilon;
      nest::double_t K;

      State_(const Parameters_&);  //!< Default initialization
      State_(const State_&);
      State_& operator=(const State_&);

      void get(DictionaryDatum&) const;
      void set(const DictionaryDatum&, const Parameters_&);
    };    

    // ---------------------------------------------------------------- 

  private:
    /**
     * Buffers of the model.
     */
    struct Buffers_ {
      Buffers_(iaf_cond_exp_bias&);                   //!<Sets buffer pointers to 0
      Buffers_(const Buffers_&, iaf_cond_exp_bias&);  //!<Sets buffer pointers to 0

      //! Logger for all analog data
      nest::UniversalDataLogger<iaf_cond_exp_bias> logger_;

      /** buffers and sums up incoming spikes/currents */
      nest::RingBuffer spike_exc_;
      nest::RingBuffer spike_inh_;
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
      double_t step_;           //!< step size in ms
      double   IntegrationStep_;//!< current integration time step, updated by GSL

      /** 
       * Input current injected by CurrentEvent.
       * This variable is used to transport the current applied into the
       * _dynamics function computing the derivative of the state vector.
       * It must be a part of Buffers_, since it is initialized once before
       * the first simulation, but not modified before later Simulate calls.
       */
      double_t I_stim_;
    };

     // ---------------------------------------------------------------- 

     /**
      * Internal variables of the model.
      */
     struct Variables_ { 
      nest::int_t    RefractoryCounts_;
     };

    // Access functions for UniversalDataLogger -------------------------------
    
    //! Read out state vector elements, used by UniversalDataLogger
    template <State_::StateVecElems elem>
    nest::double_t get_y_elem_() const { return S_.y_[elem]; }
    nest::double_t get_epsilon_() const { return S_.epsilon; }
    nest::double_t get_K_() const { return S_.K; }

    // ---------------------------------------------------------------- 

    Parameters_ P_;
    State_      S_;
    Variables_  V_;
    Buffers_    B_;

    //! Mapping of recordables names to access functions
    static nest::RecordablesMap<iaf_cond_exp_bias> recordablesMap_;
  };

  
  inline
  nest::port mynest::iaf_cond_exp_bias::check_connection(nest::Connection& c, nest::port receptor_type)
  {
    nest::SpikeEvent e;
    e.set_sender(*this);
    c.check_event(e);
    return c.get_target()->connect_sender(e, receptor_type);
  }

  inline
  nest::port mynest::iaf_cond_exp_bias::connect_sender(nest::SpikeEvent&, nest::port receptor_type)
  {
    if (receptor_type != 0)
      throw nest::UnknownReceptorType(receptor_type, get_name());
    return 0;
  }
 
  inline
  nest::port mynest::iaf_cond_exp_bias::connect_sender(nest::CurrentEvent&, nest::port receptor_type)
  {
    if (receptor_type != 0)
      throw nest::UnknownReceptorType(receptor_type, get_name());
    return 0;
  }

  inline
  nest::port mynest::iaf_cond_exp_bias::connect_sender(nest::DataLoggingRequest& dlr, 
				      nest::port receptor_type)
  {
    if (receptor_type != 0)
      throw nest::UnknownReceptorType(receptor_type, get_name());
    return B_.logger_.connect_logging_device(dlr, recordablesMap_);
  }
 
  inline
  void iaf_cond_exp_bias::get_status(DictionaryDatum &d) const
  {
    P_.get(d);
    S_.get(d);
    nest::Archiving_Node::get_status(d);

    (*d)[nest::names::recordables] = recordablesMap_.get_list();
  }

  inline
  void iaf_cond_exp_bias::set_status(const DictionaryDatum &d)
  {
    Parameters_ ptmp = P_;  // temporary copy in case of errors
    ptmp.set(d);                       // throws if BadProperty
    State_      stmp = S_;  // temporary copy in case of errors
    stmp.set(d, ptmp);                 // throws if BadProperty

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

#endif //HAVE_GSL
#endif //IAF_COND_EXP_BIAS_H
