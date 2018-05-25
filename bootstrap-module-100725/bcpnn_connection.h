/*
 *  bcpnn_connection.h
 *
 *  Written by Philip Tully
 *
 */

#ifndef BCPNN_CONNECTION_H
#define BCPNN_CONNECTION_H

/* BeginDocumentation
  Name: bcpnn_synapse - Synapse type for incremental, Bayesian spike-timing 
   dependent plasticity.

  Description:
   bcpnn_synapse is a connector to create synapses with incremental, Bayesian 
   spike timing dependent plasticity.

   tau_i	double - Primary trace presynaptic time constant
   tau_j	double - Primary trace postsynaptic time constant
   tau_e	double - Secondary trace time constant
   tau_p	double - Tertiarty trace time constant
   p_i		double - \
   p_j		double -  >- these 3 initial conditions determine weight, i.e. log(p_ij/(p_i * p_j)).
   p_ij		double - /
   K_		double - Print-now signal // Neuromodulation. Turn off learning, K = 0.
   fmax_        double - Frequency assumed as maximum firing, for match with abstract rule
   epsilon_     double - lowest possible probability of spiking, e.g. lowest assumed firing rate
   bias_        double - ANN interpretation. Only calculated here to demonstrate match to rule. 
                         Will be eliminated in future versions, where bias will be calculated postsynaptically
   gain_    double - Coefficient to scale weight as conductance, can be zero-ed out
	K_values_ vector of doubles storing the recent changes of K_ since the last pre spike occured, cleared after each send function

  Transmits: SpikeEvent
   
  References:
   [1] Wahlgren and Lansner (2001) Biological Evaluation of a Hebbian-Bayesian
       learning rule. Neurocomputing, 38-40, 433-438

   [2] Bergel, Transforming the BCPNN Learning Rule for Spiking Units to a
       Learning Rule for Non-Spiking Units (2010). KTH Masters Thesis.

  FirstVersion: November 2011
  CurrentVersion: March 2012
  Authors: Philip Tully, Bernhard Kaplan
          tully@csc.kth.se, bkaplan@kth.se
  SeeAlso: synapsedict, stdp_synapse, tsodyks_synapse, static_synapse
*/

/* for Debugging */
#include <iostream>
using namespace std;

#include "connection_het_wd.h"
#include "archiving_node.h"
#include "generic_connector.h"
#include <cmath>

namespace mynest
{
  class BCPNNConnection : public nest::ConnectionHetWD
  {
    public:
      /* Default Constructor. Sets default values for all parameters. Needed by GenericConnectorModel. */
      BCPNNConnection();

      /* Copy constructor. Needs to be defined properly in order for GenericConnector to work. */
      BCPNNConnection(const BCPNNConnection &);

      /* Default Destructor. */
      ~BCPNNConnection() {}

      void check_connection(nest::Node & s, nest::Node & r, nest::port receptor_type, nest::double_t t_lastspike);

      /* Get all properties of this connection and put them into a dictionary. */
      void get_status(DictionaryDatum & d) const;

      /* Set properties of this connection from the values given in dictionary. */
      void set_status(const DictionaryDatum & d, nest::ConnectorModel &cm);

      /* Set properties of this connection from position p in the properties array given in dictionary. */
      void set_status(const DictionaryDatum & d, nest::index p, nest::ConnectorModel &cm);

      /* Create new empty arrays for the properties of this connection in the given dictionary. It is assumed 
         that they do not exist before. */
      void initialize_property_arrays(DictionaryDatum & d) const;

      /* Append properties of this connection to the given dictionary. If the dictionary is empty, new arrays 
         are created first. */
      void append_properties(DictionaryDatum & d) const;

      /* Send an event to the receiver of this connection.  */
      void send(nest::Event& e, nest::double_t t_lastspike, const nest::CommonSynapseProperties &cp);

      /* Overloaded for all supported event types. */
      using nest::Connection::check_event;
      void check_event(nest::SpikeEvent&) {}

	  // setting the correct values for epsilon, eij, pij
      void set_initial_eps_eij_pij();

    private:
      /* data members of each connection */
      nest::double_t stp_flag_;
      nest::double_t yi_;
      nest::double_t yj_;
      nest::double_t taui_;
      nest::double_t tauj_;
      nest::double_t taue_;
      nest::double_t taup_;
      nest::double_t epsilon_;
      nest::double_t K_;
      nest::double_t bias_;
      nest::double_t fmax_;
      nest::double_t gain_;
      nest::double_t zi_;
      nest::double_t zj_;
      nest::double_t ei_;
      nest::double_t ej_;
      nest::double_t eij_;
      nest::double_t pi_;
      nest::double_t pj_;
      nest::double_t pij_;
      nest::double_t t_k_;
      std::vector<nest::double_t> times_k_changed;
      std::vector<nest::double_t> post_spiketimes;
      std::vector<nest::double_t> K_values_;
      nest::double_t U_;
      nest::double_t u_;
      nest::double_t x_;
      nest::double_t tau_rec_;
      nest::double_t tau_fac_;
  }; /* of class BCPNNConnection */

  inline 
  void BCPNNConnection::check_connection(nest::Node & s, nest::Node & r, nest::port receptor_type, nest::double_t t_lastspike)
  {
    nest::ConnectionHetWD::check_connection(s, r, receptor_type, t_lastspike);

    // For a new synapse, t_lastspike contains the point in time of the last spike.
    // So we initially read the history(t_last_spike - dendritic_delay, ...,  T_spike-dendritic_delay]
    // which increases the access counter for these entries.
    // At registration, all entries' access counters of history[0, ..., t_last_spike - dendritic_delay] will be 
    // incremented by the following call to Archiving_Node::register_stdp_connection().
    // See bug #218 for details.
    r.register_stdp_connection(t_lastspike - nest::Time(nest::Time::step(delay_)).get_ms());
  }

  /* Send an event to the receiver of this connection.
   * \param e The event to send
   * \param p The port under which this connection is stored in the Connector.
   * \param t_lastspike Time point of last spike emitted 
  
   note: every time this method is called by an outside function, a presynaptic
       event has occured and is being transmitted to the postsynaptic side. */

  inline
  void BCPNNConnection::send(nest::Event& e, nest::double_t t_lastspike, const nest::CommonSynapseProperties &)
  {
    nest::double_t t_spike = e.get_stamp().get_ms();  /* time stamp of current spike event */
    nest::double_t dendritic_delay = nest::Time(nest::Time::step(delay_)).get_ms();    /* delay from dendrite -> soma */
    nest::double_t resolution = nest::Time::get_resolution().get_ms();
    nest::int_t spike_width = nest::int_t (1. / resolution); 
    nest::double_t spike_height = 1000.0 / fmax_;     /* normalizing to match this spiking rule to abstract = 1000/FMAX (Hz)*/
    nest::int_t counter = 0;                          /* ensuring traces reverberate for duration of the spike width */

    /*nest::double_t h = e.get_stamp().get_ms() - t_lastspike;  
    nest::double_t f = std::exp(-h/tau_rec_);
    nest::double_t u_decay = (tau_fac_ < 1.0e-10) ? 0.0 : std::exp(-h/tau_fac_);*/

    /* get spike history in relevant range (t1, t2] from post-synaptic neuron */
    std::deque<nest::histentry>::iterator start;
    std::deque<nest::histentry>::iterator finish;
    target_->get_history(t_lastspike - dendritic_delay, t_spike - dendritic_delay, &start, &finish);

    while (start != finish)  {/* loop until you get to last post spike */
        post_spiketimes.push_back(start->t_);
        start++;
    }    
 
    counter = 0;
    nest::int_t number_iterations = (nest::int_t)((t_spike - t_lastspike)/resolution);
    nest::double_t K_vec_init = K_;
    if (K_values_.size() > 1) {
        K_vec_init = K_values_.front();
    }
    std::vector<nest::double_t> K_vec (number_iterations, K_vec_init);

    if (K_values_.size() > 1) {
        std::vector<nest::double_t>::iterator K_it = K_values_.end();
        std::vector<nest::double_t>::iterator time_it = times_k_changed.end();
        if (times_k_changed.back() >= t_lastspike){ 
            K_it--; 
            time_it--;
            nest::int_t idx_first = (nest::int_t) ((t_spike - t_lastspike) / resolution);
            nest::int_t idx_second;
            while (*time_it > t_lastspike){
                idx_second = (nest::int_t) ((*time_it - t_lastspike)/ resolution);
                for (nest::int_t i_k=idx_first-1; i_k >= idx_second; --i_k) {					
                    K_vec.at(i_k) = *K_it;
                } // for
                idx_first = idx_second;
                time_it--;
                K_it--;
            } // end of while
        }
        K_values_.clear();
        K_values_.push_back(K_);
        times_k_changed.clear();
        times_k_changed.push_back(*time_it);
    }
    
    /* Create a vector to represent the post spikes as a trace */
    std::vector<nest::double_t> post_active (number_iterations, 0.);
    std::vector<nest::double_t>::iterator post_it = post_spiketimes.begin(); 

    for (nest::int_t timestep = 0; timestep < number_iterations; timestep++){
        /* CASE: Default. Neither Pre nor Post spike. */
        yi_ = 0.0; 
        yj_ = 0.0;

        /* CASE: Pre without (*OR WITH post) spike - synchronous events handled automatically. */
        if(timestep == 0 && t_lastspike != 0.) {
            yi_ = spike_height * spike_width;
        }

        // if you have any post spike at all
        if (post_spiketimes.size() > 0) { 
            if (post_it != post_spiketimes.end()) { 
                if (timestep == (nest::int_t)((*post_it) - t_lastspike) / resolution){
                    yj_ = spike_height * spike_width;
                    post_it++;
                }
            }
        }

        /* Primary synaptic traces */
        zi_ += (yi_ - zi_ + epsilon_ ) * resolution / taui_;
        zj_ += (yj_ - zj_ + epsilon_ ) * resolution / tauj_;

        /* Secondary synaptic traces */
        ei_  += (zi_ - ei_) * resolution / taue_;
        ej_  += (zj_ - ej_) * resolution / taue_;
        eij_ += (zi_ * zj_ - eij_) * resolution / taue_;

        /* Tertiary synaptic traces */
        pi_  += K_vec.at(timestep) * (ei_ - pi_) * resolution / taup_;
        pj_  += K_vec.at(timestep) * (ej_ - pj_) * resolution / taup_;
        pij_ += K_vec.at(timestep) * (eij_ - pij_) * resolution / taup_;
    } /* of for */

    bias_ = std::log(pj_);
    
    if (stp_flag_ > 0.5){
        double_t h = e.get_stamp().get_ms() - t_lastspike;  
        double_t x_decay = std::exp(-h/tau_rec_);
        double_t u_decay = (tau_fac_ < 1.0e-10) ? 0.0 : std::exp(-h/tau_fac_);
        x_= 1. + (x_ -x_*u_ -1.)*x_decay; // Eq. 5 from reference [3]
        u_= U_+u_*(1.-U_)*u_decay; 
        weight_ = x_ * u_ * gain_ * (std::log(pij_ / (pi_ * pj_)));
    } else {
        weight_ = gain_ * (std::log(pij_ / (pi_ * pj_)));
    }

    /* Send the spike to the target */
    e.set_receiver(*target_);
    e.set_weight(weight_);
    e.set_delay(delay_);
    e.set_rport(rport_);
    e();
    post_spiketimes.clear();
    } //of BCPNNConnection::send
} //of namespace mynest
#endif // of #ifndef BCPNN_CONNECTION_H

