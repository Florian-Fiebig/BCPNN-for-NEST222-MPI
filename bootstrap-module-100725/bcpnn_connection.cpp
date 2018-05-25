/*
 *  bcpnn_connection.cpp
 *
 *  Written by Philip Tully, Bernhard Kaplan
 *
 */

#include "network.h"
#include "dictdatum.h"
#include "connector_model.h"
#include "common_synapse_properties.h"
#include "bcpnn_connection.h"
#include "event.h"
#include <vector>

namespace mynest
{

  BCPNNConnection::BCPNNConnection() :
    nest::ConnectionHetWD(),
    stp_flag_(0.0),
    yi_(0.0),             /* initial conditions */
    yj_(0.0),
    taui_(10.0),
    tauj_(10.0),
    taue_(100.0),
    taup_(1000.0),
    fmax_(50.0),
    K_(1.0),
    gain_(1.0),
    zi_(0.01), 		
    zj_(0.01),
    ei_(0.01),
    ej_(0.01),
    pi_(0.01),
    pj_(0.01),
    bias_(0.0),
    U_(0.25),
    u_(U_),
    x_(U_),
    tau_rec_(600.0),
    tau_fac_(0.0),
    t_k_(0.0) 
  { 
	times_k_changed.push_back(t_k_);
	K_values_.push_back(K_);
	set_initial_eps_eij_pij();
  }

  BCPNNConnection::BCPNNConnection(const BCPNNConnection &rhs) :
    nest::ConnectionHetWD(rhs)
  {
    stp_flag_ = rhs.stp_flag_;
    yi_ = rhs.yi_;
    yj_ = rhs.yj_;
    taui_ = rhs.taui_;
    tauj_ = rhs.tauj_;
    taue_ = rhs.taue_;
    taup_ = rhs.taup_;
    epsilon_ = rhs.epsilon_;
    gain_ = rhs.gain_;
    fmax_ = rhs.fmax_;
    K_ = rhs.K_;
    zi_ = rhs.zi_;
    zj_ = rhs.zj_;
    ei_ = rhs.ei_;
    ej_ = rhs.ej_;
    eij_ = rhs.eij_;
    pi_ = rhs.pi_;
    pj_ = rhs.pj_;
    pij_ = rhs.pij_;
    bias_ = rhs.bias_;
    t_k_= rhs.t_k_;
    x_ = rhs.x_;
    tau_rec_ = rhs.tau_rec_;
    tau_fac_ = rhs.tau_fac_;
    times_k_changed.push_back(rhs.t_k_);
    U_ = rhs.U_;
    u_ = rhs.u_;
    K_values_.push_back(rhs.K_);
    set_initial_eps_eij_pij();
  }

  void BCPNNConnection::set_initial_eps_eij_pij(){
    epsilon_ = 1. / (fmax_ * taup_);
    eij_ = ei_ * ej_;
    pij_ = pi_ * pj_;
  }

  void BCPNNConnection::get_status(DictionaryDatum & d) const
  {
    nest::ConnectionHetWD::get_status(d);
    def<nest::double_t>(d, "stp_flag", stp_flag_);
    def<nest::double_t>(d, "tau_i", taui_);
    def<nest::double_t>(d, "tau_j", tauj_);
    def<nest::double_t>(d, "tau_e", taue_);
    def<nest::double_t>(d, "tau_p", taup_);
    def<nest::double_t>(d, "epsilon", epsilon_);
    def<nest::double_t>(d, "fmax", fmax_);
    def<nest::double_t>(d, "bias", bias_);
    def<nest::double_t>(d, "K", K_);
    def<nest::double_t>(d, "gain", gain_);
    def<nest::double_t>(d, "p_i", pi_);
    def<nest::double_t>(d, "p_j", pj_);
    def<nest::double_t>(d, "p_ij", pij_);
    def<nest::double_t>(d, "t_k_", t_k_);
    def<nest::double_t>(d, nest::names::dU, U_);
    def<nest::double_t>(d, nest::names::u, u_);
    def<nest::double_t>(d, nest::names::tau_rec, tau_rec_);
    def<nest::double_t>(d, nest::names::tau_fac, tau_fac_);
    def<nest::double_t>(d, nest::names::x, x_);
  }

  void BCPNNConnection::set_status(const DictionaryDatum & d, nest::ConnectorModel &cm)
  {
    nest::ConnectionHetWD::set_status(d, cm);
    updateValue<nest::double_t>(d, "stp_flag", stp_flag_);
    updateValue<nest::double_t>(d, "tau_i", taui_);
    updateValue<nest::double_t>(d, "tau_j", tauj_);
    updateValue<nest::double_t>(d, "tau_e", taue_);
    updateValue<nest::double_t>(d, "tau_p", taup_);
    updateValue<nest::double_t>(d, "K", K_);
    updateValue<nest::double_t>(d, "epsilon", epsilon_);
    updateValue<nest::double_t>(d, "fmax", fmax_);
    updateValue<nest::double_t>(d, "bias", bias_);
    updateValue<nest::double_t>(d, "gain", gain_);
    updateValue<nest::double_t>(d, "p_i", pi_);
    updateValue<nest::double_t>(d, "p_j", pj_);
    updateValue<nest::double_t>(d, "p_ij", pij_);
    updateValue<nest::double_t>(d, "t_k", t_k_);
    updateValue<nest::double_t>(d, nest::names::dU, U_);
    updateValue<nest::double_t>(d, nest::names::u, u_);
    updateValue<nest::double_t>(d, nest::names::tau_rec, tau_rec_);
    updateValue<nest::double_t>(d, nest::names::tau_fac, tau_fac_);
    updateValue<nest::double_t>(d, nest::names::x, x_);
	// only update K_values_ if K is not the same it has been initialized to
	if ((t_k_) == times_k_changed.back()){
		K_values_.pop_back();
		K_values_.push_back(K_);
	} else { 
		times_k_changed.push_back(t_k_);
		K_values_.push_back(K_);
	}
  }

   /**
   * Set properties of this connection from position p in the properties
   * array given in dictionary.
   */
  void BCPNNConnection::set_status(const DictionaryDatum & d, nest::index p, nest::ConnectorModel &cm)
  {
    nest::ConnectionHetWD::set_status(d, p, cm);
    nest::set_property<nest::double_t>(d, "stp_flag", p, stp_flag_);
    nest::set_property<nest::double_t>(d, "tau_i", p, taui_);
    nest::set_property<nest::double_t>(d, "tau_j", p, tauj_);
    nest::set_property<nest::double_t>(d, "tau_e", p, taue_);
    nest::set_property<nest::double_t>(d, "tau_p", p, taup_);
    nest::set_property<nest::double_t>(d, "K", p, K_);
    nest::set_property<nest::double_t>(d, "epsilon", p, epsilon_);
    nest::set_property<nest::double_t>(d, "fmax", p, fmax_);
    nest::set_property<nest::double_t>(d, "bias", p, bias_);
    nest::set_property<nest::double_t>(d, "gain", p, gain_);
    nest::set_property<nest::double_t>(d, "p_i", p, pi_);
    nest::set_property<nest::double_t>(d, "p_j", p, pj_);
    nest::set_property<nest::double_t>(d, "p_ij", p, pij_);
    nest::set_property<nest::double_t>(d, "t_k", p, t_k_);
    nest::set_property<nest::double_t>(d, nest::names::dUs, p, U_);
    nest::set_property<nest::double_t>(d, nest::names::us, p, u_);
    nest::set_property<nest::double_t>(d, nest::names::xs, p, x_);
    nest::set_property<nest::double_t>(d, nest::names::tau_recs, p, tau_rec_);
    nest::set_property<nest::double_t>(d, nest::names::tau_facs, p, tau_fac_);
	if ((t_k_) == times_k_changed.back()){
		K_values_.pop_back();
		K_values_.push_back(K_);
	} else { 
		times_k_changed.push_back(t_k_);
		K_values_.push_back(K_);
	}
  }

  void BCPNNConnection::initialize_property_arrays(DictionaryDatum & d) const
  {
    nest::ConnectionHetWD::initialize_property_arrays(d);
    initialize_property_array(d, "stp_flag");
    initialize_property_array(d, "tau_i");
    initialize_property_array(d, "tau_j");
    initialize_property_array(d, "tau_e");
    initialize_property_array(d, "tau_p");
    initialize_property_array(d, "K");
    initialize_property_array(d, "epsilon");
    initialize_property_array(d, "fmax");
    initialize_property_array(d, "bias");
    initialize_property_array(d, "gain");
    initialize_property_array(d, "p_i");
    initialize_property_array(d, "p_j");
    initialize_property_array(d, "p_ij");
    initialize_property_array(d, "t_k");
    initialize_property_array(d, nest::names::dUs);    
    initialize_property_array(d, nest::names::us); 
    initialize_property_array(d, nest::names::tau_recs);  
    initialize_property_array(d, nest::names::tau_facs);  
    initialize_property_array(d, nest::names::xs);
  }

  /**
   * Append properties of this connection to the given dictionary. If the
   * dictionary is empty, new arrays are created first.
   */
  void BCPNNConnection::append_properties(DictionaryDatum & d) const
  {
    nest::ConnectionHetWD::append_properties(d);
    append_property<nest::double_t>(d, "stp_flag", stp_flag_);
    append_property<nest::double_t>(d, "tau_i", taui_);
    append_property<nest::double_t>(d, "tau_j", tauj_);
    append_property<nest::double_t>(d, "tau_e", taue_);
    append_property<nest::double_t>(d, "tau_p", taup_);
    append_property<nest::double_t>(d, "K", K_);
    append_property<nest::double_t>(d, "epsilon", epsilon_);
    append_property<nest::double_t>(d, "fmax", fmax_);
    append_property<nest::double_t>(d, "bias", bias_);
    append_property<nest::double_t>(d, "gain", gain_);
    append_property<nest::double_t>(d, "p_i", pi_);
    append_property<nest::double_t>(d, "p_j", pj_);
    append_property<nest::double_t>(d, "p_ij", pij_);
    append_property<nest::double_t>(d, "t_k", t_k_);
    append_property<nest::double_t>(d, nest::names::dUs, U_); 
    append_property<nest::double_t>(d, nest::names::us, u_); 
    append_property<nest::double_t>(d, nest::names::tau_recs, tau_rec_);  
    append_property<nest::double_t>(d, nest::names::tau_facs, tau_fac_);  
    append_property<nest::double_t>(d, nest::names::xs, x_);
  }
} // of namespace nest
