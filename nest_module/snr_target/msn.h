/**
 *  msn.h
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
 *  Generated from NESTML at time: 2022-09-07 11:08:01.546974
**/
#ifndef MSN
#define MSN

#include "config.h"

#ifndef HAVE_GSL
#error "The GSL library is required for neurons that require a numerical solver."
#endif

// External includes:
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

// Includes from nestkernel:
#include "archiving_node.h"
#include "connection.h"
#include "event.h"
#include "nest_types.h"
#include "ring_buffer.h"
#include "universal_data_logger.h"

// Includes from sli:
#include "dictdatum.h"

namespace nest
{
namespace msn_names
{
    const Name _V_m( "V_m" );
    const Name _U_m( "U_m" );
    const Name _I( "I" );
    const Name _g_snr__X__spikes_snr( "g_snr__X__spikes_snr" );
    const Name _g_msn__X__spikes_msn( "g_msn__X__spikes_msn" );
    const Name _a( "a" );
    const Name _b( "b" );
    const Name _c( "c" );
    const Name _d( "d" );
    const Name _I_e( "I_e" );
    const Name _V_th( "V_th" );
    const Name _V_min( "V_min" );
    const Name _Kplus( "Kplus" );
    const Name _V_r( "V_r" );
    const Name _V_T( "V_T" );
    const Name _C_m( "C_m" );
    const Name _tau_decay_msn( "tau_decay_msn" );
    const Name _tau_decay_snr( "tau_decay_snr" );
    const Name _E_rev_msn( "E_rev_msn" );
    const Name _E_rev_snr( "E_rev_snr" );
}
}



/**
 * Function computing right-hand side of ODE for GSL solver.
 * @note Must be declared here so we can befriend it in class.
 * @note Must have C-linkage for passing to GSL. Internally, it is
 *       a first-class C++ function, but cannot be a member function
 *       because of the C-linkage.
 * @note No point in declaring it inline, since it is called
 *       through a function pointer.
 * @param void* Pointer to model neuron instance.
**/
extern "C" inline int msn_dynamics( double, const double y[], double f[], void* pnode );


#include "nest_time.h"




/* BeginDocumentation
  Name: msn.

  Description:

    

  Parameters:
  The following parameters can be set in the status dictionary.
I_e [pA]  constant external input current
V_th [mV] V peak
V_r [mV] resting potential
V_T [mV] threshold potential
C_m [pF]  membrane capacitance


  Dynamic state variables:
I [pA]  input current


  Sends: nest::SpikeEvent

  Receives: Spike, Current, DataLoggingRequest
*/
class msn : public nest::ArchivingNode
{
public:
  /**
   * The constructor is only used to create the model prototype in the model manager.
  **/
  msn();

  /**
   * The copy constructor is used to create model copies and instances of the model.
   * @node The copy constructor needs to initialize the parameters and the state.
   *       Initialization of buffers and interal variables is deferred to
   *       @c init_buffers_() and @c calibrate().
  **/
  msn(const msn &);

  /**
   * Destructor.
  **/
  ~msn();

  // -------------------------------------------------------------------------
  //   Import sets of overloaded virtual functions.
  //   See: Technical Issues / Virtual Functions: Overriding, Overloading,
  //        and Hiding
  // -------------------------------------------------------------------------

  using nest::Node::handles_test_event;
  using nest::Node::handle;

  /**
   * Used to validate that we can send nest::SpikeEvent to desired target:port.
  **/
  nest::port send_test_event(nest::Node& target, nest::rport receptor_type, nest::synindex, bool);

  // -------------------------------------------------------------------------
  //   Functions handling incoming events.
  //   We tell nest that we can handle incoming events of various types by
  //   defining handle() for the given event.
  // -------------------------------------------------------------------------


  void handle(nest::SpikeEvent &);        //! accept spikes
  void handle(nest::CurrentEvent &);      //! accept input current
  void handle(nest::DataLoggingRequest &);//! allow recording with multimeter
  nest::port handles_test_event(nest::SpikeEvent&, nest::port);
  nest::port handles_test_event(nest::CurrentEvent&, nest::port);
  nest::port handles_test_event(nest::DataLoggingRequest&, nest::port);

  // -------------------------------------------------------------------------
  //   Functions for getting/setting parameters and state values.
  // -------------------------------------------------------------------------

  void get_status(DictionaryDatum &) const;
  void set_status(const DictionaryDatum &);

  // -------------------------------------------------------------------------
  //   Getters/setters for state block
  // -------------------------------------------------------------------------

  inline double get_V_m() const
  {
    return S_.ode_state[State_::V_m];
  }

  inline void set_V_m(const double __v)
  {
    S_.ode_state[State_::V_m] = __v;
  }

  inline double get_U_m() const
  {
    return S_.ode_state[State_::U_m];
  }

  inline void set_U_m(const double __v)
  {
    S_.ode_state[State_::U_m] = __v;
  }

  inline double get_I() const
  {
    return S_.ode_state[State_::I];
  }

  inline void set_I(const double __v)
  {
    S_.ode_state[State_::I] = __v;
  }

  inline double get_g_snr__X__spikes_snr() const
  {
    return S_.ode_state[State_::g_snr__X__spikes_snr];
  }

  inline void set_g_snr__X__spikes_snr(const double __v)
  {
    S_.ode_state[State_::g_snr__X__spikes_snr] = __v;
  }

  inline double get_g_msn__X__spikes_msn() const
  {
    return S_.ode_state[State_::g_msn__X__spikes_msn];
  }

  inline void set_g_msn__X__spikes_msn(const double __v)
  {
    S_.ode_state[State_::g_msn__X__spikes_msn] = __v;
  }


  // -------------------------------------------------------------------------
  //   Getters/setters for parameters
  // -------------------------------------------------------------------------

  inline double get_a() const
  {
    return P_.a;
  }

  inline void set_a(const double __v)
  {
    P_.a = __v;
  }

  inline double get_b() const
  {
    return P_.b;
  }

  inline void set_b(const double __v)
  {
    P_.b = __v;
  }

  inline double get_c() const
  {
    return P_.c;
  }

  inline void set_c(const double __v)
  {
    P_.c = __v;
  }

  inline double get_d() const
  {
    return P_.d;
  }

  inline void set_d(const double __v)
  {
    P_.d = __v;
  }

  inline double get_I_e() const
  {
    return P_.I_e;
  }

  inline void set_I_e(const double __v)
  {
    P_.I_e = __v;
  }

  inline double get_V_th() const
  {
    return P_.V_th;
  }

  inline void set_V_th(const double __v)
  {
    P_.V_th = __v;
  }

  inline double get_V_min() const
  {
    return P_.V_min;
  }

  inline void set_V_min(const double __v)
  {
    P_.V_min = __v;
  }

  inline double get_Kplus() const
  {
    return P_.Kplus;
  }

  inline void set_Kplus(const double __v)
  {
    P_.Kplus = __v;
  }

  inline double get_V_r() const
  {
    return P_.V_r;
  }

  inline void set_V_r(const double __v)
  {
    P_.V_r = __v;
  }

  inline double get_V_T() const
  {
    return P_.V_T;
  }

  inline void set_V_T(const double __v)
  {
    P_.V_T = __v;
  }

  inline double get_C_m() const
  {
    return P_.C_m;
  }

  inline void set_C_m(const double __v)
  {
    P_.C_m = __v;
  }

  inline double get_tau_decay_msn() const
  {
    return P_.tau_decay_msn;
  }

  inline void set_tau_decay_msn(const double __v)
  {
    P_.tau_decay_msn = __v;
  }

  inline double get_tau_decay_snr() const
  {
    return P_.tau_decay_snr;
  }

  inline void set_tau_decay_snr(const double __v)
  {
    P_.tau_decay_snr = __v;
  }

  inline double get_E_rev_msn() const
  {
    return P_.E_rev_msn;
  }

  inline void set_E_rev_msn(const double __v)
  {
    P_.E_rev_msn = __v;
  }

  inline double get_E_rev_snr() const
  {
    return P_.E_rev_snr;
  }

  inline void set_E_rev_snr(const double __v)
  {
    P_.E_rev_snr = __v;
  }


  // -------------------------------------------------------------------------
  //   Getters/setters for internals
  // -------------------------------------------------------------------------

  inline double get___h() const
  {
    return V_.__h;
  }

  inline void set___h(const double __v)
  {
    V_.__h = __v;
  }

  inline double get___P__g_snr__X__spikes_snr__g_snr__X__spikes_snr() const
  {
    return V_.__P__g_snr__X__spikes_snr__g_snr__X__spikes_snr;
  }

  inline void set___P__g_snr__X__spikes_snr__g_snr__X__spikes_snr(const double __v)
  {
    V_.__P__g_snr__X__spikes_snr__g_snr__X__spikes_snr = __v;
  }

  inline double get___P__g_msn__X__spikes_msn__g_msn__X__spikes_msn() const
  {
    return V_.__P__g_msn__X__spikes_msn__g_msn__X__spikes_msn;
  }

  inline void set___P__g_msn__X__spikes_msn__g_msn__X__spikes_msn(const double __v)
  {
    V_.__P__g_msn__X__spikes_msn__g_msn__X__spikes_msn = __v;
  }



protected:

private:
  void calibrate_variables(bool exclude_timestep=false);

private:
  /**
   * Synapse types to connect to
   * @note Excluded upper and lower bounds are defined as INF_, SUP_.
   *       Excluding port 0 avoids accidental connections.
  **/
  enum SynapseTypes
  {
    INF_SPIKE_RECEPTOR = 0,
      SPIKES_MSN ,
      SPIKES_SNR ,
    SUP_SPIKE_RECEPTOR
  };

#if NEST2_COMPAT
  /**
   * Reset state of neuron.
  **/
  void init_state_(const Node& proto);
#endif

  /**
   * Reset internal buffers of neuron.
  **/
  void init_buffers_();

  /**
   * Initialize auxiliary quantities, leave parameters and state untouched.
  **/
  void calibrate();

  /**
   * Take neuron through given time interval
  **/
  void update(nest::Time const &, const long, const long);

  // The next two classes need to be friends to access the State_ class/member
  friend class nest::RecordablesMap<msn>;
  friend class nest::UniversalDataLogger<msn>;

  /**
   * Free parameters of the neuron.
   *
   *
   *
   * These are the parameters that can be set by the user through @c `node.set()`.
   * They are initialized from the model prototype when the node is created.
   * Parameters do not change during calls to @c update() and are not reset by
   * @c ResetNetwork.
   *
   * @note Parameters_ need neither copy constructor nor @c operator=(), since
   *       all its members are copied properly by the default copy constructor
   *       and assignment operator. Important:
   *       - If Parameters_ contained @c Time members, you need to define the
   *         assignment operator to recalibrate all members of type @c Time . You
   *         may also want to define the assignment operator.
   *       - If Parameters_ contained members that cannot copy themselves, such
   *         as C-style arrays, you need to define the copy constructor and
   *         assignment operator to copy those members.
  **/
  struct Parameters_
  {    
    double a;
    double b;
    double c;
    double d;
    //!  constant external input current
    double I_e;
    //! V peak
    double V_th;
    double V_min;
    double Kplus;
    //! resting potential
    double V_r;
    //! threshold potential
    double V_T;
    //!  membrane capacitance
    double C_m;
    double tau_decay_msn;
    double tau_decay_snr;
    double E_rev_msn;
    double E_rev_snr;

    double __gsl_error_tol;

    /**
     * Initialize parameters to their default values.
    **/
    Parameters_();
  };

  /**
   * Dynamic state of the neuron.
   *
   *
   *
   * These are the state variables that are advanced in time by calls to
   * @c update(). In many models, some or all of them can be set by the user
   * through @c `node.set()`. The state variables are initialized from the model
   * prototype when the node is created. State variables are reset by @c ResetNetwork.
   *
   * @note State_ need neither copy constructor nor @c operator=(), since
   *       all its members are copied properly by the default copy constructor
   *       and assignment operator. Important:
   *       - If State_ contained @c Time members, you need to define the
   *         assignment operator to recalibrate all members of type @c Time . You
   *         may also want to define the assignment operator.
   *       - If State_ contained members that cannot copy themselves, such
   *         as C-style arrays, you need to define the copy constructor and
   *         assignment operator to copy those members.
  **/
  struct State_
  {
    //! Symbolic indices to the elements of the state vector y
    enum StateVecElems
    {
      V_m,
      U_m,
      g_snr__X__spikes_snr,
      g_msn__X__spikes_msn,
      // moved state variables from synapse
      I,
      STATE_VEC_SIZE
    };

    //! state vector, must be C-array for GSL solver
    double ode_state[STATE_VEC_SIZE];

    State_();
  };

  /**
   * Internal variables of the neuron.
   *
   *
   *
   * These variables must be initialized by @c calibrate, which is called before
   * the first call to @c update() upon each call to @c Simulate.
   * @node Variables_ needs neither constructor, copy constructor or assignment operator,
   *       since it is initialized by @c calibrate(). If Variables_ has members that
   *       cannot destroy themselves, Variables_ will need a destructor.
  **/
  struct Variables_
  {
    double __h;
    double __P__g_snr__X__spikes_snr__g_snr__X__spikes_snr;
    double __P__g_msn__X__spikes_msn__g_msn__X__spikes_msn;
  };

  /**
   * Buffers of the neuron.
   * Usually buffers for incoming spikes and data logged for analog recorders.
   * Buffers must be initialized by @c init_buffers_(), which is called before
   * @c calibrate() on the first call to @c Simulate after the start of NEST,
   * ResetKernel or ResetNetwork.
   * @node Buffers_ needs neither constructor, copy constructor or assignment operator,
   *       since it is initialized by @c init_nodes_(). If Buffers_ has members that
   *       cannot destroy themselves, Buffers_ will need a destructor.
  **/
  struct Buffers_
  {
    Buffers_(msn &);
    Buffers_(const Buffers_ &, msn &);

    /**
     * Logger for all analog data
    **/
    nest::UniversalDataLogger<msn> logger_;
    std::vector<long> receptor_types_;
    // -----------------------------------------------------------------------
    //   Buffers and sums of incoming spikes/currents per timestep
    // -----------------------------------------------------------------------
    std::vector< nest::RingBuffer > spike_inputs_;
    inline nest::RingBuffer& get_spikes_msn() {  return spike_inputs_[SPIKES_MSN - 1]; }
    double spikes_msn_grid_sum_;
    inline nest::RingBuffer& get_spikes_snr() {  return spike_inputs_[SPIKES_SNR - 1]; }
    double spikes_snr_grid_sum_;
    //!< Buffer for input (type: pA)
    nest::RingBuffer I_stim;
    inline nest::RingBuffer& get_I_stim() {return I_stim;}
    double I_stim_grid_sum_;

    // -----------------------------------------------------------------------
    //   GSL ODE solver data structures
    // -----------------------------------------------------------------------

    gsl_odeiv_step* __s;    //!< stepping function
    gsl_odeiv_control* __c; //!< adaptive stepsize control function
    gsl_odeiv_evolve* __e;  //!< evolution function
    gsl_odeiv_system __sys; //!< struct describing system

    // __integration_step should be reset with the neuron on ResetNetwork,
    // but remain unchanged during calibration. Since it is initialized with
    // step_, and the resolution cannot change after nodes have been created,
    // it is safe to place both here.
    double __step;             //!< step size in ms
    double __integration_step; //!< current integration time step, updated by GSL
  };

  // -------------------------------------------------------------------------
  //   Getters/setters for inline expressions
  // -------------------------------------------------------------------------
  
  inline double get_I_syn_msn() const
  {
    return get_g_msn__X__spikes_msn();
  }

  inline double get_I_syn_snr() const
  {
    return get_g_snr__X__spikes_snr();
  }

  inline double get_I_msn() const
  {
    return (get_g_msn__X__spikes_msn()) * (P_.E_rev_msn - get_V_m());
  }

  inline double get_I_snr() const
  {
    return (get_g_snr__X__spikes_snr()) * (P_.E_rev_snr - get_V_m());
  }

  inline double get_I_syn() const
  {
    return ((get_g_msn__X__spikes_msn()) * (P_.E_rev_msn - get_V_m())) + ((get_g_snr__X__spikes_snr()) * (P_.E_rev_snr - get_V_m()));
  }


  // -------------------------------------------------------------------------
  //   Getters/setters for input buffers
  // -------------------------------------------------------------------------
  
  inline nest::RingBuffer& get_spikes_msn() {return B_.get_spikes_msn();};
  inline nest::RingBuffer& get_spikes_snr() {return B_.get_spikes_snr();};
  inline nest::RingBuffer& get_I_stim() {return B_.get_I_stim();};

  // -------------------------------------------------------------------------
  //   Member variables of neuron model.
  //   Each model neuron should have precisely the following four data members,
  //   which are one instance each of the parameters, state, buffers and variables
  //   structures. Experience indicates that the state and variables member should
  //   be next to each other to achieve good efficiency (caching).
  //   Note: Devices require one additional data member, an instance of the
  //   ``Device`` child class they belong to.
  // -------------------------------------------------------------------------


  Parameters_ P_;  //!< Free parameters.
  State_      S_;  //!< Dynamic state.
  Variables_  V_;  //!< Internal Variables
  Buffers_    B_;  //!< Buffers.

  //! Mapping of recordables names to access functions
  static nest::RecordablesMap<msn> recordablesMap_;
  friend int msn_dynamics( double, const double y[], double f[], void* pnode );







}; /* neuron msn */

inline nest::port msn::send_test_event(nest::Node& target, nest::rport receptor_type, nest::synindex, bool)
{
  // You should usually not change the code in this function.
  // It confirms that the target of connection @c c accepts @c nest::SpikeEvent on
  // the given @c receptor_type.
  nest::SpikeEvent e;
  e.set_sender(*this);
  return target.handles_test_event(e, receptor_type);
}

inline nest::port msn::handles_test_event(nest::SpikeEvent&, nest::port receptor_type)
{
    assert( B_.spike_inputs_.size() == 2 );

    if ( !( INF_SPIKE_RECEPTOR < receptor_type && receptor_type < SUP_SPIKE_RECEPTOR ) )
    {
      throw nest::UnknownReceptorType( receptor_type, get_name() );
      return 0;
    }
    else
    {
      return receptor_type - 1;
    }
}

inline nest::port msn::handles_test_event(nest::CurrentEvent&, nest::port receptor_type)
{
  // You should usually not change the code in this function.
  // It confirms to the connection management system that we are able
  // to handle @c CurrentEvent on port 0. You need to extend the function
  // if you want to differentiate between input ports.
  if (receptor_type != 0)
  {
    throw nest::UnknownReceptorType(receptor_type, get_name());
  }
  return 0;
}

inline nest::port msn::handles_test_event(nest::DataLoggingRequest& dlr, nest::port receptor_type)
{
  // You should usually not change the code in this function.
  // It confirms to the connection management system that we are able
  // to handle @c DataLoggingRequest on port 0.
  // The function also tells the built-in UniversalDataLogger that this node
  // is recorded from and that it thus needs to collect data during simulation.
  if (receptor_type != 0)
  {
    throw nest::UnknownReceptorType(receptor_type, get_name());
  }

  return B_.logger_.connect_logging_device(dlr, recordablesMap_);
}

inline void msn::get_status(DictionaryDatum &__d) const
{
  // parameters
  def<double>(__d, nest::msn_names::_a, get_a());
  def<double>(__d, nest::msn_names::_b, get_b());
  def<double>(__d, nest::msn_names::_c, get_c());
  def<double>(__d, nest::msn_names::_d, get_d());
  def<double>(__d, nest::msn_names::_I_e, get_I_e());
  def<double>(__d, nest::msn_names::_V_th, get_V_th());
  def<double>(__d, nest::msn_names::_V_min, get_V_min());
  def<double>(__d, nest::msn_names::_Kplus, get_Kplus());
  def<double>(__d, nest::msn_names::_V_r, get_V_r());
  def<double>(__d, nest::msn_names::_V_T, get_V_T());
  def<double>(__d, nest::msn_names::_C_m, get_C_m());
  def<double>(__d, nest::msn_names::_tau_decay_msn, get_tau_decay_msn());
  def<double>(__d, nest::msn_names::_tau_decay_snr, get_tau_decay_snr());
  def<double>(__d, nest::msn_names::_E_rev_msn, get_E_rev_msn());
  def<double>(__d, nest::msn_names::_E_rev_snr, get_E_rev_snr());

  // initial values for state variables in ODE or kernel
  def<double>(__d, nest::msn_names::_V_m, get_V_m());
  def<double>(__d, nest::msn_names::_U_m, get_U_m());
  def<double>(__d, nest::msn_names::_I, get_I());
  def<double>(__d, nest::msn_names::_g_snr__X__spikes_snr, get_g_snr__X__spikes_snr());
  def<double>(__d, nest::msn_names::_g_msn__X__spikes_msn, get_g_msn__X__spikes_msn());

  ArchivingNode::get_status( __d );
  DictionaryDatum __receptor_type = new Dictionary();
  ( *__receptor_type )[ "SPIKES_MSN" ] = SPIKES_MSN;
  ( *__receptor_type )[ "SPIKES_SNR" ] = SPIKES_SNR;
  ( *__d )[ "receptor_types" ] = __receptor_type;

  (*__d)[nest::names::recordables] = recordablesMap_.get_list();
  def< double >(__d, nest::names::gsl_error_tol, P_.__gsl_error_tol);
  if ( P_.__gsl_error_tol <= 0. ){
    throw nest::BadProperty( "The gsl_error_tol must be strictly positive." );
  }
}

inline void msn::set_status(const DictionaryDatum &__d)
{
  // parameters
  double tmp_a = get_a();
  updateValue<double>(__d, nest::msn_names::_a, tmp_a);
  double tmp_b = get_b();
  updateValue<double>(__d, nest::msn_names::_b, tmp_b);
  double tmp_c = get_c();
  updateValue<double>(__d, nest::msn_names::_c, tmp_c);
  double tmp_d = get_d();
  updateValue<double>(__d, nest::msn_names::_d, tmp_d);
  double tmp_I_e = get_I_e();
  updateValue<double>(__d, nest::msn_names::_I_e, tmp_I_e);
  double tmp_V_th = get_V_th();
  updateValue<double>(__d, nest::msn_names::_V_th, tmp_V_th);
  double tmp_V_min = get_V_min();
  updateValue<double>(__d, nest::msn_names::_V_min, tmp_V_min);
  double tmp_Kplus = get_Kplus();
  updateValue<double>(__d, nest::msn_names::_Kplus, tmp_Kplus);
  double tmp_V_r = get_V_r();
  updateValue<double>(__d, nest::msn_names::_V_r, tmp_V_r);
  double tmp_V_T = get_V_T();
  updateValue<double>(__d, nest::msn_names::_V_T, tmp_V_T);
  double tmp_C_m = get_C_m();
  updateValue<double>(__d, nest::msn_names::_C_m, tmp_C_m);
  double tmp_tau_decay_msn = get_tau_decay_msn();
  updateValue<double>(__d, nest::msn_names::_tau_decay_msn, tmp_tau_decay_msn);
  double tmp_tau_decay_snr = get_tau_decay_snr();
  updateValue<double>(__d, nest::msn_names::_tau_decay_snr, tmp_tau_decay_snr);
  double tmp_E_rev_msn = get_E_rev_msn();
  updateValue<double>(__d, nest::msn_names::_E_rev_msn, tmp_E_rev_msn);
  double tmp_E_rev_snr = get_E_rev_snr();
  updateValue<double>(__d, nest::msn_names::_E_rev_snr, tmp_E_rev_snr);

  // initial values for state variables in ODE or kernel
  double tmp_V_m = get_V_m();
  updateValue<double>(__d, nest::msn_names::_V_m, tmp_V_m);
  double tmp_U_m = get_U_m();
  updateValue<double>(__d, nest::msn_names::_U_m, tmp_U_m);
  double tmp_I = get_I();
  updateValue<double>(__d, nest::msn_names::_I, tmp_I);
  double tmp_g_snr__X__spikes_snr = get_g_snr__X__spikes_snr();
  updateValue<double>(__d, nest::msn_names::_g_snr__X__spikes_snr, tmp_g_snr__X__spikes_snr);
  double tmp_g_msn__X__spikes_msn = get_g_msn__X__spikes_msn();
  updateValue<double>(__d, nest::msn_names::_g_msn__X__spikes_msn, tmp_g_msn__X__spikes_msn);

  // We now know that (ptmp, stmp) are consistent. We do not
  // write them back to (P_, S_) before we are also sure that
  // the properties to be set in the parent class are internally
  // consistent.
  ArchivingNode::set_status(__d);

  // if we get here, temporaries contain consistent set of properties
  set_a(tmp_a);
  set_b(tmp_b);
  set_c(tmp_c);
  set_d(tmp_d);
  set_I_e(tmp_I_e);
  set_V_th(tmp_V_th);
  set_V_min(tmp_V_min);
  set_Kplus(tmp_Kplus);
  set_V_r(tmp_V_r);
  set_V_T(tmp_V_T);
  set_C_m(tmp_C_m);
  set_tau_decay_msn(tmp_tau_decay_msn);
  set_tau_decay_snr(tmp_tau_decay_snr);
  set_E_rev_msn(tmp_E_rev_msn);
  set_E_rev_snr(tmp_E_rev_snr);
  set_V_m(tmp_V_m);
  set_U_m(tmp_U_m);
  set_I(tmp_I);
  set_g_snr__X__spikes_snr(tmp_g_snr__X__spikes_snr);
  set_g_msn__X__spikes_msn(tmp_g_msn__X__spikes_msn);


  updateValue< double >(__d, nest::names::gsl_error_tol, P_.__gsl_error_tol);
  if ( P_.__gsl_error_tol <= 0. )
  {
    throw nest::BadProperty( "The gsl_error_tol must be strictly positive." );
  }
};

#endif /* #ifndef MSN */
