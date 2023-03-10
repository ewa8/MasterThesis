/**
 *  neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml.h
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
 *  Generated from NESTML at time: 2022-11-16 08:28:56.211551
**/

#ifndef NEUROMODULATED_STDP7572C325C24847F4AFC90A8B52D50E81_NESTML__WITH_MSN7572C325C24847F4AFC90A8B52D50E81_NESTML_H
#define NEUROMODULATED_STDP7572C325C24847F4AFC90A8B52D50E81_NESTML__WITH_MSN7572C325C24847F4AFC90A8B52D50E81_NESTML_H

// C++ includes:
#include <cmath>

// Includes from nestkernel:
#include "common_synapse_properties.h"
#include "connection.h"
#include "connector_model.h"
#include "event.h"
// Includes for volume transmitter
#include "volume_transmitter.h"


// Includes from sli:
#include "dictdatum.h"
#include "dictutils.h"

/** @BeginDocumentation

**/

#define POST_NEURON_TYPE msn7572c325c24847f4afc90a8b52d50e81_nestml__with_neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml

//#define DEBUG

namespace nest
{

namespace neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml_names
{
    const Name _w( "w" );
    const Name _n( "n" );
    const Name _c( "c" );
    const Name _pre_tr( "pre_tr" );
    const Name _the_delay( "the_delay" );
    const Name _tau_tr_pre( "tau_tr_pre" );
    const Name _tau_tr_post( "tau_tr_post" );
    const Name _tau_c( "tau_c" );
    const Name _tau_n( "tau_n" );
    const Name _b( "b" );
    const Name _Wmax( "Wmax" );
    const Name _Wmin( "Wmin" );
    const Name _A_plus( "A_plus" );
    const Name _A_minus( "A_minus" );
    const Name _A_vt( "A_vt" );
}

class neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestmlCommonSynapseProperties : public CommonSynapseProperties {
public:

    neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestmlCommonSynapseProperties()
    : CommonSynapseProperties()
    {
    }

    /**
     * Get all properties and put them into a dictionary.
     */
    void get_status( DictionaryDatum& d ) const
    {
        CommonSynapseProperties::get_status( d );
    }


    /**
     * Set properties from the values given in dictionary.
     */
    void set_status( const DictionaryDatum& d, ConnectorModel& cm )
    {
      CommonSynapseProperties::set_status( d, cm );
      long vtnode_id;
      if ( updateValue< long >( d, names::vt, vtnode_id ) )
      {
        const thread tid = kernel().vp_manager.get_thread_id();
        Node* vt = kernel().node_manager.get_node_or_proxy( vtnode_id, tid );
        vt_ = dynamic_cast< volume_transmitter* >( vt );
        if ( vt_ == 0 )
        {
          throw BadProperty( "Neuromodulatory source must be volume transmitter" );
        }
      }
    }

    // N.B.: we define all parameters as public for easy reference conversion later on.
    // This may or may not benefit performance (TODO: compare with inline getters/setters)
    volume_transmitter* vt_;

    inline long get_vt_node_id() const
    {
      if ( vt_ != 0 )
      {
        return vt_->get_node_id();
      }
      else
      {
        return -1;
      }
    }
};


template < typename targetidentifierT >
class neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml : public Connection< targetidentifierT >
{
public:
  void trigger_update_weight( thread t,
    const std::vector< spikecounter >& vt_spikes,
    double t_trig,
    const neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestmlCommonSynapseProperties& cp );
private:
  double t_lastspike_;
  // time of last update, which is either time of last presyn. spike or time-driven update
  double t_last_update_;

  // vt_spikes_idx_ refers to the vt spike that has just been processed after trigger_update_weight
  // a pseudo vt spike at t_trig is stored at index 0 and vt_spikes_idx_ = 0
  index vt_spikes_idx_;

  /**
   * Dynamic state of the synapse.
   *
   * These are the state variables that are advanced in time by calls to
   * send(). In many models, some or all of them can be set by the user
   * through ``SetStatus()``.
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
  struct State_{    
    double w;
    //!  Neuromodulator concentration
    double n;
    //!  Eligibility trace
    double c;
    double pre_tr;

    State_() {};
  };

  /**
   * Free parameters of the synapse.
   *
   *
   *
   * These are the parameters that can be set by the user through @c SetStatus.
   * Parameters do not change during calls to ``send()`` and are not reset by
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
  */
  struct Parameters_{    
    //!  !!! cannot have a variable called "delay"
    double the_delay;
    //!  STDP time constant for weight changes caused by pre-before-post spike pairings.
    double tau_tr_pre;
    //!  STDP time constant for weight changes caused by post-before-pre spike pairings.
    double tau_tr_post;
    //!  Time constant of eligibility trace
    double tau_c;
    //!  Time constant of dopaminergic trace
    double tau_n;
    //!  Dopaminergic baseline concentration
    double b;
    //!  Maximal synaptic weight
    double Wmax;
    //!  Minimal synaptic weight
    double Wmin;
    //!  Multiplier applied to weight changes caused by pre-before-post spike pairings. If b (dopamine baseline concentration) is zero, then A_plus is simply the multiplier for facilitation (as in the stdp_synapse model). If b is not zero, then A_plus will be the multiplier for facilitation only if n - b is positive, where n is the instantenous dopamine concentration in the volume transmitter. If n - b is negative, A_plus will be the multiplier for depression.
    double A_plus;
    //!  Multiplier applied to weight changes caused by post-before-pre spike pairings. If b (dopamine baseline concentration) is zero, then A_minus is simply the multiplier for depression (as in the stdp_synapse model). If b is not zero, then A_minus will be the multiplier for depression only if n - b is positive, where n is the instantenous dopamine concentration in the volume transmitter. If n - b is negative, A_minus will be the multiplier for facilitation.
    double A_minus;
    //!  Multiplier applied to dopa spikes
    double A_vt;



    /** Initialize parameters to their default values. */
    Parameters_() {};
  };

  /**
   * Internal variables of the synapse.
   *
   *
   *
   * These variables must be initialized by recompute_internal_variables().
  **/
  struct Variables_
  {    
    double tau_s;    
    double __h;    
    double __P__pre_tr__pre_tr;
  };

  Parameters_ P_;  //!< Free parameters.
  State_      S_;  //!< Dynamic state.
  Variables_  V_;  //!< Internal Variables

  // -------------------------------------------------------------------------
  //   Getters/setters for state block
  // -------------------------------------------------------------------------

  inline double get_w() const
  {
    return S_.w;
  }

  inline void set_w(const double __v)
  {
    S_.w = __v;
  }

  inline double get_n() const
  {
    return S_.n;
  }

  inline void set_n(const double __v)
  {
    S_.n = __v;
  }

  inline double get_c() const
  {
    return S_.c;
  }

  inline void set_c(const double __v)
  {
    S_.c = __v;
  }

  inline double get_pre_tr() const
  {
    return S_.pre_tr;
  }

  inline void set_pre_tr(const double __v)
  {
    S_.pre_tr = __v;
  }


  // -------------------------------------------------------------------------
  //   Getters/setters for parameters
  // -------------------------------------------------------------------------

  inline double get_the_delay() const
  {
    return P_.the_delay;
  }

  inline void set_the_delay(const double __v)
  {
    P_.the_delay = __v;
  }

  inline double get_tau_tr_pre() const
  {
    return P_.tau_tr_pre;
  }

  inline void set_tau_tr_pre(const double __v)
  {
    P_.tau_tr_pre = __v;
  }

  inline double get_tau_tr_post() const
  {
    return P_.tau_tr_post;
  }

  inline void set_tau_tr_post(const double __v)
  {
    P_.tau_tr_post = __v;
  }

  inline double get_tau_c() const
  {
    return P_.tau_c;
  }

  inline void set_tau_c(const double __v)
  {
    P_.tau_c = __v;
  }

  inline double get_tau_n() const
  {
    return P_.tau_n;
  }

  inline void set_tau_n(const double __v)
  {
    P_.tau_n = __v;
  }

  inline double get_b() const
  {
    return P_.b;
  }

  inline void set_b(const double __v)
  {
    P_.b = __v;
  }

  inline double get_Wmax() const
  {
    return P_.Wmax;
  }

  inline void set_Wmax(const double __v)
  {
    P_.Wmax = __v;
  }

  inline double get_Wmin() const
  {
    return P_.Wmin;
  }

  inline void set_Wmin(const double __v)
  {
    P_.Wmin = __v;
  }

  inline double get_A_plus() const
  {
    return P_.A_plus;
  }

  inline void set_A_plus(const double __v)
  {
    P_.A_plus = __v;
  }

  inline double get_A_minus() const
  {
    return P_.A_minus;
  }

  inline void set_A_minus(const double __v)
  {
    P_.A_minus = __v;
  }

  inline double get_A_vt() const
  {
    return P_.A_vt;
  }

  inline void set_A_vt(const double __v)
  {
    P_.A_vt = __v;
  }



  // -------------------------------------------------------------------------
  //   Getters/setters for inline expressions
  // -------------------------------------------------------------------------
  

  // -------------------------------------------------------------------------
  //   Function declarations
  // -------------------------------------------------------------------------



  /**
   * Update internal state (``S_``) of the synapse according to the dynamical equations defined in the model and the statements in the ``update`` block.
  **/
  inline void
  update_internal_state_(double t_start, double timestep, const neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestmlCommonSynapseProperties& cp);

  void recompute_internal_variables();



public:
  // this line determines which common properties to use
  typedef neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestmlCommonSynapseProperties CommonPropertiesType;

  typedef Connection< targetidentifierT > ConnectionBase;

  /**
  * Default constructor.
  *
  * Sets default values for all parameters (skipping common properties).
  *
  * Needed by GenericConnectorModel.
  */
  neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml();

  /**
  * Copy constructor from a property object.
  *
  * Sets default values for all parameters (skipping common properties).
  *
  * Needs to be defined properly in order for GenericConnector to work.
  */
  neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml( const neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml& rhs );
  void process_mod_spikes_spikes_( const std::vector< spikecounter >& vt_spikes,
      double t0,
      double t1,
      const neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestmlCommonSynapseProperties& cp );

  // Explicitly declare all methods inherited from the dependent base
  // ConnectionBase. This avoids explicit name prefixes in all places these
  // functions are used. Since ConnectionBase depends on the template parameter,
  // they are not automatically found in the base class.
  using ConnectionBase::get_delay_steps;
  using ConnectionBase::set_delay_steps;
  using ConnectionBase::get_delay;
  using ConnectionBase::set_delay;
  using ConnectionBase::get_rport;
  using ConnectionBase::get_target;


  class ConnTestDummyNode : public ConnTestDummyNodeBase
  {
  public:
    // Ensure proper overriding of overloaded virtual functions.
    // Return values from functions are ignored.
    using ConnTestDummyNodeBase::handles_test_event;
    port
    handles_test_event( SpikeEvent&, rport ) override
    {
      return invalid_port_;
    }
    port
    handles_test_event( RateEvent&, rport ) override
    {
      return invalid_port_;    }
    port
    handles_test_event( DataLoggingRequest&, rport ) override
    {
      return invalid_port_;    }
    port
    handles_test_event( CurrentEvent&, rport ) override
    {
      return invalid_port_;    }
    port
    handles_test_event( ConductanceEvent&, rport ) override
    {
      return invalid_port_;    }
    port
    handles_test_event( DoubleDataEvent&, rport ) override
    {
      return invalid_port_;    }
    port
    handles_test_event( DSSpikeEvent&, rport ) override
    {
      return invalid_port_;    }
    port
    handles_test_event( DSCurrentEvent&, rport ) override
    {
      return invalid_port_;    }
  };

  /**
   *  special case for weights in NEST: only in case a NESTML state variable was decorated by @nest::weight
  **/
  inline void set_weight(double w)
  {

    // no variable was decorated by @nest::weight, so no "weight" defined from the NEST perspective
    assert(0);
  }


/*  inline double get_weight() {
  }
*/
  void
  check_connection( Node& s,
    Node& t,
    rport receptor_type,
    const CommonPropertiesType& cp )
  {
    ConnTestDummyNode dummy_target;
    ConnectionBase::check_connection_( dummy_target, s, t, receptor_type );
    try {
      dynamic_cast<msn7572c325c24847f4afc90a8b52d50e81_nestml__with_neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml&>(t);
    }
    catch (std::bad_cast &exp) {
      std::cout << "wrong type of neuron connected! Synapse 'neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml' will only work with neuron 'msn7572c325c24847f4afc90a8b52d50e81_nestml__with_neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml'.\n";
      exit(1);
    }

    if ( cp.vt_ == 0 )
    {
      throw BadProperty( "No volume transmitter has been assigned to the dopamine synapse." );
    }

    t.register_stdp_connection( t_lastspike_ - get_delay(), get_delay() );
  }

  void
  send( Event& e, const thread tid, const neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestmlCommonSynapseProperties& cp )
  {
    const double __resolution = nest::Time::get_resolution().get_ms();  // do not remove, this is necessary for the resolution() function

    auto get_thread = [tid]()
    {
        return tid;
    };

    // synapse STDP depressing/facilitation dynamics
    const double __t_spike = e.get_stamp().get_ms();
#ifdef DEBUG
    std::cout << "neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml::send(): handling pre spike at t = " << __t_spike << std::endl;
#endif
  // get history of volume transmitter spikes
  const std::vector< spikecounter >& vt_spikes = cp.vt_->deliver_spikes();
    // use accessor functions (inherited from Connection< >) to obtain delay and target
    msn7572c325c24847f4afc90a8b52d50e81_nestml__with_neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml* __target = static_cast<msn7572c325c24847f4afc90a8b52d50e81_nestml__with_neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml*>(get_target(tid));
    assert(__target);
    const double __dendritic_delay = get_delay();
    const bool pre_before_post_update = 0;
    bool pre_before_post_flag = false;

    if (t_lastspike_ < 0.)
    {
        // this is the first presynaptic spike to be processed
        t_lastspike_ = 0.;
    }
    double timestep = 0;

    {
      // get spike history in relevant range (t1, t2] from post-synaptic neuron
      std::deque< histentry__msn7572c325c24847f4afc90a8b52d50e81_nestml__with_neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml >::iterator start;
      std::deque< histentry__msn7572c325c24847f4afc90a8b52d50e81_nestml__with_neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml >::iterator finish;
      double t0 = t_last_update_;
      // For a new synapse, t_lastspike_ contains the point in time of the last
      // spike. So we initially read the
      // history(t_last_spike - dendritic_delay, ..., T_spike-dendritic_delay]
      // which increases the access counter for these entries.
      // At registration, all entries' access counters of
      // history[0, ..., t_last_spike - dendritic_delay] have been
      // incremented by Archiving_Node::register_stdp_connection(). See bug #218 for
      // details.
      __target->get_history__( t_lastspike_ - __dendritic_delay,
        __t_spike - __dendritic_delay,
        &start,
        &finish );
      // facilitation due to post-synaptic spikes since last pre-synaptic spike
      while ( start != finish )
      {
        process_mod_spikes_spikes_( vt_spikes, t0, start->t_ + __dendritic_delay, cp );
        t0 = start->t_ + __dendritic_delay;
        const double minus_dt = t_lastspike_ - ( start->t_ + __dendritic_delay );
        // get_history() should make sure that
        // start->t_ > t_lastspike_ - dendritic_delay, i.e. minus_dt < 0
        assert( minus_dt < -kernel().connection_manager.get_stdp_eps() );

        if (pre_before_post_update and start->t_ == __t_spike - __dendritic_delay)
        {
          pre_before_post_flag = true;
          break;  // this would in any case have been the last post spike to be processed
        }

#ifdef DEBUG
        std::cout << "\tprocessing post spike at t = " << start->t_ << std::endl;
#endif

        /**
         * update synapse internal state from `t_lastspike_` to `start->t_`
        **/

        update_internal_state_(t_lastspike_, (start->t_ + __dendritic_delay) - t_lastspike_, cp);

        timestep += (start->t_ + __dendritic_delay) - t_lastspike_;

        const double _tr_t = start->t_;      
      /**
       *  NESTML generated onReceive code block for postsynaptic port "post_spikes" begins here!
      **/
      S_.c += P_.A_plus * get_pre_tr();

        /**
         * internal state has now been fully updated to `start->t_ + __dendritic_delay`
        **/

        t_lastspike_ = start->t_ + __dendritic_delay;
        ++start;
      }
    }

    /**
     * update synapse internal state from `t_lastspike_` to `__t_spike`
    **/
    process_mod_spikes_spikes_( vt_spikes, t_lastspike_, __t_spike, cp );

    update_internal_state_(t_lastspike_, __t_spike - t_lastspike_, cp);

    const double _tr_t = __t_spike - __dendritic_delay;

#ifdef DEBUG
    std::cout << "\tDepressing, old w = " << S_.w << "\n";
#endif
    //std::cout << "r2 = " << get_tr_r2() << std::endl;    
    /**
     *  NESTML generated onReceive code block for presynaptic port "pre_spikes" begins here!
    **/
    S_.pre_tr += 1.0;
    S_.c -= P_.A_minus * ((POST_NEURON_TYPE*)(__target))->get_post_tr__for_neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml(_tr_t);

            set_delay( P_.the_delay );
            const long __delay_steps = nest::Time::delay_ms_to_steps( get_delay() );
            set_delay_steps(__delay_steps);
            e.set_receiver( *__target );
      e.set_weight( get_w() );
      // use accessor functions (inherited from Connection< >) to obtain delay in steps and rport
      e.set_delay_steps( get_delay_steps() );
      e.set_rport( get_rport() );
    e();
    ;

#ifdef DEBUG
    std::cout <<"\t-> new w = " << S_.w << std::endl;
#endif

    /**
     *  update all convolutions with pre spikes
    **/



    /**
     *  in case pre and post spike time coincide and pre update takes priority
    **/

    if (pre_before_post_flag)
    {      
      /**
       *  NESTML generated onReceive code block for postsynaptic port "post_spikes" begins here!
      **/
      #ifdef DEBUG
      std::cout << "\tFacilitating from c = " << S_.c << " (using trace = " << S_.pre_tr << ")";
      #endif
      S_.c += P_.A_plus * get_pre_tr();
#ifdef DEBUG
      std::cout << " to " << S_.c << std::endl;
#endif
    }

    /**
     *  synapse internal state has now been fully updated to `__t_spike`
    **/

    t_lastspike_ = __t_spike;
  }

  void get_status( DictionaryDatum& d ) const;

  void set_status( const DictionaryDatum& d, ConnectorModel& cm );
};
template < typename targetidentifierT >
void
neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml< targetidentifierT >::process_mod_spikes_spikes_( const std::vector< spikecounter >& vt_spikes,
    double t0,
    double t1,
    const neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestmlCommonSynapseProperties& cp )
{
#ifdef DEBUG
  std::cout << "\tIn process_mod_spikes_spikes_(): t0 = " << t0 << ", t1 = " << t1 << "\n";
#endif
  // process dopa spikes in (t0, t1]
  // propagate weight from t0 to t1
  if ( ( vt_spikes.size() > vt_spikes_idx_ + 1 )
    and ( t1 - vt_spikes[ vt_spikes_idx_ + 1 ].spike_time_ > -1.0 * kernel().connection_manager.get_stdp_eps() ) )
  {
    // there is at least 1 dopa spike in (t0, t1]
    // propagate up to first dopa spike
#ifdef DEBUG
  std::cout << "\t\tHandling (1) spike at t = " << vt_spikes[ vt_spikes_idx_ +1].spike_time_ << "\n";
#endif
    update_internal_state_(t0, vt_spikes[ vt_spikes_idx_ + 1 ].spike_time_ - t0, cp );
    ++vt_spikes_idx_;
#ifdef DEBUG
std::cout<<"\t\tIncrementing n_ from " << S_.n << " to " ;
#endif
    /**
     *  NESTML generated onReceive code block for volume transmitter synaptic port "mod_spikes" begins here!
    **/    
    S_.n += P_.A_vt / P_.tau_n;
#ifdef DEBUG
std::cout << S_.n << "\n";
#endif
    // process remaining dopa spikes in (t0, t1]
    double cd;
    while ( ( vt_spikes.size() > vt_spikes_idx_ + 1 )
      and ( t1 - vt_spikes[ vt_spikes_idx_ + 1 ].spike_time_ > -1.0 * kernel().connection_manager.get_stdp_eps() ) )
    {
#ifdef DEBUG
  std::cout << "\t\tHandling (2) spike at t = " << vt_spikes[ vt_spikes_idx_ ].spike_time_ << "\n";
#endif
      // propagate up to next dopa spike
      update_internal_state_(vt_spikes[ vt_spikes_idx_ ].spike_time_,
                             vt_spikes[ vt_spikes_idx_ + 1 ].spike_time_ - vt_spikes[ vt_spikes_idx_ ].spike_time_,
                             cp );
      ++vt_spikes_idx_;

      /**
       *  NESTML generated onReceive code block for volume transmitter synaptic port "mod_spikes" begins here!
      **/
#ifdef DEBUG
std::cout<<"\t\tIncrementing n_ from " << S_.n << " to " ;
#endif      
      S_.n += P_.A_vt / P_.tau_n;
#ifdef DEBUG
std::cout << S_.n << "\n";
#endif
    }
#ifdef DEBUG
  std::cout << "\t\t3\n";
#endif

    // propagate up to t1
    update_internal_state_(vt_spikes[ vt_spikes_idx_ ].spike_time_,
                           t1 - vt_spikes[ vt_spikes_idx_ ].spike_time_,
                           cp );
  }
  else
  {
#ifdef DEBUG
  std::cout << "\t\t4: updating internal state from t0 = " << t0 << " to t1 = " << t1 << "\n";
#endif

    // no dopamine spikes in (t0, t1]
    update_internal_state_( t0, t1 - t0, cp );
  }
}


template < typename targetidentifierT >
void
neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml< targetidentifierT >::get_status( DictionaryDatum& __d ) const
{
  ConnectionBase::get_status( __d );
  def< long >( __d, names::size_of, sizeof( *this ) );

  // parameters
  def< double >( __d, names::delay, P_.the_delay );  
  def<double>(__d, nest::neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml_names::_tau_tr_pre, get_tau_tr_pre());  
  def<double>(__d, nest::neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml_names::_tau_tr_post, get_tau_tr_post());  
  def<double>(__d, nest::neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml_names::_tau_c, get_tau_c());  
  def<double>(__d, nest::neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml_names::_tau_n, get_tau_n());  
  def<double>(__d, nest::neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml_names::_b, get_b());  
  def<double>(__d, nest::neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml_names::_Wmax, get_Wmax());  
  def<double>(__d, nest::neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml_names::_Wmin, get_Wmin());  
  def<double>(__d, nest::neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml_names::_A_plus, get_A_plus());  
  def<double>(__d, nest::neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml_names::_A_minus, get_A_minus());  
  def<double>(__d, nest::neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml_names::_A_vt, get_A_vt());

  // initial values for state variables in ODE or kernel  
  def<double>(__d, nest::neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml_names::_w, get_w());
  def<double>(__d, nest::neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml_names::_n, get_n());
  def<double>(__d, nest::neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml_names::_c, get_c());
  def<double>(__d, nest::neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml_names::_pre_tr, get_pre_tr());
}

template < typename targetidentifierT >
void
neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml< targetidentifierT >::set_status( const DictionaryDatum& __d,
  ConnectorModel& cm )
{
  // parameters  
  double tmp_the_delay = get_the_delay();
  updateValue<double>(__d, nest::neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml_names::_the_delay, tmp_the_delay);  
  double tmp_tau_tr_pre = get_tau_tr_pre();
  updateValue<double>(__d, nest::neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml_names::_tau_tr_pre, tmp_tau_tr_pre);  
  double tmp_tau_tr_post = get_tau_tr_post();
  updateValue<double>(__d, nest::neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml_names::_tau_tr_post, tmp_tau_tr_post);  
  double tmp_tau_c = get_tau_c();
  updateValue<double>(__d, nest::neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml_names::_tau_c, tmp_tau_c);  
  double tmp_tau_n = get_tau_n();
  updateValue<double>(__d, nest::neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml_names::_tau_n, tmp_tau_n);  
  double tmp_b = get_b();
  updateValue<double>(__d, nest::neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml_names::_b, tmp_b);  
  double tmp_Wmax = get_Wmax();
  updateValue<double>(__d, nest::neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml_names::_Wmax, tmp_Wmax);  
  double tmp_Wmin = get_Wmin();
  updateValue<double>(__d, nest::neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml_names::_Wmin, tmp_Wmin);  
  double tmp_A_plus = get_A_plus();
  updateValue<double>(__d, nest::neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml_names::_A_plus, tmp_A_plus);  
  double tmp_A_minus = get_A_minus();
  updateValue<double>(__d, nest::neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml_names::_A_minus, tmp_A_minus);  
  double tmp_A_vt = get_A_vt();
  updateValue<double>(__d, nest::neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml_names::_A_vt, tmp_A_vt);

  // initial values for state variables in ODE or kernel  
  double tmp_w = get_w();
  updateValue<double>(__d, nest::neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml_names::_w, tmp_w);
  double tmp_n = get_n();
  updateValue<double>(__d, nest::neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml_names::_n, tmp_n);
  double tmp_c = get_c();
  updateValue<double>(__d, nest::neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml_names::_c, tmp_c);
  double tmp_pre_tr = get_pre_tr();
  updateValue<double>(__d, nest::neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml_names::_pre_tr, tmp_pre_tr);


  // We now know that (ptmp, stmp) are consistent. We do not
  // write them back to (P_, S_) before we are also sure that
  // the properties to be set in the parent class are internally
  // consistent.
  ConnectionBase::set_status( __d, cm );

  // if we get here, temporaries contain consistent set of properties
  // set parameters  
  set_the_delay(tmp_the_delay);  
  set_tau_tr_pre(tmp_tau_tr_pre);  
  set_tau_tr_post(tmp_tau_tr_post);  
  set_tau_c(tmp_tau_c);  
  set_tau_n(tmp_tau_n);  
  set_b(tmp_b);  
  set_Wmax(tmp_Wmax);  
  set_Wmin(tmp_Wmin);  
  set_A_plus(tmp_A_plus);  
  set_A_minus(tmp_A_minus);  
  set_A_vt(tmp_A_vt);

  // set state  
  set_w(tmp_w);  
  set_n(tmp_n);  
  set_c(tmp_c);  
  set_pre_tr(tmp_pre_tr);

  // check invariants



  // special treatment of NEST delay
  set_delay(
get_the_delay());

  // recompute internal variables in case they are dependent on parameters or state that might have been updated in this call to set_status()
  recompute_internal_variables();
}

/**
 * NESTML internals block symbols initialisation
**/
template < typename targetidentifierT >
void neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml< targetidentifierT >::recompute_internal_variables()
{
  const double __resolution = nest::Time::get_resolution().get_ms();  // do not remove, this is necessary for the resolution() function


    
    V_.tau_s = (P_.tau_c + P_.tau_n) / (P_.tau_c * P_.tau_n); // as 1 / ms
    
    V_.__P__pre_tr__pre_tr = std::exp((-(V_.__h)) / P_.tau_tr_pre); // as real
}

/**
 * constructor
**/
template < typename targetidentifierT >
neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml< targetidentifierT >::neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml() : ConnectionBase()
{
  const double __resolution = nest::Time::get_resolution().get_ms();  // do not remove, this is necessary for the resolution() function
  
  P_.the_delay = 1; // as ms
  
  P_.tau_tr_pre = 20; // as ms
  
  P_.tau_tr_post = 20; // as ms
  
  P_.tau_c = 1000; // as ms
  
  P_.tau_n = 200; // as ms
  
  P_.b = 0.0; // as real
  
  P_.Wmax = 200.0; // as real
  
  P_.Wmin = 0.0; // as real
  
  P_.A_plus = 1.0; // as real
  
  P_.A_minus = 1.5; // as real
  
  P_.A_vt = 1.0; // as real

  V_.__h = nest::Time::get_resolution().get_ms();
  recompute_internal_variables();

  // initial values for state variables in ODE or kernel
  
  S_.w = 1.0; // as real
  
  S_.n = 0.0; // as real
  
  S_.c = 0.0; // as real
  
  S_.pre_tr = 0.0; // as real

  t_lastspike_ = 0.;
  t_last_update_ = 0.;
}

/**
 * copy constructor
**/
template < typename targetidentifierT >
neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml< targetidentifierT >::neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml( const neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml< targetidentifierT >& rhs )
: ConnectionBase( rhs )
{
    P_.the_delay = rhs.P_.the_delay;
    P_.tau_tr_pre = rhs.P_.tau_tr_pre;
    P_.tau_tr_post = rhs.P_.tau_tr_post;
    P_.tau_c = rhs.P_.tau_c;
    P_.tau_n = rhs.P_.tau_n;
    P_.b = rhs.P_.b;
    P_.Wmax = rhs.P_.Wmax;
    P_.Wmin = rhs.P_.Wmin;
    P_.A_plus = rhs.P_.A_plus;
    P_.A_minus = rhs.P_.A_minus;
    P_.A_vt = rhs.P_.A_vt;

  // state variables in ODE or kernel
    S_.w = rhs.S_.w;
    S_.n = rhs.S_.n;
    S_.c = rhs.S_.c;
    S_.pre_tr = rhs.S_.pre_tr;

    //weight_ = get_named_parameter<double>(names::weight);
    //set_weight( *rhs.weight_ );
    t_last_update_ = rhs.t_last_update_;
    t_lastspike_  = rhs.t_lastspike_;

    // special treatment of NEST delay
    set_delay(rhs.get_delay());
}

template < typename targetidentifierT >
inline void
neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml< targetidentifierT >::update_internal_state_(double t_start, double timestep, const neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestmlCommonSynapseProperties& cp)
{
  if (timestep < 1E-12)
  {
#ifdef DEBUG
    std::cout << "\tupdate_internal_state_() called with dt < 1E-12; skipping update\n" ;
#endif
    return;
  }

  const double __resolution = timestep;  // do not remove, this is necessary for the resolution() function

#ifdef DEBUG
  std::cout<< "\tUpdating internal state: t_start = " << t_start << ", dt = " << timestep << "\n";
#endif
  const double old___h = V_.__h;
  V_.__h = timestep;
  recompute_internal_variables();  
  double pre_tr__tmp = V_.__P__pre_tr__pre_tr * get_pre_tr();
  /* replace analytically solvable variables with precisely integrated values  */
  S_.pre_tr = pre_tr__tmp;
    V_.__h = old___h;
    recompute_internal_variables();  // XXX: can be skipped?
  S_.w -= get_c() * (get_n() / V_.tau_s * numerics::expm1((-(V_.tau_s)) * __resolution) - P_.b * P_.tau_c * numerics::expm1((-(__resolution)) / P_.tau_c));
  S_.w = std::max(0.0, get_w());
  S_.c = get_c() * std::exp((-(__resolution)) / P_.tau_c);
  S_.n = get_n() * std::exp((-(__resolution)) / P_.tau_n);
  t_last_update_ = t_start + timestep;
}
/**
 * Update to end of timestep ``t_trig``, while processing vt spikes and post spikes
**/
template < typename targetidentifierT >
inline void
neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml< targetidentifierT >::trigger_update_weight( thread t,
  const std::vector< spikecounter >& vt_spikes,
  const double t_trig,
  const CommonPropertiesType& cp )
{
  // propagate all state variables in the synapse to time t_trig
#ifdef DEBUG
    std::cout << "\nneuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml::trigger_update_weight(): t = " << t_trig << std::endl;
#endif
  // purely dendritic delay
  double dendritic_delay = get_delay();

  // get spike history in relevant range (t_last_update, t_trig] from postsyn. neuron
  std::deque< histentry__msn7572c325c24847f4afc90a8b52d50e81_nestml__with_neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml >::iterator start;
  std::deque< histentry__msn7572c325c24847f4afc90a8b52d50e81_nestml__with_neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml >::iterator finish;
  static_cast<msn7572c325c24847f4afc90a8b52d50e81_nestml__with_neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml*>(get_target(t))->get_history__( t_last_update_ - dendritic_delay, t_trig - dendritic_delay, &start, &finish );

  // facilitation due to postsyn. spikes since last update
  double t0 = t_last_update_;
  // double minus_dt;
  double timestep = 0;

  while ( start != finish )
  {
    process_mod_spikes_spikes_( vt_spikes, t0, start->t_ + dendritic_delay, cp );

#ifdef DEBUG
    std::cout << "\tprocessing post spike from " << t_last_update_ << " to " << start->t_ + dendritic_delay << std::endl;
#endif

    /**
     * update synapse internal state from `t_last_update_` to `start->t_`
    **/

    update_internal_state_(t_last_update_,
                           (start->t_ + dendritic_delay) - t_last_update_,
                           cp);

    const double _tr_t = start->t_;
#ifdef DEBUG
        std::cout << "\tFacilitating from c = " << S_.c << " (using trace = " << S_.pre_tr << ")";
#endif      
          /**
           *  NESTML generated onReceive code block for postsynaptic port "post_spikes" begins here!
          **/
      S_.c += P_.A_plus * get_pre_tr();

// #ifdef DEBUG
//   std::cout << "\t--> new w = " << S_.w << std::endl;
// #endif
#ifdef DEBUG
      std::cout << " to " << S_.c << std::endl;
#endif
    /**
     * internal state has now been fully updated to `start->t_ + dendritic_delay`
    **/

    t0 = start->t_ + dendritic_delay;
    // minus_dt = t_last_update_ - t0;
    t_lastspike_ = start->t_ + dendritic_delay;
    ++start;
  }

  /**
    * update synapse internal state from `t_lastspike_` to `t_trig`
  **/
  process_mod_spikes_spikes_( vt_spikes, t_lastspike_, t_trig, cp );

#ifdef DEBUG
    //std::cout << "neuromodulated_stdp7572c325c24847f4afc90a8b52d50e81_nestml__with_msn7572c325c24847f4afc90a8b52d50e81_nestml::trigger_update_weight(): \tupdating from " << t_lastspike_ << " to " << t_trig + dendritic_delay << std::endl;
#endif

 /* update_internal_state_(t_lastspike_,
                         t_trig - t_lastspike_,
                         cp);*/

  vt_spikes_idx_ = 0;
  t_lastspike_ = t_trig;
}

} // namespace

#endif /* #ifndef NEUROMODULATED_STDP7572C325C24847F4AFC90A8B52D50E81_NESTML__WITH_MSN7572C325C24847F4AFC90A8B52D50E81_NESTML_H */