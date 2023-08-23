#ifndef JIVE_IMPLICT_BFGSMODULE_H
#define JIVE_IMPLICT_BFGSMODULE_H

#include <jem/base/Flags.h>
#include <jive/implict/SolverModule.h>


JIVE_BEGIN_PACKAGE( implict )


class SolverBounds;


//-----------------------------------------------------------------------
//   class BFGSModule
//-----------------------------------------------------------------------

/** @brief 
 *  The BFGSModule class implements the L-BFGS quasi-Newton solution 
 *  technique.
 */ 

class BFGSModule : public SolverModule
{
 public:

  JEM_DECLARE_CLASS       ( BFGSModule, SolverModule );

  static const char*        TYPE_NAME;

  enum                      Option
  {
                              LINE_SEARCH = 1 << 0,
                              DELTA_CONS  = 1 << 1
  };

  typedef
    jem::Flags<Option>      Options;


  explicit                  BFGSModule

    ( const String&           name = "bfgs" );

  virtual Status            init

    ( const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat )            override;

  virtual void              shutdown

    ( const Properties&       globdat )            override;

  virtual void              configure

    ( const Properties&       props,
      const Properties&       globdat )            override;

  virtual void              getConfig

    ( const Properties&       conf,
      const Properties&       globdat )      const override;

  virtual void              advance

    ( const Properties&       globdat )            override;

  virtual void              solve

    ( const Properties&       info,
      const Properties&       globdat )            override;

  virtual void              cancel

    ( const Properties&       globdat )            override;

  virtual bool              commit

    ( const Properties&       globdat )            override;

  virtual void              setPrecision

    ( double                  eps )                override;

  virtual double            getPrecision  () const override;

  void                      setOption

    ( Option                  option,
      bool                    yesno = true );

  void                      setOptions

    ( Options                 options );

  Options                   getOptions    () const;

  virtual void              setReformIter 

    ( idx_t reformIter );

  virtual idx_t             getReformIter () const;


  virtual void              setMaxIter 

    ( idx_t maxIter );

  virtual idx_t             getMaxIter () const;


  static Ref<Module>        makeNew

    ( const String&           name,
      const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat );

  static void               declare       ();


 protected:

  virtual                  ~BFGSModule    ();


 private:

  class                     RunData_;
  class                     Work_;

  friend class              RunData_;
  friend class              Work_;


  bool                      solve_

    ( Work_&                  work,
      const Properties&       globdat );

  void                      lineSearch_

    ( Work_&                  work,
      const Properties&       globdat );

  void                      lineSearch2_

    ( Work_&                  work,
      const Properties&       globdat );

  void                      lineSearch3_  // negative line search

    ( Work_&                  work,
      const Properties&       globdat );   


 private:

  idx_t                     maxIter_;
  idx_t                     reformIter_;
  Options                   options_;

  double                    tiny_;
  double                    precision_;
  double                    lsearchTol_;
  double                    maxIncr_;
  double                    maxResIncr_;

  Ref<SolverBounds>         bounds_;
  Ref<RunData_>             rundat_;

};


JEM_DEFINE_FLAG_OPS( BFGSModule::Options )

JIVE_END_PACKAGE( implict )

#endif
