#ifndef XNAMES_H
#define XNAMES_H

//-----------------------------------------------------------------------
//   class XProps (eXtended Property Names)
//-----------------------------------------------------------------------

/** @brief 
 *  The XProps class contain property names not present in the Jive library.
 */

class XProps
{
 public:

  static const char*    ADAPT_HOW;
  static const char*    ADAPT_FROM;
  static const char*    ALLOW_TMP;
  static const char*    ARCLEN;
  static const char*    ARCLEN_0;
  static const char*    INCR_DISSIPATION;
  static const char*    DISSIPATION_FORCE;
  static const char*    FE_DISSIPATION;
  static const char*    FE_UNIT_LOAD;
  static const char*    LINE_SEARCH_TOL;
  static const char*    LOAD_SCALE_0;  
  static const char*    MAX_RES_INCR;
  static const char*    N_CONTINUES;
  static const char*    REFORM_ITER;
  static const char*    STEP_SIZE;
  static const char*    STEP_SIZE_0;

};

//-----------------------------------------------------------------------
//   class XActions (eXtended Actions)
//-----------------------------------------------------------------------

/** @brief 
 *  The XActions class contain actions not present in the Jive library.
 */

class XActions
{
 public:

  static const char*    ACCEPT;
  static const char*    ADAPT_STEP;
  static const char*    ASSEM_INCR_DISS;
  static const char*    BE_CAREFUL;
  static const char*    CHANGE_COUNT;
  static const char*    CHECK_BOUNDS;
  static const char*    CHECK_COMMIT;
  static const char*    CONTINUE;
  static const char*    CONVERGED;
  static const char*    DISCARD;
  static const char*    DO_SWITCH;
  static const char*    DONE;
  static const char*    DONE_INCR_DISS;
  static const char*    GET_ARCLEN;
  static const char*    GET_DISS_FORCE;
  static const char*    GET_STEP_SIZE;
  static const char*    REDUCED;
  static const char*    REDUCE_STEP;
  static const char*    SET_ARCLEN;
  static const char*    SET_TAYLOR_HOOD;
  static const char*    SET_STEP_SIZE;
  static const char*    STOP_ARCLEN;
  static const char*    STOP_CAREFUL;
  static const char*    STORE_LOADSCALE;
  static const char*    TEMPORARY;
  static const char*    TERMINATE;
  static const char*    TO_ARCLEN;
  static const char*    TO_BFGS;
  static const char*    TO_DISP;
  static const char*    TO_NONLIN;

};

#endif
