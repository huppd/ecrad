if(__claw_base_compiler)
  return()
endif()
set(__claw_base_compiler YES)

# For PGI compiler
set(FC_BASIC_FLAGS "-Mbyteswapio -fpic")
set(FC_FPPFLAGS "-Mpreprocess") 
set(FC_RELEASE_OPTFLAGS "-O3 -acc -ta=nvidia")
set(FC_DEBUG_OPTFLAGS "-O0")
set(FC_RELEASE_WARNFLAGS "-Minform=inform")
set(FC_DEBUG_WARNFLAGS "${FC_RELEASE_WARNFLAGS}")
set(FC_RELEASE_DEBUGFLAGS "-g")
set(FC_DEBUG_DEBUGFLAGS "${FC_RELEASE_DEBUGFLAGS}")
set(FC_OMPFLAG "-mp=nonuma")
set(FC_MODFLAG "-J")

