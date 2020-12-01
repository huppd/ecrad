file(GLOB ifsaux_SRC "*.F90" CONFIGURE_DEPENDS)

target_sources(ifsaux PUBLIC ${ifsaux_SRC})
