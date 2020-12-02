file(GLOB ifsaux_SRC CONFIGURE_DEPENDS "*.F90")

target_sources(ifsaux PUBLIC ${ifsaux_SRC})
