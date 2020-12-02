file(GLOB radiation_SRC "*.F90" CONFIGURE_DEPENDS)

target_sources(radiation PUBLIC ${radiation_SRC})

