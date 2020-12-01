file(GLOB utilities_SRC "*.F90" CONFIGURE_DEPENDS)

target_sources(utilities PUBLIC ${utilities_SRC})

