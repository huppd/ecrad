file(GLOB utilities_SRC CONFIGURE_DEPENDS "*.F90")

target_sources(utilities PUBLIC ${utilities_SRC})

