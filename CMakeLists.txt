cmake_minimum_required (VERSION 3.8)

project("permutation-decoding")

file(GLOB CPP_SOURCES CONFIGURE_DEPENDS
	"src/codec/*.cpp"
	"src/construct/*.cpp"
	"src/math/*.cpp"
	"src/running/*.cpp"
	"src/utils/*.cpp")

add_executable(Decoder "src/decoder_main.cpp")
add_executable(Builder "src/builder_main.cpp")

target_sources(Decoder PRIVATE ${CPP_SOURCES})
target_sources(Builder PRIVATE ${CPP_SOURCES})

target_include_directories(Decoder PRIVATE "src")
target_include_directories(Builder PRIVATE "src")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET Decoder PROPERTY CXX_STANDARD 20)
  set_property(TARGET Builder PROPERTY CXX_STANDARD 20)
endif()
