cmake_minimum_required(VERSION 3.26)
project(Option_Pricing)

set(CMAKE_CXX_STANDARD 17)

add_executable(Option_Pricing main.cpp
        option_pricing.cc)
