cmake_minimum_required(VERSION 3.8)
project(Floorplan)

set(CMAKE_CXX_STANDARD 14)

set(Data Data.hpp Instance.hpp)

set(Utils Utils.hpp Visualizer.hpp)

set(Algorithm
	AdaptiveSelecter.hpp
	FloorplanPacker.hpp
	BeamSearcher.hpp
	RandomLocalSearcher.hpp
)

# ���ù���Ŀ¼
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/Deploy)

add_executable(Floorplan 
	Main.cpp Config.hpp
	${Data} ${Algorithm} ${Utils}
)