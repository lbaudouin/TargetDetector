# Copyright 2014, Léo Baudouin
#
# This is free library: you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details. You should have
# received a copy of the GNU General Public License along with
# library. If not, see <http://www.gnu.org/licenses/>.

ADD_EXECUTABLE(oneblob
        oneblobtest.cpp
)
TARGET_LINK_LIBRARIES( oneblob ${PROJECT_NAME})
ADD_TEST(OneBlob oneblob ${CMAKE_SOURCE_DIR}/unitTesting/1B_8h180_8m142_bin.png)


ADD_EXECUTABLE(threeblobs
        threeblobstest.cpp
)
TARGET_LINK_LIBRARIES( threeblobs ${PROJECT_NAME})
ADD_TEST(ThreeBlobs threeblobs ${CMAKE_SOURCE_DIR}/unitTesting/3B_8h180_8m142_bin.png)


ADD_EXECUTABLE(tworings
        tworingstest.cpp
)
TARGET_LINK_LIBRARIES( tworings ${PROJECT_NAME})
ADD_TEST(TwoRings tworings ${CMAKE_SOURCE_DIR}/unitTesting/2R_8h180_8m142_bin.png)

ADD_EXECUTABLE(testtypes
        testtypes.cpp
)
TARGET_LINK_LIBRARIES( testtypes ${PROJECT_NAME})
ADD_TEST(TestType testtypes ${CMAKE_SOURCE_DIR}/unitTesting/1B_8h180_8m142_bin.png ${CMAKE_SOURCE_DIR}/unitTesting/3B_8h180_8m142_bin.png ${CMAKE_SOURCE_DIR}/unitTesting/2R_8h180_8m142_bin.png)

ADD_EXECUTABLE(autothreshold
        autothresholdtest.cpp
)
TARGET_LINK_LIBRARIES( autothreshold ${PROJECT_NAME})
ADD_TEST(AutoThreshold autothreshold ${CMAKE_SOURCE_DIR}/unitTesting/1B_8h180_8m142_bin.png)


ADD_EXECUTABLE(findtargetsgrid
        findtargetsgrid.cpp
)
TARGET_LINK_LIBRARIES( findtargetsgrid ${PROJECT_NAME})
ADD_TEST(FindTargetGrid findtargetsgrid ${CMAKE_SOURCE_DIR}/unitTesting/grid_targets_5x5_50.png ${CMAKE_SOURCE_DIR}/unitTesting/grid_targets_4x4_50_unordered.png)