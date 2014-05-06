# Copyright 2014, Léo Baudouin
#
# This is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# This software is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details. You should have
# received a copy of the GNU General Public License along with
# software. If not, see <http://www.gnu.org/licenses/>.

#Create targets detector
ADD_EXECUTABLE(sample-targets
	sample-targets.cpp
)
TARGET_LINK_LIBRARIES( sample-targets ${PROJECT_NAME})

#Create grid detector
ADD_EXECUTABLE(sample-grid
	sample-grid.cpp
)
TARGET_LINK_LIBRARIES( sample-grid ${PROJECT_NAME})