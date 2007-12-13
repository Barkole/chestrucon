#!/usr/bin/awk -f

# Copyright (C) 2007 Michael Decker <http://Inspire-Mind.de>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Authors:
# * Michael Decker (mad) <http://Inspire-Mind.de>
# * Rebecca Römer (rr)

# WARNING:
# Works ONLY with CUBIC systems
# Will NOT WORK with TRIGONAL and HEXAGONAL systems

# Tested for:
# * MgAl2O4 (Fd-3m; No. 227) (cubic) (2007-11-14)
# * Ca3N2 (Ia-3; No. 206) (cubic) (2007-11-14)

# TODO: Test for
# * Tetragonal
# * Orthorhombic
# * Moniclinic
# * Triclinic

# Version:
# * 0.1.0.1: 2007-11-20 - Add missing matrix transpose (mad)
# * 0.1.0.0: 2007-11-11 - Fixing atom vectors (usage of invert matrix) (mad), tested and proved for cubic systems by (rr)
# * 0.0.0.3: 2007-11-08 - Fixing counting (mad)
# * 0.0.0.2: 2007-11-08 - Portation to awk (mad) - needs only < 0.03 seconds instead of 0.5-2.0 seconds
# * 0.0.0.1: 2007-10-11 - First running version - bash script (mad)

# Tested:
# * 2007-11-11: Kubuntu 7.04 with
#   - GNU Awk 3.1.5
# * 2007-11-08: Kubuntu 7.04 with
#   - GNU Awk 3.1.5
# * 2007-10-11: Kubuntu 7.04 with
#   - GNU bash, version 3.2.13(1)-release (i486-pc-linux-gnu)
#   - GNU Awk 3.1.5

## Doing before anything happens

BEGIN {
	## Constants
	# Define all states
	def_state_unknown = 0;
	def_state_header = 1;
	def_state_symgrp = 2;
	def_state_struc = 3;
	def_state_struc_plat = 4;
	def_state_class = 5;
	def_state_site = 6;
	# Output accuracy
	def_print_accuracy = 8;

	## All internal variables
	current_state = def_state_unknown;
	current_atom_type_state = "";
	atom_type_counter = 0;
	atom_type_count_summary[1] = 0
	atom_type_summary[1] = "";
	struc_plat_line_counter = 1;

	## All given by input file values
	value_alat = 0;
	value_angstr = 0;
	# The Matrix:
	# 1,1 1,2 1,3
	# 2,1 2,2 2,3
	# 3,1 3,2 3,3
	value_coord_matrix[1,1] = 0;
	value_inv_coord_matrix[1,1] = 0;
	value_inv_trans_coord_matrix[1,1] = 0;

	# Outputs
	output_atomcoordinates_counter = 1;
	output_atomcoordinates[1,1] = 0.0; # First is line, second xyz
	
}

## Define state changes

/^HEADER/ {
	# Set new state
	current_state = def_state_header;
}

/^SYMGRP/ {
	# Set new state
	current_state = def_state_symgrp;
}

/^STRUC/ {
	# Set new state
	current_state = def_state_struc
}

/PLAT/ && current_state  ==  def_state_struc {
	# Set new state
	current_state = def_state_struc_plat
}


/^CLASS/ {
	# Set new state
	current_state = def_state_class
}

/^SITE/ {
	# Set new state
	current_state = def_state_site
	## Invert Matrix for recoordnate atom vectors
	# calculate determinante
	det_matrix = value_coord_matrix[1,1]*value_coord_matrix[2,2]*value_coord_matrix[3,3] + value_coord_matrix[1,2]*value_coord_matrix[2,3]*value_coord_matrix[3,1] + value_coord_matrix[1,3]*value_coord_matrix[2,1]*value_coord_matrix[3,2] - value_coord_matrix[1,3]*value_coord_matrix[2,2]*value_coord_matrix[3,1] - value_coord_matrix[1,1]*value_coord_matrix[2,3]*value_coord_matrix[3,2] - value_coord_matrix[1,2]*value_coord_matrix[2,1]*value_coord_matrix[3,3];

	# calculate invert matrix
	value_inv_coord_matrix[1,1] = (value_coord_matrix[2,2]*value_coord_matrix[3,3] - value_coord_matrix[2,3]*value_coord_matrix[3,2])/det_matrix;
	value_inv_coord_matrix[1,2] = (value_coord_matrix[1,3]*value_coord_matrix[3,2] - value_coord_matrix[1,2]*value_coord_matrix[3,3])/det_matrix;
	value_inv_coord_matrix[1,3] = (value_coord_matrix[1,2]*value_coord_matrix[2,3] - value_coord_matrix[1,3]*value_coord_matrix[2,2])/det_matrix;
	value_inv_coord_matrix[2,1] = (value_coord_matrix[2,3]*value_coord_matrix[3,1] - value_coord_matrix[2,1]*value_coord_matrix[3,3])/det_matrix;
	value_inv_coord_matrix[2,2] = (value_coord_matrix[1,1]*value_coord_matrix[3,3] - value_coord_matrix[1,3]*value_coord_matrix[3,1])/det_matrix;
	value_inv_coord_matrix[2,3] = (value_coord_matrix[1,3]*value_coord_matrix[2,1] - value_coord_matrix[1,1]*value_coord_matrix[2,3])/det_matrix;
	value_inv_coord_matrix[3,1] = (value_coord_matrix[2,1]*value_coord_matrix[3,2] - value_coord_matrix[2,2]*value_coord_matrix[3,1])/det_matrix;
	value_inv_coord_matrix[3,2] = (value_coord_matrix[1,2]*value_coord_matrix[3,1] - value_coord_matrix[1,1]*value_coord_matrix[3,2])/det_matrix;
	value_inv_coord_matrix[3,3] = (value_coord_matrix[1,1]*value_coord_matrix[2,2] - value_coord_matrix[1,2]*value_coord_matrix[2,1])/det_matrix;

	# transpose matrix
	for (i=1; i<4; i++) {
		for (j=1; j<4; j++) {
			value_inv_trans_coord_matrix[j,i] = value_inv_coord_matrix[i,j];
		}
	}
}

## Predoings on every line

{
	# Skip empty lines
	if (length($0) == 0) {
		next;
	}
	# Get line without blanks and state information
	realline = substr($0,11);
}

## Work on lines at special states

current_state  ==  def_state_header {
	# Handle line in state HEADER
	# Print only the header information
	print realline
	# Another needed hardcoded line (function?)
	printf "%." def_print_accuracy "f\n", 1;
}

current_state  ==  def_state_symgrp {
	# Handle line in state SYMGRP
	# Nothing to, only ignore this lines
}

current_state  ==  def_state_struc {
	# Handle line in state STRUC
	# Get alat value
	split(realline,splittedValues,"=")
	value_alat = splittedValues[2];
	# Calculate ANGSTR (Angstroem)
	value_angstr = value_alat*0.5291772083;
}

current_state  ==  def_state_struc_plat {
	# Handle line in substate PLAT of state STRUC
	# Get cutted matrix
	lineVector = substr(realline,6);
	# Get each value
	split(lineVector,platValue," "); # Split by space
	value_coord_matrix[struc_plat_line_counter,1] = platValue[1];
	value_coord_matrix[struc_plat_line_counter,2] = platValue[2];
	value_coord_matrix[struc_plat_line_counter,3] = platValue[3];
	printf "%." def_print_accuracy "f %." def_print_accuracy "f %." def_print_accuracy "f\n", platValue[1]*value_angstr, platValue[2]*value_angstr, platValue[3]*value_angstr;
	struc_plat_line_counter++;
}

current_state  ==  def_state_class {
	# Handle line in state CLASS
	# Nothing to, only ignore this lines
}

current_state  ==  def_state_site {
	# Handle line in state SITE
	# Split by space
	split(realline,splittedRealline," ");
	## Get atom type
	lineAtomType = substr(splittedRealline[1],6);
	# Delete all numbers from atom type
	sub(/[[:digit:]]/, "", lineAtomType);
	# New type or what?
	if (lineAtomType == current_atom_type_state) {
		# Same atom type as before -> count it
		atom_type_count_summary[atom_type_counter]++;
	} else {
		# New atom type, so add it to array and set new counter area
		atom_type_counter++;
		atom_type_count_summary[atom_type_counter] = 1;
		atom_type_summary[atom_type_counter] = lineAtomType;
		current_atom_type_state = lineAtomType;
	}
	## Get atom coordinates
	# read vector
	atom_vector[1] = substr(splittedRealline[2],5);
	atom_vector[2] = splittedRealline[3];
	atom_vector[3] = splittedRealline[4];
	# Transform old vector into new space: multiply vector with matrix
	for (i=1; i<4; i++) {
		output_atomcoordinates[output_atomcoordinates_counter,i] = 0;
		for (j=1; j<4; j++) {
			output_atomcoordinates[output_atomcoordinates_counter,i] += value_inv_trans_coord_matrix[i,j]*atom_vector[j];
		}
	}
	# Print atom coordnate at end, because I need information how many atom of similar kind available. This information is only available, AFTER parsing all "old" atom coordinates
	output_atomcoordinates_counter++;
}

current_state  ==  def_state_unknown {
	# Unknown lines
	print "Warning unknown state - ignored this line: " $0;
}

## Finishing
END {
	## Print in state SITE created information
	# Print count of atoms
	printf "%i", atom_type_count_summary[1];
	for (i = 2; i <= atom_type_counter; i++ ) {
		printf " %i", atom_type_count_summary[i];
	}
	printf "\n";
	# Hardcoded order output
	print "Direct"
	# Print all atom coordinates
	for (i = 1; i < output_atomcoordinates_counter; i++) {
		printf "%." def_print_accuracy "f %." def_print_accuracy "f %." def_print_accuracy "f\n", output_atomcoordinates[i,1], output_atomcoordinates[i,2], output_atomcoordinates[i,3];
	}
}
