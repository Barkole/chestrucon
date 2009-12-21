#!/usr/bin/awk -f

# Copyright (C) 2007 Rebecca Roemer <http://fantasia-fuoco.de>

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
# * Rebecca Roemer (rr) <http://fantasia-fuoco.de>
# * Michael Decker (mad) <http://Inspire-Mind.de>

# Version:
# * 0.0.0.4: 2008-04-17 - remove bug for angle of unit cell calculation for scaling factors other than 1.000000
# * 0.0.0.3: 2008-03-23 - generalization for varaiable number of atom sorts (mad)
# * 0.0.0.2: 2008-03-18 - include scaling factor in calculation of unit cell latice (rr)
# * 0.0.0.1: 2007-11-11 - Only working for CONTCARs containing one, two and three atom sorts and scaling factor 1.000000 for matrix (rr)

# Tested:
# * 2008-03-23: Kubuntu 7.10 with
#   - GNU Awk 3.1.5
# * 2008-03-21: Suse 10.2


BEGIN { #setting different variables
        myLineCounter = 0;
	maxLineCounter = 0;

        # cell parameters
        axis[1] = 0 ; # a-axis
        axis[2] = 0 ; # b-axis
        axis[3] = 0 ; # c-axis
        
        #norm (ger. Betrag) for a-, b, und c-vector
        betrag[1] = 0 ; # norm a-vector
        betrag[2] = 0 ; # norm b-vector
        betrag[3] = 0 ; # norm c-vector

        # values needed for calculation of cell parameters
        cosin[1] = 0 ; # cos(alpha)
        cosin[2] = 0 ; # cos(beta)
        cosin[3] = 0 ; # cos(gamma)
        
        # Scalierung der Matrixelemente
        scale[1] = 0

        # Matrix describing unit cell
        # a1 a2 a3 (line 3 of input)
        # b1 b2 b3 (line 4 of input)
        # c1 c2 c3 (line 5 of input)
        a[1] = 0 ; # a1
        a[2] = 0 ; # a2
        a[3] = 0 ; # a3 
        b[1] = 0 ; # b1
        b[2] = 0 ; # b2
        b[3] = 0 ; # b3
        c[1] = 0 ; # c1
        c[2] = 0 ; # c2 
        c[3] = 0 ; # c3

        # Unit cell angles
        alpha[1] = 0 ;
        beta[1] = 0 ;
        gamma[1] = 0 ;

        # number of atoms of different atom sorts
        atom[1] = 0 ;
        atom[2] = 0 ;
        atom[3] = 0 ;

        # total number of atoms in unit cell
        number_atoms[1] = 0;
        
        # atom sorts
        atom_sort[1] = 0;
 
	sumOfPreviousAtoms = 0;
	currentAtomCount = 0;
        }

# Functions from: http://www.gladir.com/CODER/AWK/acos.htm
function abs(a) {
   if(a<0) a=-a;
   return a;
    }

function acos(a) {
   pi=3.141592653589793 ;
   if(abs(a)==1) {
      return (1-a)*pi/2  ;  
   } else {
      return atan2(-a,sqrt(1-a*a))+2*atan2(0.5,0.5) ; 
         }
        }

{ myLineCounter++;
 }

myLineCounter == 1 {
                     for (i=1 ; i <= NF ; i++)
                        { atom_sort[i] = $i ;
                          }
                      }

myLineCounter == 3 {
                    scale[1] = $1;
                     }

myLineCounter == 4 {  
                     axis[1] = (((scale[1]*$1)^2 + (scale[1]*$2)^2 + (scale[1]*$3)^2)^(1/2)) ; 
                      } 

myLineCounter == 4 {
                     betrag[1] = ((($1)^2 + ($2)^2 + ($3)^2)^(1/2)) ;
                      }

# Choose according to your awk (** or ^ for "to the power of")
# myLineCounter == 4 {  
#                      axis[1] = (((scale[1]*$1)**2 + (scale[1]*$2)**2 + (scale[1]*$3)**2)**(1/2)) ; 
#                       }

myLineCounter == 4 {  
                      a[1] = $1 ; 
                      a[2] = $2 ;
                      a[3] = $3 ;
                        }

# Choose according to your awk (** or ^ for "to the power of")
# myLineCounter == 5 {  
#                      axis[2] = (((scale[1]*$1)**2 + (scale[1]*$2)**2 + (scale[1]*$3)**2)**(1/2)) ; 
#                       }

myLineCounter == 5 {  
                      axis[2] = (((scale[1]*$1)^2 + (scale[1]*$2)^2 + (scale[1]*$3)^2)^(1/2)) ; 
                       }

myLineCounter == 5 {
                     betrag[2] = ((($1)^2 + ($2)^2 + ($3)^2)^(1/2)) ;
                      }

myLineCounter == 5 { 
                      b[1] = $1 ;
                      b[2] = $2 ;
                      b[3] = $3 ; 
                       }

# Choose according to your awk (** or ^ for "to the power of")
# myLineCounter == 6 { 
#                       axis[3] = (((scale[1]*$1)**2 + (scale[1]*$2)**2 + (scale[1]*$3)**2)**(1/2)) ; 
#                        } 

myLineCounter == 6 { 
                      axis[3] = (((scale[1]*$1)^2 + (scale[1]*$2)^2 + (scale[1]*$3)^2)^(1/2)) ;
                       }

myLineCounter == 6 {
                     betrag[3] = ((($1)^2 + ($2)^2 + ($3)^2)^(1/2)) ;
                      }

myLineCounter == 6 { 
                      c[1] = $1 ;
                      c[2] = $2 ;
                      c[3] = $3 ; 
                       }

myLineCounter == 6 {  
                      cosin[1] = ((b[1] * c[1] + b[2] * c[2] + b[3] * c[3]) / (betrag[2] * betrag[3]));
                      cosin[2] = ((a[1] * c[1] + a[2] * c[2] + a[3] * c[3]) / (betrag[1] * betrag[3]));
                      cosin[3] = ((b[1] * a[1] + b[2] * a[2] + b[3] * a[3]) / (betrag[2] * betrag[1]));
                       }

myLineCounter == 6 { 
                      alpha[1] = ( acos(cosin[1]) * (180 / 3.141592653589793)) ;
                      beta[1] = ( acos(cosin[2]) * (180 / 3.141592653589793)) ;
                      gamma[1] = ( acos(cosin[3]) * (180 / 3.141592653589793)) ;  
                       }

myLineCounter == 7 { 
                      printf "%s %.8f %.8f %.8f %.8f %.8f %.8f\n", "CELL 0.71073 ", axis[1] , axis[2] , axis[3] , alpha[1] , beta[1] , gamma[1] ;
                       }

myLineCounter == 7 { 
                      print "LATT -1";
                      printf "%s", "SFAC ";
                       }

myLineCounter == 7 {  
                       for (i=1 ; i < NF ; i++)
                        { printf "%s %s" , atom_sort[i], "" ;
                          }
                       for (i=NF ; i == NF ; i++)
                        { printf "%s\n", atom_sort[i] ;
                          }

                      for (i=1 ; i <= NF ; i++)
                        { atom[i] = $i ;
                          }
		      max_line[0] = 8;
                      max_line[1] = ( atom[1] + 8 ) ;
                      for (i=2 ; i <= NF ; i++)
                        { max_line[i] = max_line[(i-1)] + atom[i] ;
                          } 
                        }  

(myLineCounter > max_line[maxLineCounter]) && (myLineCounter <= max_line[maxLineCounter+1]) {
			currentAtomCount++;
                      printf "%s %s %s %s %s\n", atom_sort[maxLineCounter+1] currentAtomCount, maxLineCounter+1 ,  $1,  $2,  $3 ;
                        }

# if we are on the last line of current atom, we have to move to the next atom
(myLineCounter == max_line[maxLineCounter+1]) {
				maxLineCounter++;
				currentAtomCount = 0;
			}


 END { print "END" }   
