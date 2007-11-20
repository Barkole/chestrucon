== project home ==
http://chestrucon.sourceforge.net/
https://sourceforge.net/projects/chestrucon/

== overview ==
This collection of scripts will solve the problem of converting crystal structure data files into varios other file formats. It will focus on the interchange of input-files for different quantum chemestry programms, as well as the conversion of their structure output files into SHELX INS files. INS files can be read and processed by different chemistry programms (e.g. Diamond, PLATON...).

By now the current script can convert CTRL (LMTO) files into POSCAR (VASP) files.


Now and in future this collection had to be run on difference systems. Based on linux it has to run e.g. on windows with cygwin.

The only known commercial alternative to current function is http://schmeling.ac.rwth-aachen.de/user/bernhard/wxdragon.html .

VASP is the programm using POSCAR files:
* http://cms.mpi.univie.ac.at/vasp/

LMTO is the programm using CTRL files:
* CTRL: Ist LMTO: http://www.fkf.mpg.de/andersen/LMTODOC/LMTODOC.html

The SHELX-97 homepage is:
* http://shelx.uni-ac.gwdg.de/SHELX/

Other converters are: (not converting CTRL files into POSCAR files, but handling one type)
* http://cms.mpi.univie.ac.at/odubay/p4vasp_site/download.php?list.6
* http://www.cscs.ch/~mvalle/ChemViz/doc/STM3/index.html
* http://www.fkf.mpg.de/andersen/users/kunstman/software.html

== scripts ==

=== ctrl2poscar.awk ===
./ctrl2poscar.awk CTRLFILE > POSCARFILE

Generates a POSCAR file based on CTRL fi.le

Works ONLY with CUBIC systems
Will NOT WORK with TRIGONAL and HEXAGONAL systems

Script file  contains some details.

