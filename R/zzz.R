# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
".First.lib" <- function(lib, pkg) {
    library.dynam("LatticeKrig", pkg, lib)
     packageStartupMessage(" help(LKrig) gives an overview \n
Copyright 2011, Licensed under GPL, www.gpl.org/licenses/gpl.html ")
}
