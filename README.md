# FXFEL Framework

As part of a framework, 
The FXFEL project contains the definition of a standard, intermediate format for 
accelerator (specifically, FEL) simulation codes, and scripts to  convert data 
output from various accelerator codes to and from this format. This enables an
s2e framework which is extensible to other codes without too much effort.

The standard format is in hdf5, and contains vizSchema metadata for compatibility 
with Visit. See the manual for a full spec of the format.

Currently supported codes are
  - Astra
  - Elegant
  - VSim
  - Puffin
  - Genesis
  
We gratefully acknowledge the support of STFC's ASTeC department for HPC access, using the STFC Hartree Centre,
and the John von Neumann Institute for Computing (NIC) on JUROPA at Julich Supercomputing Centre (JSC), under project HHH20

