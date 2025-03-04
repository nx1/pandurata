# New 01/03/2025
`debug.html` : notes from debugging (open in browser).

`makefile` : added a compile with debug mode (set DEBUG=1) in the makefile, this allows it to be run with `gdb ./pandurata` (see below)

`pandurata.c` :
    1. Renormalized photon momentum moved to renormalize_momentum()
    2. dt = yn[1]/20.;    --> dt = yn[1]/200.;
    3. debug prints + comments 

`/data/read_pandurata.py` : Contains function for reading pandurata output files.

`/notebooks/'Pandurata_ULX_old.ipynb` : Contains plotting functions for the output files.

`/notebooks/'Pandurata_ULX_new.ipynb` : Contains plotting functions for the output files but is on the new output data.

`scripts/plot_cyclotron_supression.py` : was unsure what the cyclotron supression function was doing.

`readme.MD` : See below 


This is the version of Pandurata that was used in Schnittman, Krolik, and Noble (2013) [ApJ 769, 156]
It uses the IDL routine read_harm3d_data2 to process the output data from Harm3d found in /rawdata, re-writing it in a Pandurata-friendly form in /data
in /data, we need to convert the ascii data to binary via 
> gcc -o write_harm3d_bin write_harm3d_bin.c -lm
> ./write_harm3d_bin

compile Pandurata in the main directory via
> make -f makepan
> ./pandurata

clean files using
> make -f makepan clean

To run in debug mode using gdb make using
> make -f makepan DEBUG=1
> gdb ./pandurata

breakpoints can be set using break main 

This will generate a bunch of data files in /data

# Python
pip install -r requirements.txt

# dat files

Inputs (converted to .bin):
bb_0624_new.dat        
gr_0624_new.dat        
ph_0624_new.dat        
rh_0624_new.dat        
ta_0624_new.dat        
te_0624_new.dat        
u0_0624_new.dat        
u1_0624_new.dat        
u2_0624_new.dat        
u3_0624_new.dat

Outputs:
dump_times.dat         
scat_cpow.0624.dat     
scat_disk.0624.dat     
scat_imag.0624.dat     
scat_ipol.0624.dat     
scat_ithp.0624.dat     
scat_line.0624.dat     
scat_spcp.0624.dat     
scat_spcr.0624.dat     
scat_spec.0624.dat     

