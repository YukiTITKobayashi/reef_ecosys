#
# Runscript for ecosys_test:
#
# ------------------ COMPILE & RUN ----------------------------
#
#gfortran -fbounds-check  mod_geochem.F90 mod_param.F90 mod_heat.F90 mod_foodweb6.F90 mod_coral3.F90 mod_macroalgae.F90 mod_seagrass.F90 mod_sedecosys.F90 mod_reef_ecosys3.F90 mod_input.F90 ecosys_test5.F90 -O2 -o ecosys_test5.exe
 gfortran -DCHAMBER_SITE5 -fbounds-check  mod_geochem.F90 mod_param.F90 mod_foodweb.F90 mod_coral.F90 mod_macroalgae.F90 mod_seagrass.F90 mod_sedecosys.F90 mod_reef_ecosys.F90 mod_input.F90 mod_output.F90 ecosys_test5.F90 -O2 -o ecosys_s5.exe
#
# ------------------ OUTPUT FILES -----------------------------
#
 ./ecosys_s5.exe
#
