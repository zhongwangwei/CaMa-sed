#==========================================================
#* Makefile for stand-alone CaMa-Flood 
#
# (C) D.Yamazaki & E. Dutra  (U-Tokyo/FCUL)  Aug 2019
#
# Licensed under the Apache License, Version 2.0 (the "License");
#   You may not use this file except in compliance with the License.
#   You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is 
#  distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
# See the License for the specific language governing permissions and limitations under the License.
#==========================================================
include   ./Mkinclude

#=============================
#*** module
OBJECTS=\
PARKIND1.o \
CMF_VARS_MOD.o \
CMF_NMLIST_MOD.o   \
CMF_UTILS_MOD.o \
CMF_SEDPAR_MOD.o \
CMF_FLDSTG_MOD.o \
CMF_OUTFLW_MOD.o \
CMF_OUTPUT_MOD.o \
CMF_INIT_MOD.o   \
CMF_MPI_MOD.o \
CMF_FORCING_MOD.o \
CMF_TIME_MOD.o \
CMF_PTHOUT_MOD.o \
CMF_STONXT_MOD.o \
CMF_DIAG_MOD.o \
CMF_LEVEE_MOD.o \
CMF_DAMOUT_MOD.o \
CMF_LAKEIN_MOD.o \
CMF_PHYSICS_MOD.o \
CMF_BOUNDARY_MOD.o \
CMF_TRACER_MOD.o \
CMF_SEDFLW_MOD.o \
CMF_DRIVER_MOD.o \
CMF_END_MOD.o \
#=============================
#*** suffix rule
.SUFFIXES : .o .F90
.F90.o:
	$(FCMP) ${FFLAGS} -c ${INC} $(MODS) ${CFLAGS} -c $<
src: $(OBJECTS)
	ar -rv libcama.a $(OBJECTS)

#=============================
#*** for main program
TARGET = MAIN_cmf
$(TARGET): $(OBJECTS) src $(TARGET).o
	$(FCMP) $(FFLAGS) $(LFLAGS) $@.o $(OBJECTS) -o $@ $(LIB)

# bug in MacOSX libcama.a does not work
#	$(FCMP) $(FFLAGS) $(LFLAGS) $@.o libcama.a -o $@ $(LIB)

#=============================
#*** general rule
all: $(TARGET)
clean:
	${RM} -rf *.o core *~ *trace temp* *.mod *.a *.s *.dSYM
cleanall:
	${RM} -rf *.o core *~ *trace temp* *.mod *.a *.s *.dSYM $(TARGET)

