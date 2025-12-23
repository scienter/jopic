EXEC = jopic
#CC = h5pcc
CC = /opt/ompi/5.0.6/bin/mpicc
OBJS = main.o parameterSetting.o findparam.o boundary.o saveFile.o fieldSolve.o loadLaser.o loadPlasma.o pml.o utility.o fieldShareX.o saveFieldHDF.o saveDensityHDF.o saveParticleHDF.o saveDumpHDF.o solveCharge.o clean.o materials.o ionization.o interpolation.o particlePush.o rearrangeParticles.o removeEdge.o movingDomain.o particleShareX.o filter.o updateCurrent.o restoreDumpHDF.o redist.o calConservation.o loadBeam.o shotLaser.o plasmaLens.o shunt_yard.o


#-------- for Beam server ----------# 
CFLAGS = -I/opt/gsl/2.8/include -I/opt/hdf5/1.14.5/include 
LDFLAGS = -L/opt/gsl/2.8/lib -L/opt/hdf5/1.14.5/lib

#-------- for PAL ----------#
#CFLAGS = -I/home/scienter/gsl2.6/include -I/usr/lib/x86_64-linux-gnu/hdf5/openmpi/include
#LDFLAGS = -L/home/scienter/gsl2.6/lib -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi

#-------- for home computer ----------#
#CFLAGS= -I/usr/include/hdf5/openmpi
#LDFLAGS= -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi




INCL = constants.h laser.h mesh.h particle.h plasma.h
LIBS = -lm -lhdf5 -lgsl -lgslcblas




$(EXEC):$(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LDFLAGS) $(LIBS) -o $(EXEC)
$(OBJS):$(INCL)
clean:
	@rm -f *.o *~ $(EXEC)
