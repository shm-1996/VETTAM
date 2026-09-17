#	Makefile for the dust module

Opacity += SemenovOpacities.o \
		  Opacity_data.o

Opacity_init.o : Opacity_data.o \
				SemenovOpacities.o

Opacity.o : SemenovOpacities.o \
			Opacity_data.o
