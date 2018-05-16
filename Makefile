CC = cc
#nrdir=$(HOME)/home/misc/numrec/recipes_c/code
nrdir=./recipes_c
#CFLAGS = -I. -I/apps/recipes_c/original -xO4 -dalign $(PAR)
#OPT= -fast -xO5
OPT= -O3
CFLAGS = -I. -I$(nrdir) $(OPT) $(PAR)
#CFLAGS = -I. -I$(nrdir) $(PAR) -g #-fnonstd
PAR=
#PAR=-xautopar -xloopinfo
LDFLAGS = -lm $(PAR) # -lsunmath -s #-L/apps/recipes_c/library -lrecipes

#CC = gcc
#CFLAGS = -I. -g
#LDFLAGS = -lm
LINTFLAGS = -lm -x -I. -I$(nrdir)

# nb VPATH forces searches in these directories...
VPATH = $(nrdir)

nrfil = nrutil.c qromb.c trapzd.c polint.c  chebev.c \
	spline.c splint.c gauleg.c
nrobs = nrutil.o qromb.o trapzd.o polint.o  chebev.o \
	spline.o splint.o gauleg.o
hfil = nr.h nrutil.h

nrlib.a: $(nrfil)
#	force nrobs to get made (default rules...)
	make $(nrobs)
	ar cr $@ $(nrobs)
	ranlib $@
	-rm $(nrobs)

vishniac: vishniac.c nrlib.a
	$(CC) $(CFLAGS) vishniac.c nrlib.a $(LDFLAGS) -o $@

lint: vish.c
	lint $(LINTFLAGS) vish.c

vish_gauss: vish.c
	$(CC) $(CFLAGS) vish.c -DDO_GAUSS  $(LDFLAGS) -o $@

vish_approx: vish.c
	$(CC) $(CFLAGS) vish.c -DUSE_S_APPROX  $(LDFLAGS) -o $@

vish_p: vish.c
	$(CC) $(CFLAGS) vish.c $(LDFLAGS) -o $@
