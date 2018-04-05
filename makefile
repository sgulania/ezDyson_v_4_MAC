CC = g++-7

CCFLAGS = -fopenmp 

LXML = -lexpat
LFLAGS = $(LXML)

CPPFLAGS =  -I /usr/local/opt/libxml++-2.6/

CODE = exedys_MAC

CCCLASSES = aobasis.C cklm.C eikr.C complexno.C klmgrid.C gauss.C orbital.C xyzgrid.C ylm.C pad.C sph.C rotnmatr.C anglegrid.C
CCMETHODS = readwrite.C simple_xml_parser.C tools.C wavetypes.C
CCMAIN = main.C dyson_main.C
CCSRC = $(CCMAIN) $(CCMETHODS) $(CCCLASSES)

CCBINOBJ = $(CCSRC:%.C=%.o)
BINOBJ = $(CCBINOBJ)

$(CODE): $(BINOBJ)
	$(CC) $^ $(LFLAGS) $(CPPFLAGS) $(CCFLAGS)  -o $(CODE)

%.o: %.C
	$(CC)  $(CCFLAGS) -c $< -o $@

clean:
	$(RM) $(CODE) *.o 
