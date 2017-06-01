#
#	Makefile for MLSS'13 practical on network modeling and information propagation
#
#	implements:
#		all (default), clean
#

include ../Makefile.config

all: generate_network generate_random_network compute_properties simulate_ic maximize_sources_ic maximize_edges_ic convert_to_gephi

# simulate_ctic maximize_sources_ctic maximize_edges_ctic

generate_network: generate_network.cpp mlss.cpp mlss.h $(CSNAP)/Snap.o
	$(CC) $(CXXFLAGS) -o generate_network generate_network.cpp mlss.cpp ../snap-adv/kronecker.cpp $(CSNAP)/Snap.o -I$(CSNAP) -I$(CSNAPADV) -I$(CGLIB) -I$(CSNAPEXP) $(LDFLAGS) $(LIBS)
	
generate_random_network: generate_random_network.cpp mlss.cpp mlss.h $(CSNAP)/Snap.o
	$(CC) $(CXXFLAGS) -o generate_random_network generate_random_network.cpp mlss.cpp ../snap-adv/kronecker.cpp $(CSNAP)/Snap.o -I$(CSNAP) -I$(CSNAPADV) -I$(CGLIB) -I$(CSNAPEXP) $(LDFLAGS) $(LIBS)

compute_properties: compute_properties.cpp mlss.cpp mlss.h $(CSNAP)/Snap.o
	$(CC) $(CXXFLAGS) -o compute_properties compute_properties.cpp mlss.cpp ../snap-adv/kronecker.cpp $(CSNAP)/Snap.o -I$(CSNAP) -I$(CSNAPADV) -I$(CGLIB) -I$(CSNAPEXP) $(LDFLAGS) $(LIBS)        

simulate_ic: simulate_ic.cpp mlss.cpp mlss.h $(CSNAP)/Snap.o
	$(CC) $(CXXFLAGS) -o simulate_ic simulate_ic.cpp mlss.cpp ../snap-adv/kronecker.cpp $(CSNAP)/Snap.o -I$(CSNAP) -I$(CSNAPADV) -I$(CGLIB) -I$(CSNAPEXP) $(LDFLAGS) $(LIBS)

maximize_sources_ic: maximize_sources_ic.cpp mlss.cpp mlss.h $(CSNAP)/Snap.o
	$(CC) $(CXXFLAGS) -o maximize_sources_ic maximize_sources_ic.cpp mlss.cpp ../snap-adv/kronecker.cpp $(CSNAP)/Snap.o -I$(CSNAP) -I$(CSNAPADV) -I$(CGLIB) -I$(CSNAPEXP) $(LDFLAGS) $(LIBS)

simulate_ctic: simulate_ctic.cpp mlss.cpp mlss.h $(CSNAP)/Snap.o
	$(CC) $(CXXFLAGS) -o simulate_ctic simulate_ctic.cpp mlss.cpp ../snap-adv/kronecker.cpp $(CSNAP)/Snap.o -I$(CSNAP) -I$(CSNAPADV) -I$(CGLIB) -I$(CSNAPEXP) $(LDFLAGS) $(LIBS)
	
maximize_sources_ctic: maximize_sources_ctic.cpp mlss.cpp mlss.h $(CSNAP)/Snap.o
	$(CC) $(CXXFLAGS) -o maximize_sources_ctic maximize_sources_ctic.cpp mlss.cpp ../snap-adv/kronecker.cpp $(CSNAP)/Snap.o -I$(CSNAP) -I$(CSNAPADV) -I$(CGLIB) -I$(CSNAPEXP) $(LDFLAGS) $(LIBS)
	
maximize_edges_ic: maximize_edges_ic.cpp mlss.cpp mlss.h $(CSNAP)/Snap.o
	$(CC) $(CXXFLAGS) -o maximize_edges_ic maximize_edges_ic.cpp mlss.cpp ../snap-adv/kronecker.cpp $(CSNAP)/Snap.o -I$(CSNAP) -I$(CSNAPADV) -I$(CGLIB) -I$(CSNAPEXP) $(LDFLAGS) $(LIBS)
	
maximize_edges_ctic: maximize_edges_ctic.cpp mlss.cpp mlss.h $(CSNAP)/Snap.o
	$(CC) $(CXXFLAGS) -o maximize_edges_ctic maximize_edges_ctic.cpp mlss.cpp ../snap-adv/kronecker.cpp $(CSNAP)/Snap.o -I$(CSNAP) -I$(CSNAPADV) -I$(CGLIB) -I$(CSNAPEXP) $(LDFLAGS) $(LIBS)
	
convert_to_gephi: convert_to_gephi.cpp mlss.cpp mlss.h $(CSNAP)/Snap.o
	$(CC) $(CXXFLAGS) -o convert_to_gephi convert_to_gephi.cpp mlss.cpp ../snap-adv/kronecker.cpp $(CSNAP)/Snap.o -I$(CSNAP) -I$(CSNAPADV) -I$(CGLIB) -I$(CSNAPEXP) $(LDFLAGS) $(LIBS)

$(CSNAP)/Snap.o:
	make -C $(CSNAP)
	
clean:
	rm -f *.o generate_network generate_random_network compute_properties simulate_ic maximize_sources_ic simulate_ctic maximize_sources_ctic maximize_edges_ic maximize_edges_ctic convert_to_gephi
	rm -rf Debug Release