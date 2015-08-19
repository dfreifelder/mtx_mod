include /usr/local/conf/ElVars

all: mtx_mod

mtx_mod: mtx_mod.cpp
	${CXX} ${EL_COMPILE_FLAGS} $< -o $@ ${EL_LINK_FLAGS} ${EL_LIBS}

clean:
	rm -rf a.out
	rm -rf mtx_mod
	rm -rf *.o
	rm -rf *~
