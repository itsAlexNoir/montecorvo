TARGET = spectral_scan

.PHONY: clean clean_src clean_bin

default: ${TARGET}


# Directives
clean: clean_src clean_bin

clean_bin:
	@rm -rf bin/*exec
clean_src:
	cd src; ${MAKE} clean

${TARGET}:
	cd src; ${MAKE}
