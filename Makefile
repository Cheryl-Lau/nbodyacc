
.PHONY: nbodyacc
nbodyacc:
	@cd build; ${MAKE} ${MAKECMDGOALS}

%::
	@cd build; ${MAKE} "${MAKECMDGOALS}"

clean:
	@cd build; ${MAKE} clean
