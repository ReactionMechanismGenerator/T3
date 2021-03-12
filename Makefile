################################################################################
#
#   Makefile for T3
#
################################################################################

install:
	bash devtools/install_arc.sh
	bash devtools/install_pyrms.sh

install-arc:
	bash devtools/install_arc.sh

install-pyrms:
	bash devtools/install_pyrms.sh

test:
	pytest -ra -vv

test-main:
	pytest tests/test_main.py -ra -vv

test-functional:
	pytest tests/test_functional.py -ra -vv
