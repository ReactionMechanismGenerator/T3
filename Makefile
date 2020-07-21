################################################################################
#
#   Makefile for T3
#
################################################################################

install-arc:
	bash devtools/install_arc.sh

install-rms:
	bash devtools/install_rms.sh

test:
	pytest -ra -vv

test-main:
	pytest tests/test_main.py -ra -vv

test-functional:
	pytest tests/test_functional.py -ra -vv
