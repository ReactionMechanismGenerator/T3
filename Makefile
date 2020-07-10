################################################################################
#
#   Makefile for T3
#
################################################################################

install-rms:
	source devtools/install_pyrms.sh

test:
	pytest -ra -vv

test-main:
	pytest tests/test_main.py -ra -vv

test-functional:
	pytest tests/test_functional.py -ra -vv
