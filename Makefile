#   Makefile for T3

install:
	bash devtools/install_all.sh

install-pyrdl:
	bash devtools/install_pyrdl.sh

test:
	pytest tests/ --cov -ra -vv

test-main:
	pytest tests/test_main.py --cov -ra -vv

test-functional:
	pytest tests/test_functional.py --cov -ra -vv
