################################################################################
#
#   Makefile for T3
#
################################################################################

test test-unittests:
	nosetests --nocapture --nologcapture --all-modules --verbose --with-coverage --cover-inclusive --cover-package=t3 --cover-erase --cover-html --cover-html-dir=testing/coverage
