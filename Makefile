.PHONY: all checks figures
all:
	python reproduce.py --all
checks:
	python reproduce.py --checks
figures:
	python reproduce.py --figures
