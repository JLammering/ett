ifeq (,$(shell sh -c 'cygpath --version 2> /dev/null'))
  # Unix
  pwd := $$(pwd)
  translate = $1
else
  # Windows mit MSys2/Cygwin
  pwd := $$(cygpath -m "$$(pwd)")
  translate = $(shell echo '$1' | sed 's/:/;/g')
endif

all: build/plotalpha.pdf

# hier Python-Skripte:
build/plotalpha.pdf: kopplung.py matplotlibrc header-matplotlib.tex | build
	TEXINPUTS="$(call translate,$(pwd):)" python kopplung.py

# hier weitere Abhängigkeiten für build/main.pdf deklarieren:
# build/main.pdf: build/plotalpha.pdf
#
# build/main.pdf: FORCE | build
# 	  TEXINPUTS="$(call translate,build:)" \
# 	  BIBINPUTS=build: \
# 	  max_print_line=1048576 \
# 	latexmk \
# 	  --lualatex \
# 	  --output-directory=build \
# 	  --interaction=nonstopmode \
# 	  --halt-on-error \
# 	main.tex

build:
	mkdir -p build

clean:
	rm -rf build

FORCE:

.PHONY: all clean
