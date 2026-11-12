SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = source
BUILDDIR      = build

.PHONY: help clean html html-zh

help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS)

clean:
	rm -rf "$(BUILDDIR)"

html:
	@$(SPHINXBUILD) -b html -W --keep-going "$(SOURCEDIR)" "$(BUILDDIR)/en" $(SPHINXOPTS)

html-zh:
	@$(SPHINXBUILD) -b html -W --keep-going zh_CN "$(BUILDDIR)/zh_CN" $(SPHINXOPTS)
