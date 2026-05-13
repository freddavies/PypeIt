# Makefile for Sphinx documentation
#

# You can set these variables from the command line.
BUILDDIR      = _build
LOCALFILES    = $(BUILDDIR)/* api/*.rst help/*.rst sphinx_warnings.out pypeit_par.rst
SPHINXAPI     = sphinx-apidoc
SPHINXAPIOPT  = members,private-members,undoc-members,show-inheritance
SPHINXBUILD   = sphinx-build
SPHINXOPTS    = -b html -d $(BUILDDIR)/doctrees -aE -w ./sphinx_warnings.out

# User-friendly check for sphinx-build
ifeq ($(shell which $(SPHINXBUILD) >/dev/null 2>&1; echo $$?), 1)
$(error The '$(SPHINXBUILD)' command was not found. Make sure you have Sphinx installed, then set the SPHINXBUILD environment variable to point to the full path of the '$(SPHINXBUILD)' executable. Alternatively you can add the directory with the executable to your PATH. If you don't have Sphinx installed, grab it from http://sphinx-doc.org/)
endif

clean:
	rm -rf $(LOCALFILES)

apirst:
	SPHINX_APIDOC_OPTIONS=$(SPHINXAPIOPT) $(SPHINXAPI) --separate -o ./api ../pypeit ../pypeit/tests/* ../pypeit/version.py ../pypeit/compiler_version*
	wget -O ./include/dev_suite_readme.rst https://raw.githubusercontent.com/pypeit/PypeIt-development-suite/main/README.rst
	python ./scripts/build_afterburn_datamodels_rst.py
	python ./scripts/build_bitmask_rst.py
	python ./scripts/build_cache_data_tbl.py
	python ./scripts/build_calib_images_rst.py
	python ./scripts/build_datacontainer_datamodels.py
	python ./scripts/build_dependency_rst.py
	python ./scripts/build_detector_table.py
	python ./scripts/build_imgproc_defaults.py
	python ./scripts/build_par_rst.py
	python ./scripts/build_specobj_rst.py
	python ./scripts/build_spectbl_rst.py
	python ./scripts/build_telluric_table.py
	python ./scripts/write_script_help.py
	python ./scripts/build_standardstbl_rst.py
	python ./scripts/patch_apidoc_for_dataclasses.py

examples:
	python ./scripts/make_example_files.py

html:
	make apirst
	make examples
	$(SPHINXBUILD) $(SPHINXOPTS) . $(BUILDDIR)/html
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

htmlnoex:
	make apirst
	$(SPHINXBUILD) $(SPHINXOPTS) . $(BUILDDIR)/html
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

htmlonly:
	$(SPHINXBUILD) $(SPHINXOPTS) . $(BUILDDIR)/html
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

picky:
	make apirst
	make examples
	$(SPHINXBUILD) $(SPHINXOPTS) -n -vvv -T -W . $(BUILDDIR)/html
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."
