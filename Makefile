src_dir = src

python_compiler = python3
requirements_file = requirements.txt

proj_init:
	mkdir $(src_dir)
	mkdir docs
	cat ".vscode" > .gitignore

py_version:
	$(python_compiler) --version > pyversion
	
pack_export:
	pip freeze > $(requirements_file)
	
pack_install: $(requirements_file)
	pip install -r $(requirements_file)
	
docs_init: docs
	cd docs && sphinx-quickstart
	
docs_gen: docs
	echo "Please add `module` into index.rst"
	sphinx-apidoc -f -o ./docs $(src_dir)
	cd docs && make html

docs_update: docs
	cd docs && make html

docs_clean: docs
	cd docs && make clean

docs_serve:
	$(python_compiler) -m http.server --directory 'docs/_build/html'

demo_run:
	cd example && streamlit run ../molfragapp/MolFragApp.py