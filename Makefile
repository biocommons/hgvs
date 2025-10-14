# Makefile for Python project

.DELETE_ON_ERROR:
.PHONY: FORCE
.PRECIOUS:
.SUFFIXES:

# TESTING sources
# UTA_DB_URL must be accessible at test time (for now)
export HGVS_CACHE_MODE=
_UTAPW=anonymous
export UTA_DB_URL=postgresql://anonymous:${_UTAPW}@uta.biocommons.org:5432/uta/uta_20241220
export HGVS_SEQREPO_URL=http://localhost:5000/seqrepo


.DEFAULT_GOAL := help
default: help

############################################################################
#= BASIC USAGE

.PHONY: help
help: ## Display help message
	@./sbin/makefile-extract-documentation ${MAKEFILE_LIST}

############################################################################
#= SETUP, INSTALLATION, PACKAGING

install: devready
.PHONY: devready
devready: ## Prepare local dev env: Create virtual env, install the pre-commit hooks
	$(call INFO_MESSAGE, Prepare local dev env: Create virtual env and install the pre-commit hooks)
	uv sync --dev
	uv run pre-commit install
	@echo '⚠️ You must activate the virtual env with `source .venv/bin/activate`'

.PHONY: build
build: ## Build package
	$(call INFO_MESSAGE, "Building package")
	rm -fr dist
	uv build

.PHONY: publish
publish: build ## publish package to PyPI
	$(call INFO_MESSAGE, "Publishing package")
	uv publish  # Requires UV_PUBLISH_TOKEN or Trusted Publishing setup

############################################################################
#= FORMATTING, TESTING, AND CODE QUALITY

.PHONY: cqa
cqa: ## Run code quality assessments
	$(call INFO_MESSAGE, "Checking lock file consistency")
	uv lock --locked
	$(call INFO_MESSAGE, "Linting and reformatting files")
	uv run pre-commit run
	$(call INFO_MESSAGE, "Checking for obsolete dependencies")
	uv run deptry src

.PHONY: test
test: ## Test the code with pytest
	@echo "🚀 Testing code: Running pytest"
	uv run pytest --cov=. --cov-report=xml

test-learn: ## add new data to test data cache
	HGVS_CACHE_MODE=learn uv run pytest -s
test-relearn: ## destroy and rebuild test data cache
	rm -fr tests/data/cache-py3.hdp tests/cassettes
	HGVS_CACHE_MODE=learn uv run pytest -s
test-relearn-iteratively: ## destroy and rebuild test data cache (biocommons/hgvs#760
	rm -fr tests/data/cache-py3.hdp tests/cassettes
	find tests/ -name 'test*.py' | HGVS_CACHE_MODE=learn xargs -tn1 -- uv run pytest --no-cov -x -s



############################################################################
#= DOCUMENTATION

.PHONY: docs-serve
docs-serve: ## Build and serve the documentation
	$(call INFO_MESSAGE, "Build and serve docs for local development")
	uv run mkdocs serve

.PHONY: docs-test
docs-test: ## Test if documentation can be built without warnings or errors
	$(call INFO_MESSAGE, "Testing whether docs can be build")
	uv run mkdocs build -s

############################################################################
#= CLEANUP

.PHONY: clean
clean:  ## Remove temporary and backup files
	$(call INFO_MESSAGE, "Remove temporary and backup files")
	find . \( -name "*~" -o -name "*.bak" \) -exec rm -frv {} +

.PHONY: cleaner
cleaner: clean  ## Remove files and directories that are easily rebuilt
	$(call INFO_MESSAGE, "Remove files and directories that are easily rebuilt")
	rm -frv .cache .DS_Store .pytest_cache .ruff_cache build coverage.xml dist docs/_build site
	find . \( -name __pycache__ -type d \) -exec rm -frv {} +
	find . \( -name "*.pyc" -o -name "*.egg-info" \) -exec rm -frv {} +
	find . \( -name "*.orig" -o -name "*.rej" \) -exec rm -frv {} +

.PHONY: cleanest
cleanest: cleaner  ## Remove all files that can be rebuilt
	$(call INFO_MESSAGE, "Remove files and directories that can be rebuilt")
	rm -frv .eggs .tox .venv venv

.PHONY: distclean
distclean: cleanest  ## Remove untracked files and other detritus
	@echo "❌ Remove untracked files and other detritus -- Too dangerous... do this yourself"
	# git clean -df
