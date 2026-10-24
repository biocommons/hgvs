# Makefile for Python project

.DELETE_ON_ERROR:
.PHONY: FORCE
.PRECIOUS:
.SUFFIXES:

# TESTING sources
# UTA_DB_URL must be accessible at test time (for now)
# export HGVS_CACHE_MODE=
_UTAPW=anonymous
export UTA_DB_URL=postgresql://anonymous:${_UTAPW}@uta.biocommons.org:5432/uta/uta_20241220


.DEFAULT_GOAL := help
default: help

define INFO_MESSAGE
	@echo "⏩\033[1;38;5;208m$(1)\033[0m"
endef

############################################################################
#= BASIC USAGE

.PHONY: help
help: ## Display help message
	@./sbin/makefile-extract-documentation ${MAKEFILE_LIST}

############################################################################
#= SETUP, INSTALLATION, PACKAGING

install: devready
.PHONY: devready
devready: ## Prepare local dev env: Create virtual env, install the git hooks (prek)
	$(call INFO_MESSAGE, Prepare local dev env: Create virtual env and install the git hooks (prek))
	uv sync --dev
	uv run prek install -f
	$(call INFO_MESSAGE, Activate the virtual env with \`source .venv/bin/activate\`)

.PHONY: build
build: ## Build package
	$(call INFO_MESSAGE, "Build package")
	rm -fr dist
	uv build

.PHONY: publish
publish: build ## Publish package to PyPI
	$(call INFO_MESSAGE, "Publish package to PyPI")
	uv publish  # Requires UV_PUBLISH_TOKEN or Trusted Publishing setup

############################################################################
#= FORMATTING, TESTING, AND CODE QUALITY

.PHONY: cqa
cqa: ## Run code quality assessments
	$(call INFO_MESSAGE, "Checking lock file consistency")
	uv lock --locked
	$(call INFO_MESSAGE, "Linting and reformatting files")
	uv run prek run
	$(call INFO_MESSAGE, "Checking for obsolete dependencies")
	uv run deptry src

.PHONY: test
test: ## Test the code with pytest
	@echo "🚀 Testing code: Running pytest"
	uv run pytest

.PHONY: test-coverage
test-coverage: ## Test the code with pytest and generate coverage report
	@echo "🚀 Testing code: Running pytest with coverage"
	uv run pytest --cov --cov-report=term-missing --cov-report=xml

.PHONY: test-matrix
test-matrix: ## Test with each supported Python version (via nox)
	@echo "🚀 Testing code: Testing supported versions with nox"
	uv run nox

test-learn: ## Add new data to test data cache
	@echo "📈 Add new data to test data cache"
	HGVS_CACHE_MODE=learn VCR_RECORD_MODE=new_episodes uv run pytest -s
test-relearn: ## Destroy and rebuild test data cache
	rm -fr tests/data/cache-py3.hdp tests/cassettes
	HGVS_CACHE_MODE=learn VCR_RECORD_MODE=new_episodes uv run pytest -s
test-relearn-iteratively: ## Destroy and rebuild test data cache (biocommons/hgvs#760
	rm -fr tests/data/cache-py3.hdp tests/cassettes
	find tests/ -name 'test*.py' | HGVS_CACHE_MODE=learn VCR_RECORD_MODE=new_episodes xargs -tn1 -- uv run pytest --no-cov -x -s

############################################################################
#= DOCUMENTATION

.PHONY: docs-serve
docs-serve: ## Build and serve the documentation
	$(call INFO_MESSAGE, "Build and serve docs for local development")
	uv run sphinx-build -b html docs docs/_build/html
	uv run python -m http.server -d docs/_build/html

.PHONY: docs-test
docs-test: ## Test if documentation can be built without errors
	$(call INFO_MESSAGE, "Testing whether docs can be build")
	uv run sphinx-build -b html docs docs/_build/html

.PHONY: docs-deploy
docs-deploy: ## Deploy documentation (no-op: Read the Docs builds hgvs docs)
	$(call INFO_MESSAGE, "Docs are built and published by Read the Docs: https://hgvs.readthedocs.io/")

############################################################################
#= UTILITIES

.PHONY: show-outdated-deps
show-outdated-deps: ## List outdated dependencies
	$(call INFO_MESSAGE, "Listing outdated dependencies")
	uv pip list --outdated

.PHONY: upgrade-deps
upgrade-deps: ## Upgrade installed dependencies
	$(call INFO_MESSAGE, "Upgrading dependencies")
	uv lock --upgrade

.PHONY: bump-deps
bump-deps: ## Bumping versions of installed dependencies
	$(call INFO_MESSAGE, "Bumping dependency versions")
	uv-bump

.PHONY: ruff-counts
ruff-counts: ## Count ruff check violations by type (ignoring per-file-ignores)
	$(call INFO_MESSAGE, "Linting code with ruff")
	uv run ruff check --no-fix --exit-zero --config pyproject.toml --config "lint.per-file-ignores={}" --output-format=json . \
		| jq --arg pwd "$$PWD/" 'map(.filename |= (ltrimstr($$pwd) | "./" + .))' \
		| jq -r 'group_by(.filename)[] | "\(.[0].filename): " + (group_by(.code) | map("\(.[0].code)×\(length)") | join(", "))'


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
cleanest: cleaner  ## Remove files and directories that can be rebuilt
	$(call INFO_MESSAGE, "Remove files and directories that can be rebuilt")
	rm -frv .eggs .nox .venv venv

.PHONY: distclean
distclean: cleanest  ## Remove untracked files and other detritus
	@echo "❌ Remove untracked files and other detritus -- Too dangerous... do this yourself"
	# git clean -df


## <LICENSE>
## Copyright 2026 Committers
##
## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at
##
##     http://www.apache.org/licenses/LICENSE-2.0
##
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.
## </LICENSE>
