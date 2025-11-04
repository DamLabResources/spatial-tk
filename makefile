.PHONY: help venv install install-dev build test test-unit test-functional test-coverage clean clean-all lint format create-test-data run

# Default target
help:
	@echo "Xenium Process Makefile"
	@echo ""
	@echo "Available targets:"
	@echo "  make venv              - Create virtual environment with conda"
	@echo "  make install           - Install package in current environment"
	@echo "  make install-dev       - Install package with dev dependencies"
	@echo "  make build             - Build distribution packages"
	@echo "  make test              - Run all tests (uses .pytest_tmp/ for temp files)"
	@echo "  make test-unit         - Run only unit tests (fast)"
	@echo "  make test-functional   - Run only functional tests (slower)"
	@echo "  make test-coverage     - Run tests with coverage report"
	@echo "  make create-test-data  - Generate subsampled test data"
	@echo "  make clean-test        - Clean up test temporary files"
	@echo "  make lint              - Run linting checks"
	@echo "  make format            - Format code with black"
	@echo "  make clean             - Remove build artifacts and caches"
	@echo "  make clean-all         - Remove build artifacts, caches, and venv"
	@echo "  make run ROOT=/path    - Run full pipeline using config.toml in ROOT directory"

# Create virtual environment
venv:
	conda create -p ./venv python=3.12 -y
	./venv/bin/pip install --upgrade pip
	./venv/bin/pip install -e ".[dev]"

# Install package
install:
	pip install -e .

# Install with development dependencies
install-dev:
	pip install -e ".[dev]"

# Build distribution packages
build:
	pip install --upgrade build
	python -m build
	@echo "Distribution packages created in dist/"

# Run all tests (with custom temp directory on larger partition)
test:
	pytest -v --basetemp=.pytest_tmp

# Run only unit tests
test-unit:
	pytest tests/unit/ -v --basetemp=.pytest_tmp

# Run only functional tests
test-functional:
	pytest tests/functional/ -v --basetemp=.pytest_tmp

# Run tests with coverage
test-coverage:
	pytest --cov=xenium_process --cov-report=html --cov-report=term --basetemp=.pytest_tmp
	@echo "Coverage report generated in htmlcov/"

# Clean up test temporary files
clean-test:
	rm -rf .pytest_tmp/
	rm -rf .pytest_cache/
	@echo "Cleaned test temporary files"

# Create test data
create-test-data:
	python scripts/create_test_data.py \
		--input-csv example.csv \
		--output-dir tests/test_data \
		--n-cells 500

# Run linting
lint:
	@echo "Running flake8..."
	-flake8 xenium_process/ tests/ --count --select=E9,F63,F7,F82 --show-source --statistics
	@echo "Running mypy..."
	-mypy xenium_process/ --ignore-missing-imports

# Format code
format:
	@echo "Formatting with black..."
	black xenium_process/ tests/ scripts/

# Clean build artifacts and caches
clean:
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	rm -rf .pytest_cache/
	rm -rf .pytest_tmp/
	rm -rf .coverage
	rm -rf htmlcov/
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete
	find . -type f -name "*.pyo" -delete
	@echo "Cleaned build artifacts and caches"

# Clean everything including venv
clean-all: clean
	rm -rf venv/
	@echo "Cleaned everything including venv"

# Development workflow shortcuts
dev-setup: venv install-dev create-test-data
	@echo "Development environment ready!"

# Quick test during development
quick-test: test-unit
	@echo "Quick unit tests passed!"

# Run full pipeline using config file
run:
	@if [ -z "$(ROOT)" ]; then \
		echo "Error: ROOT must be specified. Usage: make run ROOT=/path/to/directory"; \
		exit 1; \
	fi
	@if [ ! -d "$(ROOT)" ]; then \
		echo "Error: Directory $(ROOT) does not exist"; \
		exit 1; \
	fi
	@if [ ! -f "$(ROOT)/config.toml" ]; then \
		echo "Error: config.toml not found in $(ROOT)"; \
		exit 1; \
	fi
	@echo "Running pipeline with config: config.toml"
	@echo "Working directory: $(ROOT)"
	@echo "=========================================="
	@cd "$(ROOT)" && \
	echo "Step 1: Concatenate samples" && \
	xenium_process concat --config "config.toml" || exit 1
	@cd "$(ROOT)" && \
	echo "" && \
	echo "Step 2: Normalize data" && \
	xenium_process normalize --config "config.toml" || exit 1
	@cd "$(ROOT)" && \
	echo "" && \
	echo "Step 3: Cluster cells" && \
	xenium_process cluster --config "config.toml" || exit 1
	@cd "$(ROOT)" && \
	echo "" && \
	echo "Step 4: Annotate cell types" && \
	xenium_process annotate --config "config.toml" || exit 1
	@cd "$(ROOT)" && \
	echo "" && \
	echo "Step 5: Differential expression analysis" && \
	xenium_process differential --config "config.toml" || exit 1
	@echo ""
	@echo "=========================================="
	@echo "Pipeline completed successfully!"
