# Makefile for GIGI Python bindings

.PHONY: all build install clean test help

# Default target
all: build

# Build the extension in-place
build:
	python setup.py build_ext --inplace

# Install the package
install:
	pip install -e .

# Clean build artifacts
clean:
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	rm -f pygigi*.so
	rm -f pygigi*.pyd
	find . -name "*.pyc" -delete
	find . -name "__pycache__" -type d -exec rm -rf {} +

# Run tests
test: build
	python test_example.py

# Install dependencies
deps:
	pip install -r requirements.txt

# Build using CMake (alternative method)
cmake-build:
	mkdir -p build
	cd build && cmake .. && make -j4

# Help
help:
	@echo "Available targets:"
	@echo "  build        - Build the extension in-place using setup.py"
	@echo "  install      - Install the package in development mode"
	@echo "  clean        - Clean build artifacts"
	@echo "  test         - Run the test example"
	@echo "  deps         - Install Python dependencies"
	@echo "  cmake-build  - Build using CMake (alternative)"
	@echo "  help         - Show this help message"
