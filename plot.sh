#!/bin/bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"
git submodule update --init --recursive
cd gemmi
cmake -S. -Bbuild
cmake --build build
cd ..
matlab -nodisplay -nosplash -nodesktop -r "run('./startup.m');run('./run_tests.m');exit;"
pdflatex -output-directory=diagrams diagrams/simple_example.tex
pdflatex -output-directory=diagrams diagrams/gaussian_IMMA_test.tex
pdflatex -output-directory=diagrams diagrams/test_matmul_accuracy_berr.tex
cd diagrams
latexmk -c simple_example.tex
latexmk -c gaussian_IMMA_test.tex
latexmk -c test_matmul_accuracy_berr.tex

cd ..
