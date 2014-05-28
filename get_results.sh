#!/bin/bash

./gen_result_files_$1k.sh

python find_variants.py test10
python find_variants.py test20
python find_variants.py test30
python find_variants.py test50
python find_variants.py test100