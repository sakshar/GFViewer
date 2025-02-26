# run the following commands one-by-one to execute GFViewer on given test data sets based on your installation type

# test 1
# if installed as a conda package
gfviewer -d tests/data_test_1.xlsx -g tests/chrs_test_1-2.txt -o out_test_1 -c tests/colors_test_1.txt

# if executing from the source
gfviewer/main.py -d tests/data_test_1.xlsx -g tests/chrs_test_1-2.txt -o out_test_1 -c tests/colors_test_1.txt

# test 2
# if installed as a conda package
gfviewer -d tests/data_test_2.csv -g tests/chrs_test_1-2.txt -o out_test_2 -c tests/colors_test_2.txt -p 1 -r 1 -lpp -conc

# if executing from the source
gfviewer/main.py -d tests/data_test_2.csv -g tests/chrs_test_1-2.txt -o out_test_2 -c tests/colors_test_2.txt -p 1 -r 1 -lpp -conc

# test 3
# if installed as a conda package
gfviewer -d tests/data_test_3.tsv -g tests/chrs_test_3.fasta -o out_test_3 -cen -lpp

# if executing from the source
gfviewer/main.py -d tests/data_test_3.tsv -g tests/chrs_test_3.fasta -o out_test_3 -cen -lpp