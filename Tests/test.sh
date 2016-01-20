#This will run each of the main scripts of MetaPalette and check the results with the pre-computed values
#If in a Docker container, run with: xvfb-run ./test.sh

#Test training
python ../src/Python/Train.py -i Data/FullFileNames.txt -o TestOutput -b bcalm -r `pwd` -j jellyfish -c count_in_file -t 4 -k 5 -s 4

#Write code to check if output makes sense

#Test classifying
python ../src/Python/Classify.py -d TestOutput -o TestOutput -i Data/test-reads.fa -k default -j jellyfish -q query_per_sequence -Q C -x -t 4

#Write code to check if output makes sense

#Test plotting
python ../src/Python/Plot.py -d TestOutput -o TestOutput -p TestOutput -t genus -i Data/test-reads.fa -g Escherichia_coli

#Check if plot output makes sense