* Clean Up sequences based on CSV file

To compile: g++ cleanup_sequences.cpp -o cleanup_sequences -lpthread -I fast-cpp-csv-parser -std=c++11 -ljsoncpp -Wall
            g++ cleanup_sequences_using_sequence.cpp -o cleanup_sequences_using_sequence -lpthread -I fast-cpp-csv-parser -std=c++11 -lboost_iostreams -lmagic -Wall
            g++ cleanup_sequences_using_sequence_mrna.cpp -o cleanup_sequences_using_sequence_mrna -lpthread -I fast-cpp-csv-parser -std=c++11 -lboost_iostreams -lmagic -Wall
