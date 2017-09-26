// (C) 2016 Marcelo Soares Souza <marcelo.soares@colaborador.embrapa.br>
// This program is licensed under a LGPLv3 License.

#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <iterator>
#include <tr1/unordered_set>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <magic.h>

#include "fast-cpp-csv-parser/csv.h"
#include "cleanup_sequences_using_sequence.h"

const int THRESHOLD = 16777216;

int main(int argc,char **argv) {

  if (argc < 4)
  {
    cout << "Usage: " << argv[0] << " [CSV FILE] [FASTQ/GZ] [TYPE FASTA/FASTQ]" << "\n";

    return EXIT_FAILURE;
  }

  int count_clean = 0, count_filter = 0;
  string sequence = "";
  string csv_file = argv[1];
  string input_file = argv[2];
  string type_file = argv[3];
  string out_clean_file = input_file  + ".cleaned.result";
  string out_filter_file = input_file  + ".filtered.result";
  FastQ f;

  cout << "\nUsing " << input_file << " in " << type_file << " and " << csv_file << "\n\nProcessing..." << "\n";

  unordered_set<string> to_remove = loadCSV(csv_file);

  struct magic_set *magic = magic_open(MAGIC_MIME|MAGIC_CHECK);
  magic_load(magic, NULL);

  ifstream in_fastq(input_file, ios_base::in);
  boost::iostreams::filtering_streambuf<boost::iostreams::input> in_buf;

  if (std::strcmp("application/gzip; charset=binary", magic_file(magic, argv[2])) == 0)
    in_buf.push(boost::iostreams::gzip_decompressor());

  in_buf.push(in_fastq);
  istream in_data(&in_buf);

  ofstream out_clean(out_clean_file);
  ofstream out_filter(out_filter_file);


  std::string buffer_filter;
  buffer_filter.reserve(THRESHOLD);

  std::string buffer_clean;
  buffer_clean.reserve(THRESHOLD);

  while(!in_data.eof())
  {
    if (!getline(in_data, f.name,'\n')) break;
    if (!getline(in_data, f.sequence,'\n')) break;

    if (type_file == "fastq") {
      if (!getline(in_data, f.info,'\n')) break;
      if (!getline(in_data, f.quality,'\n')) break;
    }

    f.remove = false;

    for (auto remove : to_remove) {
      if (f.sequence.find(remove) != std::string::npos)
      {
        cout << "Found: " << remove <<  "(" << f.sequence << ")" << " in " << f.name << '\n';
        f.remove = true;
      }
    }

    if (f.remove == false) {
      if (buffer_clean.length() + 1 >= THRESHOLD) {
        out_clean << buffer_clean;
        buffer_clean.resize(0);
      }

      buffer_clean.append(f.name + "\n");
      buffer_clean.append(f.sequence + "\n");

      if (type_file == "fastq") {
        buffer_clean.append(f.info + "\n");
        buffer_clean.append(f.quality + "\n");
      }

      count_clean++;

    }
    else
    {
      if (buffer_filter.length() + 1 >= THRESHOLD) {
        out_filter << buffer_filter;
        buffer_filter.resize(0);
      }

      buffer_filter.append(f.name + "\n");
      buffer_filter.append(f.sequence + "\n");

      if (type_file == "fastq") {
        buffer_filter.append(f.info + "\n");
        buffer_filter.append(f.quality + "\n");
      }

      count_filter++;

    }
  }

  out_clean << buffer_clean;
  out_filter << buffer_filter;

  cout << "\nCheck the results in " << out_clean_file << " (" << count_clean <<  ")\n";
  cout << "\nFiltered results in " << out_filter_file << " (" << count_filter <<  ")\n";

  out_clean.close();
  out_filter.close();

  in_fastq.close();

  return EXIT_SUCCESS;
}

unordered_set<string> loadCSV(string csv_file) {
  string seq;
  unordered_set<string> to_remove;

  CSVReader<1, trim_chars<' '>, no_quote_escape<'\t'>, single_line_comment<'#'>, empty_line_comment > csv (csv_file);

  cout << "Reading CSV: " << csv_file << endl;

  while(csv.read_row(seq))
  {
    to_remove.insert(seq);
    cout << "\nMarked to Remove: " << seq << "\n";
  }

  cout << "Total Size of Sequence List: " <<  to_remove.size() << "\n";

  return to_remove;
}
