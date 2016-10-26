// (C) 2016 Marcelo Soares Souza <marcelo.soares@colaborador.embrapa.br>
// This program is licensed under a LGPLv3 License.

#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <iterator>
#include <tr1/unordered_set>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "fast-cpp-csv-parser/csv.h"
#include "cleanup_sequences_using_sequence.h"

// const int THRESHOLD = 131072;
const int THRESHOLD = 16777216;

int main(int argc,char **argv) {

  if (argc < 3)
  {
    cout << "Usage: " << argv[0] << " [CSV FILE] [FASTQ GZ]" << "\n";

    return EXIT_FAILURE;
  }

  int count_clean = 0, count_filter = 0;
  string sequence = "";
  string csv_file = argv[1];
  string fastq_file = argv[2];
  string out_clean_file = fastq_file  + ".cleaned.result";
  string out_filter_file = fastq_file  + ".filtered.result";
  FastQ f;
  unordered_set<string> to_remove = loadCSV(csv_file);

  ifstream file(fastq_file, ios_base::in | ios_base::binary);

  if(!file.is_open())
  {
    return EXIT_FAILURE;
  }

  boost::iostreams::filtering_streambuf<boost::iostreams::input> in_buf;
  in_buf.push(boost::iostreams::gzip_decompressor());
  in_buf.push(file);
  istream in_fastq(&in_buf);

  ofstream out_clean(out_clean_file);
  ofstream out_filter(out_filter_file);

  cout << "\nUsing " << fastq_file << " and " << csv_file << "\n\nProcessing..." << "\n";

  std::string buffer_filter;
  buffer_filter.reserve(THRESHOLD);

  std::string buffer_clean;
  buffer_clean.reserve(THRESHOLD);

  while(!in_fastq.eof())
  {
    if (!getline(in_fastq, f.name,'\n')) break;
    if (!getline(in_fastq, f.sequence,'\n')) break;
    if (!getline(in_fastq, f.info,'\n')) break;
    if (!getline(in_fastq, f.quality,'\n')) break;

    f.remove = false;

    for (auto remove : to_remove) {
      if (f.sequence.find(remove) != std::string::npos)
      {
        // cout << "Found: " << remove << " in " << f.name << '\n';
        f.remove = true;
      }
    }

    if (f.remove == false) {
      if (buffer_clean.length() + 1 >= THRESHOLD) {
        out_clean << buffer_clean;
        buffer_clean.resize(0);
      }

      buffer_clean.append("@" + f.name + "\n");
      buffer_clean.append(f.sequence + "\n");
      buffer_clean.append(f.info + "\n");
      buffer_clean.append(f.quality + "\n");

      count_clean++;

    }
    else
    {
      if (buffer_filter.length() + 1 >= THRESHOLD) {
        out_filter << buffer_filter;
        buffer_filter.resize(0);
      }

      buffer_filter.append("@" + f.name + "\n");
      buffer_filter.append(f.sequence + "\n");
      buffer_filter.append(f.info + "\n");
      buffer_filter.append(f.quality + "\n");

      count_filter++;

    }
  }

  out_clean << buffer_clean;
  out_filter << buffer_filter;

  cout << "\nCheck the results in " << out_clean_file << " (" << count_clean <<  ")\n";
  cout << "\nFiltered results in " << out_filter_file << " (" << count_filter <<  ")\n";

  out_clean.close();
  out_filter.close();
  file.close();

  return EXIT_SUCCESS;
}

unordered_set<string> loadCSV(string csv_file) {
  string seq;
  unordered_set<string> to_remove;

  CSVReader<1, trim_chars<' '>, no_quote_escape<'\t'>, single_line_comment<'#'> > csv (csv_file);

  while(csv.read_row(seq))
  {
    to_remove.insert(seq);
    cout << "\nMarked to Remove: " << seq << "\n";
  }

  return to_remove;
}
