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
#include <ctime>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <magic.h>

#include "fast-cpp-csv-parser/csv.h"
#include "cleanup_sequences_using_keyword.h"

// const int THRESHOLD = 131072;
const int THRESHOLD = 16777216;

int main(int argc,char **argv) {

  if (argc < 3)
  {
    cout << "Usage: " << argv[0] << " [CSV FILE] [FASTQ/GZ]" << "\n";

    return EXIT_FAILURE;
  }

  const clock_t begin_time = clock();

  int count_clean = 0, count_filter = 0;
  string sequence = "";
  string csv_file = argv[1];
  string fastq_file = argv[2];
  string out_clean_file  = fastq_file  + "_" + csv_file + ".cleaned.result";
  string out_filter_file = fastq_file  + "_" + csv_file + ".filtered.result";
  FastQ f;

  cout << "\nLoading CSV: " << csv_file << "\n";

  unordered_set<string> to_remove = loadCSV(csv_file);

  cout << "Took: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << "\n";

  struct magic_set *magic = magic_open(MAGIC_MIME|MAGIC_CHECK);
  magic_load(magic, NULL);

  ifstream in_fastq(fastq_file, ios_base::in);
  boost::iostreams::filtering_streambuf<boost::iostreams::input> in_buf;

  if (std::strcmp("application/gzip; charset=binary", magic_file(magic, argv[2])) == 0)
    in_buf.push(boost::iostreams::gzip_decompressor());

  in_buf.push(in_fastq);
  istream in_data(&in_buf);

  ofstream out_clean(out_clean_file);
  ofstream out_filter(out_filter_file);

  cout << "\nUsing " << fastq_file << " using " << csv_file << "\n\nProcessing..." << "\n";

  std::string buffer_filter;
  buffer_filter.reserve(THRESHOLD);

  std::string buffer_clean;
  buffer_clean.reserve(THRESHOLD);

  while(!in_data.eof())
  {
    if (!getline(in_data, f.name,'\n')) break;
    if (!getline(in_data, f.sequence,'\n')) break;
    if (!getline(in_data, f.info,'\n')) break;
    if (!getline(in_data, f.quality,'\n')) break;

    f.remove = false;

    for (auto remove : to_remove) {
      if (f.name.find(remove) != std::string::npos)
      {
        f.remove = true;
      }
    }

    /*
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
    else {
    */

    if (f.remove == true) {
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

  in_fastq.close();

  cout << "Total elapsed time: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC;

  return EXIT_SUCCESS;
}

unordered_set<string> loadCSV(string csv_file) {
  string seq;
  unordered_set<string> to_remove;

  CSVReader<1, trim_chars<' '>, no_quote_escape<'\t'>, single_line_comment<'#'>, empty_line_comment > csv (csv_file);

  while(csv.read_row(seq))
  {
    to_remove.insert(seq);
  }

  return to_remove;
}
