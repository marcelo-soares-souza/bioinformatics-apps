#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <iterator>
#include <tr1/unordered_set>

#include <jsoncpp/json/json.h>
#include "csv.h"
#include "cleanup_sequences.h"

// const int THRESHOLD = 131072;
const int THRESHOLD = 16777216;

int main(int argc,char **argv) {

  if (argc < 2)
  {
    cout << "Usage: " << argv[0] << " [INFO FILE]" << "\n";

    return EXIT_FAILURE;
  }

  string header = "";
  int count = 0;
  FastQ f;

  Config c = loadConfig(argv[1]);
  unordered_set<string> to_remove = loadCSV(c);

  ifstream in_fastq(c.fastq, ios::in);
  ofstream out_filter(c.filter);
  ofstream out_clean(c.clean);

  if(!in_fastq.is_open())
  {
    return EXIT_FAILURE;
  }

  cout << "\nUsing " << c.input << " (" << c.format << ") and " << c.t6 << "\n\nProcessing..." << "\n";

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

    istringstream iss(f.name);

    iss >> header;

    unordered_set<string>::const_iterator got = to_remove.find (header);

    if (got != to_remove.end())
    {
      if (buffer_filter.length() + 1 >= THRESHOLD) {
        out_filter << buffer_filter;
        buffer_filter.resize(0);
      }

      buffer_filter.append("@" + f.name + "\n");
      buffer_filter.append(f.sequence + "\n");
      buffer_filter.append(f.info + "\n");
      buffer_filter.append(f.quality + "\n");

      /*
        out_filter << "@" << f.name << "\n";
        out_filter << f.sequence << "\n";
        out_filter << f.info << "\n";
        out_filter << f.quality << "\n";
      */

      count++;
    }
    else
    {
      if (buffer_clean.length() + 1 >= THRESHOLD) {
        out_clean << buffer_clean;
        buffer_clean.resize(0);
      }

      buffer_clean.append("@" + f.name + "\n");
      buffer_clean.append(f.sequence + "\n");
      buffer_clean.append(f.info + "\n");
      buffer_clean.append(f.quality + "\n");

      /*
        out_clean << "@" << f.name << "\n";
        out_clean << f.sequence << "\n";
        out_clean << f.info << "\n";
        out_clean << f.quality << "\n";
      */
    }
  }

  out_filter << buffer_filter;
  out_clean << buffer_clean;

  cout << "\nCheck the results in " << c.clean << "\n";
  cout << "Removed Sequences (" << count << ") in " << c.filter << "\n\n";

  in_fastq.close();
  out_filter.close();
  out_clean.close();

  return EXIT_SUCCESS;
}

Config loadConfig(string filename) {
  Config c;
  Json::Value root;

  ifstream config_doc(filename, ifstream::binary);
  config_doc >> root;
  config_doc.close();

  c.fields     = root["fields"].asString();
  c.pident     = root["pident"].asFloat();
  c.qcovs      = root["qcovs"].asFloat();
  c.t6         = root["t6"].asString();
  c.input      = root["input"].asString();
  c.gi_type    = root["gi-type"].asString();
  c.blast_db   = root["blast-db"].asString();
  c.blast_type = root["blast-type"].asString();
  c.format     = c.input.substr(c.input.find_last_of(".") + 1);
  c.fastq      = c.input;
  c.input.erase(c.input.find_last_of("."));
  c.clean      = c.input  + ".nohits-" + c.blast_type + "-" + c.blast_db + "." + c.format;
  c.filter     = c.input + ".filter_pident-" + to_string(c.pident) + "_qcovs-" + to_string(c.qcovs) + "." + c.format;

  return c;
}

unordered_set<string> loadCSV(Config c) {
  string qseqid, null1, sseqid, null3, null4, null5, null6, null7, null8, null9, null10, null11, null12, sscinames;
  float pident, qcovs;
  unordered_set<string> to_remove;

  CSVReader<16, trim_chars<' '>, no_quote_escape<'\t'>, single_line_comment<'#'> > csv (c.t6);

  while(csv.read_row(qseqid, null1, sseqid, null3, null4, null5, null6, null7, null8, null9, null10, pident, null11, null12, qcovs, sscinames))
  {
    if (pident >= c.pident && qcovs >= c.qcovs )
    {
        to_remove.insert("@" + qseqid);
    }
  }

  return to_remove;
}
