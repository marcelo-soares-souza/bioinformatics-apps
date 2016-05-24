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

int main(int argc,char **argv) {

  if (argc < 2)
  {
    cout << "Usage: " << argv[0] << " [INFO FILE]" << "\n";
    return EXIT_FAILURE;
  }

  string header;
  int count = 0;
  FastQ f;

  Config c = loadConfig(argv[1]);
  unordered_set<string> to_remove = loadCSV(c);

  string filename = c.input.erase(c.input.find_last_of("."));

  ifstream in(filename + ".fastq", ios::in);
  ofstream out_filter(filename + ".filter_pident-" + to_string(c.pident) + "_qcovs-" + to_string(c.qcovs) + ".fastq");
  ofstream out_clean(filename  + ".nohits-" + c.blast_type + "-" + c.blast_db + ".fastq");

  if(!in.is_open())
  {
    return EXIT_FAILURE;
  }

  cout << "\nUsing " << c.input << " and " << c.t6 << "\nProcessing..." << "\n";

  while(!in.eof())
  {
    if(!getline(in, f.name,'\n')) break;
    if(!getline(in, f.sequence,'\n')) break;
    if(!getline(in, f.info,'\n')) break;
    if(!getline(in, f.quality,'\n')) break;

    istringstream iss(f.name);

    iss >> header;

    unordered_set<string>::const_iterator got = to_remove.find (header);

    if (got != to_remove.end())
    {
        out_filter << "@" << f.name << "\n";
        out_filter << f.sequence << "\n";
        out_filter << f.info << "\n";
        out_filter << f.quality << "\n";
        count++;
    }
    else
    {
        out_clean << "@" << f.name << "\n";
        out_clean << f.sequence << "\n";
        out_clean << f.info << "\n";
        out_clean << f.quality << "\n";
    }
  }

  cout << "Found: " << count << "\n";

  in.close();
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
