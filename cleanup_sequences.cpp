#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <iterator>
#include <tr1/unordered_set>

#include <jsoncpp/json/json.h>
#include "csv.h"

using namespace std;
using namespace io;
using namespace tr1;

typedef struct config {
  string fields;
  float pident;
  float qcovs;
  string t6;
  string input;
  string gi_type;
  string blast_db;
  string blast_type;

} Config;

typedef struct fastq {
  string name;
  string sequence;
  string info;
  string quality;
} FastQ;

Config loadConfig(string filename);
unordered_set<string> loadCSV(Config c);

int main(int argc,char **argv) {

  string header;
  int count = 0;
  FastQ f;

  Config c = loadConfig(argv[1]);
  unordered_set<string> to_remove = loadCSV(c);

  ifstream in(c.input, ios::in);

  if(!in.is_open())
  {
    return EXIT_FAILURE;
  }

  cout << "Using " << c.input << " and " << c.t6 << endl;

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
        cout << "@" << f.name << endl;
        cout << f.sequence << endl;
        cout << f.info << endl;
        cout << f.quality << endl;

        count++;
    }
  }

  in.close();

  cout << "Found: " << count << endl;

  return 0;
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
