// (C) 2016 Marcelo Soares Souza <marcelo.soares@colaborador.embrapa.br>
// This program is licensed under a LGPLv3 License.

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
  string format;
  string fastq;
  string clean;
  string filter;

} Config;

typedef struct fastq {
  string name;
  string sequence;
  string info;
  string quality;

} FastQ;

Config loadConfig(string filename);
unordered_set<string> loadCSV(Config c);
