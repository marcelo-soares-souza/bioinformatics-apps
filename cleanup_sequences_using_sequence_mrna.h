// (C) 2016 Marcelo Soares Souza <marcelo.soares@colaborador.embrapa.br>
// This program is licensed under a LGPLv3 License.

using namespace std;
using namespace io;
using namespace tr1;
using namespace boost;

typedef struct fastq {
  string name;
  string sequence;
  string info;
  string quality;
  bool remove = false;

} FastQ;

map<string,string> loadCSV(string csv_file);
