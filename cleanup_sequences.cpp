#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <iterator>
#include <tr1/unordered_set>

#include "csv.h"

using namespace std;
using namespace io;
using namespace tr1;

int main(int argc,char **argv) {

  string qseqid, null0, sseqid, null3, null4, null5, null6, null7, null8, null9, null10, null11, null12, sscinames, sequence;
  int count = 0;
  float pident, qcovs;
  unordered_set<string> to_remove;

  CSVReader<16, trim_chars<' '>, no_quote_escape<'\t'>, single_line_comment<'#'> > csv ("teste.t6");
  LineReader fastq("teste.fastq");

  while(csv.read_row(qseqid, null0, sseqid, null3, null4, null5, null6, null7, null8, null9, null10, pident, null11, null12, qcovs, sscinames))
  {
    if (pident >= 95.0) {
        to_remove.insert("@" + qseqid);
    }
  }

  string name, seq, name2, qual;

  ifstream in(argv[1],ios::in);

  if(!in.is_open()) return EXIT_FAILURE;

  cout << "Iniciando Busca" << endl;

  while(!in.eof()) {
    if(!getline(in,name,'\n')) break;
    if(!getline(in,seq,'\n')) break;
    if(!getline(in,name2,'\n')) break;
    if(!getline(in,qual,'\n')) break;

    istringstream iss(name);

    iss >> sequence;

    unordered_set<string>::const_iterator got = to_remove.find (sequence);

    if (got != to_remove.end())
    {
        // cout << "Encontrei: " << sequence << endl;
        count++;
    }
  }

  in.close();

  cout << "Encontrei: " << count << endl;

  return 0;
}
