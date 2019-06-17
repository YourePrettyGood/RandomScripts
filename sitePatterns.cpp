/****************************************************************************
 * sitePatterns.cpp                                                         *
 * Written by Patrick Reilly                                                *
 * Version 1.0 written 2017/01/16                                           *
 * Version 1.1 written 2019/05/30 Softmask fix, FOFN input, and debugging   *
 *                                                                          *
 * Description:                                                             *
 *                                                                          *
 * Syntax: sitePatterns [list of pseudoreference FASTAs]                    *
 ****************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <getopt.h>
#include <cctype>
#include <vector>
#include <map>

//Define constants for getopt:
#define no_argument 0
#define required_argument 1
#define optional_argument 2

//Version:
#define VERSION "1.1"

//Define number of bases:
#define NUM_BASES 4

//Usage/help:
#define USAGE "sitePatterns\nUsage:\n sitePatterns [list of pseudoreference FASTAs]\n"

using namespace std;

bool openFASTAs(vector<ifstream*> &input_FASTAs, vector<string> input_FASTA_paths) {
   for (auto path_iterator = input_FASTA_paths.begin(); path_iterator != input_FASTA_paths.end(); ++path_iterator) {
      ifstream *input_FASTA = new ifstream(path_iterator->c_str(), ios::in);
      if (!input_FASTA->is_open()) { //Check to make sure the input file was validly opened
         cerr << "Error opening input FASTA: " << *path_iterator << "." << endl;
         return 0;
      }
      input_FASTAs.push_back(input_FASTA);
   }
   return 1;
}

void closeFASTAs(vector<ifstream*> &input_FASTAs) {
   for (auto FASTA_iterator = input_FASTAs.begin(); FASTA_iterator != input_FASTAs.end(); ++FASTA_iterator) {
      (*FASTA_iterator)->close();
      delete *FASTA_iterator;
   }
}

bool readFASTAs(vector<ifstream*> &input_FASTAs, vector<string> &FASTA_lines) {
   unsigned long which_input_FASTA = 0;
   bool ifstream_notfail;
   for (auto FASTA_iterator = input_FASTAs.begin(); FASTA_iterator != input_FASTAs.end(); ++FASTA_iterator) {
      if ((*FASTA_iterator)->fail()) {
         cerr << "Ifstream " << which_input_FASTA+1 << " failed." << endl;
         return 0;
      }
      string FASTA_line;
      ifstream_notfail = (bool)getline(**FASTA_iterator, FASTA_line);
      if (!ifstream_notfail) {
         return ifstream_notfail;
      }
      if (FASTA_lines.size() < input_FASTAs.size()) {
         FASTA_lines.push_back(FASTA_line);
      } else {
         FASTA_lines[which_input_FASTA++] = FASTA_line;
      }
   }
   return ifstream_notfail;
}

void processScaffold(vector<string> &FASTA_headers, vector<string> &FASTA_sequences, map<string, unsigned long> &pattern_counts) {
   cerr << "Processing scaffold " << FASTA_headers[0].substr(1) << " of length " << FASTA_sequences[0].length() << endl;
   //Do all the processing for this scaffold:
   unsigned long num_sequences = FASTA_sequences.size();
   unsigned long scaffold_length = FASTA_sequences[0].length();
   for (unsigned long i = 0; i < scaffold_length; i++) {
      string site_pattern = "";
      for (unsigned long j = 0; j < num_sequences; j++) {
         switch (FASTA_sequences[j][i]) {
            case 'A':
            case 'a':
               site_pattern += "AA";
               break;
            case 'C':
            case 'c':
               site_pattern += "CC";
               break;
            case 'G':
            case 'g':
               site_pattern += "GG";
               break;
            case 'K': //G/T het site
            case 'k':
               site_pattern += "GT";
               break;
            case 'M': //A/C het site
            case 'm':
               site_pattern += "AC";
               break;
            case 'R': //A/G het site
            case 'r':
               site_pattern += "AG";
               break;
            case 'S': //C/G het site
            case 's':
               site_pattern += "CG";
               break;
            case 'T':
            case 't':
               site_pattern += "TT";
               break;
            case 'W': //A/T het site
            case 'w':
               site_pattern += "AT";
               break;
            case 'Y': //C/T het site
            case 'y':
               site_pattern += "CT";
               break;
            case 'N':
            case 'n':
            case '-':
            default: //Assume that any case not handled here is an N
               site_pattern += "NN";
               break;
         }
      }
      if (pattern_counts.count(site_pattern) > 0) {
         pattern_counts[site_pattern]++;
      } else {
         pattern_counts[site_pattern] = 1;
      }
   }
}

int main(int argc, char **argv) {
   //Variables for processing the FASTAs:
   vector<string> input_FASTA_paths;
   vector<ifstream*> input_FASTAs;
   
   //Option for debugging:
   unsigned short int debug = 0;
   //Option for input of FASTA file paths:
   string input_fofn = "";

   //Variable for storing pattern counts:
   map<string, unsigned long> pattern_counts;
   
   //Variables for getopt_long:
   int optchar;
   int structindex = 0;
   extern int optind;
   //Create the struct used for getopt:
   const struct option longoptions[] {
      {"fofn", required_argument, 0, 'f'},
      {"debug", no_argument, 0, 'd'},
      {"version", no_argument, 0, 'v'},
      {"help", no_argument, 0, 'h'}
   };
   //Read in the options:
   while ((optchar = getopt_long(argc, argv, "f:dvh", longoptions, &structindex)) > -1) {
      switch(optchar) {
         case 'f':
            cerr << "Taking input from FOFN " << optarg << endl;
            input_fofn = optarg;
            break;
         case 'd':
            cerr << "Debugging mode enabled." << endl;
            debug++;
            break;
         case 'v':
            cerr << "calculatePolymorphism version " << VERSION << endl;
            return 0;
            break;
         case 'h':
            cerr << USAGE;
            return 0;
            break;
         default:
            cerr << "Unknown option " << (unsigned char)optchar << " supplied." << endl;
            cerr << USAGE;
            return 1;
            break;
      }
   }
   //Read in the input FASTA paths:
   if (input_fofn != "") {
      string fofn_line;
      ifstream fofn;
      fofn.open(input_fofn);
      if (!fofn) {
         cerr << "Error opening file of input FASTA filenames " << input_fofn << endl;
         return 2;
      }
      while (getline(fofn, fofn_line)) {
         ifstream infile_test(fofn_line);
         if (infile_test.good()) {
            input_FASTA_paths.push_back(fofn_line);
            if (debug) {
               cerr << "Added input FASTA file " << fofn_line << " to the vector." << endl;
            }
         } else {
            cerr << "Input FASTA file " << fofn_line << " in FOFN " << input_fofn << " doesn't seem to be openable, skipping." << endl;
         }
         infile_test.close();
      }
      fofn.close();
   }
   //Read in the positional arguments:
   if (input_FASTA_paths.size() > 0 && optind < argc) {
      cerr << "Adding input FASTA files from positional arguments on top of existing set from FOFN " << input_fofn << endl;
   }
   while (optind < argc) {
      if (debug) {
         cerr << "Adding input FASTA file " << argv[optind] << " to the vector." << endl;
      }
      input_FASTA_paths.push_back(argv[optind++]);
   }
   
   //Open the input FASTAs:
   bool successfully_opened = openFASTAs(input_FASTAs, input_FASTA_paths);
   if (!successfully_opened) {
      closeFASTAs(input_FASTAs);
      return 2;
   }
   cerr << "Opened " << input_FASTAs.size() << " input FASTA files." << endl;
   
   //Set up the vector to contain each line from the n FASTA files:
   vector<string> FASTA_lines;
   FASTA_lines.reserve(input_FASTA_paths.size());
   
   //Iterate over all of the FASTAs synchronously:
   vector<string> FASTA_headers;
   vector<string> FASTA_sequences;
   while (readFASTAs(input_FASTAs, FASTA_lines)) {
      //Check if we're on a header line:
      bool all_header_lines = 1;
      bool any_header_lines = 0;
      for (auto line_iterator = FASTA_lines.begin(); line_iterator != FASTA_lines.end(); ++line_iterator) {
         all_header_lines = all_header_lines && ((*line_iterator)[0] == '>');
         any_header_lines = any_header_lines || ((*line_iterator)[0] == '>');
      }
      if (all_header_lines) {
         if (debug) {
            cerr << "Completed reading scaffold " << FASTA_lines[0].substr(1) << endl;
         }
         if (!FASTA_sequences.empty()) {
            processScaffold(FASTA_headers, FASTA_sequences, pattern_counts);
            FASTA_sequences.clear();
         }
         FASTA_headers = FASTA_lines;
         for (auto header_iterator = FASTA_headers.begin()+1; header_iterator != FASTA_headers.end(); ++header_iterator) {
            if (*header_iterator != FASTA_headers[0]) {
               cerr << "Error: FASTAs are not synchronized, headers differ." << endl;
               cerr << FASTA_headers[0] << endl;
               cerr << *header_iterator << endl;
               closeFASTAs(input_FASTAs);
               return 3;
            }
         }
      } else if (any_header_lines) {
         cerr << "Error: FASTAs are not synchronized, or not wrapped at the same length." << endl;
         closeFASTAs(input_FASTAs);
         return 4;
      } else {
         if (debug) {
            cerr << "Loading " << FASTA_lines.size() << " FASTA lines into sequences." << endl;
         }
         unsigned long sequence_index = 0;
         for (auto line_iterator = FASTA_lines.begin(); line_iterator != FASTA_lines.end(); ++line_iterator) {
            if (FASTA_sequences.size() < FASTA_lines.size()) {
               FASTA_sequences.push_back(*line_iterator);
            } else {
               FASTA_sequences[sequence_index++] += *line_iterator;
            }
         }
      }
   }
   
   //Catch any IO errors that kicked us out of the while loop:
   unsigned long which_input_FASTA = 0;
   for (auto FASTA_iterator = input_FASTAs.begin(); FASTA_iterator != input_FASTAs.end(); ++FASTA_iterator) {
      if ((*FASTA_iterator)->bad()) {
         cerr << "Error reading input FASTA: " << input_FASTA_paths[which_input_FASTA] << endl;
         cerr << "Fail bit: " << ((*FASTA_iterator)->rdstate() & ifstream::failbit) << " Bad bit: " << (*FASTA_iterator)->bad() << " EOF bit: " << (*FASTA_iterator)->eof() << endl;
         closeFASTAs(input_FASTAs);
         return 5;
      }
      which_input_FASTA++;
   }
   //If no errors kicked us out of the while loop, process the last scaffold:
   processScaffold(FASTA_headers, FASTA_sequences, pattern_counts);
   
   //Close the input FASTAs:
   closeFASTAs(input_FASTAs);
   
   //Output the pattern counts:
   for (auto pattern_iterator = pattern_counts.begin(); pattern_iterator != pattern_counts.end(); ++pattern_iterator) {
      cout << pattern_iterator->first << '\t' << pattern_iterator->second << endl;
   }
   
   return 0;
}
