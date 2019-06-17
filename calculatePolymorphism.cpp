/****************************************************************************
 * calculatePolymorphism.cpp                                                *
 * Written by Patrick Reilly                                                *
 * Version 1.0 written 2016/10/21                                           *
 * Version 1.1 written 2017/03/06 (Added base counts to debug output)       *
 * Version 1.2 written 2017/05/26 (Finally implemented inbred pi)           *
 * Version 1.3 written 2018/08/09 (Fixed bug ignoring softmasked bases)     *
 * Version 1.4 written 2018/11/08 (Omit position may output weight instead) *
 * Version 1.5 written 2019/05/06 (Option to pass FOFN instead of pos args) *
 *                                                                          *
 * Description:                                                             *
 *                                                                          *
 * Syntax: calculatePolymorphism [list of pseudoreference FASTAs]           *
 ****************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <getopt.h>
#include <cctype>
#include <vector>

//Define constants for getopt:
#define no_argument 0
#define required_argument 1
#define optional_argument 2

//Version:
#define VERSION "1.5"

//Define number of bases:
#define NUM_BASES 4

//Usage/help:
#define USAGE "calculatePolymorphism\nUsage:\n calculatePolymorphism [options] [list of pseudoreference FASTAs]\n Options:\n  --help,-h:\t\tOutput this documentation\n  --version,-v:\t\tOutput the version number\n  --fofn,-f:\t\tPass a file of filenames, rather than listing filenames\n  --segregating_sites,-s:\tOutput whether or not the site is segregating\n  --inbred,-i:\t\tAssume inbred input sequences\n  --prng_seed,-p:\t\tSet pseudo-random number generator seed for allele choice if -i is set\n  --usable_fraction,-u:\tFourth column represents fraction of unmasked bases\n  --debug,-d:\t\tOutput extra debugging info\n"

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

void processScaffold(vector<string> &FASTA_headers, vector<string> &FASTA_sequences, bool debug, bool segsites, bool inbred, bool usable) {
   cerr << "Processing scaffold " << FASTA_headers[0].substr(1) << " of length " << FASTA_sequences[0].length() << endl;
   //Do all the processing for this scaffold:
   //Polymorphism estimator: Given base frequencies at site:
   //\hat{\pi} = \(\frac{n}{n-1}\)\sum_{i=1}^{3}\sum_{j=i+1}^{4} 2\hat{p_{i}}\hat{p_{j}}
   //For biallelic sites, this reduces to the standard estimator: \(\frac{n}{n-1}\)2\hat{p}\hat{q}, since
   //p_{3} = p_{4} = 0
   unsigned long num_sequences = FASTA_sequences.size();
   unsigned long scaffold_length = FASTA_sequences[0].length();
   for (unsigned long i = 0; i < scaffold_length; i++) {
      double pi_hat = 0.0; //Accumulate the current polymorphism estimate in this variable
      unsigned long base_frequency[5] = {0, 0, 0, 0, 0}; //Store the count of A, C, G, T, N for each base
      for (unsigned long j = 0; j < num_sequences; j++) {
         switch (FASTA_sequences[j][i]) {
            case 'A':
            case 'a':
               base_frequency[0] += inbred ? 1 : 2;
               break;
            case 'C':
            case 'c':
               base_frequency[1] += inbred ? 1 : 2;
               break;
            case 'G':
            case 'g':
               base_frequency[2] += inbred ? 1 : 2;
               break;
            case 'K': //G/T het site
            case 'k':
               if (inbred) {
                  base_frequency[rand() <= (RAND_MAX-1)/2 ? 2 : 3]++;
               } else {
                  base_frequency[2]++;
                  base_frequency[3]++;
               }
               break;
            case 'M': //A/C het site
            case 'm':
               if (inbred) {
                  base_frequency[rand() <= (RAND_MAX-1)/2 ? 0 : 1]++;
               } else {
                  base_frequency[0]++;
                  base_frequency[1]++;
               }
               break;
            case 'R': //A/G het site
            case 'r':
               if (inbred) {
                  base_frequency[rand() <= (RAND_MAX-1)/2 ? 0 : 2]++;
               } else {
                  base_frequency[0]++;
                  base_frequency[2]++;
               }
               break;
            case 'S': //C/G het site
            case 's':
               if (inbred) {
                  base_frequency[rand() <= (RAND_MAX-1)/2 ? 1 : 2]++;
               } else {
                  base_frequency[1]++;
                  base_frequency[2]++;
               }
               break;
            case 'T':
            case 't':
               base_frequency[3] += inbred ? 1 : 2;
               break;
            case 'W': //A/T het site
            case 'w':
               if (inbred) {
                  base_frequency[rand() <= (RAND_MAX-1)/2 ? 0 : 3]++;
               } else {
                  base_frequency[0]++;
                  base_frequency[3]++;
               }
               break;
            case 'Y': //C/T het site
            case 'y':
               if (inbred) {
                  base_frequency[rand() <= (RAND_MAX-1)/2 ? 1 : 3]++;
               } else {
                  base_frequency[1]++;
                  base_frequency[3]++;
               }
               break;
            case 'N':
            case 'n':
               base_frequency[4] += inbred ? 1 : 2;
               break;
            case '-':
               base_frequency[4] += inbred ? 1 : 2;
               break;
            default: //Assume that any case not handled here is an N
               base_frequency[4] += inbred ? 1 : 2;
               break;
         }
      }
      double p_hat[4];
      unsigned long nonN_bases = base_frequency[0] + base_frequency[1] + base_frequency[2] + base_frequency[3];
      for (unsigned long j = 0; j < NUM_BASES; j++) {
         p_hat[j] = (double)base_frequency[j]/(double)nonN_bases;
      }
      for (unsigned long j = 0; j < NUM_BASES-1; j++) {
         for (unsigned long k = j+1; k < NUM_BASES; k++) {
            pi_hat += 2*p_hat[j]*p_hat[k];
         }
      }
      double usable_fraction = (double)nonN_bases/(double)(nonN_bases+base_frequency[4]);
      if (nonN_bases <= 1) { //The estimator doesn't work for n <= 1, so make sure this base gets ignored by the windowing script
         if (!usable) { // If we don't want to output the usable fraction, just output 1
            cout << FASTA_headers[0].substr(1) << '\t' << i+1 << '\t' << "NA" << '\t' << 1 << endl;
         } else {
            cout << FASTA_headers[0].substr(1) << '\t' << i+1 << '\t' << "NA" << '\t' << usable_fraction << endl;
         }
      } else {
         if (segsites) {
            if (!usable) { // If we don't want to output the usable fraction, just output 1
               cout << FASTA_headers[0].substr(1) << '\t' << i+1 << '\t' << (pi_hat > 0.0 ? 1 : 0) << '\t' << 0 << endl;
            } else {
               cout << FASTA_headers[0].substr(1) << '\t' << i+1 << '\t' << (pi_hat > 0.0 ? 1 : 0) << '\t' << usable_fraction << endl;
            }
         } else {
            cout << FASTA_headers[0].substr(1) << '\t' << i+1 << '\t';
            if (!usable) { // If we don't want to output the usable fraction, just output 1
               cout << (double)nonN_bases/(double)(nonN_bases-1)*pi_hat << '\t' << 0;
            } else {
               cout << (double)nonN_bases/(double)(nonN_bases-1)*pi_hat << '\t' << usable_fraction;
            }
            if (debug) {
               cout << '\t' << nonN_bases << '\t' << pi_hat;
               cout << '\t' << base_frequency[0] << '\t' << base_frequency[1] << '\t' << base_frequency[2];
               cout << '\t' << base_frequency[3] << '\t' << base_frequency[4];
            }
            cout << endl;
         }
      }
   }
}

int main(int argc, char **argv) {
   //Variables for processing the FASTAs:
   vector<string> input_FASTA_paths;
   vector<ifstream*> input_FASTAs;
   
   //Option for debugging:
   bool debug = 0;
   //Option for input of FASTA file paths:
   string input_fofn = "";
   //Option for inbred lines:
   bool inbred = 0;
   //Option to output segregating sites:
   bool segsites = 0;
   //Seed for PRNG for choosing alleles at heterozygous sites in inbred strains:
   unsigned int prng_seed = 42;
   //Option to output fraction of usable sites:
   bool usable = 0;
   
   //Variables for getopt_long:
   int optchar;
   int structindex = 0;
   extern int optind;
   //Create the struct used for getopt:
   const struct option longoptions[] {
      {"fofn", required_argument, 0, 'f'},
      {"segregating_sites", no_argument, 0, 's'},
      {"inbred", no_argument, 0, 'i'},
      {"prng_seed", required_argument, 0, 'p'},
      {"usable_fraction", no_argument, 0, 'u'},
      {"debug", no_argument, 0, 'd'},
      {"version", no_argument, 0, 'v'},
      {"help", no_argument, 0, 'h'}
   };
   //Read in the options:
   while ((optchar = getopt_long(argc, argv, "f:sip:udvh", longoptions, &structindex)) > -1) {
      switch(optchar) {
         case 'f':
            cerr << "Taking input from FOFN " << optarg << endl;
            input_fofn = optarg;
            break;
         case 's':
            cerr << "Only outputting segregating sites, not polymorphism." << endl;
            segsites = 1;
            break;
         case 'i':
            cerr << "Assuming all pseudoreferences provided are from inbred lines." << endl;
            inbred = 1;
            break;
         case 'p':
            cerr << "Using PRNG seed " << optarg << " for choosing alleles in inbred lines." << endl;
            prng_seed = atoi(optarg);
            break;
         case 'u':
            cerr << "Outputting fraction of usable sites rather than omit column" << endl;
            usable = 1;
            break;
         case 'd':
            cerr << "Debugging mode enabled." << endl;
            debug = 1;
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
         cerr << "Added input FASTA file " << argv[optind] << " to the vector." << endl;
      }
      input_FASTA_paths.push_back(argv[optind++]);
   }
   
   //Set the seed for the PRNG:
   srand(prng_seed);
   
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
            processScaffold(FASTA_headers, FASTA_sequences, debug, segsites, inbred, usable);
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
   processScaffold(FASTA_headers, FASTA_sequences, debug, segsites, inbred, usable);
   
   //Close the input FASTAs:
   closeFASTAs(input_FASTAs);
   
   return 0;
}
