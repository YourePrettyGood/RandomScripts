/*****************************************************************************
 * softmaskFromHardmask.cpp                                                  *
 * Written by Patrick Reilly                                                 *
 * Version 1.0 written 2018/06/25                                            *
 * Description:                                                              *
 *  Given an unmasked FASTA and a hardmasked FASTA, generates a softmasked   *
 *  FASTA equivalent to that produced by the -xsmall option of RepeatMasker. *
 *                                                                           *
 * Syntax: softmaskFromHardmask [options] <unmasked FASTA> <hardmasked FASTA>*
 *  unmasked FASTA:       Path to unmasked FASTA                             *
 *  hardmasked FASTA:     Path to hard-masked FASTA                          *
 *  -d,--debug            Toggle debugging output                            *
 *****************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <getopt.h>
#include <cctype>

//Define constants for getopt:
#define no_argument 0
#define required_argument 1
#define optional_argument 2

//Version:
#define version "1.0"

//Usage/help:
#define usage "softmaskFromHardmask\nUsage:\n softmaskFromHardmask [options] <unmasked FASTA> <hardmasked FASTA>\n Mandatory arguments:\n  unmasked FASTA\t\tPath to unmasked FASTA\n  hardmasked FASTA\t\t\tPath to hardmasked FASTA\n\n Options:\n  -d,--debug\t\tToggle debugging output\n\n Description:\n  Outputs a soft-masked FASTA given an unmasked and hard-masked FASTA.\n"

using namespace std;

int main(int argc, char **argv) {
   //Variables for processing the FASTAs:
   string unmasked_FASTA = "", hardmasked_FASTA = "";
   string unmasked_line = "", hardmasked_line = "";
   bool debug = 0;

   //Variables for getopt_long:
   int optchar;
   int structindex = 0;
   extern int optind;
   //Create the struct used for getopt:
   const struct option longoptions[] {
      {"version", no_argument, 0, 'v'},
      {"help", no_argument, 0, 'h'},
      {"debug", no_argument, 0, 'd'}
   };
   //Read in the options:
   while ((optchar = getopt_long(argc, argv, "dvh", longoptions, &structindex)) > -1) {
      switch(optchar) {
         case 'v':
            cerr << "softmaskFromHardmask version " << version << endl;
            return 0;
            break;
         case 'h':
            cerr << usage;
            return 0;
            break;
         case 'd':
            cerr << "Debugging turned on" << endl;
            debug = 1;
            break;
         default:
            cerr << "Unknown option " << (unsigned char)optchar << " supplied." << endl;
            cerr << usage;
            return 1;
            break;
      }
   }
   //Read in the positional arguments:
   if (optind < argc) {
      unmasked_FASTA = argv[optind++];
   }
   if (optind < argc) {
      hardmasked_FASTA = argv[optind++];
   }
   if (optind < argc) {
      cerr << "Ignoring extra positional arguments starting at " << argv[optind++] << endl;
   }
   
   //Do some error checking on the arguments:
   ifstream unmasked, hardmasked;
   unmasked.open(unmasked_FASTA);
   if (!unmasked) {
      cerr << "Error opening unmasked FASTA file." << endl;
      return 2;
   }
   getline(unmasked, unmasked_line);
   if (!unmasked) {
      cerr << "Error reading unmasked FASTA file." << endl;
      unmasked.close();
      return 3;
   }
   if (unmasked_line[0] != '>') {
      cerr << "Unmasked FASTA is not properly formatted FASTA: Does not start with >." << endl;
      unmasked.close();
      return 4;
   }
   unmasked_line = "";
   unmasked.close();
   hardmasked.open(hardmasked_FASTA);
   if (!hardmasked) {
      cerr << "Error opening hard-masked FASTA file." << endl;
      return 2;
   }
   getline(hardmasked, hardmasked_line);
   if (!hardmasked) {
      cerr << "Error reading hard-masked FASTA file." << endl;
      hardmasked.close();
      return 3;
   }
   if (hardmasked_line[0] != '>') {
      cerr << "Hard-masked FASTA is not properly formatted FASTA: Does not start with >." << endl;
      hardmasked.close();
      return 4;
   }

   hardmasked_line = "";
   hardmasked.close();
   
   //Now perform the main processing:
   unmasked.open(unmasked_FASTA);
   hardmasked.open(hardmasked_FASTA);
   getline(unmasked, unmasked_line);
   getline(hardmasked, hardmasked_line);
   while (unmasked && hardmasked) {
      if (unmasked_line.length() == 0 || hardmasked_line.length() == 0) {
         break;
      }
      if (unmasked_line[0] != '>' && hardmasked_line[0] != '>') {
         if (unmasked_line.length() != hardmasked_line.length()) {
            cerr << "Unmasked and hard-masked FASTAs do not wrap at same line length." << endl;
            if (debug) {
               cerr << unmasked_line << endl;
               cerr << hardmasked_line << endl;
            }
            unmasked.close();
            hardmasked.close();
            return 5;
         }
         unsigned long int unmaskedlength = unmasked_line.length();
         string softmasked_line = unmasked_line;
         for (unsigned long int i = 0; i < unmaskedlength; i++) {
            //Force to uppercase for comparison:
            int unmaskedbase = toupper(unmasked_line[i]);
            int hardmaskedbase = toupper(hardmasked_line[i]);
            if (unmaskedbase != 'N' && hardmaskedbase == 'N') {
               softmasked_line[i] = tolower(unmaskedbase);
            } else {
               softmasked_line[i] = unmaskedbase;
            }
         }
         cout << softmasked_line << endl;
      } else {
         cout << unmasked_line << endl;
      }
      getline(unmasked, unmasked_line);
      getline(hardmasked, hardmasked_line);
   }
   unmasked.close();
   hardmasked.close();
   
   return 0;
}
