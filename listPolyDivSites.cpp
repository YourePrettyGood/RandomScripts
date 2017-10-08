/****************************************************************************
 * listPolyDivSites.cpp                                                     *
 * Written by Patrick Reilly                                                *
 * Version 1.0 written 2016/05/03                                           *
 * Description:                                                             *
 *  Lists the differences between two IUPAC-degenerated diploid references  *
 *  either in terms of polymorphisms or divergent sites.  Sites with Ns are *
 *  not counted as polymorphism or divergence.                              *
 *                                                                          *
 * Syntax: SNPratePerBase [options] <reference FASTA> <query FASTA>         *
 *  reference FASTA:       Path to FASTA of reference diploid               *
 *  query FASTA:           Path to FASTA of query diploid                   *
 *  -p:                    Only list polymorphic sites in query             *
 *  -d:                    Only list divergent sites in query               *
 *  -n:                    Include a column that says if site is N in ref or*
 *                         query                                            *
 ****************************************************************************/

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
#define usage "listPolyDivSites\nUsage:\n listPolyDivSites [options] <reference FASTA> <query FASTA>\n Options:\n  --polymorphisms_only,-p\tOnly list polymorphic sites in query\n  --divergences_only,-d\t\tOnly list divergent sites in query\n  --list_n,-n\t\t\tInclude a column for if ref or query has an N\n\n Mandatory arguments:\n  reference FASTA\t\tPath to FASTA of reference diploid\n  query FASTA\t\t\tPath to FASTA of query diploid\n\n Description:\n  Lists the differences between two IUPAC-degenerated diploid references\n  either in terms of polymorphisms or divergent sites.\n  Sites with Ns do not count as polymorphism or divergence.\n"

using namespace std;

int main(int argc, char **argv) {
   //Variables for processing the FASTAs:
   bool polymorphisms_only = 0;
   bool divergences_only = 0;
   bool list_Ns = 0;
   string reference_FASTA = "", query_FASTA = "";
   string reference_line = "", query_line = "";
   string scaffold_name = "";
   unsigned long int scaffold_position = 1;

   //Variables for getopt_long:
   int optchar;
   int structindex = 0;
   extern int optind;
   //Create the struct used for getopt:
   const struct option longoptions[] {
      {"version", no_argument, 0, 'v'},
      {"help", no_argument, 0, 'h'},
      {"polymorphism_only", no_argument, 0, 'p'},
      {"divergences_only", no_argument, 0, 'd'},
      {"list_n", no_argument, 0, 'n'}
   };
   //Read in the options:
   while ((optchar = getopt_long(argc, argv, "pdnvh", longoptions, &structindex)) > -1) {
      switch(optchar) {
         case 'v':
            cerr << "listPolyDivSites version " << version << endl;
            return 0;
            break;
         case 'h':
            cerr << usage;
            return 0;
            break;
         case 'p':
            polymorphisms_only = 1;
            break;
         case 'd':
            divergences_only = 1;
            break;
         case 'n':
            list_Ns = 1;
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
      reference_FASTA = argv[optind++];
   }
   if (optind < argc) {
      query_FASTA = argv[optind++];
   }
   if (optind < argc) {
      cerr << "Ignoring extra positional arguments starting at " << argv[optind++] << endl;
   }
   
   //Do some error checking on the arguments:
   ifstream reference, query;
   reference.open(reference_FASTA);
   if (!reference) {
      cerr << "Error opening reference FASTA file." << endl;
      return 2;
   }
   getline(reference, reference_line);
   if (!reference) {
      cerr << "Error reading reference FASTA file." << endl;
      reference.close();
      return 3;
   }
   if (reference_line[0] != '>') {
      cerr << "Reference FASTA is not properly formatted FASTA: Does not start with >." << endl;
      reference.close();
      return 4;
   }
   reference_line = "";
   reference.close();
   query.open(query_FASTA);
   if (!query) {
      cerr << "Error opening query FASTA file." << endl;
      return 2;
   }
   getline(query, query_line);
   if (!query) {
      cerr << "Error reading query FASTA file." << endl;
      query.close();
      return 3;
   }
   if (query_line[0] != '>') {
      cerr << "Query FASTA is not properly formatted FASTA: Does not start with >." << endl;
      query.close();
      return 4;
   }

   query_line = "";
   query.close();
   
   //Now perform the main processing:
   reference.open(reference_FASTA);
   query.open(query_FASTA);
   getline(reference, reference_line);
   getline(query, query_line);
   while (reference && query) {
      if (reference_line.length() == 0 || query_line.length() == 0) {
         break;
      }
      if (reference_line[0] != '>' && query_line[0] != '>') {
         if (reference_line.length() != query_line.length()) {
            cerr << "Reference and query FASTAs do not wrap at same line length." << endl;
            reference.close();
            query.close();
            return 5;
         }
         unsigned long int reflength = reference_line.length();
         for (unsigned long int i = 0; i < reflength; i++) {
            //Ignore case:
            int refbase = toupper(reference_line[i]);
            int querybase = toupper(query_line[i]);
            //Assumption:
            //Reference does not contain degenerate bases
            
            if (refbase == 'N' || querybase == 'N') { //List masked sites as not variant
               if (list_Ns) {
                  cout << scaffold_name << '\t' << scaffold_position << '\t' << "0" << '\t' << "1" << endl;
               } else {
                  cout << scaffold_name << '\t' << scaffold_position << '\t' << "0" << endl;
               }
            } else if (polymorphisms_only) { //If the polymorphisms_only flag is set, only list query polymorphisms as variant
               if (refbase != querybase && querybase != 'A' && querybase != 'C' && querybase != 'G' && querybase != 'T') {
                  if (list_Ns) {
                     cout << scaffold_name << '\t' << scaffold_position << '\t' << "1" << '\t' << "0" << endl;
                  } else {
                     cout << scaffold_name << '\t' << scaffold_position << '\t' << "1" << endl;
                  }
               } else { //List invariant sites as not variant
                  if (list_Ns) {
                     cout << scaffold_name << '\t' << scaffold_position << '\t' << "0" << '\t' << "0" << endl;
                  } else {
                     cout << scaffold_name << '\t' << scaffold_position << '\t' << "0" << endl;
                  }
               }
            } else if (divergences_only) { //If the divergences_only flag is set, only list query divergent sites as variant
               if (refbase != querybase && (querybase == 'A' || querybase == 'C' || querybase == 'G' || querybase == 'T')) {
                  if (list_Ns) {
                     cout << scaffold_name << '\t' << scaffold_position << '\t' << "1" << '\t' << "0" << endl;
                  } else {
                     cout << scaffold_name << '\t' << scaffold_position << '\t' << "1" << endl;
                  }
               } else { //List invariant sites as not variant
                  if (list_Ns) {
                     cout << scaffold_name << '\t' << scaffold_position << '\t' << "0" << '\t' << "0" << endl;
                  } else {
                     cout << scaffold_name << '\t' << scaffold_position << '\t' << "0" << endl;
                  }
               }
            } else { //Else list all sites where query base != ref base as variant
               if (refbase != querybase) {
                  if (list_Ns) {
                     cout << scaffold_name << '\t' << scaffold_position << '\t' << "1" << '\t' << "0" << endl;
                  } else {
                     cout << scaffold_name << '\t' << scaffold_position << '\t' << "1" << endl;
                  }
               } else { //List invariant sites as not variant
                  if (list_Ns) {
                     cout << scaffold_name << '\t' << scaffold_position << '\t' << "0" << '\t' << "0" << endl;
                  } else {
                     cout << scaffold_name << '\t' << scaffold_position << '\t' << "0" << endl;
                  }
               }
            }
            scaffold_position++;
         }
      } else {
         //Get the scaffold name from the FASTA header line:
         scaffold_name = query_line.substr(1,string::npos);
         scaffold_position = 1;
      }
      getline(reference, reference_line);
      getline(query, query_line);
   }
   reference.close();
   query.close();
   
   return 0;
}