/****************************************************************************
 * nonOverlappingWindows.cpp                                                *
 * Written by Patrick Reilly                                                *
 * Version 1.0 written 2016/07/17                                           *
 * Description:                                                             *
 *  Calculates the mean of a statistic over non-overlapping windows of      *
 *  user-defined length, and can adjust the denominator of the mean based   *
 *  on a filter column.                                                     *
 *  For instance, when calculating polymorphism in a window, we must        *
 *  exclude sites with Ns from the denominator, and these Ns can be         *
 *  indicated in the filter column by a "1" if present, "0" if absent.      *
 *  Input is a TSV consisting of 3 (4 if a filter column is desired) columns*
 *  identified as scaffold, position, and statistic, respectively.          *
 *  If the scaffold length is not an integral multiple of the window size,  *
 *  the final window's mean is scaled according to its size.                *
 *                                                                          *
 * Syntax: nonOverlappingWindows [options]                                  *
 *  -i:     Path to the input TSV of per-scaffold per-site statistics       *
 *  -o:     Path to desired output TSV of per-scaffold per-window statistics*
 *  -w:     Size of the non-overlapping windows                             *
 *  -n:     Account for a fourth "filter" column to filter certain sites    *
 *          by omitting them from the numerator and denominator.            *
 ****************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <getopt.h>
#include <cctype>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <stdexcept>

//Define constants for getopt:
#define no_argument 0
#define required_argument 1
#define optional_argument 2

//Version:
#define version "1.0"

//Usage/help:
#define usage "nonOverlappingWindows\nUsage:\n nonOverlappingWindows [options]\n Options:\n  --input_tsv,-i\tPath to input TSV (default: STDIN)\n  --output_tsv,-o\tPath to output TSV (default: STDOUT)\n  --omit_n,-n\tOmit sites indicated in the filter (4th) column\n  --window_size,-w\tSize of the non-overlapping windows\n\n Description:\n  Calculates the mean of a statistic over non-overlapping windows across scaffolds in a genome.\n  Sites may be omitted from the average.\n  Input is a 3- or 4-column TSV consisting of scaffold name, position, statistic, and a filter column.\n  If the filter column is 1, the row is omitted from the average.\n  If the scaffold length is not an integral multiple of the window size,\n  the last window's average is scaled appropriately.\n"

using namespace std;

vector<string> splitString(string line_to_split, char delimiter) {
   vector<string> line_vector;
   string element;
   istringstream line_to_split_stream(line_to_split);
   while (getline(line_to_split_stream, element, delimiter)) {
      line_vector.push_back(element);
   }
   return line_vector;
}

string calcDepths(vector<double> &statistics, string scaffold, unsigned long window_size, vector<bool> &N_list) {
   string output = "";
   if (statistics.size() == 0) { //Fast skip for size 0 scaffolds
      return output;
   }
   double sum = 0.0;
   unsigned long window_start = 1;
   unsigned long last_window;
   unsigned long N_count = 0;
   unsigned long current_position = 1;
   for (auto stat_iterator = statistics.begin(); stat_iterator != statistics.end(); ++stat_iterator) {
      if (current_position % window_size == 1) {
         window_start = current_position;
         output += scaffold + '\t' + to_string(current_position) + '\t';
      }
      sum += *stat_iterator;
      //If we don't want to normalize by omitting sites, call with empty N_list vector:
      if (N_list.size() > 0 && N_list[current_position-1]) {
         N_count += N_list[current_position-1];
      }
      if (current_position % window_size == 0) {
         if (window_size-N_count > 0) {
            output += to_string(sum/(window_size-N_count)) + '\n';
         } else {
            output += "NA\n";
         }
         sum = 0.0;
         N_count = 0;
      }
      ++current_position;
   }
   //If the scaffold does not contain an integral number of windows:
   if (statistics.size() % window_size > 0) {
      last_window = statistics.size() - window_start + 1;
      output += to_string(sum/(last_window-N_count)) + '\n';
   }
   return output;
}

int main(int argc, char **argv) {
   //Variables for processing the TSV:
   bool omit_Ns = 0;
   unsigned long int window_size = 10000; //Default window size is 10kb
   string input_tsv = "", output_tsv = "";
   bool use_cin = 1, use_cout = 1; //Default to reading from cin and outputting to cout
   string input_line = "", output_line = "";
   string scaffold_name = "";

   //Variables for getopt_long:
   int optchar;
   int structindex = 0;
   extern int optind;
   //Create the struct used for getopt:
   const struct option longoptions[] {
      {"version", no_argument, 0, 'v'},
      {"help", no_argument, 0, 'h'},
      {"window_size", required_argument, 0, 'w'},
      {"input_tsv", required_argument, 0, 'i'},
      {"output_tsv", required_argument, 0, 'o'},
      {"omit_n", no_argument, 0, 'n'}
   };
   //Read in the options:
   while ((optchar = getopt_long(argc, argv, "i:o:w:nvh", longoptions, &structindex)) > -1) {
      switch(optchar) {
         case 'v':
            cerr << "nonOverlappingWindows version " << version << endl;
            return 0;
            break;
         case 'h':
            cerr << usage;
            return 0;
            break;
         case 'i':
            input_tsv = optarg;
            break;
         case 'o':
            output_tsv = optarg;
            break;
         case 'w':
            window_size = atol(optarg);
            break;
         case 'n':
            omit_Ns = 1;
            break;
         default:
            cerr << "Unknown option " << (unsigned char)optchar << " supplied." << endl;
            cerr << usage;
            return 1;
            break;
      }
   }

   //Ignore positional arguments
   if (optind < argc) {
      cerr << "Ignoring extra positional arguments starting at " << argv[optind++] << endl;
   }
   
   //Do some error checking on the arguments:
   ifstream input;
   if (input_tsv.length() > 0) { 
      input.open(input_tsv);
      if (!input) {
         cerr << "Error opening input TSV file." << endl;
         return 2;
      }
      use_cin = 0;
   }
   ofstream output;
   if (output_tsv.length() > 0) {
      output.open(output_tsv);
      if (!output) {
         cerr << "Error opening output TSV file." << endl;
         return 3;
      }
      use_cout = 0;
   }
   
   //Set initial state:
   string previous_scaffold = "";
   vector<bool> N_list;
   vector<double> scaffold_stats;
   
   //Now perform the main processing:
   while (getline(use_cin ? cin : input, input_line)) {
      vector<string> line_vector;
      line_vector = splitString(input_line, '\t');
      
      //Test that the TSV is valid and has the necessary columns:
      if (line_vector.size() < 3) { //Too few elements per line
         cerr << "Malformatted input TSV: Less than 3 columns." << endl;
         if (!use_cin) { 
            input.close();
         }
         if (!use_cout) {
            output.close();
         }
         return 4;
      } else if (line_vector.size() < 4 && omit_Ns) { //Need the filter column to perform omission
         cerr << "Malformatted input TSV: Missing filter column." << endl;
         if (!use_cin) {
            input.close();
         }
         if (!use_cout) {
            output.close();
         }
         return 5;
      }
      
      //Convert the column values from strings:
      string scaffold_name = "";
      double local_statistic;
      bool omit_position = 0;
      scaffold_name = line_vector[0];
      try {
         local_statistic = stod(line_vector[2]);
      } catch (const invalid_argument&) {
         cerr << "Unable to convert stat at " << line_vector[0] << " pos " << line_vector[1] << ": " << line_vector[2] << " to a double." << endl;
         throw;
      } catch (const out_of_range&) {
         cerr << "Stat at " << line_vector[0] << " pos " << line_vector[1] << ": " << line_vector[2] << " is out of the range of a double." << endl;
         throw;
      }
      if (line_vector.size() == 4) { 
         omit_position = stoi(line_vector[3]);
      }
      
      //If we're not on the same scaffold, output the window averages for the previous scaffold:
      if (scaffold_name != previous_scaffold) {
         if (use_cout) {
            cout << calcDepths(scaffold_stats, previous_scaffold, window_size, N_list);
         } else {
            output << calcDepths(scaffold_stats, previous_scaffold, window_size, N_list);
         }
         //Make sure to reset the window containers:
         scaffold_stats.clear();
         N_list.clear();
      }
      //Accumulate the window statistics:
      scaffold_stats.push_back(local_statistic);
      if (omit_Ns) {
         N_list.push_back(omit_position);
      }
      previous_scaffold = scaffold_name;
   }
   //Make sure to capture the last scaffold:
   if (use_cout) {
      cout << calcDepths(scaffold_stats, previous_scaffold, window_size, N_list);
   } else {
      output << calcDepths(scaffold_stats, previous_scaffold, window_size, N_list);
   }

   //Close the files if they were used instead of STDIN and STDOUT:
   if (!use_cin) {
      input.close();
   }
   if (!use_cout) {
      output.close();
   }
   
   return 0;
}
